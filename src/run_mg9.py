"""MG9 driver (spec/MG9_proposal.md draft 3, registered at 7810b8b). Deterministic only. Resumable at run level.

Stages (proposal §9 priority order):
    python src/run_mg9.py p1        A1, code check, V2, P1 (results/MG9/p1.json)
    python src/run_mg9.py sweep     pattern persistence sweep on D02, 30.50 -> 28.50 by 0.02
    python src/run_mg9.py bench     benchmark recheck of b_s (t = 600) and b_r (t = 800)
    python src/run_mg9.py grid      A3 fronts on the 31-point D02 grid
    python src/run_mg9.py refine    edge refinement, then robustness repeats at the final edges
    python src/run_mg9.py d01       mesh comparison on D01
Every integration is recorded in results/MG9/ledger.json; the run cap (230) and cutoff (10 core-h) are enforced
before each new integration (exit code 2 on a stop).
"""
import os, sys, json, time, hashlib
import numpy as np
import scipy.sparse as sp
from scipy.integrate import solve_ivp
from scipy.signal import find_peaks

sys.path.insert(0, os.path.dirname(__file__))
import landau_mg9 as L

ROOT = os.path.join(os.path.dirname(__file__), "..")
OUT = os.path.join(ROOT, "results", "MG9")
RUNS = os.path.join(OUT, "runs")
STATES = os.path.join(OUT, "states")
RUN_CAP, CUTOFF_CORE_S = 230, 10 * 3600.0
P = dict(L.P17)
LAM = L.LAM
SWEEP_B = [round(30.50 - 0.02 * i, 2) for i in range(101)]        # 30.50 ... 28.50 (start + 100 steps)
GRID_B = [round(28.50 + 0.05 * i, 2) for i in range(31)]          # 28.50 ... 30.00
STAT_TOL, ARREST_RES_TOL = 1e-6, 1e-5
BENCH = {"b_s": dict(lo=28.80, hi=29.00, t=600.0, old=(28.92859, 28.92871)),
         "b_r": dict(lo=29.30, hi=29.50, t=800.0, old=(29.3965, 29.3969))}


class Stop(Exception):
    pass


# ------------------------------------------------------------------ ledger
def _ledger_path():
    return os.path.join(OUT, "ledger.json")


def ledger():
    try:
        with open(_ledger_path()) as f:
            return json.load(f)
    except FileNotFoundError:
        return dict(integrations=[], cap=RUN_CAP, cutoff_core_s=CUTOFF_CORE_S)


def _ledger_add(key, wall):
    lg = ledger()
    lg["integrations"].append(dict(key=key, wall_s=wall, at=time.strftime("%Y-%m-%dT%H:%M:%S")))
    os.makedirs(OUT, exist_ok=True)
    with open(_ledger_path(), "w") as f:
        json.dump(lg, f, indent=1)


def _ledger_check():
    lg = ledger()
    n, used = len(lg["integrations"]), sum(r["wall_s"] for r in lg["integrations"])
    if n + 1 > RUN_CAP:
        raise Stop(f"run cap {RUN_CAP} reached")
    if used >= CUTOFF_CORE_S:
        raise Stop(f"cutoff {CUTOFF_CORE_S / 3600:.1f} core-h reached")


# ------------------------------------------------------------------ solver wrapper (§3)
def lap_periodic(N, h):
    e = np.ones(N)
    Lm = sp.diags([e[:-1], -2 * e, e[:-1]], [-1, 0, 1], format="lil")
    Lm[0, N - 1] = 1
    Lm[N - 1, 0] = 1
    return Lm.tocsr() / h ** 2


def mg9_ode(b, N, h, t_end, y0, *, seed_scale, rtol=1e-8, atol=1e-10, t_eval=None, t_start=0.0):
    """Same right-hand side, Jacobian, BDF settings and first-step retries as lattice.lattice_ode_1d, but seed_scale
    must be given explicitly and success / time reached / retry are returned. Potentials about the b-uniform state."""
    Pb = {**P, "b": b}
    X0, Y0 = L.X0, L.Y0(b)
    x = (np.arange(N) - N // 2) * h
    Lap = lap_periodic(N, h)
    g, p, q = Pb["g"], Pb["p"], Pb["q"]

    def rhs(t, y):
        pG, pX, pY = y[:N], y[N:2 * N], y[2 * N:]
        X, Y = X0 + pX, Y0 + pY
        chi = -seed_scale * np.exp(-x ** 2 / 2) * np.exp(-(t - 10) ** 2 / 18)
        return np.concatenate([Lap @ pG - q * pG + g * pX,
                               Pb["dx"] * (Lap @ pX) + p * pG - (1 + b) * pX + X ** 2 * Y - X0 ** 2 * Y0 + chi,
                               Pb["dy"] * (Lap @ pY) + b * pX - (X ** 2 * Y - X0 ** 2 * Y0)])

    def jac(t, y):
        X, Y = X0 + y[N:2 * N], Y0 + y[2 * N:]
        I = sp.identity(N, format="csr")
        return sp.bmat([[Lap - q * I, g * I, None],
                        [p * I, Pb["dx"] * Lap + sp.diags(-(1 + b) + 2 * X * Y), sp.diags(X ** 2)],
                        [None, sp.diags(b - 2 * X * Y), Pb["dy"] * Lap - sp.diags(X ** 2)]], format="csc")

    sol, retry = None, None
    for retry in (None, 1e-3, 1e-2):
        sol = solve_ivp(rhs, (t_start, t_end), np.asarray(y0, float), method="BDF", jac=jac, rtol=rtol, atol=atol,
                        max_step=0.5, first_step=retry, t_eval=t_eval)
        if sol.success:
            break
    t_reached = float(sol.t[-1]) if sol.t.size else t_start
    complete = bool(sol.success) and abs(t_reached - t_end) <= 1e-9 * max(1.0, abs(t_end))
    yf = sol.y[:, -1]
    return dict(x=x, t=sol.t, y=sol.y, y_final=yf, residual=float(np.abs(rhs(t_end, yf)).max()),
                success=bool(sol.success), t_reached=t_reached, complete=complete, retry=retry)


def integrate(key, b, N, h, t_end, y0, *, seed_scale, **kw):
    """Ledger-checked integration (one count per call)."""
    _ledger_check()
    t0 = time.time()
    r = mg9_ode(b, N, h, t_end, y0, seed_scale=seed_scale, **kw)
    _ledger_add(key, time.time() - t0)
    return r


# ------------------------------------------------------------------ field helpers
def carry_concentration(y, N, b_old, b_new):
    """§3: preserve concentrations; only phi_Y changes, by Y0(b_old) - Y0(b_new)."""
    y = np.array(y, float)
    y[2 * N:] += L.Y0(b_old) - L.Y0(b_new)
    return y


def fourier_resample(y, N_old, N_new):
    """Zero-padding Fourier interpolation of each field (both grids start at x = -L/2)."""
    out = []
    for c in range(3):
        f = np.fft.rfft(y[c * N_old:(c + 1) * N_old])
        F = np.zeros(N_new // 2 + 1, complex)
        m = min(len(f), len(F))
        F[:m] = f[:m]
        if N_old % 2 == 0 and m == N_old // 2 + 1 and N_new > N_old:
            F[m - 1] *= 0.5                                      # split the Nyquist term
        out.append(np.fft.irfft(F, n=N_new) * (N_new / N_old))
    return np.concatenate(out)


def fourier_shift(y, N, h, s):
    """Translate each field by distance s (continuous; exact for the represented modes)."""
    k = 2 * np.pi * np.fft.rfftfreq(N, d=h)
    out = []
    for c in range(3):
        f = np.fft.rfft(y[c * N:(c + 1) * N]) * np.exp(-1j * k * s)
        out.append(np.fft.irfft(f, n=N))
    return np.concatenate(out)


def kernel_weights(h, lam=LAM):
    """§6.2: boxcar of physical width exactly lam; node j covers [j-1/2, j+1/2]; weights sum to lam/h."""
    half = lam / 2 / h
    n = int(np.floor(half)) + 1
    w = {}
    for j in range(-n, n + 1):
        lo, hi = max(j - 0.5, -half), min(j + 0.5, half)
        if hi > lo:
            w[j] = hi - lo
    return w


def boxcar(f, h, lam=LAM):
    w = kernel_weights(h, lam)
    out = np.zeros_like(f, dtype=float)
    for j, wt in w.items():
        out += wt * np.roll(f, -j)
    return out / sum(w.values())


def envelope_rms(phiX, h, lam=LAM):
    """Draft 3 §6.2 as registered (local RMS with one lam-boxcar). It FAILED its registered synthetic validation
    (phase-dependent edge error up to 0.27 lam); kept for the record and the test that documents the failure."""
    m = boxcar(phiX, h, lam)
    return np.sqrt(2 * boxcar((phiX - m) ** 2, h, lam))


def envelope(phiX, h, lam=LAM):
    """Amendment 1 (spec/MG9_amendment_1.md): complex demodulation at k = 2 pi / lam with the lam-boxcar applied
    twice as the low-pass. m = K2[phi_X]; z = (phi_X - m) e^{-ikx}; E = 2 |K2[z]|. x = (j - N/2) h."""
    N = len(phiX)
    x = (np.arange(N) - N // 2) * h
    K2 = lambda f: boxcar(boxcar(f, h, lam), h, lam)
    m = K2(phiX)
    z = (phiX - m) * np.exp(-2j * np.pi * x / lam)
    return 2 * np.abs(K2(z.real) + 1j * K2(z.imag))


def cell_count(phiX, h, E_pat, lam=LAM):
    """Minima of phi_X with prominence >= E_pat/2 and separation >= lam/2 (periodic)."""
    N = len(phiX)
    ext = np.concatenate([phiX, phiX, phiX])
    pk, _ = find_peaks(-ext, prominence=0.5 * E_pat, distance=max(1, int(np.ceil(lam / 2 / h))))
    return int(np.sum((pk >= N) & (pk < 2 * N)))


def components(mask, h, lam=LAM):
    """§6.2: periodic components of a boolean mask; gaps < lam merged, pieces < lam dropped.
    Returns ('fill', None) if one component covers the ring with no gap >= lam, else list of index arrays."""
    N = len(mask)
    m = np.asarray(mask, bool).copy()
    if not m.any():
        return []
    if m.all():
        return "fill"
    runs = _runs(m)                                               # list of (start, length, value) on the ring
    for s, ln, v in runs:
        if not v and ln * h < lam:
            m[(s + np.arange(ln)) % N] = True
    if m.all():
        return "fill"
    comps = [((s + np.arange(ln)) % N) for s, ln, v in _runs(m) if v and ln * h >= lam]
    return comps


def _runs(m):
    N = len(m)
    if m.all() or (~m).all():
        return [(0, N, bool(m[0]))]
    start = int(np.argmax(m != m[-1]))                            # first index where value changes from m[-1]
    out, s = [], start
    for i in range(1, N + 1):
        j = (start + i) % N
        if i == N or m[j] != m[(start + i - 1) % N]:
            out.append((s % N, (start + i) - s, bool(m[s % N])))
            s = start + i
    return out


def splice_patch(state, N, h, half_width):
    """§6.1: state (potentials) times w(x) = 1/2[tanh((x+a)/l) - tanh((x-a)/l)], l = lam/2; zero elsewhere."""
    x = (np.arange(N) - N // 2) * h
    w = 0.5 * (np.tanh((x + half_width) / (LAM / 2)) - np.tanh((x - half_width) / (LAM / 2)))
    return np.concatenate([w * state[c * N:(c + 1) * N] for c in range(3)])


def phase_convention(state, N, h):
    """Translate by whole nodes so the cell minimum of phi_X nearest x = 0 sits at node N//2 (x = 0)."""
    pX = state[N:2 * N]
    E_pat = float(envelope(pX, h).mean())
    pk, _ = find_peaks(-np.concatenate([pX, pX, pX]), prominence=0.5 * E_pat,
                       distance=max(1, int(np.ceil(LAM / 2 / h))))
    pk = pk[(pk >= N) & (pk < 2 * N)] - N
    j = int(pk[np.argmin(np.abs(pk - N // 2))])
    return np.concatenate([np.roll(state[c * N:(c + 1) * N], N // 2 - j) for c in range(3)])


# ------------------------------------------------------------------ front tracking and classification (§6.2-6.3)
class Tracker:
    def __init__(self, N, h, E_pat, lam=LAM):
        self.N, self.h, self.E_pat, self.lam = N, h, E_pat, lam
        self.prev = None
        self.xL = self.xR = None

    def sample(self, phiX):
        E = envelope(phiX, self.h, self.lam)
        comps = components(E >= 0.5 * self.E_pat, self.h, self.lam)
        n = cell_count(phiX, self.h, self.E_pat, self.lam)
        rec = dict(n=n, status="ok", breakup=False, remote=False, xL=None, xR=None)
        if comps == "fill":
            rec["status"] = "fill"
            self.prev = np.arange(self.N)
            return rec
        if not comps:
            rec["status"] = "loss"
            self.prev = None
            return rec
        if self.prev is None:                                     # t = 0, or after loss
            ctr = self.N // 2
            hit = [c for c in comps if ctr in set(c.tolist())]
            if not hit:
                rec["status"] = "untracked"
                return rec
            tracked = hit[0]
        else:
            prev = set(self.prev.tolist())
            ov = [len(prev.intersection(c.tolist())) for c in comps]
            if max(ov) == 0:
                rec["status"] = "untracked"
                rec["remote"] = True
                return rec
            i = int(np.argmax(ov))
            tracked = comps[i]
            big = [o for o in ov if o * self.h >= self.lam]
            rec["breakup"] = len(big) >= 2
            rec["remote"] = any(o == 0 for o in ov)
        self.prev = tracked
        xl, xr = self._ends(E, tracked)
        rec["xL"], rec["xR"] = self._unwrap(xl, xr)
        return rec

    def _ends(self, E, comp):
        N, h, th = self.N, self.h, 0.5 * self.E_pat
        s, e = int(comp[0]), int(comp[-1])                        # comp is ordered along the ring
        x = lambda j: (j - N // 2) * h
        def cross(i_in, i_out, side):
            a, c = E[i_in % N] - th, E[i_out % N] - th
            f = a / (a - c) if a != c else 0.5
            return x(i_in) + side * f * h
        return cross(s, s - 1, -1), cross(e, e + 1, +1)

    def _unwrap(self, xl, xr):
        Lx = self.N * self.h
        if self.xL is None:
            if xl > xr:                                           # the patch wraps through x = +/- L/2
                xl -= Lx
            self.xL, self.xR = xl, xr
        else:
            xl += Lx * np.round((self.xL - xl) / Lx)
            xr += Lx * np.round((self.xR - xr) / Lx)
            self.xL, self.xR = xl, xr
        return self.xL, self.xR


def has_reversal(x, d):
    """Times t1 < t2 < t3 with opposite moves of at least d each."""
    x = np.asarray(x, float)
    if x.size < 3:
        return False
    mn_b, mx_b = np.minimum.accumulate(x), np.maximum.accumulate(x)
    mn_a, mx_a = np.minimum.accumulate(x[::-1])[::-1], np.maximum.accumulate(x[::-1])[::-1]
    up_down = (x - mn_b >= d) & (x - mn_a >= d)
    down_up = (mx_b - x >= d) & (mx_a - x >= d)
    return bool(np.any(up_down | down_up))


def classify(track, T, residual, complete, lam=LAM, t_win=100.0):
    """§6.3 ordered procedure. track: list of dict(t, status, xL, xR, n, breakup, remote) every 5 time units."""
    if not complete or not track or track[0]["status"] != "ok":
        return "UNRESOLVED-NUMERICAL", {}
    st = [r["status"] for r in track]
    first = next((i for i, s in enumerate(st) if s in ("fill", "loss")), None)
    if first is not None:
        ev = st[first]
        holds = all(s == ev for s in st[first:])
        info = dict(event=ev, t_event=track[first]["t"])
        if holds:
            return ("INVADES" if ev == "fill" else "RETREATS"), info
        return "UNRESOLVED-TEMPORARY", info
    if any(s == "untracked" for s in st):
        return "UNRESOLVED-INTERNAL", dict(reason="untracked")
    w = [r for r in track if r["t"] >= t_win - 1e-9]
    if any(r["breakup"] or r["remote"] for r in w):
        return "UNRESOLVED-INTERNAL", dict(reason="breakup/remote")
    t = np.array([r["t"] for r in w])
    oL = -(np.array([r["xL"] for r in w]) - w[0]["xL"])          # outward displacement, left
    oR = np.array([r["xR"] for r in w]) - w[0]["xR"]
    if has_reversal(oL, lam / 2) or has_reversal(oR, lam / 2):
        return "UNRESOLVED-REVERSAL", {}
    dL, dR = float(oL[-1]), float(oR[-1])
    info = dict(dL=dL, dR=dR)
    if dL >= lam and dR >= lam:
        return "INVADES", info
    if dL <= -lam and dR <= -lam:
        return "RETREATS", info
    if max(abs(dL), abs(dR)) >= lam:
        return "UNRESOLVED-ASYMMETRIC", info
    late = t >= t[-1] - 200.0 - 1e-9
    vL = float(np.polyfit(t[late], oL[late], 1)[0])
    vR = float(np.polyfit(t[late], oR[late], 1)[0])
    n = np.array([r["n"] for r in w])
    info.update(vL=vL, vR=vR, n_range=int(n.max() - n.min()), residual=residual)
    if (abs(dL) <= lam / 2 and abs(dR) <= lam / 2 and abs(vL) <= lam / 2000 and abs(vR) <= lam / 2000
            and residual <= ARREST_RES_TOL and n.max() - n.min() <= 1):
        margin = abs(dL) <= lam / 4 and abs(dR) <= lam / 4 and abs(vL) <= lam / 4000 and abs(vR) <= lam / 4000
        info["margin"] = bool(margin)
        return "ARRESTED", info
    return "SLOW", info


def needs_extension(label, info):
    return label == "SLOW" or (label == "ARRESTED" and not info.get("margin", False))


# ------------------------------------------------------------------ records
def _key(*parts):
    return "_".join(str(p) for p in parts)


def rec_path(key):
    return os.path.join(RUNS, key + ".json")


def load_rec(key):
    try:
        with open(rec_path(key)) as f:
            return json.load(f)
    except FileNotFoundError:
        return None


def save_rec(key, rec):
    os.makedirs(RUNS, exist_ok=True)
    with open(rec_path(key), "w") as f:
        json.dump(rec, f, indent=1, default=float)


def state_path(key):
    return os.path.join(STATES, key + ".npy")


def save_state(key, y):
    os.makedirs(STATES, exist_ok=True)
    np.save(state_path(key), np.asarray(y, float))


def load_state(key):
    return np.load(state_path(key))


def sha(y):
    return hashlib.sha256(np.ascontiguousarray(y, dtype=np.float64).tobytes()).hexdigest()


# ------------------------------------------------------------------ P1 stage
def mode_field(b, N, h, A0, with_slaved=True):
    """A0 cos(kx) phi + A0^2 [w0 + w2 cos(2kx)] on the D02-type mesh (problem (i) vectors at the given b)."""
    k = L.K34
    out = L.cubic_coefficient(L.M(b, L.kappa(k, h)), L.M(b, 0.0), L.M(b, L.kappa(2 * k, h)), L.Q_mg(b), L.C_mg)
    x = (np.arange(N) - N // 2) * h
    c1, c2 = np.cos(k * x), np.cos(2 * k * x)
    f = [A0 * c1 * out["phi"][i] + (A0 ** 2 * (out["w0"][i] + out["w2"][i] * c2) if with_slaved else 0) for i in range(3)]
    return np.concatenate(f), out


def project34(y, N, h, psi):
    x = (np.arange(N) - N // 2) * h
    c = np.cos(L.K34 * x)
    u = np.array([2.0 / N * np.dot(y[i * N:(i + 1) * N], c) for i in range(3)])
    return float(psi @ u)


def v2_run(A0, b, sig_unused=None, rtol=1e-8, atol=1e-10, tag=""):
    d = L.DOMAINS["D02"]
    N, h = d["N"], d["h"]
    y0, cc = mode_field(b, N, h, A0)
    key = _key("V2", f"A{A0:g}", tag) if tag else _key("V2", f"A{A0:g}")
    if load_rec(key):
        return load_rec(key)
    _ledger_check()
    t0w = time.time()
    ts, As, y, t, stop_t, ok = [0.0], [project34(y0, N, h, cc["psi"])], y0, 0.0, None, True
    while t < 3000.0 - 1e-9 and stop_t is None:
        t1 = min(t + 100.0, 3000.0)
        te = np.arange(t + 1.0, t1 + 1e-9, 1.0)
        r = mg9_ode(b, N, h, t1, y, seed_scale=0.0, rtol=rtol, atol=atol, t_eval=te, t_start=t)
        if not r["complete"]:
            ok = False
            break
        for j, tj in enumerate(r["t"]):
            a = project34(r["y"][:, j], N, h, cc["psi"])
            ts.append(float(tj)); As.append(a)
            if abs(a) >= 2 * A0 or abs(a) <= A0 / 2:
                stop_t = float(tj)
                break
        y, t = r["y_final"], t1
    _ledger_add(key, time.time() - t0w)
    rec = dict(key=key, A0=A0, b=b, rtol=rtol, atol=atol, t=ts, A=As, stop_t=stop_t if stop_t is not None else t,
               complete=ok)
    save_rec(key, rec)
    return rec


def stage_p1():
    out = {}
    bc, kc = L.bc_continuous()
    out["A1"] = dict(continuous=dict(bc=bc, kappa_c=kc),
                     boxes={d: L.box_thresholds(d) for d in L.DOMAINS},
                     alpha=dict(D02_mode34=L.alpha(L.K34, 0.2), D01_mode34=L.alpha(L.K34, 0.1),
                                h01_mode33=L.alpha(2 * np.pi * 33 / 100, 0.1)))
    bs = np.round(np.arange(28.5, 30.5 + 1e-9, 0.05), 2)
    out["A1"]["sigma_scan"] = dict(b=bs.tolist(),
                                   continuous=[L.sigma_max_continuous(b)[0] for b in bs],
                                   D02_mode34=[L.sigma(b, L.kappa(L.K34, 0.2)) for b in bs],
                                   D01_mode34=[L.sigma(b, L.kappa(L.K34, 0.1)) for b in bs])
    sh = {c: L.l3_swift_hohenberg(c) for c in (0.0, 0.5, 1.0)}
    code_ok = all(abs(sh[c] - (19 * c * c / 18 - 0.75)) <= 1e-8 for c in sh)
    l3 = {p: L.l3_modelg(p) for p in ("i", "ii", "iii", "iv")}
    out["code_check"] = dict(values=sh, passed=code_ok)
    out["L3"] = {p: dict(L3=v["L3"], terms=v["terms"], b=v["b"], k=v["k"], kappa=v["kappa"], kappa2=v["kappa2"],
                         phi=v["phi"].tolist(), psi=v["psi"].tolist(), null_residual=v["null_residual"])
                 for p, v in l3.items()}
    b_i = l3["i"]["b"]
    sig = L.sigma(b_i, L.kappa(L.K34, 0.2))
    out["sigma34_at_b_i"] = sig
    S = [0.05, 0.1, 0.2]
    runs = [v2_run(a, b_i) for a in S]
    dec = L.v2_decision([dict(r) for r in runs], l3["i"]["L3"], sig)
    doubled = False
    if dec["n_identifiable"] == 0:
        S, doubled = [0.1, 0.2, 0.4], True
        runs = [v2_run(a, b_i) for a in S]
    mid = S[1]
    rep = v2_run(mid, b_i, rtol=1e-10, atol=1e-12, tag="tight")
    rr = [dict(r) for r in runs]
    base_fit = L.fit_l3(next(r for r in runs if r["A0"] == mid)["t"], next(r for r in runs if r["A0"] == mid)["A"], sig)
    rep_fit = L.fit_l3(rep["t"], rep["A"], sig) if rep["stop_t"] >= 40 else dict(L3=np.nan)
    dec = L.v2_decision(rr, l3["i"]["L3"], sig, repeat=(base_fit["L3"], rep_fit["L3"]))
    if not all(r["complete"] for r in runs + [rep]):
        dec["passed"] = False
    out["V2"] = dict(amplitudes=S, doubled=doubled, decision=dec,
                     fits={str(r["A0"]): dict(status=r.get("status"), fit=r.get("fit"), stop_t=r["stop_t"])
                           for r in rr},
                     repeat=dict(A0=mid, base=base_fit["L3"], tight=rep_fit["L3"]))
    out["P1"] = L.p1_outcome(code_ok, dec["passed"], {p: v["L3"] for p, v in l3.items()}, l3["i"]["terms"])
    os.makedirs(OUT, exist_ok=True)
    with open(os.path.join(OUT, "p1.json"), "w") as f:
        json.dump(out, f, indent=1, default=float)
    print("P1:", out["P1"], flush=True)
    return out


# ------------------------------------------------------------------ preparation (sweep and settles)
def prep_classify(r, N, h, E_ref):
    if not r["complete"] or r["residual"] > STAT_TOL:
        return "PREP-UNRESOLVED"
    pX = r["y_final"][N:2 * N]
    if np.abs(pX).max() < 0.1:
        return "PREP-LOST"
    E_pat = float(envelope(pX, h).mean())
    n = cell_count(pX, h, E_pat)
    return "PERSISTS" if n == int(round(N * h / LAM)) else "PREP-CHANGED"


def settle(key, b, N, h, y0, t_settle=1000.0, t_ext=3000.0):
    """Integrate to t_settle, extend once (same integration count) to t_ext if not stationary."""
    _ledger_check()
    t0w = time.time()
    r = mg9_ode(b, N, h, t_settle, y0, seed_scale=0.0)
    ext = False
    if r["complete"] and r["residual"] > STAT_TOL:
        r = mg9_ode(b, N, h, t_ext, r["y_final"], seed_scale=0.0, t_start=t_settle)
        ext = True
    _ledger_add(key, time.time() - t0w)
    status = prep_classify(r, N, h, None)
    return r, status, ext


def stage_sweep():
    d = L.DOMAINS["D02"]
    N, h = d["N"], d["h"]
    prev_b, prev_y = None, None
    for i, b in enumerate(SWEEP_B):
        key = _key("sweep", f"{b:.2f}")
        rec = load_rec(key)
        if rec is None:
            if i == 0:
                y0 = _mode_field_unstable(b, N, h, 0.5)
            else:
                y0 = carry_concentration(prev_y, N, prev_b, b)
            r, status, ext = settle(key, b, N, h, y0)
            rec = dict(key=key, b=b, status=status, extended=ext, complete=r["complete"], success=r["success"],
                       t_reached=r["t_reached"], residual=r["residual"], retry=r["retry"], from_b=prev_b,
                       carry="concentration", seed_scale=0.0)
            if r["complete"]:
                pX = r["y_final"][N:2 * N]
                E_pat = float(envelope(pX, h).mean())
                rec.update(E_pat=E_pat, cells=cell_count(pX, h, E_pat) if np.abs(pX).max() >= 0.1 else 0,
                           sha256=sha(r["y_final"]))
                save_state(key, r["y_final"])
            save_rec(key, rec)
            print(key, status, flush=True)
            if i == 0 and status != "PERSISTS":
                print("sweep start failed: PREP-UNRESOLVED for P2", flush=True)
                return
        # continuation (§5): from the changed / uniform / last complete stationary state
        if rec["status"] in ("PERSISTS", "PREP-CHANGED", "PREP-LOST"):
            prev_y = load_state(key) if rec["status"] != "PREP-LOST" else np.zeros(3 * N)
            prev_b = b
        else:
            pass                                                  # PREP-UNRESOLVED: keep the last stationary state


def _mode_field_unstable(b, N, h, A0):
    """Sweep start: A0 cos(kx) times the X-normalized leading eigenvector of M(kappa_34) at b (mode 34 unstable)."""
    w, V = np.linalg.eig(L.M(b, L.kappa(L.K34, h)))
    v = np.real(V[:, int(np.argmax(w.real))])
    v = v / v[1]
    x = (np.arange(N) - N // 2) * h
    c = np.cos(L.K34 * x)
    return np.concatenate([A0 * v[i] * c for i in range(3)])


def is_sweep_b(b):
    return any(abs(b - s) < 1e-9 for s in SWEEP_B)


def prepared_state(b):
    """§6.1: returns (status, state or None, provenance) for D02 at b, settling (and recording) if needed."""
    d = L.DOMAINS["D02"]
    N, h = d["N"], d["h"]
    b = round(b, 4)
    if is_sweep_b(b):
        key = _key("sweep", f"{b:.2f}")
        rec = load_rec(key)
        if rec is None:
            return "NOT-RUN", None, key
        st = "PERSISTS" if rec["status"] == "PERSISTS" else rec["status"]
        return st, (load_state(key) if st == "PERSISTS" else None), key
    above = min(s for s in SWEEP_B if s > b + 1e-9)
    src = _key("sweep", f"{above:.2f}")
    srec = load_rec(src)
    key = _key("prep", f"{b:.4f}")
    rec = load_rec(key)
    if rec is None:
        if srec is None:
            return "NOT-RUN", None, key                           # nothing recorded until the sweep has run
        if srec["status"] != "PERSISTS":
            rec = dict(key=key, b=b, status=srec["status"], source=src,
                       note="source sweep state not PERSISTS; its status is propagated, nothing settled")
            save_rec(key, rec)
        else:
            y0 = carry_concentration(load_state(src), N, above, b)
            r, status, ext = settle(key, b, N, h, y0)
            rec = dict(key=key, b=b, status=status, extended=ext, complete=r["complete"], success=r["success"],
                       t_reached=r["t_reached"], residual=r["residual"], retry=r["retry"], source=src,
                       carry="concentration", seed_scale=0.0)
            if status == "PERSISTS":
                save_state(key, r["y_final"])
                rec["sha256"] = sha(r["y_final"])
            save_rec(key, rec)
    st = rec["status"]
    return st, (load_state(key) if st == "PERSISTS" else None), key


# ------------------------------------------------------------------ front runs (§6.1-6.3)
def front_run(key, b, N, h, y_patch, E_pat):
    """Integrate a spliced patch to 600, extend to 1800 if needed (one integration count); classify."""
    rec = load_rec(key)
    if rec is not None:
        return rec
    _ledger_check()
    t0w = time.time()
    tr = Tracker(N, h, E_pat)
    track = [dict(t=0.0, **tr.sample(y_patch[N:2 * N]))]
    def run(t_a, t_b, y):
        te = np.arange(t_a + 5.0, t_b + 1e-9, 5.0)
        r = mg9_ode(b, N, h, t_b, y, seed_scale=0.0, t_eval=te, t_start=t_a)
        for j, tj in enumerate(r["t"]):
            track.append(dict(t=float(tj), **tr.sample(r["y"][N:2 * N, j])))
        return r
    r = run(0.0, 600.0, y_patch)
    label, info = classify(track, 600.0, r["residual"], r["complete"])
    T, extended = 600.0, False
    if r["complete"] and needs_extension(label, info):
        r = run(600.0, 1800.0, r["y_final"])
        T, extended = 1800.0, True
        label, info = classify(track, 1800.0, r["residual"], r["complete"])
        if label == "SLOW":
            label = "UNRESOLVED-SLOW"
    elif label == "SLOW":
        label = "UNRESOLVED-SLOW"
    _ledger_add(key, time.time() - t0w)
    rec = dict(key=key, b=b, N=N, h=h, E_pat=E_pat, T=T, extended=extended, label=label, info=info,
               complete=r["complete"], success=r["success"], t_reached=r["t_reached"], residual=r["residual"],
               retry=r["retry"], seed_scale=0.0, track=track)
    save_rec(key, rec)
    print(key, label, flush=True)
    return rec


def d02_front(b, tag="grid"):
    d = L.DOMAINS["D02"]
    N, h = d["N"], d["h"]
    st, y, prov = prepared_state(b)
    key = _key("front", tag, f"{b:.4f}")
    if st == "NOT-RUN":
        raise RuntimeError(f"sweep state for b = {b} not run yet; run the sweep stage first")
    if st != "PERSISTS":
        rec = load_rec(key) or dict(key=key, b=b, label=("UNRESOLVED-PREP" if st in ("PREP-UNRESOLVED", "NOT-RUN")
                                                         else "UNAVAILABLE-" + st.replace("PREP-", "")),
                                    prep=prov)
        save_rec(key, rec)
        return rec
    E_pat = float(envelope(y[N:2 * N], h).mean())
    yp = splice_patch(phase_convention(y, N, h), N, h, 25.0)
    rec = front_run(key, b, N, h, yp, E_pat)
    rec["prep"] = prov
    save_rec(key, rec)
    return rec


def front_class(label):
    if label == "ARRESTED":
        return "A"
    if label in ("INVADES", "RETREATS"):
        return "M"
    if label.startswith("UNAVAILABLE"):
        return "N"
    return "U"


def stage_grid():
    for b in GRID_B:
        d02_front(b)


# ------------------------------------------------------------------ edges, refinement, robustness (§6.4)
def grid_classes(tag_points):
    """tag_points: dict b -> label. Returns sorted list of (b, class)."""
    return sorted((b, front_class(l)) for b, l in tag_points.items())


def arrest_edges(seq):
    """seq: sorted list of (b, class). Returns dict with interval, edges (status, inner b, outer b) or a reason."""
    A = [i for i, (_, c) in enumerate(seq) if c == "A"]
    if not A:
        return dict(intervals=0)
    runs, cur = [], [A[0]]
    for i in A[1:]:
        if i == cur[-1] + 1:
            cur.append(i)
        else:
            runs.append(cur); cur = [i]
    runs.append(cur)
    if len(runs) > 1:
        return dict(intervals=len(runs))
    lo, hi = runs[0][0], runs[0][-1]
    def edge(i_in, i_out):
        if i_out < 0 or i_out >= len(seq):
            return dict(status="unresolved", inner=seq[i_in][0], outer=None, reason="scan end")
        c = seq[i_out][1]
        st = {"M": "resolved", "N": "censored", "U": "unresolved"}[c]
        return dict(status=st, inner=seq[i_in][0], outer=seq[i_out][0])
    return dict(intervals=1, lower=edge(lo, lo - 1), upper=edge(hi, hi + 1),
                interval=(seq[lo][0], seq[hi][0]))


def all_d02_fronts():
    out = {}
    for fn in sorted(os.listdir(RUNS)) if os.path.isdir(RUNS) else []:
        if fn.startswith("front_grid_") or fn.startswith("front_refine_"):
            r = load_rec(fn[:-5])
            out[round(r["b"], 4)] = r["label"]
    return out


def stage_refine():
    for side in ("lower", "upper"):
        for _ in range(2):
            ed = arrest_edges(grid_classes(all_d02_fronts()))
            if ed.get("intervals") != 1 or ed[side]["status"] != "resolved":
                break
            mid = round(0.5 * (ed[side]["inner"] + ed[side]["outer"]), 4)
            if mid in all_d02_fronts():
                break
            d02_front(mid, tag="refine")
    ed = arrest_edges(grid_classes(all_d02_fronts()))
    if ed.get("intervals") == 1:
        for side in ("lower", "upper"):
            robustness(ed[side]["inner"])
    print(json.dumps(final_edges(), indent=1, default=float), flush=True)


def robustness(b):
    """§6.4 (a) D02L tiled, settled 300, 100-unit patch; (b) D02 with lam/2 Fourier shift, patch on a maximum."""
    d, dl = L.DOMAINS["D02"], L.DOMAINS["D02L"]
    N, h, NL = d["N"], d["h"], dl["N"]
    st, y, prov = prepared_state(b)
    res = {}
    keyL = _key("robust", "D02L", f"{b:.4f}")
    if load_rec(keyL) is None:
        tiled = np.concatenate([np.tile(y[c * N:(c + 1) * N], 2) for c in range(3)])
        skey = _key("prepL", f"{b:.4f}")
        srec = load_rec(skey)
        if srec is None:
            _ledger_check(); t0w = time.time()
            r = mg9_ode(b, NL, h, 300.0, tiled, seed_scale=0.0)
            _ledger_add(skey, time.time() - t0w)
            ok = r["complete"] and r["residual"] <= STAT_TOL
            n = cell_count(r["y_final"][NL:2 * NL], h, float(envelope(r["y_final"][NL:2 * NL], h).mean())) if ok else None
            srec = dict(key=skey, b=b, complete=r["complete"], residual=r["residual"], cells=n,
                        status="PERSISTS" if ok and n == 68 else "PREP-UNRESOLVED" if not ok else "PREP-CHANGED")
            save_rec(skey, srec)
            if srec["status"] == "PERSISTS":
                save_state(skey, r["y_final"])
        if srec["status"] == "PERSISTS":
            yL = load_state(skey)
            E_pat = float(envelope(yL[NL:2 * NL], h).mean())
            front_run(keyL, b, NL, h, splice_patch(phase_convention(yL, NL, h), NL, h, 50.0), E_pat)
        else:
            save_rec(keyL, dict(key=keyL, b=b, label="UNRESOLVED-PREP", prep=skey))
    res["D02L"] = load_rec(keyL)["label"]
    keyS = _key("robust", "shift", f"{b:.4f}")
    if load_rec(keyS) is None:
        E_pat = float(envelope(y[N:2 * N], h).mean())
        ys = fourier_shift(phase_convention(y, N, h), N, h, LAM / 2)
        front_run(keyS, b, N, h, splice_patch(ys, N, h, 25.0), E_pat)
    res["shift"] = load_rec(keyS)["label"]
    return res


def final_edges():
    """Edges after refinement with robustness applied (a failed repeat turns the point into U)."""
    pts = all_d02_fronts()
    for b in list(pts):
        if pts[b] == "ARRESTED":
            rl, rs = load_rec(_key("robust", "D02L", f"{b:.4f}")), load_rec(_key("robust", "shift", f"{b:.4f}"))
            if rl is not None and rs is not None and not (rl["label"] == "ARRESTED" and rs["label"] == "ARRESTED"):
                pts[b] = "UNRESOLVED-ROBUSTNESS"
    ed = arrest_edges(grid_classes(pts))
    if ed.get("intervals") == 1:
        for side in ("lower", "upper"):
            b = ed[side]["inner"]
            if load_rec(_key("robust", "D02L", f"{b:.4f}")) is None:
                ed[side]["robustness"] = "NOT RUN"
    ed["points"] = pts
    return ed


# ------------------------------------------------------------------ D01 mesh comparison (§6.4)
def d01_front(b):
    d, d1 = L.DOMAINS["D02"], L.DOMAINS["D01"]
    N, h, N1, h1 = d["N"], d["h"], d1["N"], d1["h"]
    key = _key("front", "D01", f"{b:.4f}")
    if load_rec(key):
        return load_rec(key)
    st, y, prov = prepared_state(b)
    if st == "NOT-RUN":
        raise RuntimeError(f"sweep state for b = {b} not run yet; run the sweep stage first")
    if st != "PERSISTS":
        rec = dict(key=key, b=b, label=("UNRESOLVED-PREP" if st in ("PREP-UNRESOLVED", "NOT-RUN")
                                        else "UNAVAILABLE-" + st.replace("PREP-", "")), prep=prov)
        save_rec(key, rec)
        return rec
    skey = _key("prep01", f"{b:.4f}")
    srec = load_rec(skey)
    if srec is None:
        y1 = fourier_resample(y, N, N1)
        r, status, ext = settle(skey, b, N1, h1, y1)
        srec = dict(key=skey, b=b, status=status, extended=ext, complete=r["complete"], residual=r["residual"],
                    source=prov)
        if status == "PERSISTS":
            save_state(skey, r["y_final"])
        save_rec(skey, srec)
    if srec["status"] != "PERSISTS":
        rec = dict(key=key, b=b, label=("UNRESOLVED-PREP" if srec["status"] == "PREP-UNRESOLVED"
                                        else "UNAVAILABLE-" + srec["status"].replace("PREP-", "")), prep=skey)
        save_rec(key, rec)
        return rec
    y1 = load_state(skey)
    E_pat = float(envelope(y1[N1:2 * N1], h1).mean())
    rec = front_run(key, b, N1, h1, splice_patch(phase_convention(y1, N1, h1), N1, h1, 25.0), E_pat)
    rec["prep"] = skey
    save_rec(key, rec)
    return rec


def stage_d01():
    ed = final_edges()
    if ed.get("intervals") != 1:
        print("D01: no single arrest interval at h = 0.2; MESH-SENSITIVE NOT EVALUATED", flush=True)
        return
    extra = 0
    for side, sgn in (("lower", -1), ("upper", +1)):
        e = ed[side]
        if e["status"] != "resolved":
            continue
        b_in, b_out = e["inner"], e["outer"]
        b1, b2 = sorted((b_in, b_out))
        pts = [round(b1 - 0.05, 4), b1, b2, round(b2 + 0.05, 4)]
        labs = {b: d01_front(b)["label"] for b in pts}
        cls = [front_class(labs[b]) for b in sorted(labs)]
        if not ("A" in cls and "M" in cls) and extra < 2:
            nb = round((min(pts) - 0.05) if sgn < 0 else (max(pts) + 0.05), 4)
            d01_front(nb)
            extra += 1


# ------------------------------------------------------------------ benchmark recheck (§6.5)
def stage_bench():
    d = L.DOMAINS["D02"]
    N, h = d["N"], d["h"]
    rkey = "bench_ref_b29"
    if load_rec(rkey) is None:
        r = integrate(rkey, 29.0, N, h, 400.0, np.zeros(3 * N), seed_scale=1.0)
        save_rec(rkey, dict(key=rkey, complete=r["complete"], residual=r["residual"], sha256=sha(r["y_final"])))
        save_state(rkey, r["y_final"])
    yref = load_state(rkey)

    def exists_bs(b):
        key = _key("bench", "bs", f"{b:.6f}")
        rec = load_rec(key)
        if rec is None:
            r = integrate(key, b, N, h, BENCH["b_s"]["t"], yref, seed_scale=0.0)   # potentials carried unchanged
            rec = dict(key=key, b=b, complete=r["complete"], phiX0=float(r["y_final"][N + N // 2]),
                       residual=r["residual"], carry="potential (MG003 convention)", seed_scale=0.0)
            save_rec(key, rec)
        return None if not rec["complete"] else rec["phiX0"] < -1.0

    def single_br(b):
        key = _key("bench", "br", f"{b:.6f}")
        rec = load_rec(key)
        if rec is None:
            r = integrate(key, b, N, h, BENCH["b_r"]["t"], np.zeros(3 * N), seed_scale=1.0)
            pX = r["y_final"][N:2 * N]
            nmin = int(np.sum((pX[1:-1] < pX[:-2]) & (pX[1:-1] < pX[2:]) & (pX[1:-1] < -4)))
            rec = dict(key=key, b=b, complete=r["complete"], nmin=nmin, residual=r["residual"], seed_scale=1.0)
            save_rec(key, rec)
        return None if not rec["complete"] else rec["nmin"] == 1

    out = {}
    for name, test, want_lo, want_hi in (("b_s", exists_bs, False, True), ("b_r", single_br, True, False)):
        lo, hi = BENCH[name]["lo"], BENCH[name]["hi"]
        vlo, vhi = test(lo), test(hi)
        if vlo is not want_lo or vhi is not want_hi:
            out[name] = dict(status="BENCHMARK-UNRESOLVED", endpoints=(vlo, vhi))
            continue
        ok = True
        for _ in range(8):
            m = round(0.5 * (lo + hi), 6)
            v = test(m)
            if v is None:
                ok = False
                break
            if v == want_lo:
                lo = m
            else:
                hi = m
        mid_old = 0.5 * sum(BENCH[name]["old"])
        out[name] = dict(status="OK" if ok else "BENCHMARK-UNRESOLVED", bracket=(lo, hi),
                         shift=0.5 * (lo + hi) - mid_old, sensitive=abs(0.5 * (lo + hi) - mid_old) > 0.05)
    with open(os.path.join(OUT, "benchmarks.json"), "w") as f:
        json.dump(out, f, indent=1, default=float)
    print(json.dumps(out, indent=1, default=float), flush=True)
    return out


STAGES = dict(p1=stage_p1, sweep=stage_sweep, bench=stage_bench, grid=stage_grid, refine=stage_refine, d01=stage_d01)

if __name__ == "__main__":
    stage = sys.argv[1] if len(sys.argv) > 1 else ""
    if stage not in STAGES:
        print(__doc__)
        sys.exit(1)
    try:
        STAGES[stage]()
    except Stop as e:
        print("STOP:", e, flush=True)
        sys.exit(2)

"""MG9d driver (spec/MG9d_proposal.md draft 3). Deterministic only. Resumable at run level.

    python src/run_mg9d.py floor     static measurement-floor validation (no integration) -> results/MG9d/floor.json
    python src/run_mg9d.py checks    state-transfer checks (no integration)               -> results/MG9d/transfer.json
    python src/run_mg9d.py runs      R0, R1, R3, R2a, R2b, R2c, R4 (ledger-checked; cap 10, cutoff 1 core-h)
    python src/run_mg9d.py pack      pack snapshots into results/MG9d/snapshots.npz with a manifest

`runs` refuses to start unless results/MG9d/REGISTERED exists (written only after PI registration).
Frozen MG9 machinery (envelope, Tracker, fourier_shift, fourier_resample, cell_count) is imported from run_mg9.
"""
import os, sys, json, time, hashlib
import numpy as np
import scipy.sparse as sp
from scipy.integrate import solve_ivp
from scipy.optimize import minimize_scalar
from scipy.signal import find_peaks

sys.path.insert(0, os.path.dirname(__file__))
import landau_mg9 as L
import run_mg9 as R

ROOT = os.path.join(os.path.dirname(__file__), "..")
OUT = os.path.join(ROOT, "results", "MG9d")
RUNS = os.path.join(OUT, "runs")
SNAP = os.path.join(OUT, "snapshots")
SRC_NPZ = os.path.join(ROOT, "results", "MG9c", "snapshots.npz")
SRC_MANIFEST = os.path.join(ROOT, "results", "MG9c", "snapshots_manifest.json")
SRC_KEY = "front_D02_29.30000_t30000"
SRC_REC = os.path.join(ROOT, "results", "MG9c", "runs", "front_D02_29.30000.json")
REGISTERED_FLAG = os.path.join(OUT, "REGISTERED")

B = 29.30
LAM = R.LAM
N02, H02 = 500, 0.2
N01, H01 = 1000, 0.1
T_END, SEG, DT_DIAG, DT_TRACK, DT_SNAP = 10000.0, 1000.0, 50.0, 5.0, 1000.0
RUN_CAP, CUTOFF_CORE_S = 10, 3600.0
SOLVER_REG = dict(rtol=1e-8, atol=1e-10, max_step=0.5)
SOLVER_TIGHT = dict(rtol=1e-10, atol=1e-12, max_step=0.1)
ORDER = ("R0", "R1", "R3", "R2a", "R2b", "R2c", "R4")
SHIFTS = dict(R2a=0.25 * H02, R2b=0.5 * H02, R2c=0.75 * H02)
FLOOR_SHIFTS = (0.001, 0.004, 0.01, 0.25 * H02, 0.5 * H02, 0.75 * H02)
WIN_A, WIN_B, WIN_V = (2000.0, 6000.0), (6000.0, 10000.0), (2000.0, 10000.0)

# §3 validity and §4 bands (provisional operational choices, frozen at registration)
AMBIG_TOL = 1e-3
S_MAX = 1e-3
CONS_REL, CONS_ABS = 0.10, 5e-4      # PROPOSED (pre-run spec item P2): draft 3 had 2e-4, below C_edge's static error
VALID_REF_T = 2000.0                 # PROPOSED (pre-run spec item P1): draft 3 measured items 5-6 from t' = 0
AGREE_REL, DIFFER_REL, DIFFER_VF = 0.05, 0.20, 3.0
SYM_REL = 0.01
RESOLVED_VF = 3.0
VF_MIN = 5e-9
Y0_MAX_CHANGE, ROUNDTRIP_TOL, R4_BACK_TOL = 2e-5, 1e-12, 1e-12


class Stop(Exception):
    pass


def e_pat():
    """§3: the exact recorded E_pat of the MG9c source run (rounded decimals elsewhere are display values)."""
    with open(SRC_REC) as f:
        return float(json.load(f)["E_pat"])


def sha(y):
    return hashlib.sha256(np.ascontiguousarray(y, dtype=np.float64).tobytes()).hexdigest()


def platform():
    import platform as pf, scipy
    return dict(python=pf.python_version(), numpy=np.__version__, scipy=scipy.__version__, system=pf.platform())


# ------------------------------------------------------------------ solver (§2): mg9_ode with max_step exposed
def make_rhs_jac(b, N, h, seed_scale=0.0):
    """The right-hand side and Jacobian of run_mg9.mg9_ode, unchanged, factored out so the initial RHS can be
    evaluated without stepping."""
    Pb = {**R.P, "b": b}
    X0, Y0 = L.X0, L.Y0(b)
    x = (np.arange(N) - N // 2) * h
    Lap = R.lap_periodic(N, h)
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
    return x, rhs, jac


def mg9d_ode(b, N, h, t_end, y0, *, seed_scale, rtol=1e-8, atol=1e-10, max_step=0.5, t_eval=None, t_start=0.0):
    """run_mg9.mg9_ode with max_step exposed. At max_step 0.5 it must reproduce mg9_ode bit for bit (tested)."""
    x, rhs, jac = make_rhs_jac(b, N, h, seed_scale)
    sol, retry = None, None
    for retry in (None, 1e-3, 1e-2):
        sol = solve_ivp(rhs, (t_start, t_end), np.asarray(y0, float), method="BDF", jac=jac, rtol=rtol, atol=atol,
                        max_step=max_step, first_step=retry, t_eval=t_eval)
        if sol.success:
            break
    t_reached = float(sol.t[-1]) if sol.t.size else t_start
    complete = bool(sol.success) and abs(t_reached - t_end) <= 1e-9 * max(1.0, abs(t_end))
    yf = sol.y[:, -1]
    return dict(x=x, t=sol.t, y=sol.y, y_final=yf, residual=float(np.abs(rhs(t_end, yf)).max()),
                success=bool(sol.success), t_reached=t_reached, complete=complete, retry=retry)


def rhs_max(y, N, h, b=B):
    _, rhs, _ = make_rhs_jac(b, N, h, 0.0)
    return float(np.abs(rhs(0.0, np.asarray(y, float))).max())


# ------------------------------------------------------------------ start states (§2)
def source_state():
    with np.load(SRC_NPZ) as z:
        y = np.array(z[SRC_KEY], float)
    with open(SRC_MANIFEST) as f:
        want = json.load(f)[SRC_KEY]
    want = want["sha256"] if isinstance(want, dict) else want
    if sha(y) != want:
        raise Stop(f"source hash mismatch: {sha(y)} != {want}")
    return y


def nyquist_filter(y, N):
    out = []
    for c in range(3):
        f = np.fft.rfft(y[c * N:(c + 1) * N])
        if N % 2 == 0:
            f[-1] = 0.0
        out.append(np.fft.irfft(f, n=N))
    return np.concatenate(out)


def reflect(y, N):
    """Reflection about the node x = 0 (index N//2): j -> (N - j) mod N, each field. Exact."""
    idx = (N - np.arange(N)) % N
    return np.concatenate([y[c * N:(c + 1) * N][idx] for c in range(3)])


def to_d01(y):
    return R.fourier_resample(y, N02, N01)


def from_d01_nodes(y1):
    """Sample a D01 state back at the D02 nodes (every second node; both grids start at x = -L/2)."""
    return np.concatenate([y1[c * N01:(c + 1) * N01][::2] for c in range(3)])


def common_start():
    return nyquist_filter(source_state(), N02)


def start_states(y0=None):
    """dict run -> dict(y, N, h, solver)."""
    y0 = common_start() if y0 is None else y0
    st = dict(R0=dict(y=y0, N=N02, h=H02, solver=SOLVER_REG),
              R1=dict(y=y0, N=N02, h=H02, solver=SOLVER_TIGHT),
              R3=dict(y=reflect(y0, N02), N=N02, h=H02, solver=SOLVER_REG),
              R4=dict(y=to_d01(y0), N=N01, h=H01, solver=SOLVER_REG))
    for k, s in SHIFTS.items():
        st[k] = dict(y=R.fourier_shift(y0, N02, H02, s), N=N02, h=H02, solver=SOLVER_REG, shift=s)
    return st


def relax_estimate(y, N, h, rel_cut=1e-6):
    """Static (no stepping) size of the shape adjustment a start state needs: the pseudo-inverse correction
    dy = -J^+ RHS(y), dropping singular values below rel_cut * max (the near-zero translation mode). Linearised,
    so an order-of-magnitude estimate only. Returns max |dphi_X| and ||dphi_X|| / ||phi_X||."""
    _, rhs, jac = make_rhs_jac(B, N, h, 0.0)
    r = rhs(0.0, y)
    U, s, Vt = np.linalg.svd(jac(0.0, y).toarray())
    k = s > rel_cut * s[0]
    dy = Vt[k].T @ ((U[:, k].T @ (-r)) / s[k])
    f = y[N:2 * N]
    return dict(max_dphiX=float(np.abs(dy[N:2 * N]).max()),
                rel_norm_dphiX=float(np.linalg.norm(dy[N:2 * N]) / np.linalg.norm(f)))


def transfer_checks(with_relax=False):
    """§2 state-transfer checks (no stepping). Returns dict with ok and the values."""
    src = source_state()                                          # raises on a hash mismatch
    y0 = nyquist_filter(src, N02)
    st = start_states(y0)
    out = dict(source_sha256=sha(src), y0_sha256=sha(y0), y0_max_change=float(np.abs(y0 - src).max()))
    rt = {}
    for k, s in SHIFTS.items():
        rt[k] = float(np.abs(R.fourier_shift(st[k]["y"], N02, H02, -s) - y0).max())
    rt["R3"] = float(np.abs(reflect(st["R3"]["y"], N02) - y0).max())
    out["roundtrip"] = rt
    out["R4_back"] = float(np.abs(from_d01_nodes(st["R4"]["y"]) - y0).max())
    out["rhs_max"] = dict(source=rhs_max(src, N02, H02), **{k: rhs_max(v["y"], v["N"], v["h"]) for k, v in st.items()})
    out["start_sha256"] = {k: sha(v["y"]) for k, v in st.items()}
    if with_relax:
        out["relax_estimate"] = {k: relax_estimate(st[k]["y"], st[k]["N"], st[k]["h"]) for k in ("R0", "R2a", "R2b",
                                                                                                  "R4")}
    out["ok"] = bool(out["y0_max_change"] <= Y0_MAX_CHANGE and max(rt.values()) <= ROUNDTRIP_TOL
                     and out["R4_back"] <= R4_BACK_TOL)
    return out


# ------------------------------------------------------------------ diagnostics (§3), any N and h
def xgrid(N, h):
    return (np.arange(N) - N // 2) * h


def shift_field(f, h, s):
    N = len(f)
    k = 2 * np.pi * np.fft.rfftfreq(N, d=h)
    return np.fft.irfft(np.fft.rfft(f) * np.exp(-1j * k * s), n=N)


def centres(phiX, h):
    N = len(phiX)
    x = xgrid(N, h)
    E = R.envelope(phiX, h)
    wE, w2 = E ** 2, phiX ** 2
    return float(np.sum(x * wE) / np.sum(wE)), float(np.sum(x * w2) / np.sum(w2))


def minima(phiX, h, E_pat):
    """Minima of phi_X with prominence >= E_pat/2 (as cell_count), sub-grid by a 3-point parabola, sorted."""
    N = len(phiX)
    x = xgrid(N, h)
    ext = np.concatenate([phiX, phiX, phiX])
    pk, _ = find_peaks(-ext, prominence=0.5 * E_pat, distance=max(1, int(np.ceil(LAM / 2 / h))))
    pk = pk[(pk >= N) & (pk < 2 * N)] - N
    xs = []
    for j in pk:
        a, b, c = phiX[(j - 1) % N], phiX[j], phiX[(j + 1) % N]
        d = a - 2 * b + c
        off = 0.5 * (a - c) / d if d != 0 else 0.0
        xs.append(x[j] + off * h)
    return np.sort(np.array(xs))


def match_minima(ref, m):
    """Each reference minimum -> nearest current minimum if within lam/4, else None."""
    out = []
    for x0 in ref:
        if len(m) == 0:
            out.append(None)
            continue
        j = int(np.argmin(np.abs(m - x0)))
        out.append(float(m[j]) if abs(m[j] - x0) < LAM / 4 else None)
    return out


def align(f0, f, h, centre):
    """d minimising ||f - T_d f0||, bounded to centre +- lam/4; ambiguous within AMBIG_TOL of a bound."""
    obj = lambda s: float(np.sum((f - shift_field(f0, h, s)) ** 2))
    lo, hi = centre - LAM / 4, centre + LAM / 4
    r = minimize_scalar(obj, bounds=(lo, hi), method="bounded", options=dict(xatol=1e-10))
    d = float(r.x)
    S = float(np.sqrt(obj(d)) / np.linalg.norm(f0))
    return d, S, bool(min(d - lo, hi - d) < AMBIG_TOL)


class Diag:
    """Per-run diagnostics, measured from the run's own t' = 0 state."""

    def __init__(self, phiX0, h, E_pat, edges0):
        self.h, self.E_pat, self.f0 = h, E_pat, np.array(phiX0, float)
        m0 = minima(self.f0, h, E_pat)
        self.ref = m0[1:-1]                                       # interior minima (all but the outermost each side)
        self.c0 = centres(self.f0, h)
        self.e0 = edges0
        self.d_prev = 0.0
        self.f_ref = self.d_ref = None

    def sample(self, t, phiX, edges):
        d, S, amb = align(self.f0, phiX, self.h, self.d_prev)
        self.d_prev = d
        if self.f_ref is None and t >= VALID_REF_T - 1e-9:
            self.f_ref, self.d_ref = np.array(phiX, float), d
        S_ref = None
        if self.f_ref is not None:
            S_ref = float(np.linalg.norm(phiX - shift_field(self.f_ref, self.h, d - self.d_ref))
                          / np.linalg.norm(self.f_ref))
        cA, c2 = centres(phiX, self.h)
        xl, xr = edges if edges is not None else (None, None)
        ok_e = xl is not None and self.e0[0] is not None
        return dict(t=float(t), d=d, S=S, S_ref=S_ref, ambiguous=amb,
                    C_edge=((xl + xr) / 2 - (self.e0[0] + self.e0[1]) / 2) if ok_e else None,
                    W=(xr - xl) if xl is not None else None,
                    C_amp=cA - self.c0[0], C_phi2=c2 - self.c0[1],
                    minima=match_minima(self.ref, minima(phiX, self.h, self.E_pat)))


def minima_series(diag, ref):
    """m(t'): mean shift over ONE subset, the reference minima matched at every sample."""
    keep = [i for i in range(len(ref)) if all(s["minima"][i] is not None for s in diag)]
    m = [float(np.mean([s["minima"][i] - ref[i] for i in keep])) if keep else None for s in diag]
    return m, keep


# ------------------------------------------------------------------ static measurement floor (§3)
def recover(f0, f, h, E_pat, N):
    """All five measures of the displacement of f relative to f0 (one-shot, as at a single sample)."""
    tr0, tr1 = R.Tracker(N, h, E_pat), R.Tracker(N, h, E_pat)
    r0, r1 = tr0.sample(f0), tr1.sample(f)
    dg = Diag(f0, h, E_pat, (r0["xL"], r0["xR"]))
    s = dg.sample(0.0, f, (r1["xL"], r1["xR"]))
    m, keep = minima_series([s], dg.ref)
    return dict(d=s["d"], C_edge=s["C_edge"], C_amp=s["C_amp"], C_phi2=s["C_phi2"], m=m[0], S=s["S"],
                ambiguous=s["ambiguous"], matched=len(keep), of=len(dg.ref))


def floor_validation():
    y0 = common_start()
    E = e_pat()
    out = dict(E_pat=E, cases={})
    worst_d, per_measure = 0.0, {k: 0.0 for k in ("d", "C_edge", "C_amp", "C_phi2", "m")}
    for grid, (y, N, h) in dict(D02=(y0, N02, H02), D01=(to_d01(y0), N01, H01)).items():
        f0 = y[N:2 * N]
        for s0 in FLOOR_SHIFTS:
            for s in (s0, -s0):
                f = shift_field(f0, h, s)
                r = recover(f0, f, h, E, N)
                err = {k: (abs(r[k] - s) if r[k] is not None else None) for k in per_measure}
                out["cases"][f"{grid}_{s:+.4f}"] = dict(true=s, **r, err=err)
                for k, e in err.items():
                    per_measure[k] = max(per_measure[k], e if e is not None else np.inf)
                worst_d = max(worst_d, err["d"])
    out["delta"] = worst_d
    out["max_err_per_measure"] = per_measure
    out["v_f"] = vf_from_delta(worst_d)
    out["resolved_threshold"] = RESOLVED_VF * out["v_f"]
    out["consistency_band_at_d0"] = CONS_ABS
    return out


def vf_from_delta(delta):
    return max(VF_MIN, 4.0 * delta / 8000.0)


# ------------------------------------------------------------------ ledger
def _ledger_path():
    return os.path.join(OUT, "ledger.json")


def ledger():
    try:
        with open(_ledger_path()) as f:
            return json.load(f)
    except FileNotFoundError:
        return dict(integrations=[], cap=RUN_CAP, cutoff_core_s=CUTOFF_CORE_S)


def _ledger_write(lg):
    os.makedirs(OUT, exist_ok=True)
    with open(_ledger_path(), "w") as f:
        json.dump(lg, f, indent=1)


def _ledger_start(key):
    """Counted at start, so an interrupted integration still counts toward the cap."""
    lg = ledger()
    n = len(lg["integrations"])
    used = sum(r["wall_s"] or 0.0 for r in lg["integrations"])
    if n + 1 > RUN_CAP:
        raise Stop(f"run cap {RUN_CAP} reached")
    if used >= CUTOFF_CORE_S:
        raise Stop(f"cutoff {CUTOFF_CORE_S / 3600:.1f} core-h reached")
    lg["integrations"].append(dict(key=key, wall_s=None, started=time.strftime("%Y-%m-%dT%H:%M:%S")))
    _ledger_write(lg)
    return len(lg["integrations"]) - 1


def _ledger_finish(i, wall):
    lg = ledger()
    lg["integrations"][i]["wall_s"] = wall
    _ledger_write(lg)


# ------------------------------------------------------------------ records
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
        json.dump(rec, f, separators=(",", ":"), default=float)


def save_snap(name, y):
    os.makedirs(SNAP, exist_ok=True)
    np.save(os.path.join(SNAP, name + ".npy"), np.asarray(y, float))


# ------------------------------------------------------------------ one run (§2-§3)
def one_run(key, st, E_pat, ode=None):
    rec = load_rec(key)
    if rec is not None:
        return rec
    ode = ode or mg9d_ode
    N, h, y = st["N"], st["h"], np.array(st["y"], float)
    i_led = _ledger_start(key)
    t0w = time.time()
    tr = R.Tracker(N, h, E_pat)
    r0 = tr.sample(y[N:2 * N])
    track = [dict(t=0.0, **r0)]
    dg = Diag(y[N:2 * N], h, E_pat, (r0["xL"], r0["xR"]))
    diag = [dg.sample(0.0, y[N:2 * N], (r0["xL"], r0["xR"]))]
    snaps = [dict(name=f"{key}_t0", t=0.0, sha256=sha(y))]
    save_snap(f"{key}_t0", y)
    r, ta = None, 0.0
    while ta < T_END - 1e-9:
        tb = ta + SEG
        te = np.arange(ta + DT_TRACK, tb + 1e-9, DT_TRACK)
        r = ode(B, N, h, tb, y, seed_scale=0.0, t_eval=te, t_start=ta, **st["solver"])
        for j, tj in enumerate(r["t"]):
            f = r["y"][N:2 * N, j]
            s = tr.sample(f)
            track.append(dict(t=float(tj), **s))
            if abs(tj / DT_DIAG - round(tj / DT_DIAG)) < 1e-9:
                diag.append(dg.sample(tj, f, (s["xL"], s["xR"])))
        if not r["complete"]:
            break
        y = r["y_final"]
        save_snap(f"{key}_t{int(round(tb))}", y)
        snaps.append(dict(name=f"{key}_t{int(round(tb))}", t=tb, sha256=sha(y)))
        ta = tb
    wall = time.time() - t0w
    _ledger_finish(i_led, wall)
    rec = dict(key=key, b=B, N=N, h=h, solver=st["solver"], shift=st.get("shift"), E_pat=E_pat,
               init_sha256=sha(st["y"]), complete=bool(r["complete"]), success=bool(r["success"]),
               t_reached=float(r["t_reached"]), residual=float(r["residual"]), retry=r["retry"], seed_scale=0.0,
               wall_s=wall, ref_minima=[float(v) for v in dg.ref], snapshots=snaps, track=track, diag=diag,
               platform=platform())
    save_rec(key, rec)
    print(key, "complete" if rec["complete"] else "INCOMPLETE", f"{wall:.0f}s", flush=True)
    return rec


# ------------------------------------------------------------------ validity and speeds (§3)
def _at(diag, t):
    for s in diag:
        if abs(s["t"] - t) < 1e-6:
            return s
    return None


def validity(rec, t_ref=None):
    """§3 items 1-6. Items 5-6 are measured over [t_ref, T] (t_ref = VALID_REF_T as proposed; 0 = draft 3 as
    written). Returns dict(valid, failed=[...], values)."""
    t_ref = VALID_REF_T if t_ref is None else t_ref
    failed, val = [], {}
    diag, track = rec["diag"], rec["track"]
    n_exp = int(round(T_END / DT_DIAG)) + 1
    if not rec["complete"] or len(diag) != n_exp:
        failed.append(1)
    bad = [s["t"] for s in track if s["status"] != "ok" or s["breakup"] or s["remote"]]
    val["tracker_bad_samples"] = len(bad)
    if bad:
        failed.append(2)
    n0 = track[0]["n"]
    ns = [s["n"] for s in track]
    val["cells"] = [n0, min(ns), max(ns)]
    W0 = diag[0]["W"]
    WT = diag[-1]["W"] if diag else None
    val["dW"] = (WT - W0) if (W0 is not None and WT is not None) else None
    if min(ns) < n0 - 1 or max(ns) > n0 + 1 or val["dW"] is None or abs(val["dW"]) >= LAM / 2:
        failed.append(3)
    val["ambiguous_samples"] = sum(1 for s in diag if s["ambiguous"])
    if val["ambiguous_samples"]:
        failed.append(4)
    val["S_T_from_0"] = diag[-1]["S"]
    val["S_T"] = diag[-1]["S"] if t_ref == 0 else diag[-1].get("S_ref")
    val["t_ref"] = t_ref
    if not (val["S_T"] is not None and val["S_T"] <= S_MAX):
        failed.append(5)
    m, keep = minima_series(diag, rec["ref_minima"])
    val["minima_matched"] = [len(keep), len(rec["ref_minima"])]
    i_ref = next((i for i, x in enumerate(diag) if abs(x["t"] - t_ref) < 1e-6), None)
    if i_ref is None:
        failed.append(6)
    else:
        A, Z = diag[i_ref], diag[-1]
        dd = Z["d"] - A["d"]
        band = CONS_REL * abs(dd) + CONS_ABS
        diff = lambda k: (Z[k] - A[k]) if (Z[k] is not None and A[k] is not None) else None
        cons = dict(C_edge=diff("C_edge"), C_amp=diff("C_amp"), C_phi2=diff("C_phi2"),
                    m=(m[-1] - m[i_ref]) if (m[-1] is not None and m[i_ref] is not None) else None)
        val["consistency"] = dict(d=dd, band=band, **cons)
        if any(v is None or abs(v - dd) > band for v in cons.values()):
            failed.append(6)
    return dict(valid=not failed, failed=sorted(set(failed)), values=val)


def speeds(diag):
    """Endpoint slopes of d: a over [2000, 6000], b over [6000, 10000], v over [2000, 10000]."""
    def sl(w):
        A, Bq = _at(diag, w[0]), _at(diag, w[1])
        return (Bq["d"] - A["d"]) / (w[1] - w[0]) if (A and Bq) else None
    return dict(a=sl(WIN_A), b=sl(WIN_B), v=sl(WIN_V))


# ------------------------------------------------------------------ decision rules (§4)
def resolved(v, vf):
    return v is not None and abs(v) >= RESOLVED_VF * vf


def compare(A, Bref, vf):
    d = abs(A - Bref)
    if d <= AGREE_REL * abs(Bref) + vf:
        return "AGREE"
    if d > DIFFER_REL * abs(Bref) + DIFFER_VF * vf:
        return "DIFFER"
    return "INCONCLUSIVE"


def c1(sp0, sp1, valid0, valid1, vf):
    if not (valid0 and valid1):
        return "INCONCLUSIVE"
    cmp = [compare(sp1[k], sp0[k], vf) for k in ("a", "b", "v")]
    if all(c == "AGREE" for c in cmp):
        return "NUMERICS-AGREE"
    if any(c == "DIFFER" for c in cmp):
        return "NUMERICS-SENSITIVE"
    return "INCONCLUSIVE"


def reflect_adjust(v_raw):
    """§3: applied once, in the evaluator only."""
    return None if v_raw is None else -v_raw


def c3(v0, v3_raw, valid0, valid3, vf):
    if not (valid0 and valid3):
        return "NOT EVALUATED"
    v3 = reflect_adjust(v3_raw)
    return "SYMMETRIC" if abs(v3 - v0) <= SYM_REL * abs(v0) + vf else "ASYMMETRY-FLAG"


def c2(v0, valid0, arms, vf):
    """arms: list of (v, valid) for R2a-c."""
    if not valid0:
        return "INCONCLUSIVE"
    if any(ok and compare(v, v0, vf) == "DIFFER" for v, ok in arms):
        return "ALIGNMENT-DEPENDENT"
    if all(ok for _, ok in arms) and all(compare(v, v0, vf) == "AGREE" for v, _ in arms):
        return "ALIGNMENT-INDEPENDENT"
    return "INCONCLUSIVE"


def c4_label(a, b, vf):
    ra, rb = resolved(a, vf), resolved(b, vf)
    if not ra and not rb:
        return "NO RESOLVED MOTION"
    if not ra and rb:
        return "ONSET"
    if ra and rb and np.sign(a) != np.sign(b):
        return "REVERSAL"
    if ra and not rb:
        return "DECREASED BELOW RESOLUTION"
    rho = b / a
    if rho <= 0.90:
        return "DECREASING"
    if 0.95 <= rho <= 1.05:
        return "PERSISTENT"
    if rho >= 1.10:
        return "INCREASING"
    return "INCONCLUSIVE"


def c4(sp0, sp1, valid0, valid1, c1_label, vf):
    if not valid0:
        return "INCONCLUSIVE"
    lab = c4_label(sp0["a"], sp0["b"], vf)
    if c1_label == "NUMERICS-AGREE":
        if not valid1 or c4_label(sp1["a"], sp1["b"], vf) != lab:
            return "INCONCLUSIVE"
    return lab


def c5(v0, v4, valid4, vf):
    if not valid4:
        return "R4 INVALID / LOST-IDENTITY"
    if not resolved(v0, vf):
        return "NOT EVALUATED"
    if resolved(v4, vf) and np.sign(v4) != np.sign(v0):
        return "REVERSED"
    if abs(v4 - v0) <= 0.20 * abs(v0) + vf:
        return "SIMILAR"
    if abs(v4) < 0.5 * abs(v0):
        return "REDUCED"
    return "OTHER"


def overall(lab, sp0, sp1, valid0, valid1, vf):
    """lab: dict C1..C5. Ordered rules 1-7 of §4."""
    if not valid0 or lab["C3"] in ("ASYMMETRY-FLAG", "NOT EVALUATED") or lab["C1"] == "INCONCLUSIVE":
        return "INCONCLUSIVE"
    if lab["C1"] == "NUMERICS-SENSITIVE":
        return "NUMERICAL ACCURACY UNRESOLVED"
    if valid1 and all(not resolved(s[k], vf) for s in (sp0, sp1) for k in ("a", "b", "v")):
        return "NO RESOLVED DRIFT"
    if lab["C2"] == "ALIGNMENT-DEPENDENT" or lab["C5"] in ("REDUCED", "REVERSED"):
        return "MESH EFFECT SUPPORTED"
    unavailable = lab["C5"] in ("R4 INVALID / LOST-IDENTITY", "NOT EVALUATED")
    qual = " (on D02; the finer-mesh comparison is unavailable)" if unavailable else ""
    if lab["C2"] == "ALIGNMENT-INDEPENDENT" and (lab["C5"] == "SIMILAR" or unavailable):
        if lab["C4"] in ("DECREASING", "DECREASED BELOW RESOLUTION"):
            return f"RELAXATION SUPPORTED [C4 {lab['C4']}]" + qual
        if lab["C4"] == "PERSISTENT":
            return "MOVING STATE SUPPORTED" + qual
    return f"INCONCLUSIVE (C4 {lab['C4']})"


def evaluate(recs, vf):
    """recs: dict run -> record (or None). Returns the full evaluation."""
    V = {k: (validity(r) if r else dict(valid=False, failed=["not run"], values={})) for k, r in recs.items()}
    S = {k: (speeds(r["diag"]) if r else dict(a=None, b=None, v=None)) for k, r in recs.items()}
    ok = {k: V[k]["valid"] for k in recs}
    lab = {}
    lab["C1"] = c1(S["R0"], S["R1"], ok["R0"], ok["R1"], vf)
    lab["C3"] = c3(S["R0"]["v"], S["R3"]["v"], ok["R0"], ok["R3"], vf)
    lab["C2"] = c2(S["R0"]["v"], ok["R0"], [(S[k]["v"], ok[k]) for k in ("R2a", "R2b", "R2c")], vf)
    lab["C4"] = c4(S["R0"], S["R1"], ok["R0"], ok["R1"], lab["C1"], vf)
    lab["C5"] = c5(S["R0"]["v"], S["R4"]["v"], ok["R4"], vf)
    reading = overall(lab, S["R0"], S["R1"], ok["R0"], ok["R1"], vf)
    S_adj = dict(S)
    S_adj["R3_reflection_adjusted"] = {k: reflect_adjust(v) for k, v in S["R3"].items()}
    return dict(reading=reading, components=lab, validity=V, speeds_raw=S, speeds_R3_adjusted=S_adj[
        "R3_reflection_adjusted"], v_f=vf, resolved_threshold=RESOLVED_VF * vf)


# ------------------------------------------------------------------ stages
def stage_floor():
    out = floor_validation()
    out["platform"] = platform()
    os.makedirs(OUT, exist_ok=True)
    with open(os.path.join(OUT, "floor.json"), "w") as f:
        json.dump(out, f, indent=1, default=float)
    print(json.dumps(dict(delta=out["delta"], v_f=out["v_f"], max_err_per_measure=out["max_err_per_measure"]),
                     indent=1, default=float))
    return out


def stage_checks(with_relax=True):
    out = transfer_checks(with_relax)
    out["platform"] = platform()
    os.makedirs(OUT, exist_ok=True)
    with open(os.path.join(OUT, "transfer.json"), "w") as f:
        json.dump(out, f, indent=1, default=float)
    print(json.dumps({k: v for k, v in out.items() if k not in ("start_sha256",)}, indent=1, default=float))
    return out


def stage_runs():
    if not os.path.exists(REGISTERED_FLAG):
        raise Stop("MG9d is not registered: results/MG9d/REGISTERED is missing")
    chk = stage_checks(with_relax=False)
    if not chk["ok"]:
        raise Stop("state-transfer checks failed; nothing integrated")
    with open(os.path.join(OUT, "floor.json")) as f:
        fl = json.load(f)
    E = e_pat()
    st = start_states()
    for k in ORDER:
        one_run(k, st[k], E)
    recs = {k: load_rec(k) for k in ORDER}
    ev = evaluate(recs, fl["v_f"])
    lg = ledger()
    ev["integrations"] = len(lg["integrations"])
    ev["core_h"] = sum(r["wall_s"] or 0.0 for r in lg["integrations"]) / 3600.0
    ev["platform"] = platform()
    with open(os.path.join(OUT, "evaluation.json"), "w") as f:
        json.dump(ev, f, indent=1, default=float)
    print(json.dumps(dict(reading=ev["reading"], components=ev["components"]), indent=1))


def stage_pack():
    names = sorted(n[:-4] for n in os.listdir(SNAP) if n.endswith(".npy"))
    arrs = {n: np.load(os.path.join(SNAP, n + ".npy")) for n in names}
    np.savez_compressed(os.path.join(OUT, "snapshots.npz"), **arrs)
    with open(os.path.join(OUT, "snapshots_manifest.json"), "w") as f:
        json.dump({n: sha(a) for n, a in arrs.items()}, f, indent=1)
    print(len(names), "snapshots packed")


STAGES = dict(floor=stage_floor, checks=stage_checks, runs=stage_runs, pack=stage_pack)

if __name__ == "__main__":
    try:
        STAGES[sys.argv[1]]()
    except Stop as e:
        print("STOP:", e)
        sys.exit(2)

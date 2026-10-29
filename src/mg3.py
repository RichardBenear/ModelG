"""MG3 library (protocol v1.3.4; GATES.md MG3 incl. the pre-run implementation specification of 2026-10-28).

Contents
- sim_run: one homogeneous-start trajectory, simulated in 10-time-unit segments with online detection, tracking and
  field-level confirmation (protocol §5.4) of every persistent birth. Stops 50 time units after the run's first
  confirmed event (PI decision 2026-10-28), else at T = 500. Resumable from a saved state (until=...).
- first-passage statistics (protocol §6): exposure, MLE, Garwood CI, stationary-hazard check, conditional-binomial
  rate ratio, activated-scaling fit.
- Stage-A outcome counting and Stage-B cell selection (GATES.md MG3 with Amendment 4 §3).
Stream keys (protocol §7): seed_for(3, b_idx, om_idx, stage, replicate, segment); segment = int(t0 / SEG).
"""
import os, sys, time, math
import numpy as np
from scipy import stats
sys.path.insert(0, os.path.dirname(__file__))
from lattice import tau_leap_euler_1d, P17, steady_state
from measure import smooth_periodic, detect, track, fourier_shift, circ
from confirm import confirm_field
from runner import reference, seed_for

H, L, N = 0.2, 100.0, 500
SEG = 10.0                    # segment length (time units)
T_END, T_B, T_C = 500.0, 50.0, 480.0
PERSIST = 20.0
POST = 50.0                   # continue this long after the first confirmed event, then stop
PROF_WIN = (30.0, 50.0)       # C2(b)/(c) profile window after detection
PROF_K = 4                    # blocks in the profile window (5 time units each)
PROF_MIN = 7                  # eligible snapshots required per block (of 10)
WIN = 10.0                    # |x| < WIN profile window (as MG2)


def _x():
    return (np.arange(N) - N // 2) * H


# ------------------------------------------------------------------------------------------------ simulation
def new_state(b, Om, tau, key, b_det=29.0):
    """key = (gate, b_idx, om_idx, stage, replicate). Homogeneous Poisson start at the b steady state."""
    S0 = np.array(steady_state({**P17, "b": b}))
    rs = np.random.default_rng(seed_for(*key, segment=0) ^ 0x5A5A5A5A)
    n = [rs.poisson(Om * S0[k], N).astype(np.int64) for k in range(3)]
    return dict(key=tuple(key), b=b, b_det=b_det, Om=float(Om), tau=tau, t=0.0, nG=n[0], nX=n[1], nY=n[2],
                phi=[], times=[], dets=[], births=[], conf_cache={}, first=None, stop=T_END, done=False,
                wall=0.0, clips=0, steps=0)


def sim_run(st, until=T_END):
    """Advance state st to min(until, st['stop']); returns st. Pure function of (st, until) given the streams."""
    if st["done"]:
        return st
    ref = reference(st["b_det"], H)
    x = ref["x"]
    S0 = np.array(steady_state({**P17, "b": st["b"]}))
    P = {**P17, "b": st["b"]}
    tau = st["tau"]
    rec = int(round(0.5 / tau))
    w0 = time.time()
    while st["t"] < min(until, st["stop"]) - 1e-9:
        seg = int(round(st["t"] / SEG))
        nsteps = int(round(SEG / tau))
        rG, rX, rY, rt, clips = tau_leap_euler_1d(st["nG"], st["nX"], st["nY"], st["Om"], P["a"], P["b"], P["g"],
                                                  P["p"], P["q"], P["dx"], P["dy"], H, tau, nsteps, st["t"], 0.0, x,
                                                  rec, seed_for(*st["key"], segment=seg))
        st["clips"] += int(clips); st["steps"] += nsteps
        st["t"] = (seg + 1) * SEG                     # exact segment boundary (avoids float drift)
        tt = np.round(st["t"] - SEG + 0.5 * np.arange(1, rt.size + 1), 6)
        phi = np.stack([rG / st["Om"] - S0[0], rX / st["Om"] - S0[1], rY / st["Om"] - S0[2]], axis=1).astype(np.float32)
        ps = smooth_periodic(phi[:, 1, :].astype(float), H)
        for s in range(rt.size):
            st["phi"].append(phi[s]); st["times"].append(float(tt[s]))
            st["dets"].append(detect(ps[s], x, H, ref["thr"], ref["lam"]))
        _process_births(st, ref)
    st["wall"] += time.time() - w0
    if st["t"] >= st["stop"] - 1e-9 or st["t"] >= T_END - 1e-9:
        st["done"] = True
    return st


def _process_births(st, ref):
    if st["first"] is not None:      # births after the first confirmed event are not used (C2a counts up to it)
        return
    times = np.array(st["times"])
    tr = track(st["dets"], times, L, x0=None)
    seen = {(b["t_start"], b["track"]) for b in st["births"]}
    cand = sorted((q for q in tr["tracks"] if q["t_last"] - q["t_start"] >= PERSIST - 1e-9),
                  key=lambda q: (q["t_start"], q["id"]))
    for q in cand:
        k = (float(q["t_start"]), int(q["id"]))
        if k in seen:
            continue
        s0 = q["hits"][0][0]
        sc = int(np.argmin(np.abs(times - (q["t_start"] + PERSIST))))
        pos = float(st["dets"][s0][q["hits"][0][1]])
        if sc not in st["conf_cache"]:
            c = confirm_field(st["phi"][sc].astype(float), st["b"], H, ref["A"], ref["lam"], det_pos=pos)
            st["conf_cache"][sc] = dict(cls=c["cls"], n_persistent=c["n_persistent"], n_deep_end=c["n_deep_end"],
                                        min_phiX_end=c["min_phiX_end"], maxdt_end=c["maxdt_end"])
            shared = False
        else:
            shared = True                                    # outcome shared (counted once, protocol §5.4)
        cc = st["conf_cache"][sc]
        near = None
        if not shared:
            near = bool(c["near_detection"])
        b = dict(t_start=k[0], track=k[1], pos=pos, s_conf=sc, t_conf=float(times[sc]), cls=cc["cls"],
                 shared=shared, near=near)
        st["births"].append(b)
        if cc["cls"] == "persistent":
            st["first"] = dict(b)
            st["first"]["field"] = st["phi"][sc].copy()
            st["stop"] = min(k[0] + POST, T_END)
            return


def summarize(st):
    """Compact per-run record (no snapshot arrays except the first-event confirmation field)."""
    nd = np.array([len(d) for d in st["dets"]])
    times = np.array(st["times"])
    fe = st["first"]
    rec = dict(key=st["key"], b=st["b"], b_det=st["b_det"], Om=st["Om"], tau=st["tau"], t_end=st["t"],
               wall=st["wall"], clip_frac=st["clips"] / (3 * N * max(st["steps"], 1)),
               mean_dets_t_lt_50=float(nd[times < T_B].mean()) if np.any(times < T_B) else float("nan"),
               n_dets=nd.astype(np.int16), births=[{k: v for k, v in b.items()} for b in st["births"]],
               first=None if fe is None else {k: v for k, v in fe.items()})
    rec.update(first_passage(None if fe is None else fe["t_start"], st["t"]))
    rec["profile"] = nucleated_profile(st) if fe is not None else None
    return rec


def first_passage(t_first, t_reached, t_b=T_B, t_c=T_C):
    """Protocol §6 for one run. t_first: detection time of the first confirmed event (None if none).
    Returns dict(startup, event (counted, in [t_b, t_c]), exposure, expo_w1 ([t_b,200)), expo_w2 ([200,t_c]),
    in_w1, complete (the run reached the time needed to classify it))."""
    if t_first is not None and t_first < t_b:
        return dict(startup=True, event=False, exposure=0.0, expo_w1=0.0, expo_w2=0.0, in_w1=False, complete=True)
    if t_first is not None and t_first <= t_c:
        e = t_first
        return dict(startup=False, event=True, exposure=e - t_b, expo_w1=min(e, 200.0) - t_b,
                    expo_w2=max(0.0, e - 200.0), in_w1=e < 200.0, complete=True)
    complete = t_reached >= T_END - 1e-9 or t_first is not None
    return dict(startup=False, event=False, exposure=t_c - t_b, expo_w1=200.0 - t_b, expo_w2=t_c - 200.0,
                in_w1=False, complete=complete)


def nucleated_profile(st):
    """C2(b)/(c): aligned profile of the first nucleated particle over [t_e+30, t_e+50] in PROF_K blocks.
    Eligible snapshot: exactly one detection and it belongs to the nucleated track. Run eligible if every block
    has >= PROF_MIN eligible snapshots (survival-conditioned)."""
    fe = st["first"]
    times = np.array(st["times"])
    tr = track(st["dets"], times, L, x0=None)
    q = next(q for q in tr["tracks"] if q["id"] == fe["track"] and abs(q["t_start"] - fe["t_start"]) < 1e-9)
    hit = {s: j for s, j in q["hits"]}
    x = _x(); win = np.abs(x) < WIN
    edges = fe["t_start"] + np.linspace(PROF_WIN[0], PROF_WIN[1], PROF_K + 1)
    blocks = np.zeros((PROF_K, 3, int(win.sum()))); full = np.zeros((3, N)); nb = np.zeros(PROF_K, int)
    for s, t in enumerate(times):
        if not (edges[0] - 1e-9 <= t < edges[-1] - 1e-9 or abs(t - edges[-1]) < 1e-9):
            continue
        if len(st["dets"][s]) != 1 or s not in hit:
            continue
        kb = min(int(np.searchsorted(edges, t + 1e-9, side="right") - 1), PROF_K - 1)
        al = fourier_shift(st["phi"][s].astype(float), H, float(st["dets"][s][0]))
        blocks[kb] += al[:, win]; full += al; nb[kb] += 1
    ok = bool(np.all(nb >= PROF_MIN))
    for k in range(PROF_K):
        if nb[k]:
            blocks[k] /= nb[k]
    if nb.sum():
        full /= nb.sum()
    return dict(eligible=ok, n_per_block=nb.tolist(), blocks=blocks.astype(np.float32), full=full.astype(np.float32))


def wavelength_core(x, phiX, lo=1.0, hi=7.0):
    """Pre-stated noise-robust wavelength: 2 x mean spacing of zero crossings of the x-symmetrised profile in
    lo < x < hi. Applied identically to the deterministic soliton and to nucleated mean profiles."""
    i0 = x.size // 2
    f = 0.5 * (phiX[i0:] + phiX[i0::-1][: x.size - i0])
    r = x[i0:]
    sg = np.sign(f)
    idx = np.where(sg[:-1] * sg[1:] < 0)[0]
    z = r[idx] - f[idx] * (r[idx + 1] - r[idx]) / (f[idx + 1] - f[idx])
    z = z[(z > lo) & (z < hi)]
    return 2.0 * float(np.mean(np.diff(z))) if z.size >= 2 else float("nan")


# ------------------------------------------------------------------------------------------------ statistics
def garwood(n, level=0.95):
    a = 1 - level
    lo = 0.0 if n == 0 else stats.chi2.ppf(a / 2, 2 * n) / 2
    hi = stats.chi2.ppf(1 - a / 2, 2 * n + 2) / 2
    return lo, hi


def cell_rate(runs, L_=L):
    """Constant-hazard MLE per unit length and time from first-passage records (complete runs only)."""
    rr = [r for r in runs if r["complete"]]
    n = sum(r["event"] for r in rr)
    E = sum(r["exposure"] for r in rr)
    lo, hi = garwood(n)
    k = n / (L_ * E) if E > 0 else float("nan")
    out = dict(n=int(n), exposure=float(E), k=k, k_lo=lo / (L_ * E) if E > 0 else float("nan"),
               k_hi=hi / (L_ * E) if E > 0 else float("nan"),
               k_upper95_onesided=2.995732 / (L_ * E) if (n == 0 and E > 0) else None,
               n_startup=int(sum(r["startup"] for r in rr)), n_runs=len(rr), n_incomplete=len(runs) - len(rr))
    out.update(stationary_check(rr))
    return out


def stationary_check(runs):
    """Conditional binomial test: events in [t_b,200) vs [200,t_c] given exposures (two-sided exact)."""
    n1 = sum(r["event"] and r["in_w1"] for r in runs)
    n2 = sum(r["event"] and not r["in_w1"] for r in runs)
    E1 = sum(r["expo_w1"] for r in runs); E2 = sum(r["expo_w2"] for r in runs)
    if n1 + n2 == 0 or E1 + E2 == 0:
        return dict(stat_n1=n1, stat_n2=n2, stat_p=float("nan"), stationary_ok=True)
    p = stats.binomtest(n1, n1 + n2, E1 / (E1 + E2)).pvalue
    return dict(stat_n1=int(n1), stat_n2=int(n2), stat_E1=E1, stat_E2=E2, stat_p=float(p), stationary_ok=bool(p >= 0.01))


def rate_ratio(n1, E1, n2, E2, level=0.95):
    """Ratio k1/k2 with the exact conditional-binomial (Clopper-Pearson) CI."""
    n = n1 + n2
    if n == 0:
        return dict(ratio=float("nan"), lo=0.0, hi=float("inf"))
    ci = stats.binomtest(n1, n).proportion_ci(confidence_level=level, method="exact")
    f = lambda p: (p / (1 - p)) * (E2 / E1) if p < 1 else float("inf")
    ratio = (n1 / E1) / (n2 / E2) if n2 > 0 else float("inf")
    return dict(ratio=ratio, lo=f(ci.low), hi=f(ci.high))


def activated_fit(cells):
    """C3: weighted LS of ln k vs Omega, weights = n (inverse Poisson variance of ln k). cells: list of
    (Omega, n, k). Slope CI from the known-variance SE (normal); t-based CI (df m-2) reported as sensitivity."""
    Om = np.array([c[0] for c in cells], float); n = np.array([c[1] for c in cells], float)
    y = np.log(np.array([c[2] for c in cells], float)); w = n
    xb = np.sum(w * Om) / w.sum(); yb = np.sum(w * y) / w.sum()
    Sxx = np.sum(w * (Om - xb) ** 2)
    slope = np.sum(w * (Om - xb) * (y - yb)) / Sxx
    icpt = yb - slope * xb
    res = y - (icpt + slope * Om)
    r2 = 1 - np.sum(w * res**2) / np.sum(w * (y - yb) ** 2)
    se = math.sqrt(1.0 / Sxx)
    m = len(cells)
    s2 = np.sum(w * res**2) / (m - 2) if m > 2 else float("nan")
    se_t = math.sqrt(s2 / Sxx) if m > 2 else float("nan")
    tq = stats.t.ppf(0.975, m - 2) if m > 2 else float("nan")
    return dict(S=-slope, slope=slope, intercept=icpt, r2=float(r2), ci=(slope - 1.96 * se, slope + 1.96 * se),
                ci_t=(slope - tq * se_t, slope + tq * se_t), m=m)


# ------------------------------------------------------------------------------------------------ adaptive design
def stage_a_count(runs, lo=T_B, hi=150.0):
    """Stage-A statistic: runs whose first confirmed event (post burn-in) is detected in [50, 150]."""
    return sum(1 for r in runs if r["event"] and r["first"] is not None and lo <= r["first"]["t_start"] <= hi)


def round2sf(v):
    return float(f"{v:.2g}")


def select_stage_b(screened, max_cells=4):
    """screened: list of dict(Om, count, excluded) in Stage-A order (descending Omega, midpoints inserted at their
    Omega). Returns (selected Omegas, note). GATES.md MG3 + Amendment 4 §3."""
    ok = sorted([c for c in screened if not c["excluded"]], key=lambda c: -c["Om"])
    if not ok:
        return [], "no screened cell (outcome e)"
    if all(c["count"] == 0 for c in ok) and not any(c.get("any_event") for c in ok):
        return [], "no confirmed event in any screened cell (outcome a)"
    if ok[0]["count"] == 4 and ok[0]["Om"] == 300:
        return [300], "4/4 at the top of the grid (outcome d)"
    S = [c["Om"] for c in ok if 1 <= c["count"] <= 3]
    oms = [c["Om"] for c in ok]
    if len(S) > max_cells:
        S = sorted(S)[:max_cells]
        return sorted(S, reverse=True), "more than 4 partial cells: the 4 lowest Omega kept"
    if len(S) < 3:
        # anchor: S, or the 0/4 -> 4/4 transition (the last 4/4 cell and its upper neighbour)
        four = [c["Om"] for c in ok if c["count"] == 4]
        if not S and four:            # anchor at the 0/4 -> 4/4 transition
            i = oms.index(four[0]); S = [oms[i]] + ([oms[i - 1]] if i > 0 else [])
        elif not S:                   # events only outside [50, 150]: anchor at the lowest-Omega cell with an event
            S = [min(c["Om"] for c in ok if c.get("any_event"))]
        while len(S) < 3:
            cand = [o for o in oms if o not in S]
            if not cand:
                break
            idx = [oms.index(o) for o in S] or [0]
            dist = lambda o: min(abs(oms.index(o) - i) for i in idx)
            cnt = {c["Om"]: c["count"] for c in ok}
            cand.sort(key=lambda o: (dist(o), -cnt[o], o))
            S.append(cand[0])
    return sorted(S, reverse=True), "selected"

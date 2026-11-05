"""MG9b library (GATES.md MG9b; spec/MG9b_proposal.md draft 3, registered 2026-11-05, frozen at 185aceb).

Pure functions shared by the driver (run_mg9b.py), the evaluator (evaluate_mg9b.py) and the synthetic calibration
(calib_mg9b.py):
- design rules: Part A placement (§2.1), Part B candidate order (§2.2), the stop20 batching rule (MG3b), the
  structural identifiability rule (§2.2);
- fits (§4): per-b Poisson GLM log k = alpha_b - S_b * Omega with exposure offset, profile-likelihood CI, deviance
  lack of fit; the trend model S_b = S0 + S1 (b - 29.175) with free alpha_b, profile CI for S1 with fixed brackets
  (§4.4); the 1-df likelihood-ratio check of the trend against separate slopes (§4.2);
- decision rules (§5): primary label, exploratory rule, qualifiers.

Conventions: k is per unit length and time (L = 100); cells are (Omega, n, E) with n post-burn-in events and E the
summed post-burn-in exposure in time units (startup runs contribute zero exposure). S is positive when k falls with
Omega.
"""
import math
import numpy as np
from scipy import stats

L = 100.0
TB, TC = 50.0, 480.0
U = TC - TB
TAU0 = 6.666666666666667e-4
B_CENTRE = 29.175
B_NEW = (29.20, 29.35)
B_HIST = 29.00
B_IDX = {29.00: 0, 29.20: 1, 29.35: 2}
GATE = 9
STAGE = dict(PARTA=0, CLIP=1, SCREEN=2, PROD=3, VB=4)
GRID = (250.0, 200.0, 170.0, 150.0, 140.0, 130.0, 120.0, 110.0, 100.0, 90.0, 80.0)
NA = 4                                   # Part A runs per Omega
BATCH, TARGET, CAP = 16, 20, 64
CANDIDATE_OFFSETS = (-15.0, -30.0, 15.0, 30.0)
S1_BRACKET = (-2.0, 2.0)                 # §4.4: profile search bracket for S1 (per unit Omega per unit b)
SB_BRACKET = (-1.0, 1.0)                 # §4.4: profile search bracket for each S_b (per unit Omega)
CLIP_MAX = 1e-4
CROWD_MAX = 2.0
CUTOFF_CORE_H = 50.0
VB2_LATE = np.arange(1900.0, 2000.0 + 1e-9, 10.0)
VB3_TOL = dict(rtol=1e-10, atol=1e-12)
BRANCH_MIN = 400                         # §4.3: retained replicates for branch support
COV_MIN, SIDE_MAX = 0.93, 0.035


def om_idx(Om):
    """Stream index for an Omega (all MG9b Omega values are whole numbers)."""
    return int(round(float(Om)))


# ------------------------------------------------------------------------------------------------ design rules
def parta_next(counts):
    """Part A (§2.1). counts: list of dicts in grid order, each {Om, screen_ok, count} (count = runs of 4 with a
    confirmed event in [50, 480]; None if the cell was not run because its screen failed before the runs).
    Returns dict(stop, reason, next_Om, Om_T). Called after every cell; next_Om is the next grid Omega to scan."""
    for i, c in enumerate(counts):
        if not c["screen_ok"]:
            return _parta_end(counts[:i], f"screen failed at Omega {c['Om']:g}")
        if c["count"] >= 3:
            return _parta_end(counts[:i + 1], f"3 of 4 reached at Omega {c['Om']:g}")
    if len(counts) >= len(GRID):
        return _parta_end(counts, "end of grid")
    return dict(stop=False, reason=None, next_Om=GRID[len(counts)], Om_T=None)


def _parta_end(valid, reason):
    tops = [c["Om"] for c in valid if c["count"] is not None and c["count"] >= 2]
    return dict(stop=True, reason=reason, next_Om=None, Om_T=max(tops) if tops else None)


def candidates(Om_T):
    """Ordered unique Part B candidate list (§2.2)."""
    return [Om_T + d for d in CANDIDATE_OFFSETS]


def partb_cells(Om_T, screen_results):
    """screen_results: list of bools for the candidates tried so far, in candidate order.
    Returns (accepted cells including Om_T, next candidate to screen or None when done)."""
    cand = candidates(Om_T)
    acc = [Om_T] + [cand[i] for i, ok in enumerate(screen_results) if ok]
    if len(acc) >= 3 or len(screen_results) >= len(cand):
        return acc[:3], None
    return acc, cand[len(screen_results)]


def stop20_need(events_by_rep):
    """MG3b batching. events_by_rep: list of 0/1 post-burn-in events for finished reps 0..m-1 (in rep order).
    Returns (done, runs_needed): done when, at the end of a batch, events >= 20, or at the cap."""
    for b in range(1, CAP // BATCH + 1):
        n = b * BATCH
        if len(events_by_rep) < n:
            return False, n
        if sum(events_by_rep[:n]) >= TARGET or n >= CAP:
            return True, n
    return True, CAP


def identifiable(cells):
    """§2.2 registered slope-eligibility rule (a conservative sufficient condition for a finite slope estimate, not a
    necessary one; MG034 r3): at least two distinct Omega with >= 1 post-burn-in event. Zero-event cells with
    positive exposure stay in an eligible fit."""
    return len({float(Om) for Om, n, E in cells if n >= 1 and E > 0}) >= 2


# ------------------------------------------------------------------------------------------------ likelihood pieces
def _prep(cells):
    Om = np.array([c[0] for c in cells], float)
    n = np.array([c[1] for c in cells], float)
    E = np.array([c[2] for c in cells], float)
    keep = E > 0
    Om, n, E = Om[keep], n[keep], E[keep]
    x = Om - (np.sum(n * Om) / n.sum() if n.sum() > 0 else Om.mean())
    return dict(Om=Om, n=n, E=E, x=x, N=n.sum(), sx=float(np.sum(n * x)))


def _lb(p, S):
    """Profile log likelihood of one b's slope S (alpha_b maximised; constant dropped), with gradient and Hessian.
    log mu_i = alpha - S x_i + log(L E_i)."""
    z = -S * p["x"]
    zm = z.max()
    w = L * p["E"] * np.exp(z - zm)
    W = w.sum()
    m1 = np.sum(w * p["x"]) / W
    m2 = np.sum(w * p["x"] ** 2) / W
    ll = -p["N"] * (math.log(W) + zm) - S * p["sx"]
    g = p["N"] * m1 - p["sx"]
    h = -p["N"] * (m2 - m1 * m1)
    return ll, g, h


def _alpha(p, S):
    """MLE of alpha_b (on the centred x) given S."""
    return math.log(p["N"]) - math.log(np.sum(L * p["E"] * np.exp(-S * p["x"]))) if p["N"] > 0 else -math.inf


def _newton1(f, x0, lo=-5.0, hi=5.0, it=200, tol=1e-12):
    """Maximise a concave 1D function f(x) -> (ll, g, h). Returns (x, ok)."""
    x = x0
    for _ in range(it):
        ll, g, h = f(x)
        if not (math.isfinite(ll) and math.isfinite(g)):
            return x, False
        step = -g / h if h < 0 else math.copysign(0.5, g)
        step = max(-0.5, min(0.5, step))
        for _ in range(60):
            if f(x + step)[0] >= ll - 1e-13:
                break
            step *= 0.5
        x = min(hi, max(lo, x + step))
        if abs(step) < tol:
            return x, True
    return x, abs(f(x)[1]) < 1e-6 * max(1.0, abs(f(x)[0]))


def _bisect(f, a, b, it=80):
    """Root of f on [a, b] with f(a) < 0 < f(b) or f(a) > 0 > f(b)."""
    fa = f(a)
    for _ in range(it):
        m = 0.5 * (a + b)
        fm = f(m)
        if (fm > 0) == (fa > 0):
            a, fa = m, fm
        else:
            b = m
    return 0.5 * (a + b)


def _profile_ci(prof, xhat, lhat, bracket, level=0.95):
    """CI {x : 2 (lhat - prof(x)) <= chi2_1(level)} within the bracket. A side whose boundary lies beyond the bracket
    is unbounded (+-inf), not a failure (§4.4)."""
    c = stats.chi2.ppf(level, 1) / 2
    f = lambda x: (lhat - prof(x)) - c                       # < 0 inside
    lo_b, hi_b = bracket
    lo = -math.inf if f(lo_b) <= 0 else _bisect(f, lo_b, xhat)
    hi = math.inf if f(hi_b) <= 0 else _bisect(f, xhat, hi_b)
    return lo, hi


NOT_LOCATED = "not located within the registered search range; conservatively extended to infinity"


def ci_note(ci):
    """Wording for a bound outside the bracket (MG034 r3): not a proven unbounded interval, not a failure."""
    notes = [f"{side} bound {NOT_LOCATED}" for side, v in (("lower", ci[0]), ("upper", ci[1])) if math.isinf(v)]
    return "; ".join(notes) or None


# ------------------------------------------------------------------------------------------------ fits
def fit_b(cells, with_ci=True):
    """Per-b fit (§4.1). cells: list of (Omega, n, E). Returns dict(S, ci, alpha_at, deviance, df, p_lof, ok)."""
    if not identifiable(cells):
        return dict(identifiable=False)
    p = _prep(cells)
    S, ok = _newton1(lambda s: _lb(p, s), 0.05, *SB_BRACKET)
    if not ok or not math.isfinite(S):
        return dict(identifiable=True, ok=False, reason="optimizer did not converge")
    ll, _, h = _lb(p, S)
    a = _alpha(p, S)
    mu = L * p["E"] * np.exp(a - S * p["x"])
    dev = float(2 * np.sum(np.where(p["n"] > 0, p["n"] * np.log(np.where(p["n"] > 0, p["n"], 1) / mu), 0.0)
                           - (p["n"] - mu)))
    df = int(np.unique(p["Om"]).size - 2)
    out = dict(identifiable=True, ok=True, S=float(S), se=float(math.sqrt(-1.0 / h)), deviance=dev, df=df,
               p_lof=float(stats.chi2.sf(dev, df)) if df > 0 else None, ll=float(ll), n_cells=int(p["Om"].size),
               k_fit={f"{o:g}": float(m / (L * e)) for o, m, e in zip(p["Om"], mu, p["E"])})
    if with_ci:
        out["ci"] = _profile_ci(lambda s: _lb(p, s)[0], S, ll, SB_BRACKET)
        out["ci_note"] = ci_note(out["ci"])
    return out


def _trend_parts(by_b):
    """by_b: {b: cells}. Uses only identifiable b values."""
    ps = {b: _prep(c) for b, c in by_b.items() if identifiable(c)}
    d = {b: b - B_CENTRE for b in ps}
    return ps, d


def _ll_trend(ps, d, S0, S1):
    ll = g0 = g1 = h00 = h01 = h11 = 0.0
    for b, p in ps.items():
        l, g, h = _lb(p, S0 + S1 * d[b])
        ll += l; g0 += g; g1 += g * d[b]; h00 += h; h01 += h * d[b]; h11 += h * d[b] ** 2
    return ll, np.array([g0, g1]), np.array([[h00, h01], [h01, h11]])


def _max_S0(ps, d, S1, S0=0.05):
    def f(s0):
        ll, g, h = _ll_trend(ps, d, s0, S1)
        return ll, g[0], h[0, 0]
    return _newton1(f, S0, *SB_BRACKET)


def fit_trend(by_b, with_ci=True):
    """Trend model (§4.2). by_b: {b: cells}. Returns dict with S0, S1, ci (profile, bracketed), LR check (3 b only),
    status 'ok' | 'insufficient' | 'numerical'."""
    ps, d = _trend_parts(by_b)
    bs = sorted(ps)
    if len(bs) < 2:
        return dict(status="insufficient", b_used=bs)
    # 2D Newton with step halving on the concave log likelihood
    th = np.array([0.05, 0.0])
    ok = False
    for _ in range(300):
        ll, g, H = _ll_trend(ps, d, *th)
        if not np.all(np.isfinite(g)):
            break
        try:
            step = -np.linalg.solve(H, g)
        except np.linalg.LinAlgError:
            break
        if not np.all(np.isfinite(step)):
            break
        sc = min(1.0, 0.5 / max(1e-300, np.max(np.abs(step))))
        step = step * sc
        for _ in range(60):
            if _ll_trend(ps, d, *(th + step))[0] >= ll - 1e-13:
                break
            step = step * 0.5
        th = th + step
        if np.max(np.abs(step)) < 1e-12:
            ok = True
            break
    S0, S1 = float(th[0]), float(th[1])
    if not ok or not (math.isfinite(S0) and math.isfinite(S1)):
        return dict(status="numerical", b_used=bs, reason="trend optimizer did not converge")
    ll_hat = _ll_trend(ps, d, S0, S1)[0]
    out = dict(status="ok", b_used=bs, S0=S0, S1=S1, ll=ll_hat)
    if with_ci:
        def prof(s1):
            s0, okk = _max_S0(ps, d, s1, S0)
            if not okk:
                raise FloatingPointError("profile inner optimisation failed")
            return _ll_trend(ps, d, s0, s1)[0]
        try:
            out["ci"] = _profile_ci(prof, S1, ll_hat, S1_BRACKET)
            out["ci_note"] = ci_note(out["ci"])
        except FloatingPointError as e:
            return dict(status="numerical", b_used=bs, reason=str(e), S0=S0, S1=S1)
    if len(bs) == 3:
        ll_sep = 0.0
        for b, p in ps.items():
            s, okk = _newton1(lambda s: _lb(p, s), S0 + S1 * d[b], *SB_BRACKET)
            if not okk:
                return dict(status="numerical", b_used=bs, reason="separate-slope fit failed", S0=S0, S1=S1)
            ll_sep += _lb(p, s)[0]
        lr = max(0.0, 2 * (ll_sep - ll_hat))
        out["lr"] = dict(stat=float(lr), df=1, p=float(stats.chi2.sf(lr, 1)))
    return out


# ------------------------------------------------------------------------------------------------ branches and rules
def branch_of(b_used):
    s = set(round(b, 2) for b in b_used)
    if {29.0, 29.2, 29.35} <= s:
        return "B3"
    if s == {29.0, 29.2}:
        return "B2-20"
    if s == {29.0, 29.35}:
        return "B2-35"
    if len(s) >= 2:
        return "B2-new"                            # {29.20, 29.35} without b = 29.00: unsupported (MG034 r3)
    return "B-none"


def direction(ci):
    lo, hi = ci
    if hi < 0:
        return "NEGATIVE"
    if lo > 0:
        return "POSITIVE"
    return None


def primary(trend, branch_calibrated, budget_truncated):
    """§5 primary outcome (with amendment 1). trend: fit_trend output. branch_calibrated: bool for the branch the
    study ends in. budget_truncated: the study-level flag: the budget cutoff truncated the planned rate-data
    collection at any stage (placement, screening, a whole b, an unstarted accepted cell or a production cell)."""
    if trend["status"] == "insufficient":
        return dict(label="NOT RESOLVED (insufficient)", confirmatory="NOT RESOLVED", exploratory=False)
    if trend["status"] == "numerical":
        return dict(label="NOT RESOLVED (NUMERICAL)", confirmatory="NOT RESOLVED", exploratory=False)
    dr = direction(trend["ci"])
    if dr is None:
        return dict(label="NOT RESOLVED", confirmatory="NOT RESOLVED", exploratory=False)
    lab = f"{dr} LINEAR TREND"
    if (not branch_calibrated) or budget_truncated:
        why = [w for w, c in (("branch not calibrated", not branch_calibrated),
                              ("BUDGET-TRUNCATED", budget_truncated)) if c]
        return dict(label=f"{lab} (EXPLORATORY)", confirmatory="NOT RESOLVED", exploratory=True, why=why)
    return dict(label=lab, confirmatory=lab, exploratory=False)


def nonmonotone(fits):
    """Descriptive qualifier: middle S_b outside the range of the other two by more than its CI half-width."""
    bs = sorted(b for b, f in fits.items() if f.get("ok"))
    if len(bs) != 3:
        return None
    lo_, mid, hi_ = (fits[b] for b in bs)
    rng = (min(lo_["S"], hi_["S"]), max(lo_["S"], hi_["S"]))
    half = 0.5 * (mid["ci"][1] - mid["ci"][0])
    out = max(rng[0] - mid["S"], mid["S"] - rng[1])
    return bool(out > half)

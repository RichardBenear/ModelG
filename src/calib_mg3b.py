"""MG3b Part 3: synthetic calibration (registered design: spec/MG3b_followup_proposal.md, Part 3; GATES.md MG3b).

No lattice simulation. Each synthetic run reproduces MG3b's first-passage bookkeeping:
  - startup event (first confirmed event detected at t < 50) with probability p_s -> no exposure;
  - otherwise a first event from hazard h(t) on t >= 50; an event counts if detected at t <= 480, else the run is
    right-censored at 480; exposure = min(t, 480) - 50.
p_s is varied independently of h(t). Designs: 'fixed64', 'fixed96' (Part 1 arms / the Omega = 120 cell) and
'stop20' (batches of 16 until >= 20 events at a batch end, cap 64).

Calibrated (candidates frozen in the proposal):
  - per-cell rate CI: (a) Garwood, (b) exponential profile likelihood, (c) parametric bootstrap percentile (B = 999);
  - Part 1 rate-ratio CI: (a) conditional binomial Clopper-Pearson, (b) conditional binomial mid-p,
    (c) Wald on log ratio with 0.5 added to each count; Part 1 operating characteristics;
  - activated slope CI over the full decision tree: (a) Poisson GLM Wald, (b) GLM profile likelihood,
    (c) parametric bootstrap from the fitted GLM (B = 999) re-running the batching/stopping rule;
  - stationary-hazard check: split (a) t = 200, (b) the cell's exposure median; size and power.

Run:  python src/calib_mg3b.py bench     (timing of 200 replicates per candidate and configuration)
      python src/calib_mg3b.py run       (full calibration; reads the benchmark to choose R for the bootstrap)
Outputs: results/MG3b/calibration/ (bench.json, part3_results.json, decisions.json), results/MG3b/budget_ledger.csv.
Workers: MODELG_WORKERS (default 2).
"""
import os, sys, json, time, platform, csv
import numpy as np
from scipy import stats, special
from concurrent.futures import ProcessPoolExecutor

L = 100.0
TB, TC = 50.0, 480.0
U = TC - TB                       # 430: maximum post-burn-in exposure per run
BATCH, TARGET, CAP = 16, 20, 64
ENTROPY = 20261029
R_MAIN = 4000
R_BOOT_FALLBACK = 1000
B_BOOT = 999
BENCH_R = 200
CEILING_CORE_H = 8.0
OUT = os.path.join(os.path.dirname(__file__), "..", "results", "MG3b", "calibration")
LEDGER = os.path.join(os.path.dirname(__file__), "..", "results", "MG3b", "budget_ledger.csv")

OMEGAS = np.array([130.0, 120.0, 110.0, 100.0])
K120 = 4.6e-5                     # per unit length and time (MG3: 4.633e-5)
PS_CELLS = (0.15, 0.35, 0.33, 0.31)
MARGIN = (0.5, 2.0)

# eligibility thresholds (proposal, Part 3)
COV95, TAIL95 = 0.94, 0.035
COV90, TAIL90 = 0.89, 0.06
EQ_ERR = 0.05
STAT_ALPHA, STAT_SIZE = 0.01, 0.012
BRANCH_MIN = 1000


# ------------------------------------------------------------------------------------------------ hazards
SHAPES = ("const", "weib0.7", "weib1.5", "step2", "step0.5", "early2")


def _hazard_params(shape, lam):
    """Hazards on u = t - 50 in [0, U], scaled so that H(U) = lam * U (same mean hazard as the constant case)."""
    HU = lam * U
    if shape == "const":
        return dict(kind="pw", knots=np.array([0.0, U]), rates=np.array([lam]))
    if shape.startswith("weib"):
        return dict(kind="weib", k=float(shape[4:]), HU=HU)
    if shape.startswith("step"):
        c, u0 = float(shape[4:]), 150.0                     # step at t = 200
        a = HU / (u0 + c * (U - u0))
        return dict(kind="pw", knots=np.array([0.0, u0, U]), rates=np.array([a, c * a]))
    if shape == "early2":
        a = HU / (U + 50.0)                                 # x2 on t in [50, 100)
        return dict(kind="pw", knots=np.array([0.0, 50.0, U]), rates=np.array([2 * a, a]))
    raise ValueError(shape)


def cum_hazard(shape, lam, u):
    p = _hazard_params(shape, lam)
    if p["kind"] == "weib":
        return p["HU"] * (np.asarray(u) / U) ** p["k"]
    knots, rates = p["knots"], p["rates"]
    H = np.zeros_like(np.asarray(u, float))
    for j in range(len(rates)):
        H = H + rates[j] * np.clip(np.asarray(u) - knots[j], 0.0, knots[j + 1] - knots[j])
    return H


def sample_u(rng, shape, lam, size):
    """Post-burn-in first-event times u (np.inf when no event in [0, U]). lam may be an array broadcastable to size."""
    E = rng.exponential(1.0, size)
    lam = np.broadcast_to(np.asarray(lam, float), size)
    if shape == "const":
        u = E / lam
        return np.where(u <= U, u, np.inf)
    out = np.full(size, np.inf)
    if shape.startswith("weib"):
        k, HU = float(shape[4:]), lam * U
        ok = E <= HU
        out[ok] = U * (E[ok] / HU[ok]) ** (1.0 / k)
        return out
    # piecewise-constant: invert segment by segment (lam scalar per call for these shapes)
    p = _hazard_params(shape, float(np.mean(lam)))
    assert np.allclose(lam, lam.flat[0]), "piecewise shapes need a scalar rate"
    knots, rates = p["knots"], p["rates"]
    Hk = np.concatenate([[0.0], np.cumsum(rates * np.diff(knots))])
    for j in range(len(rates)):
        m = (E >= Hk[j]) & (E < Hk[j + 1])
        out[m] = knots[j] + (E[m] - Hk[j]) / rates[j]
    return out


# ------------------------------------------------------------------------------------------------ designs
def sim_cell(rng, design, lam, ps, R, shape="const"):
    """Simulate R replicates of one cell. lam = per-run post-burn-in hazard k*L (scalar or shape (R, 1)).
    Returns dict of arrays: u (R, nmax; inf = censored), st (startup), used (run in the sample), n, E, nruns, nst."""
    nmax = {"fixed64": 64, "fixed96": 96, "stop20": CAP}[design]
    lam_a = np.broadcast_to(np.asarray(lam, float).reshape(-1, 1) if np.ndim(lam) else np.asarray(lam, float),
                            (R, nmax))
    st = rng.random((R, nmax)) < np.broadcast_to(np.asarray(ps, float).reshape(-1, 1) if np.ndim(ps) else ps,
                                                 (R, nmax))
    u = sample_u(rng, shape, lam_a, (R, nmax))
    ev = (~st) & np.isfinite(u)
    ex = np.where(st, 0.0, np.minimum(u, U))
    if design == "stop20":
        nb = nmax // BATCH
        evb = ev.reshape(R, nb, BATCH).sum(2).cumsum(1)
        hit = evb >= TARGET
        j = np.where(hit.any(1), np.argmax(hit, 1), nb - 1)
        used = np.arange(nmax)[None, :] < ((j + 1) * BATCH)[:, None]
    else:
        used = np.ones((R, nmax), bool)
    n = (ev & used).sum(1)
    E = (ex * used).sum(1)
    return dict(u=u, st=st, ev=ev, ex=ex, used=used, n=n, E=E, nruns=used.sum(1), nst=(st & used).sum(1))


# ------------------------------------------------------------------------------------------------ per-cell rate CIs
def ci_garwood(n, E, level=0.95):
    a = 1 - level
    lo = np.where(n > 0, stats.chi2.ppf(a / 2, 2 * np.maximum(n, 1)) / 2, 0.0)
    hi = stats.chi2.ppf(1 - a / 2, 2 * n + 2) / 2
    return lo / (L * E), hi / (L * E)


def ci_profile_rate(n, E, level=0.95):
    """Exponential profile likelihood: {k : 2 [l(khat) - l(k)] <= chi2_1}. Closed form with Lambert W."""
    c = stats.chi2.ppf(level, 1) / 2
    n = np.asarray(n, float)
    hi = np.empty_like(n); lo = np.zeros_like(n)
    z = n == 0
    hi[z] = c                                                    # x = k L E solves x = c
    nz = ~z
    arg = -np.exp(-1.0 - c / n[nz])
    lo[nz] = -n[nz] * np.real(special.lambertw(arg, 0))
    hi[nz] = -n[nz] * np.real(special.lambertw(arg, -1))
    return lo / (L * E), hi / (L * E)


def ci_boot_rate(rng, design, n, E, nruns, nst, level=0.95, B=B_BOOT, chunk=40):
    """Parametric bootstrap percentile: k* from re-simulating the design at (khat, p_s hat), constant hazard."""
    R = len(n)
    lo = np.zeros(R); hi = np.zeros(R)
    kh = np.where(E > 0, n / (L * np.maximum(E, 1e-300)), 0.0)
    psh = nst / nruns
    a = 1 - level
    for i0 in range(0, R, chunk):
        i1 = min(R, i0 + chunk)
        m = i1 - i0
        lam = np.repeat(kh[i0:i1] * L, B)
        ps = np.repeat(psh[i0:i1], B)
        c = sim_cell(rng, design, lam, ps, m * B)
        ks = np.where(c["E"] > 0, c["n"] / (L * np.maximum(c["E"], 1e-300)), 0.0).reshape(m, B)
        lo[i0:i1] = np.quantile(ks, a / 2, axis=1)
        hi[i0:i1] = np.quantile(ks, 1 - a / 2, axis=1)
    return lo, hi


# ------------------------------------------------------------------------------------------------ rate ratio (Part 1)
def _ratio_from_p(p, E1, E2):
    """r = k2/k1 from p = P(event in arm 2 | total) = k2 E2 / (k1 E1 + k2 E2)."""
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.where(p >= 1, np.inf, p / (1 - p) * (E1 / E2))


def ci_ratio_cp(n1, E1, n2, E2, level):
    a = 1 - level
    n = n1 + n2
    plo = np.where(n2 > 0, stats.beta.ppf(a / 2, np.maximum(n2, 1), n - n2 + 1), 0.0)
    phi = np.where(n2 < n, stats.beta.ppf(1 - a / 2, n2 + 1, np.maximum(n - n2, 1)), 1.0)
    return _ratio_from_p(plo, E1, E2), _ratio_from_p(phi, E1, E2)


def _bisect(f, lo, hi, it=60):
    for _ in range(it):
        mid = 0.5 * (lo + hi)
        g = f(mid) > 0
        lo = np.where(g, mid, lo); hi = np.where(g, hi, mid)
    return 0.5 * (lo + hi)


def ci_ratio_midp(n1, E1, n2, E2, level):
    a = 1 - level
    n = (n1 + n2).astype(float); x = n2.astype(float)
    # lower: P(X > x) + 0.5 P(X = x) = a/2 (increasing in p); upper: P(X < x) + 0.5 P(X = x) = a/2 (decreasing in p)
    f_lo = lambda p: stats.binom.sf(x, n, p) + 0.5 * stats.binom.pmf(x, n, p) - a / 2
    f_hi = lambda p: a / 2 - (stats.binom.cdf(x - 1, n, p) + 0.5 * stats.binom.pmf(x, n, p))
    z = np.zeros_like(n); o = np.ones_like(n)
    plo = np.where(x > 0, _bisect(lambda p: -f_lo(p), z, o), 0.0)
    phi = np.where(x < n, _bisect(lambda p: -f_hi(p), z, o), 1.0)
    return _ratio_from_p(plo, E1, E2), _ratio_from_p(phi, E1, E2)


def ci_ratio_wald(n1, E1, n2, E2, level):
    zq = stats.norm.ppf(1 - (1 - level) / 2)
    lr = np.log((n2 + 0.5) / E2) - np.log((n1 + 0.5) / E1)
    se = np.sqrt(1 / (n1 + 0.5) + 1 / (n2 + 0.5))
    return np.exp(lr - zq * se), np.exp(lr + zq * se)


RATIO_METHODS = {"a": ci_ratio_cp, "b": ci_ratio_midp, "c": ci_ratio_wald}


def part1_decision(method, n1, E1, n2, E2):
    """0 = EQUIVALENT, 1 = DIFFERENT, 2 = INCONCLUSIVE (method None = no eligible method -> always INCONCLUSIVE)."""
    if method is None:
        return np.full(len(n1), 2)
    f = RATIO_METHODS[method]
    lo90, hi90 = f(n1, E1, n2, E2, 0.90)
    lo95, hi95 = f(n1, E1, n2, E2, 0.95)
    eq = (lo90 >= MARGIN[0]) & (hi90 <= MARGIN[1])
    diff = ~eq & ((lo95 > 1) | (hi95 < 1))
    return np.where(eq, 0, np.where(diff, 1, 2))


# ------------------------------------------------------------------------------------------------ slope (Poisson GLM)
XC = OMEGAS - OMEGAS.mean()


def glm_fit(n, E):
    """log mu_i = a + s x_i + log(L E_i), x = Omega - mean. Vectorised Newton. Returns s, se(s), deviance."""
    n = n.astype(float); LE = L * E
    s = np.zeros(n.shape[0])
    for _ in range(50):
        w = LE * np.exp(s[:, None] * XC)
        a = np.log(n.sum(1)) - np.log(w.sum(1))
        mu = w * np.exp(a)[:, None]
        g = ((n - mu) * XC).sum(1)
        h00 = mu.sum(1); h01 = (mu * XC).sum(1); h11 = (mu * XC * XC).sum(1)
        info = h11 - h01 ** 2 / h00                              # profile information for s
        step = np.clip(g / info, -0.5, 0.5)
        # step halving: the profile log likelihood is concave in s, so a halved step that increases it exists
        l0 = _lprof(s, n, E)
        for _ in range(40):
            worse = ~(_lprof(s + step, n, E) >= l0 - 1e-12)
            if not worse.any():
                break
            step = np.where(worse, 0.5 * step, step)
        s = s + step
        if np.max(np.abs(step)) < 1e-12:
            break
    w = LE * np.exp(s[:, None] * XC)
    a = np.log(n.sum(1)) - np.log(w.sum(1))
    mu = w * np.exp(a)[:, None]
    h00 = mu.sum(1); h01 = (mu * XC).sum(1); h11 = (mu * XC * XC).sum(1)
    se = np.sqrt(1.0 / (h11 - h01 ** 2 / h00))
    with np.errstate(divide="ignore", invalid="ignore"):
        dev = 2 * np.where(n > 0, n * np.log(np.where(n > 0, n, 1) / np.where(mu > 0, mu, 1)), 0.0) - 2 * (n - mu)
    return s, se, dev.sum(1), a


def _lprof(s, n, E):
    """Profile log likelihood of s (up to a constant)."""
    N = n.sum(1)
    return -N * np.log((L * E * np.exp(s[:, None] * XC)).sum(1)) + s * (n * XC).sum(1)


def ci_slope_profile(n, E, s_hat, se, level=0.95):
    c = stats.chi2.ppf(level, 1) / 2
    n = n.astype(float)
    l0 = _lprof(s_hat, n, E)
    f = lambda s: (l0 - _lprof(s, n, E)) - c                    # < 0 inside the interval
    w = 20 * se
    hi = _bisect(lambda s: -f(s), s_hat, s_hat + w)             # f increases away from s_hat
    lo = _bisect(lambda s: f(s), s_hat - w, s_hat)
    return lo, hi


def ci_slope_wald(s, se, level=0.95):
    z = stats.norm.ppf(1 - (1 - level) / 2)
    return s - z * se, s + z * se


def ci_slope_boot(rng, n, E, nst, nruns, n120, level=0.95, B=B_BOOT, chunk=10):
    """Parametric bootstrap from the fitted GLM: cells 130/110/100 re-run with the stop20 rule, the 120 cell as a
    fixed sample of its size; per-cell p_s hat for startup loss. Percentile CI of s*."""
    s, se, dev, a = glm_fit(n, E)
    R = n.shape[0]
    lo = np.zeros(R); hi = np.zeros(R)
    al = 1 - level
    psh = nst / nruns
    for i0 in range(0, R, chunk):
        i1 = min(R, i0 + chunk); m = i1 - i0
        kcell = np.exp(a[i0:i1, None] + s[i0:i1, None] * XC)            # per-cell rate per unit length
        nn = np.zeros((m * B, 4)); EE = np.zeros((m * B, 4))
        for j in range(4):
            lam = np.repeat(kcell[:, j] * L, B)
            ps = np.repeat(psh[i0:i1, j], B)
            if j == 1:
                for nr in np.unique(n120[i0:i1]):
                    sel = np.repeat(n120[i0:i1] == nr, B)
                    c = sim_cell(rng, "fixed64" if nr == 64 else "fixed96", lam[sel], ps[sel], int(sel.sum()))
                    nn[sel, j] = c["n"]; EE[sel, j] = c["E"]
            else:
                c = sim_cell(rng, "stop20", lam, ps, m * B)
                nn[:, j] = c["n"]; EE[:, j] = c["E"]
        ss, _, _, _ = glm_fit(nn, EE)
        ss = ss.reshape(m, B)
        lo[i0:i1] = np.quantile(ss, al / 2, axis=1); hi[i0:i1] = np.quantile(ss, 1 - al / 2, axis=1)
    return lo, hi


# ------------------------------------------------------------------------------------------------ stationary check
def binom_two_sided_p(k, n, p):
    """Exact two-sided binomial p-value (scipy's 'minlike' rule), vectorised over replicates."""
    k = np.asarray(k); n = np.asarray(n); p = np.asarray(p, float)
    nmax = int(n.max()) if n.size else 0
    j = np.arange(nmax + 1)[None, :]
    pm = stats.binom.pmf(j, n[:, None], p[:, None])
    pm = np.where(j <= n[:, None], pm, 0.0)
    d = stats.binom.pmf(k, n, p)
    out = (pm * (pm <= d[:, None] * (1 + 1e-7))).sum(1)
    return np.where(n == 0, 1.0, np.minimum(out, 1.0))


def stationary_p(c, split):
    """split: 'fixed' (t = 200) or 'median' (the cell's exposure median). Conditional binomial on events in
    [50, split) vs [split, 480] given window exposures."""
    used = c["used"] & ~c["st"]
    ex = np.where(used, c["ex"], 0.0)
    if split == "fixed":
        s = np.full(ex.shape[0], 150.0)
    else:
        s = _exposure_median(ex, used)
    E1 = np.minimum(ex, s[:, None]).sum(1)
    E2 = np.maximum(ex - s[:, None], 0.0).sum(1)
    ev = c["ev"] & used
    n1 = (ev & (c["u"] < s[:, None])).sum(1)
    n2 = (ev & (c["u"] >= s[:, None])).sum(1)
    with np.errstate(invalid="ignore", divide="ignore"):
        p0 = np.where(E1 + E2 > 0, E1 / (E1 + E2), 0.5)
    return binom_two_sided_p(n1, n1 + n2, p0)


def _exposure_median(ex, used):
    """Per replicate: s such that sum_i min(ex_i, s) = 0.5 * sum_i ex_i (piecewise linear in s)."""
    R, n = ex.shape
    xs = np.sort(np.where(used, ex, 0.0), axis=1)
    tot = xs.sum(1)
    cum = np.cumsum(xs, axis=1)
    # f(s) at s = xs[:, j]: cum[:, j] + xs[:, j] * (n - j - 1)
    f = cum + xs * (n - np.arange(n) - 1)[None, :]
    half = 0.5 * tot
    j = np.argmax(f >= half[:, None], axis=1)
    prev_x = np.where(j > 0, xs[np.arange(R), np.maximum(j - 1, 0)], 0.0)
    prev_f = np.where(j > 0, f[np.arange(R), np.maximum(j - 1, 0)], 0.0)
    slope = (n - j).astype(float)
    return np.where(tot > 0, prev_x + (half - prev_f) / slope, 150.0)


# ------------------------------------------------------------------------------------------------ configurations
def cell_configs():
    out = []
    for design in ("fixed64", "fixed96", "stop20"):
        for ps in (0.0, 0.2, 0.4, 0.6):
            for kl in (5e-4, 1.5e-3, 5e-3, 1.5e-2):
                for shape in SHAPES:
                    out.append(dict(kind="cell", design=design, ps=ps, kl=kl, shape=shape))
        out.append(dict(kind="cell", design=design, ps=0.35, kl=K120 * L, shape="const", central=True))
    return out


def part1_configs():
    out = []
    for nr in (64, 96):
        for pair in ((0.35, 0.05), (0.35, 0.35)):
            for r in (0.5, 0.71, 0.8, 1.0, 1.25, 1.41, 2.0):
                out.append(dict(kind="part1", nr=nr, pair=pair, r=r))
    return out


def tree_configs():
    base = [dict(S=S, pref=1.0, ps=PS_CELLS, label=f"S{S}") for S in (0.05, 0.073, 0.10)]
    capped = [dict(S=0.073, pref=1 / 3, ps=PS_CELLS, label="cap_pref/3"),
              dict(S=0.073, pref=1.0, ps=(0.6,) * 4, label="cap_ps0.6"),
              dict(S=0.073, pref=1 / 3, ps=(0.6,) * 4, label="cap_both")]
    trees = [(1.0, 1.0), (0.71, 1.0), (0.5, 1.0), (0.5, 0.71)]
    out = []
    for g in base + capped:
        for r, r2 in trees:
            for nr in (64, 96):
                out.append(dict(kind="tree", nr=nr, r=r, r2=r2, **g))
    out.append(dict(kind="tree", nr=64, r=1.0, r2=1.0, S=0.073, pref=1.0, ps=PS_CELLS, label="lof_130/3",
                    lof=True))
    return out


def rates_for(cfg):
    k = cfg["pref"] * K120 * np.exp(-cfg["S"] * (OMEGAS - 120.0))
    if cfg.get("lof"):
        k = k.copy(); k[0] /= 3.0
    return k


# ------------------------------------------------------------------------------------------------ jobs
def _cov(lo, hi, true):
    lo = np.asarray(lo); hi = np.asarray(hi)
    return dict(cov=float(np.mean((lo <= true) & (true <= hi))), miss_lo=float(np.mean(true < lo)),
                miss_hi=float(np.mean(true > hi)))


def job_cell(cfg, seed, R, boot):
    rng = np.random.default_rng(seed)
    t0 = time.process_time()
    c = sim_cell(rng, cfg["design"], cfg["kl"], cfg["ps"], R, cfg["shape"])
    n, E = c["n"], c["E"]
    out = dict(cfg=cfg, R=R, mean_n=float(n.mean()), mean_runs=float(c["nruns"].mean()),
               p_capped=float(np.mean((c["nruns"] == CAP) & (n < TARGET))) if cfg["design"] == "stop20" else None)
    k_true = cfg["kl"] / L
    kh = np.where(E > 0, n / (L * np.maximum(E, 1e-300)), 0.0)
    if cfg["shape"] == "const":
        res = {}
        pos = n > 0                                     # log-width is defined where the lower limit is > 0
        for m, (lo, hi) in (("a", ci_garwood(n, E)), ("b", ci_profile_rate(n, E))):
            res[m] = _cov(lo, hi, k_true)
            res[m]["mlogw"] = float(np.mean(np.log(hi[pos] / lo[pos]))) if np.any(pos) else None
        if boot is not None:
            Rb = min(R, boot)
            lo, hi = ci_boot_rate(rng, cfg["design"], n[:Rb], E[:Rb], c["nruns"][:Rb], c["nst"][:Rb])
            res["c"] = _cov(lo, hi, k_true); res["c"]["R"] = Rb
            ok = lo > 0
            res["c"]["mlogw"] = float(np.mean(np.log(hi[ok] / lo[ok]))) if np.any(ok) else None
        out["rate_ci"] = res
        out["mean_khat_over_k"] = float(kh.mean() / k_true)
    else:
        # bias of khat relative to the exposure-weighted mean hazard (expected events per unit exposure)
        used = c["used"] & ~c["st"]
        Hsum = (np.where(used, cum_hazard(cfg["shape"], cfg["kl"], np.where(used, c["ex"], 0.0)), 0.0)).sum(1)
        target = Hsum.sum() / (L * E.sum())
        out["khat_bias_vs_exposure_weighted_hazard"] = float(kh.mean() / target - 1)
    out["stat"] = {sp: float(np.mean(stationary_p(c, sp) < STAT_ALPHA)) for sp in ("fixed", "median")}
    out["cpu_s"] = time.process_time() - t0
    return out


def job_part1(cfg, seed, R):
    rng = np.random.default_rng(seed)
    t0 = time.process_time()
    d = "fixed64" if cfg["nr"] == 64 else "fixed96"
    k1 = K120 * L
    a1 = sim_cell(rng, d, k1, cfg["pair"][0], R)
    a2 = sim_cell(rng, d, cfg["r"] * k1, cfg["pair"][1], R)
    n1, E1, n2, E2 = a1["n"], a1["E"], a2["n"], a2["E"]
    res = {}
    for m, f in RATIO_METHODS.items():
        lo95, hi95 = f(n1, E1, n2, E2, 0.95)
        lo90, hi90 = f(n1, E1, n2, E2, 0.90)
        dec = part1_decision(m, n1, E1, n2, E2)
        res[m] = dict(ci95=_cov(lo95, hi95, cfg["r"]), ci90=_cov(lo90, hi90, cfg["r"]),
                      p_equiv=float(np.mean(dec == 0)), p_diff=float(np.mean(dec == 1)),
                      p_inconc=float(np.mean(dec == 2)),
                      mlogw95=float(np.mean(np.log(hi95 / lo95)[np.isfinite(hi95) & (lo95 > 0)])))
    return dict(cfg=cfg, R=R, mean_n1=float(n1.mean()), mean_n2=float(n2.mean()), res=res,
                cpu_s=time.process_time() - t0)


def job_tree(cfg, seed, R, ratio_method, boot):
    """Full decision tree: Part 1 (+1b) -> production tau -> Part 2 cells -> slope CIs."""
    rng = np.random.default_rng(seed)
    t0 = time.process_time()
    k = rates_for(cfg)
    d = "fixed64" if cfg["nr"] == 64 else "fixed96"
    p120 = cfg["ps"][1]
    a1 = sim_cell(rng, d, k[1] * L, p120, R)                         # tau arm
    a2 = sim_cell(rng, d, cfg["r"] * k[1] * L, p120, R)              # tau/2 arm
    dec = part1_decision(ratio_method, a1["n"], a1["E"], a2["n"], a2["E"])
    # tau/4 arm (used only where DIFFERENT). Part 1b always has 64 runs [corrected in addendum 1; the main run used d]
    a4 = sim_cell(rng, "fixed64", cfg["r"] * cfg["r2"] * k[1] * L, p120, R)
    dec1b = part1_decision(ratio_method, a2["n"], a2["E"], a4["n"], a4["E"])
    branch = np.where(dec == 0, "tau", np.where(dec == 2, "tau2_inc", np.where(dec1b == 0, "tau2_1b", "stop")))
    fac = np.where(branch == "tau", 1.0, cfg["r"])                   # production-tau rate factor (constant in Omega)
    n = np.zeros((R, 4), int); E = np.zeros((R, 4)); nst = np.zeros((R, 4)); nruns = np.zeros((R, 4))
    for j in (0, 2, 3):
        c = sim_cell(rng, "stop20", (fac * k[j] * L)[:, None], cfg["ps"][j], R)
        n[:, j], E[:, j], nst[:, j], nruns[:, j] = c["n"], c["E"], c["nst"], c["nruns"]
        if j == 0:
            capped130 = (c["nruns"] == CAP) & (c["n"] < TARGET)
    sel120 = branch == "tau"
    for key, arr in (("n", n), ("E", E), ("nst", nst), ("nruns", nruns)):
        arr[:, 1] = np.where(sel120, a1[key], a2[key])
    live = branch != "stop"
    s, se, dev, _ = glm_fit(n[live], E[live])
    true_s = -cfg["S"]
    br = branch[live]
    cis = {"a": ci_slope_wald(s, se), "b": ci_slope_profile(n[live], E[live], s, se)}
    if boot is not None:
        Rb = min(int(live.sum()), boot)
        lo, hi = ci_slope_boot(rng, n[live][:Rb], E[live][:Rb], nst[live][:Rb], nruns[live][:Rb],
                               np.full(Rb, cfg["nr"]))
        cis["c"] = (lo, hi)
    res = {}
    for m, (lo, hi) in cis.items():
        m_br = br[:len(lo)]
        per = {}
        for b in ("tau", "tau2_inc", "tau2_1b"):
            sel = m_br == b
            if sel.sum() == 0:
                continue
            cv = _cov(lo[sel], hi[sel], true_s) if not cfg.get("lof") else {}
            cv["n_rep"] = int(sel.sum())
            per[b] = cv
        allc = _cov(lo, hi, true_s) if not cfg.get("lof") else {}
        res[m] = dict(by_branch=per, all=allc, mlogw=float(np.mean(np.log(hi - lo))), n_rep=len(lo))
    out = dict(cfg=cfg, R=R, ratio_method=ratio_method,
               branch_counts={b: int(np.sum(branch == b)) for b in ("tau", "tau2_inc", "tau2_1b", "stop")},
               p_capped130=float(capped130.mean()), mean_events=n.mean(0).tolist(), slope=res,
               dev_reject_05=float(np.mean(dev > stats.chi2.ppf(0.95, 2))), mean_S_hat=float(-s.mean()),
               cpu_s=time.process_time() - t0)
    return out


# ------------------------------------------------------------------------------------------------ driver
def _seeds(phase, n):
    root = np.random.SeedSequence(ENTROPY)
    ph = root.spawn(2)[{"bench": 0, "run": 1}[phase]]
    return ph.spawn(n)


def _platform():
    import numpy, scipy
    return dict(machine=platform.node(), os=platform.platform(), python=platform.python_version(),
                numpy=numpy.__version__, scipy=scipy.__version__, numba="n/a (not used)",
                workers=int(os.environ.get("MODELG_WORKERS", "2")))


def _ledger(item, kind, wall_s, extra=""):
    os.makedirs(os.path.dirname(LEDGER), exist_ok=True)
    new = not os.path.exists(LEDGER)
    p = _platform()
    with open(LEDGER, "a", newline="") as f:
        w = csv.writer(f)
        if new:
            w.writerow(["item", "key", "kind", "wall_s", "workers", "machine", "os", "python", "numpy", "scipy",
                        "numba", "utc_time"])
        w.writerow([item, extra, kind, f"{wall_s:.1f}", p["workers"], p["machine"], p["os"], p["python"],
                    p["numpy"], p["scipy"], p["numba"], time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())])


def _run_jobs(jobs):
    nw = int(os.environ.get("MODELG_WORKERS", "2"))
    with ProcessPoolExecutor(nw) as ex:
        futs = [ex.submit(*j) for j in jobs]
        return [f.result() for f in futs]


def all_jobs(phase, R, boot_cell, boot_tree, ratio_method):
    cc, pc, tc = cell_configs(), part1_configs(), tree_configs()
    seeds = _seeds(phase, len(cc) + len(pc) + len(tc))
    jobs = []
    for i, c in enumerate(cc):
        jobs.append((job_cell, c, seeds[i], R, boot_cell if c["shape"] == "const" else None))
    for i, c in enumerate(pc):
        jobs.append((job_part1, c, seeds[len(cc) + i], R))
    for i, c in enumerate(tc):
        jobs.append((job_tree, c, seeds[len(cc) + len(pc) + i], R, ratio_method, boot_tree))
    return jobs


def bench():
    os.makedirs(OUT, exist_ok=True)
    w0 = time.time()
    jobs = all_jobs("bench", BENCH_R, BENCH_R, BENCH_R, "a")
    res = _run_jobs(jobs)
    wall = time.time() - w0
    cpu = sum(r["cpu_s"] for r in res)
    # projection: analytic parts scale with R_MAIN / BENCH_R; bootstrap parts dominate and scale the same way
    proj_core_h = cpu * (R_MAIN / BENCH_R) / 3600
    out = dict(bench_R=BENCH_R, wall_s=wall, cpu_s=cpu, projected_core_h_at_R4000=proj_core_h,
               ceiling_core_h=CEILING_CORE_H, platform=_platform(),
               by_kind={k: sum(r["cpu_s"] for r in res if r["cfg"]["kind"] == k) for k in ("cell", "part1", "tree")})
    json.dump(out, open(os.path.join(OUT, "bench.json"), "w"), indent=1)
    _ledger("Part 3 calibration", "calibration", cpu, "benchmark (cpu-s summed over workers)")
    print(json.dumps(out, indent=1))


def _eligible_rate(cells):
    out = {}
    for m in ("a", "b", "c"):
        rs = [c["rate_ci"][m] for c in cells if c["cfg"]["shape"] == "const" and m in c["rate_ci"]]
        if not rs:
            continue
        mincov = min(r["cov"] for r in rs)
        ok = all(r["cov"] >= COV95 and r["miss_lo"] <= TAIL95 and r["miss_hi"] <= TAIL95 for r in rs)
        cen = [c["rate_ci"][m]["mlogw"] for c in cells if c["cfg"].get("central") and c["cfg"]["design"] == "fixed64"]
        out[m] = dict(eligible=ok, min_cov=mincov, n_configs=len(rs), central_mlogw=cen[0] if cen else None)
    return out


def _select(elig, key="central_mlogw"):
    """Smallest mean log-width at the central configuration among eligible; widths within 1% tie -> earlier letter."""
    ok = [m for m in sorted(elig) if elig[m]["eligible"]]
    if not ok:
        best = max(sorted(elig), key=lambda m: elig[m]["min_cov"])
        return dict(method=best, approximate=True, reason="no candidate eligible: highest minimum coverage")
    wmin = min(elig[m][key] for m in ok)
    for m in ok:                                             # alphabetical = order listed
        if elig[m][key] <= wmin + 0.01 * abs(wmin):
            return dict(method=m, approximate=False, reason="eligible; smallest central width (1% tie rule)")


def _eligible_ratio(p1):
    out = {}
    for m in ("a", "b", "c"):
        rows = [p["res"][m] for p in p1]
        ok95 = all(r["ci95"]["cov"] >= COV95 and r["ci95"]["miss_lo"] <= TAIL95 and r["ci95"]["miss_hi"] <= TAIL95
                   for r in rows)
        ok90 = all(r["ci90"]["cov"] >= COV90 and r["ci90"]["miss_lo"] <= TAIL90 and r["ci90"]["miss_hi"] <= TAIL90
                   for r in rows)
        okeq = all(p["res"][m]["p_equiv"] <= EQ_ERR for p in p1 if p["cfg"]["r"] in (0.5, 2.0))
        cen = [p["res"][m]["mlogw95"] for p in p1 if p["cfg"]["nr"] == 64 and p["cfg"]["pair"] == (0.35, 0.05)
               and p["cfg"]["r"] == 1.0]
        out[m] = dict(eligible=ok95 and ok90 and okeq, ok95=ok95, ok90=ok90, ok_equiv_error=okeq,
                      min_cov=min(r["ci95"]["cov"] for r in rows),
                      central_mlogw=cen[0])
    return out


def _eligible_slope(trees, R_boot=R_MAIN):
    out = {}
    for m in ("a", "b", "c"):
        rows = []
        for t in trees:
            if t["cfg"].get("lof") or m not in t["slope"]:
                continue
            # a branch decides eligibility if it has >= 1000 replicates; for the bootstrap at the fallback R = 1000 the
            # threshold is scaled pro rata (250). [Corrected 2026-10-29 after the main run: the first version scaled by
            # the bootstrap's live-replicate count, giving ~995 instead of 1000 at R_boot = 4000. See `recheck`.]
            thr = BRANCH_MIN * (R_boot / R_MAIN if m == "c" else 1.0)
            for b, cv in t["slope"][m]["by_branch"].items():
                if cv["n_rep"] >= thr:
                    rows.append((t["cfg"]["label"], t["cfg"]["nr"], t["cfg"]["r"], t["cfg"]["r2"], b, cv))
        if not rows:
            continue
        ok = all(cv["cov"] >= COV95 and cv["miss_lo"] <= TAIL95 and cv["miss_hi"] <= TAIL95 for *_, cv in rows)
        cen = [t["slope"][m]["mlogw"] for t in trees if t["cfg"]["label"] == "S0.073" and t["cfg"]["nr"] == 64
               and t["cfg"]["r"] == 1.0 and not t["cfg"].get("lof")]
        out[m] = dict(eligible=ok, min_cov=min(cv["cov"] for *_, cv in rows), n_branch_cells=len(rows),
                      central_mlogw=cen[0] if cen else None,
                      failing=[dict(label=a, nr=b_, r=c, r2=d, branch=e, **cv) for a, b_, c, d, e, cv in rows
                               if not (cv["cov"] >= COV95 and cv["miss_lo"] <= TAIL95 and cv["miss_hi"] <= TAIL95)])
    return out


def _stat_choice(cells):
    out = {}
    for sp in ("fixed", "median"):
        size = [c["stat"][sp] for c in cells if c["cfg"]["shape"] == "const"]
        powr = [c["stat"][sp] for c in cells if c["cfg"]["shape"] in ("step2", "step0.5")]
        out[sp] = dict(max_size=max(size), eligible=max(size) <= STAT_SIZE, mean_power_steps=float(np.mean(powr)))
    ok = [sp for sp in ("fixed", "median") if out[sp]["eligible"]]
    if not ok:
        choice = None
    elif len(ok) == 1:
        choice = ok[0]
    else:
        pf, pm = out["fixed"]["mean_power_steps"], out["median"]["mean_power_steps"]
        choice = "median" if pm > pf + 0.01 else "fixed"     # tie (within 0.01) -> (a) fixed
    return out, choice


def run():
    os.makedirs(OUT, exist_ok=True)
    b = json.load(open(os.path.join(OUT, "bench.json")))
    R_boot = R_MAIN if b["projected_core_h_at_R4000"] <= CEILING_CORE_H - b["cpu_s"] / 3600 else R_BOOT_FALLBACK
    w0 = time.time()
    # phase A: per-cell and Part 1 (independent of the ratio method)
    cc, pc, tc = cell_configs(), part1_configs(), tree_configs()
    seeds = _seeds("run", len(cc) + len(pc) + len(tc))
    jobsA = [(job_cell, c, seeds[i], R_MAIN, R_boot if c["shape"] == "const" else None) for i, c in enumerate(cc)]
    jobsA += [(job_part1, c, seeds[len(cc) + i], R_MAIN) for i, c in enumerate(pc)]
    resA = _run_jobs(jobsA)
    cells = [r for r in resA if r["cfg"]["kind"] == "cell"]
    p1 = [r for r in resA if r["cfg"]["kind"] == "part1"]
    el_rate = _eligible_rate(cells); sel_rate = _select(el_rate)
    el_ratio = _eligible_ratio(p1)
    ok_ratio = [m for m in sorted(el_ratio) if el_ratio[m]["eligible"]]
    ratio_method = _select(el_ratio)["method"] if ok_ratio else None
    if ratio_method:
        pe = [p["res"][ratio_method]["p_equiv"] for p in p1 if p["cfg"]["nr"] == 64
              and p["cfg"]["pair"] == (0.35, 0.05) and p["cfg"]["r"] == 1.0][0]
        n_arm = 64 if pe >= 0.80 else 96
    else:
        pe, n_arm = None, 64
    stat_rows, stat_split = _stat_choice(cells)
    # phase B: decision tree with the selected ratio method
    jobsB = [(job_tree, c, seeds[len(cc) + len(pc) + i], R_MAIN, ratio_method, R_boot) for i, c in enumerate(tc)]
    trees = _run_jobs(jobsB)
    el_slope = _eligible_slope(trees, R_boot); sel_slope = _select(el_slope)
    wall = time.time() - w0
    cpu = sum(r["cpu_s"] for r in resA + trees)
    lof = [t for t in trees if t["cfg"].get("lof")][0]
    act = [t for t in trees if not t["cfg"].get("lof")]
    decisions = dict(
        registered_design="spec/MG3b_followup_proposal.md at 202c00f; GATES.md MG3b",
        R=R_MAIN, R_bootstrap=R_boot, B=B_BOOT,
        mc_se_at_0p95=dict(R4000=float(np.sqrt(0.95 * 0.05 / R_MAIN)), R_boot=float(np.sqrt(0.95 * 0.05 / R_boot))),
        n_configurations=dict(cell=len(cc), cell_constant=sum(c["shape"] == "const" for c in cc), part1=len(pc),
                              tree=len(tc)),
        per_cell_rate_ci=sel_rate,
        rate_ratio_ci=dict(method=ratio_method, fallback_part1_always_inconclusive=ratio_method is None),
        part1_runs_per_arm=n_arm, p_equiv_at_r1_64=pe,
        slope_ci=sel_slope, stationary_split=stat_split,
        stationary_fallback_uncalibrated=stat_split is None,
        platform=_platform(), cpu_core_h=cpu / 3600, wall_s=wall)
    results = dict(eligibility=dict(rate=el_rate, ratio=el_ratio, slope=el_slope, stationary=stat_rows),
                   cells=cells, part1=p1, trees=trees,
                   lack_of_fit=dict(dev_power_05=lof["dev_reject_05"], mean_S_hat=lof["mean_S_hat"],
                                    size_05_activated=dict(min=min(t["dev_reject_05"] for t in act),
                                                           max=max(t["dev_reject_05"] for t in act))))
    json.dump(results, open(os.path.join(OUT, "part3_results.json"), "w"), indent=1, default=str)
    json.dump(decisions, open(os.path.join(OUT, "decisions.json"), "w"), indent=1, default=str)
    _ledger("Part 3 calibration", "calibration", cpu, f"main run R={R_MAIN} R_boot={R_boot} (cpu-s summed)")
    print(json.dumps(decisions, indent=1, default=str))


# ------------------------------------------------------------------------------------------------ addendum 1
R_ADD = 20000


def addendum_configs():
    out = []
    for pair in ((0.35, 0.35, 0.35), (0.35, 0.05, 0.05)):
        for r in (0.5, 0.71, 0.8, 1.0, 1.25, 1.41, 2.0):
            for r2 in (0.5, 0.71, 1.0, 1.41, 2.0):
                out.append(dict(kind="addendum1", r=r, r2=r2, pair=pair))
    return out


def job_addendum(cfg, seed, R, method="a"):
    """(A) Garwood coverage of the reused Omega = 120 arm's rate, by production branch; (B) Part 1b triggered and
    whole-study equivalence error. Part 1 and 1b only, 64 runs per arm."""
    rng = np.random.default_rng(seed)
    t0 = time.process_time()
    k1 = K120 * L
    a1 = sim_cell(rng, "fixed64", k1, cfg["pair"][0], R)
    a2 = sim_cell(rng, "fixed64", cfg["r"] * k1, cfg["pair"][1], R)
    a4 = sim_cell(rng, "fixed64", cfg["r"] * cfg["r2"] * k1, cfg["pair"][2], R)
    dec = part1_decision(method, a1["n"], a1["E"], a2["n"], a2["E"])
    dec1b = part1_decision(method, a2["n"], a2["E"], a4["n"], a4["E"])
    branch = np.where(dec == 0, "tau", np.where(dec == 2, "tau2_inc", np.where(dec1b == 0, "tau2_1b", "stop")))
    rate = {}
    for b, arm, true in (("tau", a1, K120), ("tau2_inc", a2, cfg["r"] * K120), ("tau2_1b", a2, cfg["r"] * K120)):
        sel = branch == b
        if sel.sum() == 0:
            continue
        lo, hi = ci_garwood(arm["n"][sel], arm["E"][sel])
        cv = _cov(lo, hi, true); cv["n_rep"] = int(sel.sum())
        rate[b] = cv
    trig = dec == 1
    eq1b = trig & (dec1b == 0)
    out = dict(cfg=cfg, R=R, rate_by_branch=rate,
               part1b=dict(p_triggered=float(trig.mean()), p_equiv_given_triggered=float(eq1b.sum() / max(trig.sum(), 1)),
                           p_triggered_and_equiv=float(eq1b.mean()), n_triggered=int(trig.sum())),
               branch_counts={b: int(np.sum(branch == b)) for b in ("tau", "tau2_inc", "tau2_1b", "stop")},
               cpu_s=time.process_time() - t0)
    return out


def addendum():
    cfgs = addendum_configs()
    root = np.random.SeedSequence(ENTROPY, spawn_key=(2,))
    seeds = root.spawn(len(cfgs))
    w0 = time.time()
    res = _run_jobs([(job_addendum, c, seeds[i], R_ADD) for i, c in enumerate(cfgs)])
    # (A) per branch: qualifying configurations and the observed minimum
    A = {}
    for b in ("tau", "tau2_inc", "tau2_1b"):
        rows = [(r["cfg"], r["rate_by_branch"][b]) for r in res if b in r["rate_by_branch"]
                and r["rate_by_branch"][b]["n_rep"] >= BRANCH_MIN]
        fail = [dict(r=c["r"], r2=c["r2"], pair=c["pair"], **cv) for c, cv in rows
                if not (cv["cov"] >= COV95 and cv["miss_lo"] <= TAIL95 and cv["miss_hi"] <= TAIL95)]
        A[b] = dict(n_qualifying=len(rows), min_cov=min((cv["cov"] for _, cv in rows), default=None),
                    max_miss_lo=max((cv["miss_lo"] for _, cv in rows), default=None),
                    max_miss_hi=max((cv["miss_hi"] for _, cv in rows), default=None),
                    label_approximate=bool(fail), n_failing=len(fail), failing=fail)
    # (B) boundaries r2 in {0.5, 2}
    bnd = [r for r in res if r["cfg"]["r2"] in (0.5, 2.0)]
    B = dict(max_whole_study=max(r["part1b"]["p_triggered_and_equiv"] for r in bnd),
             max_conditional=max(r["part1b"]["p_equiv_given_triggered"] for r in bnd if r["part1b"]["n_triggered"] >= 100),
             rows=[dict(r=r["cfg"]["r"], r2=r["cfg"]["r2"], pair=r["cfg"]["pair"], **r["part1b"]) for r in bnd])
    B["label_approximate"] = B["max_whole_study"] > 0.05
    cpu = sum(r["cpu_s"] for r in res)
    out = dict(R=R_ADD, n_configs=len(cfgs), A_reused_rate=A, B_part1b=B, platform=_platform(), cpu_core_h=cpu / 3600,
               wall_s=time.time() - w0, results=res)
    d = os.path.join(OUT, "addendum1"); os.makedirs(d, exist_ok=True)
    json.dump(out, open(os.path.join(d, "addendum1_results.json"), "w"), indent=1, default=str)
    _ledger("Part 3 calibration", "calibration", cpu, f"addendum 1 R={R_ADD} (cpu-s summed)")
    print(json.dumps({k: v for k, v in out.items() if k != "results"}, indent=1, default=str)[:6000])


def recheck():
    """Recompute the slope eligibility and selection from part3_results.json with the corrected branch threshold
    (no new simulation). Writes eligibility_recheck.json; decisions.json is left as the run wrote it."""
    r = json.load(open(os.path.join(OUT, "part3_results.json")))
    d = json.load(open(os.path.join(OUT, "decisions.json")))
    el = _eligible_slope(r["trees"], d["R_bootstrap"])
    sel = _select(el)
    out = dict(slope_eligibility=el, slope_selection=sel, run_selection=d["slope_ci"],
               selection_unchanged=sel == d["slope_ci"])
    json.dump(out, open(os.path.join(OUT, "eligibility_recheck.json"), "w"), indent=1, default=str)
    print(json.dumps({k: v for k, v in out.items() if k != "slope_eligibility"}, indent=1),
          {m: (v["eligible"], round(v["min_cov"], 4), v["n_branch_cells"]) for m, v in el.items()})


if __name__ == "__main__":
    {"bench": bench, "run": run, "recheck": recheck, "addendum": addendum}[sys.argv[1]]()

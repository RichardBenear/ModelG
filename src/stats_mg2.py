"""MG003b: inference layer for profile-bias statistics (MG2 C1/C2/C6, MG3 C2b, V2, V4).

Point estimate (reported): the U-statistic for the squared relative bias,
    Ehat2 = sum_{i!=j} <d_i, d_j> / (M(M-1) |phi*|^2),  d_i = m_i - phi*  (unbiased; may be < 0).
Decision bounds (used for PASS/FAIL), by the triangle inequality around the run-mean m_bar:
    | |mu-phi*| - |m_bar-phi*| | <= |m_bar - mu|,
and |m_bar - mu| is bounded by k * sqrt(tr(Sigma_hat)/M), with k = t_{0.975, M-1}.
This is exact in the rank-1 Gaussian case and conservative otherwise; coverage is CALIBRATED in
src/calibrate_bounds.py across Gaussian, low-rank, two-point and heavy-tailed run distributions.
Three-way decision against a tolerance eps:
    PASS if upper < eps;  FAIL if lower > eps;  INCONCLUSIVE otherwise.
"""
import numpy as np
from scipy import stats


def u_stat_bias2(runs, ref):
    """runs: (M, d) run-level profiles; ref: (d,). Unbiased estimate of |mu-ref|^2/|ref|^2."""
    D = runs - ref
    M = D.shape[0]
    s = D.sum(0)
    tot = s @ s - np.einsum("ij,ij->", D, D)
    return tot / (M * (M - 1)) / (ref @ ref)


def bias_bounds(runs, ref, level=0.95):
    """Relative (lower, upper) bounds on |mu - ref|/|ref| (two-sided level), plus point values."""
    M = runs.shape[0]
    mbar = runs.mean(0)
    nref = np.sqrt(ref @ ref)
    tot = np.sqrt(((mbar - ref) @ (mbar - ref))) / nref
    trS = np.sum(runs.var(0, ddof=1))
    k = stats.t.ppf(0.5 + level / 2, M - 1)
    half = k * np.sqrt(trS / M) / nref
    e2 = u_stat_bias2(runs, ref)
    return dict(E_total=tot, Eb_point=np.sqrt(max(e2, 0.0)), Eb2_ustat=e2,
                lower=max(tot - half, 0.0), upper=tot + half, halfwidth=half, k=k, M=M)


def decide(bounds, eps):
    if bounds["upper"] < eps:
        return "PASS"
    if bounds["lower"] > eps:
        return "FAIL"
    return "INCONCLUSIVE"


def diff_bounds(runsA, runsB, ref_norm, level=0.95):
    """Upper bound on |muA - muB|/ref_norm (Welch-type k from the smaller sample)."""
    mA, mB = runsA.mean(0), runsB.mean(0)
    d = np.sqrt((mA - mB) @ (mA - mB)) / ref_norm
    v = np.sum(runsA.var(0, ddof=1)) / runsA.shape[0] + np.sum(runsB.var(0, ddof=1)) / runsB.shape[0]
    k = stats.t.ppf(0.5 + level / 2, min(runsA.shape[0], runsB.shape[0]) - 1)
    return dict(diff=d, upper=d + k * np.sqrt(v) / ref_norm, lower=max(d - k * np.sqrt(v) / ref_norm, 0.0))


def alpha_fit(omegas, e2):
    """Least-squares fit of log(Ehat2) = c - 2 alpha log(Omega) over cells with Ehat2 > 0."""
    om, e = np.asarray(omegas, float), np.asarray(e2, float)
    ok = e > 0
    if ok.sum() < 3:
        return np.nan
    A = np.vstack([np.ones(ok.sum()), -2 * np.log(om[ok])]).T
    coef, *_ = np.linalg.lstsq(A, np.log(e[ok]), rcond=None)
    return coef[1]


def alpha_ci(cells, omegas, B=2000, level=0.95, widen=1.0, rng=None):
    """Run-level nonparametric bootstrap CI for alpha: resample runs within each Omega cell,
    recompute the U-statistics, refit. widen: calibrated multiplicative half-width factor."""
    rng = np.random.default_rng(rng)
    ref = cells["ref"]
    a_hat = alpha_fit(omegas, [u_stat_bias2(r, ref) for r in cells["runs"]])
    boots = []
    for _ in range(B):
        e2 = []
        for r in cells["runs"]:
            idx = rng.integers(0, r.shape[0], r.shape[0])
            if len(np.unique(idx)) < 2:
                idx[0] = (idx[0] + 1) % r.shape[0]
            e2.append(u_stat_bias2(r[idx], ref))
        boots.append(alpha_fit(omegas, e2))
    boots = np.array(boots)
    frac_valid = np.mean(np.isfinite(boots))
    boots = boots[np.isfinite(boots)]
    if len(boots) < 0.8 * B:
        return a_hat, (np.nan, np.nan), frac_valid
    lo, hi = np.quantile(boots, [(1 - level) / 2, (1 + level) / 2])
    lo, hi = a_hat - widen * (a_hat - lo), a_hat + widen * (hi - a_hat)
    return a_hat, (lo, hi), frac_valid

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


def u_stat_bias2_boot(runs, idx, ref):
    """U-statistic on a bootstrap resample, using only pairs of DISTINCT original runs (MG003b fix:
    pairing a duplicated run with itself adds |d|^2 terms and biases the bootstrap upward)."""
    D = runs - ref
    c = np.bincount(idx, minlength=runs.shape[0]).astype(float)
    S = (c[:, None] * D).sum(0)
    self_terms = np.sum(c**2 * np.einsum("ij,ij->i", D, D))
    npairs = idx.size**2 - np.sum(c**2)
    if npairs <= 0:
        return np.nan
    return (S @ S - self_terms) / npairs / (ref @ ref)


def select_cells(cells, omegas, level=0.95):
    """Pre-registered cell selection for the alpha fit: Omega >= 100 cells whose bias lower bound > 0
    (significant bias), decided once on the original data."""
    keep = []
    for o, r in zip(omegas, cells["runs"]):
        if o >= 100 and bias_bounds(r, cells["ref"], level)["lower"] > 0:
            keep.append(True)
        else:
            keep.append(False)
    return np.array(keep)


def alpha_fit(omegas, e2):
    """LS fit of log(Ehat2) = c - 2 alpha log(Omega); requires >= 3 cells, all Ehat2 > 0."""
    om, e = np.asarray(omegas, float), np.asarray(e2, float)
    if len(e) < 3 or np.any(~np.isfinite(e)) or np.any(e <= 0):
        return np.nan
    A = np.vstack([np.ones(len(e)), -2 * np.log(om)]).T
    coef, *_ = np.linalg.lstsq(A, np.log(e), rcond=None)
    return coef[1]


def alpha_ci(cells, omegas, B=2000, level=0.95, widen=1.0, rng=None):
    """alpha point estimate and run-level bootstrap percentile CI on the pre-selected cells.
    Returns (alpha_hat, (lo, hi), fraction of valid bootstrap fits, selected-cell mask)."""
    rng = np.random.default_rng(rng)
    ref = cells["ref"]
    om = np.asarray(omegas, float)
    keep = select_cells(cells, om, level)
    if keep.sum() < 3:
        return np.nan, (np.nan, np.nan), 0.0, keep
    runs = [r for r, k in zip(cells["runs"], keep) if k]
    a_hat = alpha_fit(om[keep], [u_stat_bias2(r, ref) for r in runs])
    boots = np.empty(B)
    for b in range(B):
        e2 = []
        for r in runs:
            idx = rng.integers(0, r.shape[0], r.shape[0])
            e2.append(u_stat_bias2_boot(r, idx, ref))
        boots[b] = alpha_fit(om[keep], e2)
    ok = np.isfinite(boots)
    if ok.mean() < 0.8 or not np.isfinite(a_hat):
        return a_hat, (np.nan, np.nan), ok.mean(), keep
    lo, hi = np.quantile(boots[ok], [(1 - level) / 2, (1 + level) / 2])
    return a_hat, (a_hat - widen * (a_hat - lo), a_hat + widen * (hi - a_hat)), ok.mean(), keep

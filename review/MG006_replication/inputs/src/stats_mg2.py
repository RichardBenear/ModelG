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


# ============================================================================ MG003c additions
# Assumption A1 (stated in protocol v1.3 §5.3): run-level aligned profiles of eligible runs are i.i.d.
# with approximately Gaussian fluctuations, i.e. no rare discrete run-level modes beyond the tracked
# events (loss / replication / merge), which are removed by survival conditioning. The confidence
# claims of bias_bounds / diff_bounds are conditional on A1. a1_diagnostics() checks A1 on the actual
# decision data; a cell that fails the check is INCONCLUSIVE.

def welch_diff_bounds(runsA, runsB, ref_norm, level=0.95):
    """Two-sample bounds on |muA - muB| / ref_norm with Welch-Satterthwaite df on the trace variances.
    ref_norm must be a DETERMINISTIC norm (e.g. |phi*| of the lattice-ODE soliton)."""
    MA, MB = runsA.shape[0], runsB.shape[0]
    mA, mB = runsA.mean(0), runsB.mean(0)
    d = np.sqrt((mA - mB) @ (mA - mB)) / ref_norm
    vA = np.sum(runsA.var(0, ddof=1)) / MA
    vB = np.sum(runsB.var(0, ddof=1)) / MB
    df = (vA + vB) ** 2 / (vA**2 / (MA - 1) + vB**2 / (MB - 1)) if (vA + vB) > 0 else min(MA, MB) - 1
    k = stats.t.ppf(0.5 + level / 2, df)
    half = k * np.sqrt(vA + vB) / ref_norm
    return dict(diff=d, upper=d + half, lower=max(d - half, 0.0), df=df)


def _robust_z(x):
    med = np.median(x)
    mad = np.median(np.abs(x - med)) * 1.4826
    if mad > 0:
        return (x - med) / mad
    # MAD = 0: a majority of identical values; any differing value is an infinite outlier
    return np.where(np.abs(x - med) > 0, np.inf, 0.0)


def _studentized_loo(x):
    """Externally studentised leave-one-out residuals; ~ t_{M-2} under i.i.d. Gaussian.
    Zero leave-one-out spread with a differing value -> infinite (flagged)."""
    M = x.size
    out = np.empty(M)
    for i in range(M):
        r = np.delete(x, i)
        s = r.std(ddof=1)
        dev = x[i] - r.mean()
        out[i] = dev / (s * np.sqrt(1 + 1 / (M - 1))) if s > 0 else (np.inf if dev != 0 else 0.0)
    return out


def a1_diagnostics(block_profiles, ref, z_block=5.0, alpha_run=0.01, p_min=1e-3):
    """block_profiles: (M runs, K blocks, d). Checks A1 on the decision data:
      (i)  within-run block deviations (block mean minus run mean), projected on the run-mean bias
           direction and the first two principal components: Shapiro-Wilk p >= p_min and no
           robust |z| > z_block;
      (ii) run-level means projected on the same directions: no externally studentised leave-one-out
           residual beyond the Bonferroni t_{M-2} critical value at overall level alpha_run
           (MG003c: replaces a MAD-based robust z, whose null false-alarm rate at M = 12 was ~21%).
    Returns dict(ok=bool, details=...)."""
    M, K, d = block_profiles.shape
    runs = block_profiles.mean(1)
    mbar = runs.mean(0)
    dirs = []
    b = mbar - ref
    if np.linalg.norm(b) > 0:
        dirs.append(b / np.linalg.norm(b))
    W = (block_profiles - runs[:, None, :]).reshape(M * K, d)
    _, _, Vt = np.linalg.svd(W - W.mean(0), full_matrices=False)
    dirs += [Vt[0], Vt[1]]
    det = {"shapiro_p": [], "max_block_z": [], "max_run_t": []}
    t_crit = stats.t.ppf(1 - alpha_run / (2 * M * len(dirs)), M - 2)
    det["run_t_crit"] = float(t_crit)
    ok = True
    for v in dirs:
        pb = W @ v
        p = stats.shapiro(pb).pvalue if np.std(pb) > 0 else 1.0
        zb = np.max(np.abs(_robust_z(pb)))
        tr = np.max(np.abs(_studentized_loo((runs - mbar) @ v)))
        det["shapiro_p"].append(float(p)); det["max_block_z"].append(float(zb)); det["max_run_t"].append(float(tr))
        ok &= (p >= p_min) and (zb <= z_block) and (tr <= t_crit)
    return dict(ok=bool(ok), details=det)


def decide_cell(get_runs, ref, eps, level_stage=0.95, m_min=12, m_max_total=24, ext=12):
    """Full pre-registered single-cell procedure (MG2-C1/C6, V4):
      1. obtain eligible runs (survival-conditioned) until m_min eligible or m_max_total attempted;
         fewer than m_min eligible -> INCONCLUSIVE ('insufficient eligible runs');
      2. A1 diagnostics fail -> INCONCLUSIVE ('A1 not supported');
      3. three-way decision at level_stage; if INCONCLUSIVE, ONE extension of `ext` eligible runs and
         re-decide on all runs at level_stage (calibrated two-look error rates, protocol §5.3).
    get_runs(n) must return (block_profiles (n_eligible, K, d), n_attempted)."""
    blocks, att = get_runs(m_min)
    if blocks.shape[0] < m_min:
        return "INCONCLUSIVE", dict(reason="insufficient eligible runs", attempted=att)
    if not a1_diagnostics(blocks, ref)["ok"]:
        return "INCONCLUSIVE", dict(reason="A1 not supported")
    dcs = decide(bias_bounds(blocks.mean(1), ref, level_stage), eps)
    if dcs != "INCONCLUSIVE":
        return dcs, dict(stage=1)
    more, att2 = get_runs(ext)
    allb = np.concatenate([blocks, more], 0)
    if more.shape[0] < ext:
        return "INCONCLUSIVE", dict(reason="insufficient eligible runs in extension")
    if not a1_diagnostics(allb, ref)["ok"]:
        return "INCONCLUSIVE", dict(reason="A1 not supported (extension)")
    return decide(bias_bounds(allb.mean(1), ref, level_stage), eps), dict(stage=2)


# ---------------------------------------------------------------- MG005: corrected gate-level audit
def audit_within_cell(cell_projections, alpha=0.01):
    """Corrected gate-level A1 audit (MG005; replaces the registered pooled audit, which standardised each
    cell before pooling and therefore could never exceed its critical value: max |t| <= ~3.43 < 4.65 for
    13 cells of 12 runs).
    cell_projections: list over (cell, field, direction) of 1-D arrays of run-level projections.
    Each value is externally studentised against the other runs of its own cell (~ t_{M-2} under i.i.d.
    Gaussian runs) and Bonferroni-corrected over every run x field x direction test at overall level alpha.
    Returns dict(max_stat, crit, n_tests, flagged=list of (index, run, stat))."""
    from scipy import stats as _st
    n_tests = int(sum(p.size for p in cell_projections))
    out, mx, crit_by_m = [], 0.0, {}
    for i, p in enumerate(cell_projections):
        M = p.size
        if M not in crit_by_m:
            crit_by_m[M] = _st.t.ppf(1 - alpha / (2 * n_tests), M - 2)
        t = np.abs(_studentized_loo(p))
        mx = max(mx, float(t.max()))
        for r in np.where(t > crit_by_m[M])[0]:
            out.append((i, int(r), float(t[r])))
    crit = min(crit_by_m.values()) if crit_by_m else np.nan
    return dict(max_stat=mx, crit=float(crit), n_tests=n_tests, flagged=out)

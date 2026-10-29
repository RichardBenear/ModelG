"""Tests for src/calib_mg3b.py (MG3b Part 3 calibration helpers). Run: python tests/test_calib_mg3b.py"""
import os, sys
import numpy as np
from scipy import stats, optimize
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
import calib_mg3b as C
import mg3


def test_garwood_matches_mg3():
    for n in (0, 1, 5, 10, 40):
        lo, hi = C.ci_garwood(np.array([n]), np.array([1.0 / C.L]))      # L E = 1 -> limits are on the count scale
        glo, ghi = mg3.garwood(n)
        assert abs(lo[0] - glo) < 1e-6 and abs(hi[0] - ghi) < 1e-6


def test_profile_rate_against_root_finding():
    c = stats.chi2.ppf(0.95, 1) / 2
    for n in (1, 3, 20):
        ll = lambda x: n * np.log(x) - x
        f = lambda x: ll(n) - ll(x) - c
        lo_r = optimize.brentq(f, 1e-9, n); hi_r = optimize.brentq(f, n, 50 * n + 50)
        lo, hi = C.ci_profile_rate(np.array([n]), np.array([1.0 / C.L]))
        assert abs(lo[0] - lo_r) < 1e-8 and abs(hi[0] - hi_r) < 1e-8
    lo, hi = C.ci_profile_rate(np.array([0]), np.array([1.0 / C.L]))
    assert lo[0] == 0 and abs(hi[0] - c) < 1e-12


def test_ratio_cp_matches_mg3():
    for n1, E1, n2, E2 in ((10, 2158.0, 14, 2500.0), (3, 100.0, 0, 120.0), (5, 50.0, 7, 80.0)):
        # mg3.rate_ratio gives k1/k2 for arguments (n1, E1, n2, E2); C gives k2/k1 -> call mg3 swapped
        ref = mg3.rate_ratio(n2, E2, n1, E1)
        lo, hi = C.ci_ratio_cp(np.array([n1]), np.array([E1]), np.array([n2]), np.array([E2]), 0.95)
        assert abs(lo[0] - ref["lo"]) < 1e-6 * max(1, ref["lo"])
        assert (np.isinf(hi[0]) and np.isinf(ref["hi"])) or abs(hi[0] - ref["hi"]) < 1e-6 * max(1, ref["hi"])


def test_midp_inside_cp():
    n1 = np.array([10, 30, 4]); n2 = np.array([14, 25, 9]); E1 = np.array([2000.0, 3000, 500]); E2 = E1 * 1.1
    a = C.ci_ratio_cp(n1, E1, n2, E2, 0.95); b = C.ci_ratio_midp(n1, E1, n2, E2, 0.95)
    assert np.all(b[0] >= a[0] - 1e-12) and np.all(b[1] <= a[1] + 1e-12)


def test_samplers_mean_hazard():
    rng = np.random.default_rng(1)
    lam = 3e-3
    for shape in C.SHAPES:
        u = C.sample_u(rng, shape, lam, 400000)
        # P(no event in [0, U]) = exp(-H(U)) = exp(-lam U) for every shape (same mean hazard)
        p_cens = np.mean(np.isinf(u))
        assert abs(p_cens - np.exp(-lam * C.U)) < 4e-3, (shape, p_cens)
        # empirical CDF at 150 matches 1 - exp(-H(150))
        F = np.mean(u <= 150.0)
        assert abs(F - (1 - np.exp(-C.cum_hazard(shape, lam, 150.0)))) < 4e-3, shape


def test_stop_rule():
    rng = np.random.default_rng(2)
    c = C.sim_cell(rng, "stop20", 0.02, 0.0, 2000)               # fast rate: 16 runs give ~16 events
    assert set(np.unique(c["nruns"])) <= {16, 32, 48, 64}
    hit = c["n"] >= C.TARGET
    # whenever the run count is below the cap, the target was reached at that batch end
    assert np.all(hit[c["nruns"] < C.CAP])
    # and not reached at the previous batch end
    ev = c["ev"]
    prev = np.array([ev[i, :c["nruns"][i] - 16].sum() for i in range(2000)])
    assert np.all(prev[c["nruns"] > 16] < C.TARGET)


def test_glm_recovers_slope_and_matches_profile_center():
    rng = np.random.default_rng(3)
    S = 0.073
    k = C.K120 * np.exp(-S * (C.OMEGAS - 120))
    E = np.full((1, 4), 1e5)
    n = rng.poisson(C.L * E * k)
    s, se, dev, a = C.glm_fit(n, E)
    # MLE: score zero
    mu = C.L * E * np.exp(a[:, None] + s[:, None] * C.XC)
    assert abs(((n - mu) * C.XC).sum()) < 1e-6 and abs((n - mu).sum()) < 1e-6
    lo, hi = C.ci_slope_profile(n, E, s, se)
    assert lo[0] < s[0] < hi[0] and abs((hi[0] - lo[0]) / (2 * 1.96 * se[0]) - 1) < 0.05


def test_two_sided_binom_matches_scipy():
    rng = np.random.default_rng(4)
    for _ in range(50):
        n = int(rng.integers(0, 40)); k = int(rng.integers(0, n + 1)) if n else 0; p = float(rng.uniform(0.05, 0.95))
        ours = C.binom_two_sided_p(np.array([k]), np.array([n]), np.array([p]))[0]
        ref = stats.binomtest(k, n, p).pvalue if n > 0 else 1.0
        assert abs(ours - ref) < 1e-9, (k, n, p, ours, ref)


def test_exposure_median():
    ex = np.array([[0.0, 100.0, 430.0, 50.0]]); used = np.ones_like(ex, bool)
    s = C._exposure_median(ex, used)[0]
    assert abs(np.minimum(ex, s).sum() - 0.5 * ex.sum()) < 1e-9


def test_decision_codes():
    big = np.array([400]); E = np.array([1e5])
    assert C.part1_decision("a", big, E, big, E)[0] == 0                  # equal, precise -> EQUIVALENT
    assert C.part1_decision("a", big, E, big * 3, E)[0] == 1              # ratio 3 -> DIFFERENT
    assert C.part1_decision("a", np.array([3]), E, np.array([4]), E)[0] == 2
    assert C.part1_decision(None, big, E, big, E)[0] == 2                 # no eligible method -> INCONCLUSIVE


if __name__ == "__main__":
    tests = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    for t in tests:
        t(); print("ok", t.__name__)
    print(f"{len(tests)} passed")

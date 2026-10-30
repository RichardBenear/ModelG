"""Independent recomputation of the MG3b headline numbers from the raw run records (CLAUDE.md: recompute before
stating a result). Shares no code with src/evaluate_mg3b.py, calib_mg3b.py or mg3.py: exposures are rebuilt from the
first confirmed event's detection time and t_end, CIs come from scipy.stats directly, and the GLM is fitted by
scipy.optimize on the full (alpha, S) likelihood. Run: python -I results/MG3b/independent/recompute.py
"""
import glob, os, pickle, json
import numpy as np
from scipy import stats, optimize

D = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
L, TB, TC = 100.0, 50.0, 480.0


def cell(name):
    n = E = st = 0
    times = []
    for f in glob.glob(f"{D}/runs/{name}/*.pkl"):
        r = pickle.load(open(f, "rb"))
        t = None if r["first"] is None else r["first"]["t_start"]
        if t is not None and t < TB:
            st += 1
        elif t is not None and t <= TC:
            n += 1; E += t - TB; times.append(t)
        else:
            assert r["t_end"] >= 500 - 1e-6 or t is not None, (name, f)
            E += TC - TB
    return dict(runs=len(glob.glob(f"{D}/runs/{name}/*.pkl")), n=n, E=E, startup=st, times=times)


def garwood(n, E):
    lo = stats.gamma.ppf(0.025, n) if n else 0.0
    return lo / (L * E), stats.gamma.ppf(0.975, n + 1) / (L * E)


def ratio_cp(n1, E1, n2, E2, level):
    ci = stats.binomtest(n2, n1 + n2).proportion_ci(level, "exact")
    f = lambda p: p / (1 - p) * E1 / E2
    return f(ci.low), f(ci.high)


out = {}
c = {k: cell(k) for k in ("P1_tau", "P1_tau2", "P2_130", "P2_110", "P2_100")}
a, b = c["P1_tau"], c["P1_tau2"]
r90, r95 = ratio_cp(a["n"], a["E"], b["n"], b["E"], 0.90), ratio_cp(a["n"], a["E"], b["n"], b["E"], 0.95)
eq = 0.5 <= r90[0] and r90[1] <= 2.0
out["part1"] = dict(n=(a["n"], b["n"]), E=(a["E"], b["E"]), startup=(a["startup"], b["startup"]),
                    r=(b["n"] / b["E"]) / (a["n"] / a["E"]), ci90=r90, ci95=r95,
                    decision="EQUIVALENT" if eq else ("DIFFERENT" if r95[0] > 1 or r95[1] < 1 else "INCONCLUSIVE"))
Om = np.array([130.0, 120.0, 110.0, 100.0])
cs = [c["P2_130"], c["P1_tau"], c["P2_110"], c["P2_100"]]
n = np.array([x["n"] for x in cs], float); E = np.array([x["E"] for x in cs])
out["cells"] = {f"{o:g}": dict(runs=x["runs"], n=x["n"], E=x["E"], startup=x["startup"], k=x["n"] / (L * x["E"]),
                               ci=garwood(x["n"], x["E"])) for o, x in zip(Om, cs)}


def nll(p, Om, n, E):
    mu = L * E * np.exp(p[0] - p[1] * (Om - 115.0))
    return float(np.sum(mu - n * np.log(mu)))


def fit(Om, n, E):
    res = optimize.minimize(nll, [np.log(n.sum() / (L * E.sum())), 0.07], args=(Om, n, E), method="Nelder-Mead",
                            options=dict(xatol=1e-10, fatol=1e-12, maxiter=20000))
    S = res.x[1]; m0 = res.fun
    prof = lambda s: optimize.minimize_scalar(lambda al: nll([al, s], Om, n, E), bounds=(-30, 0), method="bounded",
                                              options=dict(xatol=1e-12)).fun - m0 - stats.chi2.ppf(0.95, 1) / 2
    lo = optimize.brentq(prof, S - 0.1, S); hi = optimize.brentq(prof, S, S + 0.1)
    mu = L * E * np.exp(res.x[0] - S * (Om - 115.0))
    dev = 2 * np.sum(np.where(n > 0, n * np.log(n / mu), 0) - (n - mu))
    return dict(S=S, profile95=(lo, hi), deviance=dev, p=stats.chi2.sf(dev, len(n) - 2))


out["slope"] = fit(Om, n, E)
out["slope_without_100"] = fit(Om[:3], n[:3], E[:3])
print(json.dumps(out, indent=1, default=float))

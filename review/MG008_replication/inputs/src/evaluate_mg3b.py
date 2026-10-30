"""MG3b evaluator (GATES.md MG3b; spec/MG3b_followup_proposal.md Parts 1 and 2; methods from
results/MG3b/calibration/decisions.json and the Part 3 addendum 1 labels).

Reads results/MG3b/runs/<cell>/<rep>.pkl, plan.json and budget_ledger.csv. Writes results/MG3b/evaluation.json and
prints a summary. Everything is recomputed from the run records; the driver's own decisions are checked against it.
Run: python src/evaluate_mg3b.py
"""
import os, sys, json, glob, pickle
import numpy as np
from scipy import stats
sys.path.insert(0, os.path.dirname(__file__))
import mg3b
from calib_mg3b import ci_garwood, ci_ratio_cp, part1_decision, stationary_p, MARGIN

ROOT = os.path.join(os.path.dirname(__file__), "..")
OUT = os.environ.get("MG3B_OUT", os.path.join(ROOT, "results", "MG3b"))
L, TB, TC, U = 100.0, 50.0, 480.0, 430.0
CELLS = {130.0: "P2_130", 110.0: "P2_110", 100.0: "P2_100"}
POWER_R, POWER_SEED, ALPHA = 2000, 20261030, 0.01
# Registered labels (GATES.md MG3b: Part 3 result and addendum 1)
SLOPE_LABEL = "approximate (observed minimum coverage 0.937 over qualifying calibration cells; not a guaranteed floor)"
OM120_LABEL = {"tau": 0.943, "tau2_inc": 0.875, "tau2_1b": 0.841}


def load(cell):
    fs = sorted(glob.glob(f"{OUT}/runs/{cell}/*.pkl"), key=lambda f: int(os.path.basename(f)[:-4]))
    return [pickle.load(open(f, "rb")) for f in fs]


def t_first(r):
    return None if r["first"] is None else float(r["first"]["t_start"])


def cell_stats(runs):
    """n, exposure, startup, Garwood k-hat (decisions.json: per-cell rate CI (a)), risk sets, landmarks."""
    assert all(r["complete"] for r in runs), "incomplete run in a finished cell"
    n = int(sum(r["event"] for r in runs)); E = float(sum(r["exposure"] for r in runs))
    lo, hi = ci_garwood(np.array([n]), np.array([E]))
    tf = [t_first(r) for r in runs]
    alive = lambda t: sum(1 for r, x in zip(runs, tf) if not r["startup"] and (x is None or x >= t))
    out = dict(n_runs=len(runs), n=n, exposure=E, n_startup=int(sum(r["startup"] for r in runs)),
               k=n / (L * E), k_lo=float(lo[0]), k_hi=float(hi[0]),
               risk_set={str(t): alive(t) for t in (50, 100, 150)},
               below_target_precision=bool(n < 20 and len(runs) >= 64))
    lm = {}
    for t0 in (100.0, 150.0):          # left truncation at the landmark: runs alive at t0, exposure from t0
        rs = [(r, x) for r, x in zip(runs, tf) if not r["startup"] and (x is None or x >= t0)]
        nn = sum(1 for r, x in rs if r["event"])
        EE = sum((x if r["event"] else TC) - t0 for r, x in rs)
        a, b = ci_garwood(np.array([nn]), np.array([float(EE)]))
        lm[str(int(t0))] = dict(risk_set=len(rs), n=nn, exposure=EE, k=nn / (L * EE) if EE > 0 else float("nan"),
                                k_lo=float(a[0]), k_hi=float(b[0]))
    # direct comparison: k on [100, 480] / k on [50, 100), disjoint exposure, Clopper-Pearson (rate-ratio method (a))
    ok = [(r, x) for r, x in zip(runs, tf) if not r["startup"]]
    n1 = sum(1 for r, x in ok if r["event"] and x < 100.0)
    E1 = sum(min(x if r["event"] else TC, 100.0) - TB for r, x in ok)
    n2 = sum(1 for r, x in ok if r["event"] and x >= 100.0)
    E2 = sum(max((x if r["event"] else TC) - 100.0, 0.0) for r, x in ok)
    rl, rh = ci_ratio_cp(np.array([n1]), np.array([E1]), np.array([n2]), np.array([E2]), 0.95)
    lm["ratio_late_over_early"] = dict(n_early=n1, E_early=E1, n_late=n2, E_late=E2,
                                       ratio=(n2 / E2) / (n1 / E1) if n1 > 0 and E2 > 0 else float("nan"),
                                       lo=float(rl[0]), hi=float(rh[0]))
    out["landmarks"] = lm
    return out


def as_arrays(runs):
    """One replicate in calib_mg3b's array form: u = post-burn-in event time (inf if none), st, ev, ex, used."""
    u = np.array([[t_first(r) - TB if r["event"] else np.inf for r in runs]])
    st = np.array([[bool(r["startup"]) for r in runs]])
    ex = np.array([[float(r["exposure"]) for r in runs]])
    return dict(u=u, st=st, ev=np.isfinite(u) & ~st, ex=ex, used=np.ones_like(st))


def stationary(runs, design, seed):
    """Registered split (b): the cell's exposure median. Power by simulating the cell's procedure (batching, stopping,
    cap, p_s hat, median recomputed per replicate) with hazard k-hat before the observed split and c k-hat after."""
    c = as_arrays(runs)
    from calib_mg3b import _exposure_median
    s_obs = float(_exposure_median(np.where(c["used"] & ~c["st"], c["ex"], 0.0), c["used"] & ~c["st"])[0])
    p_obs = float(stationary_p(c, "median")[0])
    n = c["ev"].sum(); E = c["ex"].sum(); lam = n / E                # per-run hazard (k L)
    ps = float(c["st"].mean())
    rng = np.random.default_rng(seed)
    pw = {}
    for fac in (2.0, 0.5):
        pw[str(fac)] = float((stationary_p(sim_design(rng, design, lam, ps, s_obs, fac), "median") < ALPHA).mean())
    return dict(split=s_obs + TB, p=p_obs, power_x2=pw["2.0"], power_x0_5=pw["0.5"],
                label="informative" if min(pw.values()) >= 0.5 else "uninformative")


def sim_design(rng, design, lam, ps, s, fac, R=POWER_R):
    """Piecewise-constant hazard lam on [0, s), fac*lam on [s, U] (post-burn-in time); fixed64 or stop20."""
    nmax = 64
    st = rng.random((R, nmax)) < ps
    e1 = rng.exponential(1.0, (R, nmax)) / lam
    e2 = s + rng.exponential(1.0, (R, nmax)) / (fac * lam)
    u = np.where(e1 < s, e1, e2)
    u = np.where(u <= U, u, np.inf)
    ev = (~st) & np.isfinite(u)
    ex = np.where(st, 0.0, np.minimum(u, U))
    if design == "stop20":
        evb = ev.reshape(R, 4, 16).sum(2).cumsum(1)
        hit = evb >= 20
        j = np.where(hit.any(1), np.argmax(hit, 1), 3)
        used = np.arange(nmax)[None, :] < ((j + 1) * 16)[:, None]
    else:
        used = np.ones((R, nmax), bool)
    return dict(u=u, st=st, ev=ev, ex=ex, used=used)


# ------------------------------------------------------------------------------------------------ Poisson GLM
def glm(Om, n, E):
    """log k = alpha - S Omega; Newton on (alpha, s) with s = -S; profile-likelihood and Wald 95% CIs for S,
    deviance lack of fit. Omega centred at its mean."""
    x = Om - Om.mean(); off = np.log(L * E)
    beta = np.array([np.log(n.sum() / (L * E).sum()), 0.0])
    X = np.column_stack([np.ones_like(x), x])
    for _ in range(100):
        mu = np.exp(X @ beta + off)
        step = np.linalg.solve(X.T @ (mu[:, None] * X), X.T @ (n - mu))
        beta += step
        if np.max(np.abs(step)) < 1e-13:
            break
    mu = np.exp(X @ beta + off)
    cov = np.linalg.inv(X.T @ (mu[:, None] * X))
    s, se = beta[1], float(np.sqrt(cov[1, 1]))
    ll = lambda b: float(np.sum(n * (X @ b + off) - np.exp(X @ b + off)))

    def prof(sv):                                   # alpha profiled out in closed form
        a = np.log(n.sum() / np.sum(np.exp(sv * x + off)))
        return ll(np.array([a, sv]))
    l0 = prof(s); c = stats.chi2.ppf(0.95, 1) / 2
    from scipy.optimize import brentq
    hi = brentq(lambda v: l0 - prof(v) - c, s, s + 30 * se)
    lo = brentq(lambda v: l0 - prof(v) - c, s - 30 * se, s)
    with np.errstate(divide="ignore", invalid="ignore"):
        dev = float(2 * np.sum(np.where(n > 0, n * np.log(n / mu), 0.0) - (n - mu)))
    df = len(n) - 2
    z = stats.norm.ppf(0.975)
    return dict(S=-s, se=se, S_profile_lo=-hi, S_profile_hi=-lo, S_wald_lo=-s - z * se, S_wald_hi=-s + z * se,
                alpha_at_mean=float(beta[0]), Omega_mean=float(Om.mean()), deviance=dev, df=df,
                deviance_p=float(stats.chi2.sf(dev, df)) if df > 0 else float("nan"),
                fitted_k={f"{o:g}": float(m / (L * e)) for o, m, e in zip(Om, mu, E)})


# ------------------------------------------------------------------------------------------------ Part 1 startup
def newcombe(x1, n1, x2, n2, level=0.95):
    """Newcombe (1998) hybrid score CI for p1 - p2 (Wilson intervals)."""
    z = stats.norm.ppf(1 - (1 - level) / 2)

    def wilson(x, n):
        p = x / n; d = 1 + z * z / n
        c = (p + z * z / (2 * n)) / d; h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
        return c - h, c + h
    p1, p2 = x1 / n1, x2 / n2
    l1, u1 = wilson(x1, n1); l2, u2 = wilson(x2, n2)
    d = p1 - p2
    return d, d - np.sqrt((p1 - l1) ** 2 + (u2 - p2) ** 2), d + np.sqrt((u1 - p1) ** 2 + (p2 - l2) ** 2)


def part1(plan, dec):
    a, b = load("P1_tau"), load("P1_tau2")
    ca, cb = cell_stats(a), cell_stats(b)
    n1, E1, n2, E2 = (np.array([v]) for v in (ca["n"], ca["exposure"], cb["n"], cb["exposure"]))
    l90, h90 = ci_ratio_cp(n1, E1, n2, E2, 0.90); l95, h95 = ci_ratio_cp(n1, E1, n2, E2, 0.95)
    d = {0: "EQUIVALENT", 1: "DIFFERENT", 2: "INCONCLUSIVE"}[int(part1_decision(dec["rate_ratio_ci"]["method"],
                                                                                  n1, E1, n2, E2)[0])]
    xa, xb = ca["n_startup"], cb["n_startup"]
    diff, dl, dh = newcombe(xa, len(a), xb, len(b))
    fisher = float(stats.fisher_exact([[xa, len(a) - xa], [xb, len(b) - xb]]).pvalue)
    surv = {c: sorted(t_first(r) for r in rr if r["first"] is not None) for c, rr in (("P1_tau", a), ("P1_tau2", b))}
    return dict(tau=ca, tau2=cb, ratio=float((n2 / E2)[0] / (n1 / E1)[0]),
                ci90=[float(l90[0]), float(h90[0])], ci95=[float(l95[0]), float(h95[0])], margin=list(MARGIN),
                decision=d, driver_decision=plan["part1"], agrees_with_driver=d == plan["part1"],
                startup=dict(frac_tau=xa / len(a), frac_tau2=xb / len(b), diff=float(diff),
                             newcombe95=[float(dl), float(dh)], fisher_p=fisher,
                             wording="a difference was detected at the 0.05 level" if fisher < 0.05
                             else "no difference was detected at the 0.05 level"),
                first_event_times=surv)


def main():
    plan = json.load(open(f"{OUT}/plan.json"))
    dec = json.load(open(os.path.join(OUT, "calibration", "decisions.json")))
    res = dict(plan_branch=plan["branch"], production_tau=plan["production_tau"], platform=plan["platform"])
    res["part1"] = part1(plan, dec)
    if plan["branch"] == "stop":
        res["part2"] = None
    else:
        cells = {120.0: plan["omega120_cell"], **CELLS}
        p2 = {}
        for i, (Om, cell) in enumerate(sorted(cells.items(), reverse=True)):
            runs = load(cell)
            cs = cell_stats(runs)
            if cell.startswith("P2_"):
                done, need = mg3b.stop20_done(runs)
                cs["stop20_ok"] = bool(done and need == len(runs))
                design = "stop20"
            else:
                design = "fixed64"
                cs["ci_label"] = f"approximate (observed minimum branch-conditional coverage {OM120_LABEL[plan['branch']]})"
            cs["stationary"] = stationary(runs, design, POWER_SEED + i)
            cs["cell"] = cell
            p2[f"{Om:g}"] = cs
        Om = np.array(sorted(cells, reverse=True))
        n = np.array([p2[f"{o:g}"]["n"] for o in Om], float); E = np.array([p2[f"{o:g}"]["exposure"] for o in Om])
        fit = glm(Om, n, E); fit["label"] = SLOPE_LABEL
        sens = glm(Om[Om != 100.0], n[Om != 100.0], E[Om != 100.0]); sens["label"] = SLOPE_LABEL
        res["part2"] = dict(cells=p2, slope=fit, slope_without_100=sens,
                            conditional_on="production tau " + str(plan["production_tau"]))
    rows = mg3b.ledger_rows(f"{OUT}/budget_ledger.csv")
    by = {}
    for r in rows:
        by[r["item"]] = by.get(r["item"], 0.0) + float(r["wall_s"]) / 3600
    res["ledger_core_h"] = dict(by_item=by, total=sum(by.values()),
                                production_total=sum(v for k, v in by.items() if "calibration" not in k))
    json.dump(res, open(f"{OUT}/evaluation.json", "w"), indent=1, default=float)
    p = res["part1"]
    print(f"Part 1: n {p['tau']['n']}/{p['tau2']['n']}, E {p['tau']['exposure']:.1f}/{p['tau2']['exposure']:.1f}, "
          f"r = {p['ratio']:.3f}, 90% CI {p['ci90'][0]:.3f}-{p['ci90'][1]:.3f}, 95% CI {p['ci95'][0]:.3f}-"
          f"{p['ci95'][1]:.3f} -> {p['decision']} (driver {p['driver_decision']}); startup {p['tau']['n_startup']}/"
          f"{p['tau2']['n_startup']}, Fisher p {p['startup']['fisher_p']:.3g}")
    if res["part2"]:
        for k, c in res["part2"]["cells"].items():
            st = c["stationary"]
            print(f"Omega {k}: runs {c['n_runs']}, n {c['n']}, startup {c['n_startup']}, E {c['exposure']:.0f}, "
                  f"k {c['k']:.3e} [{c['k_lo']:.3e}, {c['k_hi']:.3e}], risk {c['risk_set']}, stationary p "
                  f"{st['p']:.3g} power {st['power_x2']:.2f}/{st['power_x0_5']:.2f} {st['label']}")
        for nm in ("slope", "slope_without_100"):
            f = res["part2"][nm]
            print(f"{nm}: S {f['S']:.4f} profile [{f['S_profile_lo']:.4f}, {f['S_profile_hi']:.4f}] Wald "
                  f"[{f['S_wald_lo']:.4f}, {f['S_wald_hi']:.4f}], deviance {f['deviance']:.2f} df {f['df']} "
                  f"p {f['deviance_p']:.3g}")
    print("ledger core-h:", {k: round(v, 3) for k, v in by.items()}, "total", round(res["ledger_core_h"]["total"], 3))


if __name__ == "__main__":
    main()

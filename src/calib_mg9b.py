"""MG9b synthetic calibration (spec/MG9b_proposal.md draft 3 §4.3; registered 2026-11-05, frozen at 185aceb).

Each replicate:
1. Historical b = 29.00 with MG3b's selection: Part 1 tau and tau/2 arms at Omega = 120 (64 runs each, equal true
   rates, p_s 0.14 / 0.16), MG3b's registered Part 1 rule (ratio method a, 90% CI inside [0.5, 2]). Only EQUIVALENT
   replicates are retained. Part 2 cells 130, 110, 100 by stop20 (p_s 0.00, 0.40, 0.34); the tau arm is the 120 cell.
2. Each new b (29.20, 29.35) through the §2 algorithm: Part A descent with the 3-of-4 stop and 2-of-4 top rule;
   Part B candidates in order; stop20 production. One screen draw per screening decision fails with prob. pi_f.
3. Fits (mg9b.fit_trend), branch, label inputs.
First-passage model: startup with prob. p_s (zero exposure, no event); otherwise constant hazard k L after t = 50,
censored at 480. Rates k_b(Omega) = 3.41e-5 exp(-S_b (Omega - Omega0_b)).

Usage: python src/calib_mg9b.py [--reps R] [--cases O-0.10,...]   (MODELG_WORKERS sets processes)
Output: results/MG9b/calibration/{replicates_<case>.npz, calibration.json, ledger.json}
"""
import os, sys, json, time, platform, argparse
import numpy as np
from concurrent.futures import ProcessPoolExecutor
sys.path.insert(0, os.path.dirname(__file__))
import mg9b
from mg9b import L, U, GRID
from calib_mg3b import part1_decision

K_REF = 3.41e-5
S_HIST = 0.064
OM0_HIST = 120.0
PS_HIST = {130.0: 0.00, 110.0: 0.40, 100.0: 0.34}
PS_P1 = (0.14, 0.16)                       # tau, tau/2 arms
ENTROPY = 20261105
REPS = 4000
CORE_H_CAP = 4.0
OUT = os.path.join(os.path.dirname(__file__), "..", "results", "MG9b", "calibration")
WORKERS = int(os.environ.get("MODELG_WORKERS", "2"))


def _case(name, S_b, Om0, ps, pif, S1):
    return dict(name=name, S_b=S_b, Om0=Om0, ps=ps, pi_f=pif, S1=S1)


def cases():
    """§4.3 scenario table. S_b keyed by b; S1 = None for the non-monotone case (no unique true S1)."""
    out = []
    lin = lambda s1: {29.0: S_HIST, 29.2: S_HIST + s1 * 0.20, 29.35: S_HIST + s1 * 0.35}
    for s1 in (-0.10, -0.05, 0.0, 0.10):
        out.append(_case(f"O{s1:+.2f}", lin(s1), {29.2: 140.0, 29.35: 160.0}, 0.25, 0.0, s1))
    for s1 in (-0.10, 0.0):
        out.append(_case(f"P{s1:+.2f}", lin(s1), {29.2: 140.0, 29.35: 160.0}, 0.50, 0.0, s1))
    for s1 in (-0.10, 0.0):
        out.append(_case(f"Q{s1:+.2f}", lin(s1), {29.2: 140.0, 29.35: 75.0}, 0.25, 0.2, s1))
    out.append(_case("N", {29.0: 0.064, 29.2: 0.080, 29.35: 0.020}, {29.2: 140.0, 29.35: 160.0}, 0.25, 0.0, None))
    return out


# ------------------------------------------------------------------------------------------------ generators
def runs(rng, lam, ps, m):
    """m runs: (events 0/1, exposures). Startup runs: no event, zero exposure."""
    st = rng.random(m) < ps
    u = rng.exponential(1.0, m) / lam
    ev = (~st) & (u <= U)
    ex = np.where(st, 0.0, np.minimum(u, U))
    return ev.astype(int), ex


def stop20_cell(rng, lam, ps):
    ev, ex = runs(rng, lam, ps, mg9b.CAP)
    done, need = mg9b.stop20_need(list(ev))
    return int(ev[:need].sum()), float(ex[:need].sum()), need, bool(need >= mg9b.CAP and ev[:need].sum() < mg9b.TARGET)


def historical(rng):
    """Returns (retained, cells) for b = 29.00 under MG3b's Part 1 selection and Part 2 rule."""
    lam = K_REF * L
    ev1, ex1 = runs(rng, lam, PS_P1[0], 64)
    ev2, ex2 = runs(rng, lam, PS_P1[1], 64)
    d = part1_decision("a", np.array([ev1.sum()]), np.array([ex1.sum()]), np.array([ev2.sum()]),
                       np.array([ex2.sum()]))[0]
    if d != 0:
        return False, None
    cells = [(120.0, int(ev1.sum()), float(ex1.sum()))]
    for Om, ps in PS_HIST.items():
        n, E, _, _ = stop20_cell(rng, K_REF * L * np.exp(-S_HIST * (Om - OM0_HIST)), ps)
        cells.append((Om, n, E))
    return True, cells


def new_b(rng, S, Om0, ps, pif):
    """Part A placement, Part B screens and production for one new b. Returns dict(cells, Om_T, runs, capped)."""
    k = lambda Om: K_REF * np.exp(-S * (Om - Om0))
    counts, nruns = [], 0
    while True:
        st = mg9b.parta_next(counts)
        if st["stop"]:
            break
        Om = st["next_Om"]
        if pif > 0 and rng.random() < pif:
            counts.append(dict(Om=Om, screen_ok=False, count=None))
            continue
        ev, _ = runs(rng, k(Om) * L, ps, mg9b.NA)
        nruns += mg9b.NA
        counts.append(dict(Om=Om, screen_ok=True, count=int(ev.sum())))
    Om_T = st["Om_T"]
    if Om_T is None:
        return dict(cells=[], Om_T=None, runs=nruns, capped=0, placed=False)
    res = []
    while True:
        acc, nxt = mg9b.partb_cells(Om_T, res)
        if nxt is None:
            break
        res.append(not (pif > 0 and rng.random() < pif))
    cells, capped = [], 0
    for Om in acc:
        n, E, need, cap = stop20_cell(rng, k(Om) * L, ps)
        cells.append((Om, n, E)); nruns += need; capped += cap
    return dict(cells=cells, Om_T=Om_T, runs=nruns, capped=capped, placed=True)


def replicate(case, rep):
    rng = np.random.default_rng(np.random.SeedSequence(entropy=ENTROPY, spawn_key=(case["idx"], rep)))
    kept, hcells = historical(rng)
    if not kept:
        return dict(retained=False)
    by_b, info = {29.0: hcells}, {}
    for b in mg9b.B_NEW:
        r = new_b(rng, case["S_b"][b], case["Om0"][b], case["ps"], case["pi_f"])
        by_b[b] = r["cells"]; info[b] = r
    tr = mg9b.fit_trend({b: c for b, c in by_b.items() if c})
    br = mg9b.branch_of(tr.get("b_used", []))
    out = dict(retained=True, branch=br, status=tr["status"], S1=tr.get("S1", np.nan),
               lo=tr.get("ci", (np.nan, np.nan))[0], hi=tr.get("ci", (np.nan, np.nan))[1],
               lr_p=tr.get("lr", {}).get("p", np.nan),
               ncells={b: len(info[b]["cells"]) for b in mg9b.B_NEW},
               placed={b: info[b]["placed"] for b in mg9b.B_NEW},
               runs=sum(info[b]["runs"] for b in mg9b.B_NEW), capped=sum(info[b]["capped"] for b in mg9b.B_NEW))
    return out


def _chunk(args):
    case, r0, r1 = args
    t0 = time.process_time()
    res = [replicate(case, r) for r in range(r0, r1)]
    return res, time.process_time() - t0


# ------------------------------------------------------------------------------------------------ summaries
def _prop(k, n):
    if n == 0:
        return dict(p=None, se=None, k=int(k), n=int(n))
    p = k / n
    return dict(p=float(p), se=float(np.sqrt(p * (1 - p) / n)), k=int(k), n=int(n))


def summarize(case, res):
    R = [r for r in res if r["retained"]]
    S1 = case["S1"]
    out = dict(case=case["name"], S1_true=S1, generated=len(res), retained=len(R),
               branch_freq={b: sum(r["branch"] == b for r in R) for b in ("B3", "B2-20", "B2-35", "B2-new", "B-none")})

    def block(rs, n_den):
        ok = [r for r in rs if r["status"] == "ok"]
        num = sum(r["status"] == "numerical" for r in rs)
        neg = sum(1 for r in ok if r["hi"] < 0)
        pos = sum(1 for r in ok if r["lo"] > 0)
        b = dict(n=len(rs), numerical=num, neg=_prop(neg, n_den), pos=_prop(pos, n_den))
        if S1 is not None:
            ident = [r for r in rs if r["status"] in ("ok", "numerical")]
            cov = sum(1 for r in ok if r["lo"] <= S1 <= r["hi"])      # numerical failures count as not covering
            b["coverage"] = _prop(cov, len(ident))
        lr = [r for r in ok if np.isfinite(r["lr_p"])]
        if lr:
            b["lr_reject_005"] = _prop(sum(r["lr_p"] < 0.05 for r in lr), len(lr))
        return b

    out["pooled"] = block(R, len(R))                      # B-none counts as no directional call
    out["branches"] = {br: block([r for r in R if r["branch"] == br], sum(r["branch"] == br for r in R))
                       for br in ("B3", "B2-20", "B2-35", "B2-new")}
    out["no_trend_fraction"] = _prop(sum(r["status"] == "insufficient" for r in R), len(R))
    out["cells_per_b"] = {str(b): {str(k): sum(r["ncells"][b] == k for r in R) for k in range(4)} for b in mg9b.B_NEW}
    out["not_placed"] = {str(b): sum(not r["placed"][b] for r in R) for b in mg9b.B_NEW}
    out["runs_new_b"] = dict(mean=float(np.mean([r["runs"] for r in R])) if R else None,
                             p95=float(np.percentile([r["runs"] for r in R], 95)) if R else None)
    out["capped_cells_mean"] = float(np.mean([r["capped"] for r in R])) if R else None
    return out


def labels(summ):
    """§4.3 per-branch calibration labels from the eight linear cases."""
    lin = [s for s in summ if s["S1_true"] is not None]
    lab = {}
    lab["B2-new"] = dict(label="not calibrated", reason="unsupported branch (b = 29.00 absent; MG034 r3)")
    for br in ("B3", "B2-20", "B2-35"):
        sup0 = any(s["S1_true"] == 0 and s["branches"][br]["n"] >= mg9b.BRANCH_MIN for s in lin)
        sup1 = any(s["S1_true"] != 0 and s["branches"][br]["n"] >= mg9b.BRANCH_MIN for s in lin)
        if not (sup0 and sup1):
            lab[br] = dict(label="not calibrated", reason="not supported (fewer than 400 retained replicates)")
            continue
        fails = []
        for s in lin:
            b = s["branches"][br]
            if b["n"] < mg9b.BRANCH_MIN:
                continue
            if b["coverage"]["p"] is None or b["coverage"]["p"] < mg9b.COV_MIN:
                fails.append(f"{s['case']}: coverage {b['coverage']['p']}")
            if s["S1_true"] == 0:
                for side in ("neg", "pos"):
                    if b[side]["p"] > mg9b.SIDE_MAX:
                        fails.append(f"{s['case']}: {side} rate {b[side]['p']:.4f}")
        lab[br] = dict(label="calibrated" if not fails else "not calibrated", fails=fails)
    return lab


def run(reps=REPS, names=None):
    os.makedirs(OUT, exist_ok=True)
    cs = cases()
    for i, c in enumerate(cs):
        c["idx"] = i
    if names:
        cs = [c for c in cs if c["name"] in names]
    t0 = time.time()
    summ, cpu = [], 0.0
    chunk = 100
    with ProcessPoolExecutor(WORKERS) as ex:
        for c in cs:
            jobs = [(c, r0, min(reps, r0 + chunk)) for r0 in range(0, reps, chunk)]
            res = []
            for rr, ct in ex.map(_chunk, jobs):
                res += rr; cpu += ct
            s = summarize(c, res)
            summ.append(s)
            print(f"[{c['name']}] retained {s['retained']}/{reps}; branches {s['branch_freq']}; "
                  f"pooled {s['pooled'].get('coverage', {}).get('p')}; cpu so far {cpu / 3600:.2f} core-h", flush=True)
            if cpu / 3600 > CORE_H_CAP:
                print(f"[BUDGET] calibration exceeded {CORE_H_CAP} core-h; stopping", flush=True)
                break
    out = dict(spec="spec/MG9b_proposal.md draft 3 (frozen 185aceb), §4.3", reps_generated=reps, entropy=ENTROPY,
               cases=summ, labels=labels(summ) if len(summ) == len(cases()) else None,
               platform=dict(python=platform.python_version(), numpy=np.__version__, os=platform.platform()),
               workers=WORKERS, wall_s=time.time() - t0, cpu_core_h=cpu / 3600)
    json.dump(out, open(os.path.join(OUT, "calibration.json"), "w"), indent=1, default=float)
    return out


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--reps", type=int, default=REPS)
    ap.add_argument("--cases", default=None)
    a = ap.parse_args()
    o = run(a.reps, a.cases.split(",") if a.cases else None)
    print(json.dumps(o["labels"], indent=1))

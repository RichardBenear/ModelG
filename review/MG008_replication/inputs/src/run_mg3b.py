"""MG3b driver (GATES.md MG3b, with its pre-run implementation specification and Part 3 addendum 1).

Order: decisions.json gate -> tau/4 timing runs -> clip pilots for the new Omega (130, 110) -> Part 1 (tau and tau/2
arms, interleaved by batch of 16) -> Part 1 decision -> Part 1b (tau/4) if DIFFERENT -> production branch ->
Part 2 cells 130, 110, 100 (stop20) at the production tau. The Omega = 120 cell is the selected Part 1 arm.

Resumable at run level: finished runs -> <OUT>/runs/<cell>/<rep>.pkl; unfinished runs checkpoint every 100 time
units to <OUT>/partial/ (gitignored). Every decision is recomputed from the saved runs. Budget: binding run caps and
per-item core-hour ceilings, checked before every batch (a started batch always completes); every finished run is
written to the ledger. A budget stop prints counts and costs only.
Usage: MODELG_WORKERS=9 python src/run_mg3b.py      (MG3B_OUT overrides the output folder, e.g. for a scratch test)
"""
import os, sys, json, pickle, time, glob
import numpy as np
from multiprocessing import Pool
sys.path.insert(0, os.path.dirname(__file__))
import mg3
import mg3b
from mg3b import TAU0, GATE, STAGE, BATCH, om_idx
from runner import reference

ROOT = os.path.join(os.path.dirname(__file__), "..")
OUT = os.environ.get("MG3B_OUT", os.path.join(ROOT, "results", "MG3b"))
DECISIONS = os.environ.get("MG3B_DECISIONS", os.path.join(ROOT, "results", "MG3b", "calibration", "decisions.json"))
RUNS, PART = f"{OUT}/runs", f"{OUT}/partial"
LEDGER = f"{OUT}/budget_ledger.csv"
WORKERS = int(os.environ.get("MODELG_WORKERS", "2"))
ARM_TAU = dict(P1_tau=TAU0, P1_tau2=TAU0 / 2, P1_tau4=TAU0 / 4)


def log(s):
    print(s, flush=True)
    os.makedirs(OUT, exist_ok=True)
    with open(f"{OUT}/run_log.txt", "a") as f:
        f.write(s + "\n")


def load_decisions(path=DECISIONS):
    """Enforcement rule 5: the driver refuses to start without the committed Part 3 decisions."""
    if not os.path.exists(path):
        log(f"[decisions] {path} not found: MG3b refuses to start (GATES.md MG3b, enforcement rule 5)")
        raise SystemExit(5)
    d = json.load(open(path))
    for k in ("per_cell_rate_ci", "rate_ratio_ci", "part1_runs_per_arm", "slope_ci"):
        if k not in d:
            raise SystemExit(f"decisions.json lacks {k}")
    return d


# ------------------------------------------------------------------------------------------------ run level
def run_job(a):
    """a: dict(cell, rep, Om, tau, key, until). Returns the finished record (new=True if completed now)."""
    fin = f"{RUNS}/{a['cell']}/{a['rep']}.pkl"
    if os.path.exists(fin):
        rec = pickle.load(open(fin, "rb")); rec["new"] = False
        return rec
    part = f"{PART}/{a['cell']}/{a['rep']}.pkl"
    st = pickle.load(open(part, "rb")) if os.path.exists(part) else mg3.new_state(29.0, a["Om"], a["tau"], a["key"])
    os.makedirs(os.path.dirname(part), exist_ok=True)
    until = a.get("until", mg3.T_END)
    while not st["done"] and st["t"] < until - 1e-9:
        st = mg3.sim_run(st, until=min(until, st["t"] + 100.0))
        if not st["done"] and st["t"] < until - 1e-9:
            pickle.dump(st, open(part + ".tmp", "wb")); os.replace(part + ".tmp", part)
    rec = mg3b.summarize_b(st); rec.update(rep=a["rep"], cell=a["cell"], item=a["item"])
    os.makedirs(os.path.dirname(fin), exist_ok=True)
    pickle.dump(rec, open(fin + ".tmp", "wb")); os.replace(fin + ".tmp", fin)
    if os.path.exists(part):
        os.remove(part)
    rec["new"] = True
    return rec


def finished(cell):
    fs = sorted(glob.glob(f"{RUNS}/{cell}/*.pkl"), key=lambda f: int(os.path.basename(f)[:-4]))
    return [pickle.load(open(f, "rb")) for f in fs]


def partial_wall(cells):
    return sum(pickle.load(open(f, "rb"))["wall"] for c in cells for f in glob.glob(f"{PART}/{c}/*.pkl"))


def batch_started(cell, reps):
    return any(os.path.exists(f"{RUNS}/{cell}/{r}.pkl") or os.path.exists(f"{PART}/{cell}/{r}.pkl") for r in reps)


def run_batch(pool, cell, item, Om, tau, stage, reps, ctx):
    """One batch of runs, after the per-batch budget check (skipped for a batch already started)."""
    reps = list(reps)
    if all(os.path.exists(f"{RUNS}/{cell}/{r}.pkl") for r in reps):
        return
    if not batch_started(cell, reps):
        cells = ctx["item_cells"][item]
        used = mg3b.ledger_used_s(LEDGER, item) + partial_wall(cells)
        walls = [r["wall"] for c in cells for r in finished(c)]
        ok, info = mg3b.batch_allowed(used, walls, len(reps), item, ctx["n_arm"], ctx["tau4_factor"])
        if not ok:
            ncell = {c: len(finished(c)) for c in cells}
            log(f"[BUDGET] {item}: used {info['used_h']:.2f} + projected {info['projected_h']:.2f} core-h > ceiling "
                f"{info['ceiling_h']:.2f}; batch {cell} reps {reps[0]}-{reps[-1]} not submitted. Runs so far {ncell}. "
                f"MG3b stops here; this item is INCONCLUSIVE by the budget rule unless the PI records an extension.")
            raise SystemExit(2)
    jobs = [dict(cell=cell, item=item, rep=r, Om=float(Om), tau=tau, key=(GATE, 0, om_idx(Om), stage, r)) for r in reps]
    for rec in pool.imap_unordered(run_job, jobs, chunksize=1):
        if rec["new"]:
            mg3b.ledger_append(LEDGER, item, f"{cell}/{rec['rep']}", "run", rec["wall"], WORKERS)
        log(f"  [{cell} rep {rec['rep']}] t_end {rec['t_end']:.0f}, births {len(rec['births'])}, "
            f"{'startup' if rec['startup'] else 'event' if rec['event'] else 'censored'}, wall {rec['wall']:.0f}s")


# ------------------------------------------------------------------------------------------------ preliminaries
def timing(pool):
    """Pre-run spec: 2 runs of 100 time units at tau/4, Omega = 120, timing stream; rescales the Part 1b ceiling."""
    fn = f"{OUT}/timing_tau4.json"
    if os.path.exists(fn):
        return json.load(open(fn))
    jobs = [dict(cell="TIMING_tau4", item="timing", rep=r, Om=120.0, tau=TAU0 / 4, until=100.0,
                 key=(GATE, 0, om_idx(120.0), STAGE["TIMING"], r)) for r in range(2)]
    recs = list(pool.imap(run_job, jobs))
    for r in recs:
        if r["new"]:
            mg3b.ledger_append(LEDGER, "timing", f"TIMING_tau4/{r['rep']}", "timing", r["wall"], WORKERS)
    per = float(np.mean([r["wall"] / r["t_end"] for r in recs]))
    out = dict(core_s_per_time_unit=per, factor=per / mg3b.TAU4_ASSUMED, runs=[(r["t_end"], r["wall"]) for r in recs],
               platform=mg3b._versions(), workers=WORKERS)
    json.dump(out, open(fn, "w"), indent=1)
    log(f"[timing] tau/4: {per:.2f} core-s per time unit (assumed {mg3b.TAU4_ASSUMED}); Part 1b ceiling factor "
        f"{out['factor']:.2f}")
    return out


def clip_pilots():
    """V3-procedure clip pilot at tau for the Omega values not on MG3's grid (130, 110): T = 20, clip fraction
    <= 1e-4. The registered tau is fixed, so a failing pilot stops MG3b for a PI decision."""
    fn = f"{OUT}/clip_pilots.json"
    if os.path.exists(fn):
        res = json.load(open(fn))
    else:
        res = {}
        for Om in (130.0, 110.0):
            st = mg3.new_state(29.0, Om, TAU0, (GATE, 0, om_idx(Om), STAGE["PILOT"], 0))
            st = mg3.sim_run(st, until=20.0)
            res[f"{Om:g}"] = dict(clip_frac=st["clips"] / (3 * mg3.N * st["steps"]), wall=st["wall"])
            mg3b.ledger_append(LEDGER, "pilot", f"PILOT/{Om:g}", "pilot", st["wall"], 1)
        json.dump(res, open(fn, "w"), indent=1)
    bad = {k: v for k, v in res.items() if v["clip_frac"] > 1e-4}
    log(f"[pilot] clip fractions at tau: {{{', '.join(f'{k}: {v['clip_frac']:.1e}' for k, v in res.items())}}}")
    if bad:
        log("[pilot] clip fraction > 1e-4: the registered tau is not usable there; MG3b stops for a PI decision")
        raise SystemExit(4)


# ------------------------------------------------------------------------------------------------ parts
def fixed_arm(pool, cell, n_arm, ctx, Om=120.0):
    for b in range(n_arm // BATCH):
        run_batch(pool, cell, cell, Om, ARM_TAU[cell], STAGE[cell], range(b * BATCH, (b + 1) * BATCH), ctx)
    return finished(cell)


def part2_cell(pool, Om, item, tau, ctx):
    cell = f"P2_{Om:g}"
    while True:
        runs = finished(cell)
        done, need = mg3b.stop20_done(runs)
        if done and len(runs) >= need:
            break
        run_batch(pool, cell, item, Om, tau, STAGE["P2"], range(need - BATCH, need), ctx)
    r = mg3.cell_rate(finished(cell))
    log(f"[cell] {cell}: {r['n_runs']} runs, events {r['n']}, startup {r['n_startup']}, exposure {r['exposure']:.0f}")


def main():
    t0 = time.time()
    dec = load_decisions()
    n_arm, method = int(dec["part1_runs_per_arm"]), dec["rate_ratio_ci"]["method"]
    log(f"=== MG3b driver start {time.strftime('%Y-%m-%d %H:%M:%S')}; {WORKERS} workers; Part 1 arms {n_arm} runs; "
        f"ratio method {method}; output {OUT}")
    reference(29.0, 0.2)                                   # main process, before the pool (shared cache)
    with Pool(WORKERS) as pool:
        tm = timing(pool)
        clip_pilots()
        ctx = dict(n_arm=n_arm, tau4_factor=tm["factor"],
                   item_cells=dict(P1_tau=["P1_tau"], P1_tau2=["P1_tau2"], P1_tau4=["P1_tau4"],
                                   P2_tau=["P2_130", "P2_110", "P2_100"], P2_tau2=["P2_130", "P2_110", "P2_100"]))
        # Part 1: arms interleaved by batch
        for b in range(n_arm // BATCH):
            for cell in ("P1_tau", "P1_tau2"):
                run_batch(pool, cell, cell, 120.0, ARM_TAU[cell], STAGE[cell], range(b * BATCH, (b + 1) * BATCH), ctx)
        out1, c1 = mg3b.compare(method, finished("P1_tau"), finished("P1_tau2"))
        log(f"[Part 1] k(tau/2)/k(tau): {out1} ({c1})")
        out1b, c1b = None, None
        branch, tau, cell120 = mg3b.production_branch(out1)
        if branch == "part1b_needed":
            fixed_arm(pool, "P1_tau4", 64, ctx)            # Part 1b always 64 runs
            out1b, c1b = mg3b.compare(method, finished("P1_tau2"), finished("P1_tau4"))
            log(f"[Part 1b] k(tau/4)/k(tau/2): {out1b} ({c1b})")
            branch, tau, cell120 = mg3b.production_branch(out1, out1b)
        plan = dict(part1=out1, part1_counts=c1, part1b=out1b, part1b_counts=c1b, branch=branch, production_tau=tau,
                    omega120_cell=cell120, decisions=DECISIONS, timing=tm, platform=mg3b._versions(), workers=WORKERS)
        json.dump(plan, open(f"{OUT}/plan.json", "w"), indent=1, default=str)
        log(f"[branch] {branch}; production tau {tau}; Omega = 120 cell {cell120}")
        if branch == "stop":
            log("=== MG3b stops: rate results INCONCLUSIVE on time-step grounds (Part 1b not EQUIVALENT)")
            return
        item = "P2_tau" if branch == "tau" else "P2_tau2"
        for Om in (130.0, 110.0, 100.0):
            part2_cell(pool, Om, item, tau, ctx)
    log(f"=== MG3b driver complete {time.strftime('%Y-%m-%d %H:%M:%S')} (session wall {time.time() - t0:.0f}s); "
        f"ledger total {mg3b.ledger_used_s(LEDGER) / 3600:.2f} core-h")


if __name__ == "__main__":
    main()

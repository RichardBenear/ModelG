"""MG3 gate driver (protocol v1.3.4; GATES.md MG3 with A3, A4 and the pre-run implementation specification).

Resumable at run level: every finished run -> results/MG3/runs/<cell>/<rep>.pkl; unfinished runs checkpoint every
100 time units to results/MG3/partial/ (not committed). All adaptive decisions are recomputed from the saved runs.
Order: Stage A (with crowding screen, midpoints) -> Stage-B selection -> Stage B (+ one extension) -> V5 -> C4 ->
C2(b) reference ensembles -> N1. Budget: 18 core-hours (PI decision 2026-10-28), summed per-run wall time.
Usage: python src/run_mg3.py
"""
import os, sys, json, pickle, time, glob
import numpy as np
from multiprocessing import Pool
sys.path.insert(0, os.path.dirname(__file__))
import mg3
from runner import one_run, seed_for, reference

GATE = 3
CFG = json.load(open(os.path.join(os.path.dirname(__file__), "..", "configs", "MG3_main.json")))
GRID = CFG["Omega_grid"]
B_LIST = [29.0, 29.16, 29.3, 28.5]           # b_idx 0..3 (C4_b list, then the N1 b)
TAU0 = 6.666666666666667e-4
BUDGET = 18 * 3600.0
OUT = os.environ.get("MG3_OUT", "results/MG3")
RUNS, PART = f"{OUT}/runs", f"{OUT}/partial"
STAGE = dict(A=0, B=1, C4=2, N1=3, V5=4, REF=5, PILOT=6)
for d in (RUNS, PART):
    os.makedirs(d, exist_ok=True)


def log(s):
    print(s, flush=True)
    with open(f"{OUT}/run_log.txt", "a") as f:
        f.write(s + "\n")


def om_idx(Om, mids):
    return GRID.index(Om) if Om in GRID else 9 + mids.index(Om)


def used():
    tot = 0.0
    for fn in glob.glob(f"{RUNS}/*/*.pkl") + glob.glob(f"{PART}/*/*.pkl"):
        tot += pickle.load(open(fn, "rb"))["wall"]
    for fn in glob.glob(f"{OUT}/refs/*.pkl") + glob.glob(f"{OUT}/pilots/*.pkl"):
        tot += pickle.load(open(fn, "rb"))["wall"]
    return tot


def started(cell):
    return bool(glob.glob(f"{RUNS}/{cell}/*.pkl") or glob.glob(f"{PART}/{cell}/*.pkl"))


def check_budget(where, cell=None):
    """Budget is checked when a cell (or extension) starts, as in MG2. A cell already started before a forced restart
    (2-hour background limit) is not re-checked, so restarts do not change the stopping rule."""
    if cell is not None and started(cell):
        return used()
    u = used()
    if u > BUDGET:
        log(f"[BUDGET] {u/3600:.2f} core-hours used > cap before {where}; MG3 stops here (INCONCLUSIVE by the budget rule)")
        raise SystemExit(2)
    return u


# ------------------------------------------------------------------------------------------------ run level
def run_job(a):
    """a: dict(cell, rep, b, b_det, Om, tau, key, until)."""
    fin = f"{RUNS}/{a['cell']}/{a['rep']}.pkl"
    if os.path.exists(fin):
        return pickle.load(open(fin, "rb"))
    part = f"{PART}/{a['cell']}/{a['rep']}.pkl"
    st = pickle.load(open(part, "rb")) if os.path.exists(part) else mg3.new_state(a["b"], a["Om"], a["tau"], a["key"], a["b_det"])
    os.makedirs(os.path.dirname(part), exist_ok=True)
    while not st["done"] and st["t"] < a["until"] - 1e-9:
        st = mg3.sim_run(st, until=min(a["until"], st["t"] + 100.0))
        if not st["done"]:
            pickle.dump(st, open(part + ".tmp", "wb")); os.replace(part + ".tmp", part)
    if st["done"]:
        rec = mg3.summarize(st); rec["rep"] = a["rep"]; rec["cell"] = a["cell"]
        os.makedirs(os.path.dirname(fin), exist_ok=True)
        pickle.dump(rec, open(fin + ".tmp", "wb")); os.replace(fin + ".tmp", fin)
        if os.path.exists(part):
            os.remove(part)
        return rec
    return dict(partial=True, cell=a["cell"], rep=a["rep"], t=st["t"], wall=st["wall"],
                n_dets=np.array([len(d) for d in st["dets"]]), times=np.array(st["times"]))


def run_cell(pool, cell, b, Om, tau, stage, b_idx, reps, mids, until=mg3.T_END, b_det=None):
    jobs = [dict(cell=cell, rep=r, b=b, b_det=b if b_det is None else b_det, Om=float(Om), tau=tau,
                 key=(GATE, b_idx, om_idx(Om, mids), stage, r), until=until) for r in reps]
    out = []
    for rec in pool.imap_unordered(run_job, jobs, chunksize=1):
        out.append(rec)
        if not rec.get("partial"):
            fe = rec["first"]
            log(f"  [{cell} rep {rec['rep']}] t_end {rec['t_end']:.0f}, births {len(rec['births'])}, first "
                f"{'none' if fe is None else round(fe['t_start'], 1)} ({'startup' if rec['startup'] else 'event' if rec['event'] else 'censored'}), "
                f"wall {rec['wall']:.0f}s")
    return sorted(out, key=lambda r: r["rep"])


def finished(cell):
    return [pickle.load(open(f, "rb")) for f in sorted(glob.glob(f"{RUNS}/{cell}/*.pkl"))]


# ------------------------------------------------------------------------------------------------ pilots
def clip_tau(pool, b, Om, mids, b_idx):
    """V3-procedure clip pilot for a new (b, Omega): one homogeneous run, T = 20, clip fraction <= 1e-4."""
    os.makedirs(f"{OUT}/pilots", exist_ok=True)
    fn = f"{OUT}/pilots/b{b}_Om{Om:g}.pkl"
    if os.path.exists(fn):
        return pickle.load(open(fn, "rb"))["tau"]
    tau, wall, rows = TAU0, 0.0, []
    for halv in range(3):
        st = mg3.new_state(b, float(Om), tau, (GATE, b_idx, om_idx(Om, mids), STAGE["PILOT"], halv), b_det=29.0)
        st = mg3.sim_run(st, until=20.0)            # T = 20 < persistence: no births, so no confirmations
        frac = st["clips"] / (3 * mg3.N * st["steps"])
        wall += st["wall"]; rows.append((tau, st["clips"], frac))
        if frac <= 1e-4:
            break
        tau /= 2
    pickle.dump(dict(b=b, Om=Om, tau=tau, rows=rows, wall=wall), open(fn, "wb"))
    log(f"[pilot] b={b} Omega={Om:g}: tau={tau:.4e} (rows {rows})")
    return tau


def tau_for(pool, b, Om, mids):
    if b == 29.0 and Om in GRID:
        return TAU0                       # V3 / v3_newgrid: no clips at any grid Omega
    return clip_tau(pool, b, Om, mids, B_LIST.index(b))


# ------------------------------------------------------------------------------------------------ stage A
def stage_a_cell(pool, Om, mids):
    cell = f"A_{Om:g}"
    fin = finished(cell)
    if len(fin) == 4:
        return cell_summary_a(cell, Om, fin)
    check_budget(cell)
    tau = tau_for(pool, 29.0, Om, mids)
    ph1 = run_cell(pool, cell, 29.0, Om, tau, STAGE["A"], 0, range(4), mids, until=mg3.T_B)
    md = []
    for r in ph1:
        if r.get("partial"):
            md.append(float(r["n_dets"][r["times"] < mg3.T_B].mean()))
        else:
            md.append(r["mean_dets_t_lt_50"])
    crowd = float(np.mean(md))
    if crowd > CFG["crowding_max_mean_dets_t_lt_50"]:
        json.dump(dict(cell=cell, Om=Om, excluded=True, mean_dets_t_lt_50=crowd, per_run=md),
                  open(f"{OUT}/excluded_{cell}.json", "w"))
        log(f"[A] {cell}: EXCLUDED by crowding (mean detections/snapshot in t<50 = {crowd:.2f} > 2); descent stops")
        return dict(Om=Om, count=0, excluded=True, any_event=False, crowd=crowd)
    fin = run_cell(pool, cell, 29.0, Om, tau, STAGE["A"], 0, range(4), mids)
    return cell_summary_a(cell, Om, fin)


def cell_summary_a(cell, Om, runs):
    c = mg3.stage_a_count(runs)
    crowd = float(np.mean([r["mean_dets_t_lt_50"] for r in runs]))
    s = dict(Om=Om, count=c, excluded=False, any_event=any(r["first"] is not None for r in runs), crowd=crowd,
             n_startup=sum(r["startup"] for r in runs), n_event=sum(r["event"] for r in runs))
    log(f"[A] {cell}: {c}/4 with a confirmed event in [50,150]; post-burn-in events {s['n_event']}, startup {s['n_startup']}; "
        f"mean dets/snapshot t<50 {crowd:.2f}; budget used {used()/3600:.2f} h")
    return s


def stage_a(pool):
    screened, mids = [], []
    for Om in GRID:
        s = stage_a_cell(pool, Om, mids)
        screened.append(s)
        if s["excluded"] or s["count"] == 4:
            break
    # Amendment 4 §3(c): midpoints between adjacent screened, non-excluded cells bracketing 0/4 -> 4/4
    ok = [c for c in screened if not c["excluded"]]
    if len(ok) >= 2 and ok[-1]["count"] == 4 and ok[-2]["count"] == 0:
        hi, lo = ok[-2]["Om"], ok[-1]["Om"]
        for _ in range(CFG["max_midpoints"]):
            m = mg3.round2sf(np.sqrt(hi * lo))
            if m in (hi, lo) or m in [c["Om"] for c in screened]:
                break
            mids.append(m)
            log(f"[A] midpoint inserted: Omega = {m:g} between {hi:g} and {lo:g}")
            s = stage_a_cell(pool, m, mids)
            screened.append(s)
            if s["excluded"] or 1 <= s["count"] <= 3:
                break
            if s["count"] == 0:
                hi = m
            else:
                lo = m
    return screened, mids


# ------------------------------------------------------------------------------------------------ later stages
def counted_events(runs):
    return sum(r["event"] for r in runs if r["complete"])


def grow_cell(pool, cell, b, Om, tau, stage, b_idx, mids, b_det=None):
    """16 runs; one extension of 16 if fewer than 5 counted events."""
    runs = finished(cell)
    if len(runs) < 16:
        check_budget(cell, cell)
        runs = run_cell(pool, cell, b, Om, tau, stage, b_idx, range(16), mids, b_det=b_det)
    n = counted_events(runs)
    if n < CFG["stageB_min_events"] and len(runs) < 32:
        log(f"[ext] {cell}: {n} events < 5; one extension of 16 runs")
        check_budget(cell + " extension")
        runs = run_cell(pool, cell, b, Om, tau, stage, b_idx, range(32), mids, b_det=b_det)
    r = mg3.cell_rate(runs)
    log(f"[cell] {cell}: {len(runs)} runs, events {r['n']}, startup {r['n_startup']}, exposure {r['exposure']:.0f}, "
        f"k = {r['k']:.3e} [{r['k_lo']:.3e}, {r['k_hi']:.3e}], stationary p = {r['stat_p']}; budget used {used()/3600:.2f} h")
    return runs


def middle(sel):
    s = sorted(sel, reverse=True)
    return s[len(s) // 2]


def v5(pool, Om, mids):
    cell = f"V5_{Om:g}"
    base = finished(f"B_{Om:g}")
    runs = finished(cell)
    if len(runs) < 16:
        check_budget(cell)
        runs = run_cell(pool, cell, 29.0, Om, TAU0 / 2, STAGE["V5"], 0, range(16), mids)
    rb, rh = mg3.cell_rate(base), mg3.cell_rate(runs)
    rr = mg3.rate_ratio(rh["n"], rh["exposure"], rb["n"], rb["exposure"])
    ok = (rr["lo"] <= 1 <= rr["hi"]) and (0.67 <= rr["ratio"] <= 1.5)
    log(f"[V5] Omega={Om:g}: k(tau/2)/k(tau) = {rr['ratio']:.3f} [{rr['lo']:.3f}, {rr['hi']:.3f}] "
        f"(events {rh['n']} vs {rb['n']}) -> {'PASS' if ok else 'FAIL'}")
    return ok


def refs(pool, Om, mids):
    os.makedirs(f"{OUT}/refs", exist_ok=True)
    fn = f"{OUT}/refs/R_{Om:g}.pkl"
    if os.path.exists(fn):
        return
    check_budget(f"R_{Om:g}")
    jobs = [dict(b=29.0, Om=float(Om), h=0.2, tau=TAU0, T=100.0, start="soliton", kind="mg2",
                 rng=seed_for(GATE, 0, om_idx(Om, mids), STAGE["REF"], r)) for r in range(CFG["C2_reference_runs"])]
    res = list(pool.imap(one_run, jobs, chunksize=1))
    keep = ("rng", "eligible", "events", "blocks", "n_per_block", "farvar", "wall", "clip_frac", "n_elig_snap")
    runs = [{k: r[k] for k in keep if k in r} for r in res]
    pickle.dump(dict(Om=Om, runs=runs, wall=sum(r["wall"] for r in runs)), open(fn, "wb"))
    log(f"[ref] R_{Om:g}: {sum(r['eligible'] for r in runs)} eligible of {len(runs)} soliton-start runs (C2(b) needs >= 8)")


def main():
    t0 = time.time()
    log(f"=== MG3 driver start {time.strftime('%Y-%m-%d %H:%M:%S')}; budget used {used()/3600:.2f} h")
    for b in (29.16, 29.3):
        reference(b, 0.2)                 # computed once in the main process (shared cache file)
    with Pool(int(os.environ.get("MODELG_WORKERS", "2"))) as pool:
        screened, mids = stage_a(pool)
        sel, note = mg3.select_stage_b(screened)
        json.dump(dict(screened=screened, midpoints=mids, selected=sel, note=note), open(f"{OUT}/plan.json", "w"), indent=1)
        log(f"[select] Stage-B cells {sel} ({note})")
        if not sel:
            log("=== MG3 Stage A terminal outcome; no Stage B"); return
        for Om in sel:
            grow_cell(pool, f"B_{Om:g}", 29.0, Om, tau_for(pool, 29.0, Om, mids), STAGE["B"], 0, mids)
        mid = middle(sel)
        if not v5(pool, mid, mids):
            log("[V5] FAIL: protocol requires repeating all MG3 runs at tau/2; stopping for a PI budget decision")
            raise SystemExit(3)
        for b in (29.16, 29.3):
            grow_cell(pool, f"C4_b{b}_{mid:g}", b, mid, tau_for(pool, b, mid, mids), STAGE["C4"], B_LIST.index(b), mids)
        for Om in sel:
            refs(pool, Om, mids)
        nOm = min(sel)
        cell = f"N1_b28.5_{nOm:g}"
        if len(finished(cell)) < CFG["N1_runs"]:
            check_budget(cell, cell)
            run_cell(pool, cell, 28.5, nOm, tau_for(pool, 28.5, nOm, mids), STAGE["N1"], 3, range(CFG["N1_runs"]), mids, b_det=29.0)
        log(f"[cell] {cell}: done; budget used {used()/3600:.2f} h")
    log(f"=== MG3 driver complete {time.strftime('%Y-%m-%d %H:%M:%S')} (session wall {time.time()-t0:.0f}s); budget used {used()/3600:.2f} h")


if __name__ == "__main__":
    main()

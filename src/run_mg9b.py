"""MG9b driver (GATES.md MG9b; spec/MG9b_proposal.md draft 3, registered 2026-11-05, frozen at 185aceb).

Order (§6): V-b1 for both b; then at 29.35 Part A, Part B screens, Part B production, V-b3 (Omega_T cell), V-b2;
then the same at 29.20. The synthetic calibration (calib_mg9b.py) runs separately in the cloud, before this.

Stages:  python src/run_mg9b.py vb1     (deterministic reference checks; cheap; allowed before the go)
         python src/run_mg9b.py runs    (everything else; refuses unless results/MG9b/GO exists, written at the PI go)
Resumable at run level. Finished runs: results/MG9b/runs/<cell>/<rep>.json (compact, committed) and
results/MG9b/fields/<cell>/<rep>.npz (confirmation fields, gitignored except the V-b2 sample, copied to vb2_fields.npz).
Unfinished runs checkpoint to results/MG9b/partial/ (gitignored). Every design decision is recomputed from the saved
records, so a restart reproduces it. Budget: one ledger, cutoff 50 core-h, checked before every pilot, screen, batch and
control (cooperative: runs in flight finish and count). Workers: MODELG_WORKERS (9 on the PI's PC).
"""
import os, sys, json, glob, time, pickle, platform, subprocess
import numpy as np
from multiprocessing import Pool
sys.path.insert(0, os.path.dirname(__file__))
import mg3
import mg9b
from mg9b import TAU0, GATE, STAGE, B_IDX, om_idx
from runner import reference
from confirm import confirm_field
from lattice import lattice_ode_1d, P17

ROOT = os.path.join(os.path.dirname(__file__), "..")
OUT = os.environ.get("MG9B_OUT", os.path.join(ROOT, "results", "MG9b"))
WORKERS = int(os.environ.get("MODELG_WORKERS", "2"))
H = mg3.H
PRIOR_RUN_S = 300.0                      # prior per-run cost before any run has finished (MG3b: ~280 core-s on the PC)


class BudgetStop(Exception):
    pass


def P(*a):
    return os.path.join(OUT, *a)


def log(s):
    print(s, flush=True)
    os.makedirs(OUT, exist_ok=True)
    with open(P("run_log.txt"), "a") as f:
        f.write(s + "\n")


def platform_info():
    try:
        commit = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip()
    except Exception:
        commit = None
    return dict(python=platform.python_version(), numpy=np.__version__, os=platform.platform(),
                machine=platform.node(), code_commit=commit, workers=WORKERS)


def _dump(obj, fn):
    os.makedirs(os.path.dirname(fn), exist_ok=True)
    json.dump(obj, open(fn + ".tmp", "w"), indent=1, default=_js)
    os.replace(fn + ".tmp", fn)


def _js(o):
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, (np.floating,)):
        return float(o)
    if isinstance(o, np.ndarray):
        return o.tolist()
    if isinstance(o, tuple):
        return list(o)
    raise TypeError(type(o))


# ------------------------------------------------------------------------------------------------ ledger
def ledger():
    fn = P("ledger.json")
    return json.load(open(fn)) if os.path.exists(fn) else dict(entries=[], budget_stop=None)


def ledger_add(item, cell, rep, wall):
    lg = ledger()
    lg["entries"].append(dict(item=item, cell=cell, rep=rep, wall_s=float(wall), time=time.strftime("%Y-%m-%d %H:%M:%S")))
    _dump(lg, P("ledger.json"))


def used_s():
    part = sum(pickle.load(open(f, "rb"))["wall"] for f in glob.glob(P("partial", "*", "*.pkl")))
    return sum(e["wall_s"] for e in ledger()["entries"]) + part


def mean_run_s():
    w = [e["wall_s"] for e in ledger()["entries"] if e["item"] in ("A", "B")]
    return float(np.mean(w)) if w else PRIOR_RUN_S


def check_budget(n_runs, what, cost_each=None):
    c = mean_run_s() if cost_each is None else cost_each
    used = used_s()
    if used + n_runs * c > mg9b.CUTOFF_CORE_H * 3600:
        lg = ledger()
        lg["budget_stop"] = dict(at=what, used_h=used / 3600, projected_h=n_runs * c / 3600,
                                 time=time.strftime("%Y-%m-%d %H:%M:%S"))
        _dump(lg, P("ledger.json"))
        log(f"[BUDGET] {what}: used {used / 3600:.2f} + projected {n_runs * c / 3600:.2f} core-h > cutoff "
            f"{mg9b.CUTOFF_CORE_H}; MG9b stops here for the PI (partial analysis per §6)")
        raise BudgetStop(what)


# ------------------------------------------------------------------------------------------------ runs
def compact(st):
    """Committed per-run record: first passage, births, confirmations (no fields)."""
    r = mg3.summarize(st)
    times = np.array(st["times"])
    conf = [dict(sc=int(sc), t_conf=float(times[sc]), **cc) for sc, cc in sorted(st["conf_cache"].items())]
    fe = r["first"]
    return dict(key=list(st["key"]), b=st["b"], b_det=st["b_det"], Om=st["Om"], tau=st["tau"], t_end=st["t"],
                wall=st["wall"], clip_frac=r["clip_frac"], mean_dets_t_lt_50=r["mean_dets_t_lt_50"],
                startup=r["startup"], event=r["event"], exposure=r["exposure"], complete=r["complete"],
                t_first=None if fe is None else fe["t_start"],
                first=None if fe is None else {k: v for k, v in fe.items() if k != "field"},
                births=r["births"], confirmations=conf)


def simulate(a):
    """One run from scratch or from its checkpoint, to a['until']. Returns (record, fields dict)."""
    part = P("partial", a["cell"], f"{a['rep']}.pkl")
    if os.path.exists(part):
        st = pickle.load(open(part, "rb"))
    else:
        st = mg3.new_state(a["b"], a["Om"], TAU0, a["key"], b_det=a["b"])
    os.makedirs(os.path.dirname(part), exist_ok=True)
    until = a["until"]
    while not st["done"] and st["t"] < until - 1e-9:
        st = mg3.sim_run(st, until=min(until, st["t"] + 100.0))
        if not st["done"] and st["t"] < until - 1e-9:
            pickle.dump(st, open(part + ".tmp", "wb")); os.replace(part + ".tmp", part)
    rec = compact(st)
    fields = {str(int(sc)): np.asarray(st["phi"][sc], np.float32) for sc in st["conf_cache"]}
    if os.path.exists(part):
        os.remove(part)
    return rec, fields


def run_job(a):
    fin = P("runs", a["cell"], f"{a['rep']}.json")
    if os.path.exists(fin):
        rec = json.load(open(fin)); rec["new"] = False
        return rec
    rec, fields = simulate(a)
    assert rec["b_det"] == rec["b"] == a["b"], "b_det must equal b (draft 3 §2)"
    rec.update(cell=a["cell"], rep=a["rep"], item=a["item"])
    if fields:
        fn = P("fields", a["cell"], f"{a['rep']}.npz")
        os.makedirs(os.path.dirname(fn), exist_ok=True)
        np.savez_compressed(fn, **fields)
    _dump(rec, fin)
    rec["new"] = True
    return rec


def finished(cell):
    fs = glob.glob(P("runs", cell, "*.json"))
    return sorted((json.load(open(f)) for f in fs), key=lambda r: r["rep"])


def cell_name(b, kind, Om):
    return f"b{b:.2f}_{kind}_{Om:g}"


def run_set(pool, b, kind, Om, reps, until, item):
    """Run (or reload) a set of replicates of one cell. Budget is checked for any not-yet-finished rep."""
    cell = cell_name(b, kind, Om)
    stage = dict(A=STAGE["PARTA"], CLIP=STAGE["CLIP"], SCR=STAGE["SCREEN"], B=STAGE["PROD"])[kind]
    todo = [r for r in reps if not os.path.exists(P("runs", cell, f"{r}.json"))]
    if todo:
        cost = mean_run_s() * (until / mg3.T_END) if kind in ("CLIP", "SCR") else None
        check_budget(len(todo), f"{cell} reps {todo[0]}-{todo[-1]}", cost)
    jobs = [dict(b=b, Om=float(Om), cell=cell, rep=r, item=item, until=until,
                 key=(GATE, B_IDX[b], om_idx(Om), stage, r)) for r in reps]
    out = []
    for rec in (pool.imap(run_job, jobs) if pool is not None else map(run_job, jobs)):
        if rec["new"]:
            ledger_add(item, cell, rec["rep"], rec["wall"])
        out.append(rec)
    return out


# ------------------------------------------------------------------------------------------------ Part A / Part B
def clip_ok(pool, b, Om, kind, rep):
    """Clip pilot at tau (no halving): one run to T = 20; clip fraction <= 1e-4."""
    cell = cell_name(b, kind, Om)
    stage = STAGE["CLIP"] if kind == "CLIP" else STAGE["SCREEN"]
    fin = P("runs", cell, f"{rep}.json")
    if not os.path.exists(fin):
        check_budget(1, f"clip pilot {cell}", mean_run_s() * 20 / mg3.T_END)
    rec = run_job(dict(b=b, Om=float(Om), cell=cell, rep=rep, item="pilot", until=20.0,
                       key=(GATE, B_IDX[b], om_idx(Om), stage, rep)))
    if rec["new"]:
        ledger_add("pilot", cell, rep, rec["wall"])
    return rec["clip_frac"] <= mg9b.CLIP_MAX, rec["clip_frac"]


def part_a(pool, b):
    counts, rows = [], []
    while True:
        st = mg9b.parta_next(counts)
        if st["stop"]:
            break
        Om = st["next_Om"]
        ok, cf = clip_ok(pool, b, Om, "CLIP", 0)
        if not ok:
            counts.append(dict(Om=Om, screen_ok=False, count=None))
            rows.append(dict(Om=Om, clip_frac=cf, screen="clip failed"))
            log(f"[A {b:.2f}] Omega {Om:g}: clip fraction {cf:.1e} > 1e-4; unavailable, descent stops")
            continue
        recs = run_set(pool, b, "A", Om, range(mg9b.NA), mg3.T_END, "A")
        crowd = float(np.mean([r["mean_dets_t_lt_50"] for r in recs]))
        cnt = sum(1 for r in recs if r["event"] and r["complete"])
        ok = crowd <= mg9b.CROWD_MAX
        counts.append(dict(Om=Om, screen_ok=ok, count=cnt if ok else None))
        rows.append(dict(Om=Om, clip_frac=cf, crowd=crowd, count=cnt, screen="ok" if ok else "crowded",
                         startup=sum(r["startup"] for r in recs)))
        log(f"[A {b:.2f}] Omega {Om:g}: {cnt}/4 events in [50, 480], crowding {crowd:.2f}"
            + ("" if ok else " > 2: excluded, descent stops"))
    res = dict(b=b, cells=rows, stop_reason=st["reason"], Om_T=st["Om_T"])
    _dump(res, P(f"partA_{b:.2f}.json"))
    log(f"[A {b:.2f}] stop: {st['reason']}; Omega_T = {st['Om_T']}")
    return res


def screen_ok(pool, b, Om):
    """Pre-production screen of a Part B candidate: clip pilot (rep 0) + 4 crowding runs to t = 50 (reps 1-4)."""
    ok, cf = clip_ok(pool, b, Om, "SCR", 0)
    if not ok:
        return False, dict(Om=Om, clip_frac=cf, screen="clip failed")
    recs = run_set(pool, b, "SCR", Om, range(1, 5), 50.0, "screen")
    crowd = float(np.mean([r["mean_dets_t_lt_50"] for r in recs]))
    return crowd <= mg9b.CROWD_MAX, dict(Om=Om, clip_frac=cf, crowd=crowd,
                                          screen="ok" if crowd <= mg9b.CROWD_MAX else "crowded")


def production(pool, b, Om):
    """stop20 batching; a started batch always completes. Returns cell status."""
    cell = cell_name(b, "B", Om)
    while True:
        recs = finished(cell)
        ev = [int(r["event"] and r["complete"]) for r in recs]
        done, need = mg9b.stop20_need(ev)
        if done and len(recs) >= need:
            capped = need >= mg9b.CAP and sum(ev[:need]) < mg9b.TARGET
            return dict(Om=Om, runs=need, events=sum(ev[:need]), capped=capped, budget_short=False)
        start = (len(recs) // mg9b.BATCH) * mg9b.BATCH
        run_set(pool, b, "B", Om, range(start, need), mg3.T_END, "B")


def part_b(pool, b, Om_T):
    fn = P(f"partB_{b:.2f}.json")
    res = json.load(open(fn)) if os.path.exists(fn) else dict(b=b, Om_T=Om_T, screens=[], cells=[])
    while True:
        acc, nxt = mg9b.partb_cells(Om_T, [s["screen"] == "ok" for s in res["screens"]])
        if nxt is None:
            break
        ok, row = screen_ok(pool, b, nxt)
        res["screens"].append(row); _dump(res, fn)
        log(f"[B {b:.2f}] candidate {nxt:g}: {row['screen']}")
    res["accepted"] = acc
    _dump(res, fn)
    cells = []
    for Om in acc:
        try:
            c = production(pool, b, Om)
        except BudgetStop:
            recs = finished(cell_name(b, "B", Om))
            cells.append(dict(Om=Om, runs=len(recs), events=sum(int(r["event"] and r["complete"]) for r in recs),
                              capped=False, budget_short=True))
            res["cells"] = cells; _dump(res, fn)
            raise
        cells.append(c)
        res["cells"] = cells; _dump(res, fn)
        log(f"[B {b:.2f}] cell {Om:g}: {c['runs']} runs, {c['events']} events" + (" (capped)" if c["capped"] else ""))
    return res


# ------------------------------------------------------------------------------------------------ controls
def deep_minima(f, thr):
    l, r = np.roll(f, 1), np.roll(f, -1)
    return np.where((f < l) & (f <= r) & (f < thr))[0]


def parabolic_centre(fX, x):
    i = int(np.argmin(fX)); N = fX.size
    ym, y0, yp = fX[(i - 1) % N], fX[i], fX[(i + 1) % N]
    den = ym - 2 * y0 + yp
    off = 0.5 * (ym - yp) / den if den != 0 else 0.0
    return float(x[i] + off * H)


def vb1(b):
    """§3 V-b1 reference validity: reference continued to t = 2000; checks 1900-2000."""
    fn = P(f"vb1_{b:.2f}.json")
    if os.path.exists(fn):
        return json.load(open(fn))
    os.makedirs(OUT, exist_ok=True)
    ref = reference(b, H)
    N = ref["x"].size; Lx = N * H; x = ref["x"]
    y0 = np.concatenate(ref["phi"])
    t0 = time.time()
    xx, ts, Y, mdt = lattice_ode_1d({**P17, "b": b}, N, H, 2000.0, seed_scale=0.0, y0=y0, t_eval=mg9b.VB2_LATE)
    wall = time.time() - t0
    complete = len(ts) == len(mg9b.VB2_LATE)
    thr = -0.5 * abs(ref["A"])
    circ = lambda a, c: float(min(abs(a - c) % Lx, Lx - abs(a - c) % Lx))
    out = dict(b=b, complete=bool(complete), lam=ref["lam"], A=ref["A"], wall_s=wall, maxdt_end=mdt,
               rhs0=float(lattice_ode_1d({**P17, "b": b}, N, H, 1e-6, seed_scale=0.0, y0=y0)[-1]))
    if complete:
        mins = [deep_minima(Y[N:2 * N, k], thr) for k in range(Y.shape[1])]
        one_each = all(m.size == 1 for m in mins)
        chain = one_each and all(circ(x[mins[k][0]], x[mins[k - 1][0]]) <= ref["lam"] for k in range(1, len(mins)))
        c0 = parabolic_centre(ref["phi"][1], x)
        c1 = parabolic_centre(Y[N:2 * N, -1], x)
        drift = circ(c1, c0)
        depth = float(Y[N:2 * N, -1].min())
        out.update(n_minima=[int(m.size) for m in mins], one_chain=bool(chain), c0=c0, c2000=c1, drift=drift,
                   drift_limit=ref["lam"] / 4, depth_2000=depth, depth_rel=abs(depth - ref["A"]) / abs(ref["A"]),
                   pass_=bool(chain and drift <= ref["lam"] / 4 and abs(depth - ref["A"]) <= 0.05 * abs(ref["A"])))
        np.savez_compressed(P(f"vb1_{b:.2f}_final.npz"), phi=Y[:, -1].reshape(3, N).astype(np.float32))
    else:
        out["pass_"] = False
    out["result"] = "PASS" if out["pass_"] else "REFERENCE-INVALID"
    _dump(out, fn)
    ledger_add("vb", f"vb1_{b:.2f}", 0, wall)
    log(f"[V-b1 {b:.2f}] {out['result']}")
    return out


def _confirm_job(a):
    ref = reference(a["b"], H)
    f = np.load(P("fields", a["cell"], f"{a['rep']}.npz"))[str(a["sc"])].astype(float)
    t0 = time.time()
    c = confirm_field(f, a["b"], H, ref["A"], ref["lam"], **a["kw"])
    c["wall"] = time.time() - t0
    return {**{k: a[k] for k in ("cell", "rep", "sc")}, **c}


def confirmation_order(b, parta, partb):
    """Unique confirmation fields at b in the deterministic §3 order."""
    out = []
    for row in sorted(parta["cells"], key=lambda r: -r["Om"]):
        if row.get("screen") == "clip failed":
            continue
        for r in finished(cell_name(b, "A", row["Om"])):
            out += [(r["cell"], r["rep"], c["sc"], c["t_conf"], c) for c in sorted(r["confirmations"], key=lambda c: c["t_conf"])]
    for Om in partb.get("accepted", []):
        for r in finished(cell_name(b, "B", Om)):
            out += [(r["cell"], r["rep"], c["sc"], c["t_conf"], c) for c in sorted(r["confirmations"], key=lambda c: c["t_conf"])]
    return out


def vb2(pool, b, parta, partb):
    fn = P(f"vb2_{b:.2f}.json")
    if os.path.exists(fn):
        return json.load(open(fn))
    sample = confirmation_order(b, parta, partb)[:10]
    if len(sample) < 10:
        res = dict(b=b, n_fields=len(sample), result="CONTROL-INSUFFICIENT")
        _dump(res, fn); log(f"[V-b2 {b:.2f}] CONTROL-INSUFFICIENT ({len(sample)} fields)"); return res
    check_budget(10, f"V-b2 {b:.2f}", 50.0)
    jobs = [dict(b=b, cell=c, rep=r, sc=sc, kw=dict(checks=mg9b.VB2_LATE)) for c, r, sc, t, _ in sample]
    late = list(pool.imap(_confirm_job, jobs)) if pool is not None else list(map(_confirm_job, jobs))
    pairs = []
    for (c, r, sc, t, std), lt in zip(sample, late):
        ledger_add("vb", f"vb2_{c}", r, lt["wall"])
        pairs.append(dict(cell=c, rep=r, sc=sc, t_conf=t, std_cls=std["cls"], std_complete=std.get("complete", True),
                          late_cls=lt["cls"], late_complete=lt["complete"],
                          agree=(std["cls"] == "persistent") == (lt["cls"] == "persistent")))
    all_complete = all(p["std_complete"] and p["late_complete"] for p in pairs)
    n_agree = sum(p["agree"] for p in pairs)
    result = ("CONTROL-NUMERICAL-FAILURE" if not all_complete else
              "PASS" if n_agree >= 9 else "CONFIRMATION-SENSITIVE")
    res = dict(b=b, n_fields=10, pairs=pairs, all_complete=all_complete, n_agree=n_agree,
               three_class_agree=sum(p["std_cls"] == p["late_cls"] for p in pairs), result=result)
    os.makedirs(OUT, exist_ok=True)
    np.savez_compressed(P(f"vb2_{b:.2f}_fields.npz"),
                        **{f"{c}__{r}__{sc}": np.load(P("fields", c, f"{r}.npz"))[str(sc)] for c, r, sc, _, _ in sample})
    _dump(res, fn); log(f"[V-b2 {b:.2f}] {result} ({n_agree}/10 agree; complete {all_complete})")
    return res


def vb3(pool, b, Om_T):
    fn = P(f"vb3_{b:.2f}.json")
    if os.path.exists(fn):
        return json.load(open(fn))
    cell = cell_name(b, "B", Om_T)
    items = [(r["rep"], c) for r in finished(cell) for c in r["confirmations"]]
    check_budget(len(items), f"V-b3 {b:.2f}", 10.0)
    jobs = [dict(b=b, cell=cell, rep=rep, sc=c["sc"], kw=dict(mg9b.VB3_TOL)) for rep, c in items]
    tight = list(pool.imap(_confirm_job, jobs)) if pool is not None else list(map(_confirm_job, jobs))
    pairs = []
    for (rep, std), tg in zip(items, tight):
        ledger_add("vb", f"vb3_{cell}", rep, tg["wall"])
        pairs.append(dict(rep=rep, sc=std["sc"], std_cls=std["cls"], tight_cls=tg["cls"], std_np=std["n_persistent"],
                          tight_np=tg["n_persistent"], tight_complete=tg["complete"],
                          std_complete=std.get("complete", True),
                          same=(std["cls"] == "persistent") == (tg["cls"] == "persistent")))
    ok = all(p["same"] and p["tight_complete"] and p["std_complete"] for p in pairs)
    res = dict(b=b, cell=cell, n=len(pairs), pairs=pairs, result="PASS" if ok else "CONFIRMATION-TOLERANCE")
    _dump(res, fn); log(f"[V-b3 {b:.2f}] {res['result']} ({len(pairs)} confirmations)")
    return res


# ------------------------------------------------------------------------------------------------ main
def stage_vb1():
    for b in (29.35, 29.20):
        reference(b, H)                                   # main process, before any pool (shared cache)
    return {b: vb1(b) for b in (29.35, 29.20)}


def stage_runs():
    if not os.path.exists(P("GO")):
        log(f"[guard] {P('GO')} not found: production runs need the PI go after ChatGPT's implementation review")
        raise SystemExit(3)
    if ledger().get("budget_stop"):
        log("[guard] a budget stop is recorded; continuing needs a PI decision (clear it explicitly)")
        raise SystemExit(2)
    log(f"=== MG9b runs start {time.strftime('%Y-%m-%d %H:%M:%S')}; {WORKERS} workers; {platform_info()}")
    v1 = stage_vb1()
    with Pool(WORKERS) as pool:
        try:
            for b in (29.35, 29.20):
                if not v1[b]["pass_"]:
                    log(f"[{b:.2f}] REFERENCE-INVALID: Part A/B not run (§3)")
                    continue
                pa = part_a(pool, b)
                if pa["Om_T"] is None:
                    log(f"[{b:.2f}] NOT PLACED")
                    continue
                pb = part_b(pool, b, pa["Om_T"])
                vb3(pool, b, pa["Om_T"])
                vb2(pool, b, pa, pb)
        except BudgetStop:
            log("=== MG9b stopped by the budget rule")
            return
    log(f"=== MG9b runs complete {time.strftime('%Y-%m-%d %H:%M:%S')}; ledger {used_s() / 3600:.2f} core-h")


if __name__ == "__main__":
    st = sys.argv[1] if len(sys.argv) > 1 else "runs"
    {"vb1": stage_vb1, "runs": stage_runs}[st]()

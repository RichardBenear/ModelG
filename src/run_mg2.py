"""MG005: MG2 gate driver (protocol v1.3.3, frozen; criteria in GATES.md).

Resumable: every completed cell is stored in results/MG2/cells/<name>.pkl and skipped on restart.
Budget: 12 core-hours (sum of per-run wall times; raised from 8 by PI decision 2026-10-28; P6 deferred as D3). The order runs the decisive cells first:
  P1  C1 cells   Omega = 1e8, b in {29, 29.16} (with the single extension decided inline)
  P2  C6         h = 0.1, b = 29, Omega = 1e8
  P3  ladder     Omega in {1e4, 3e3, 1e3, 300, 100, 1e5, 1e6}, both b   (C2, C3, C4)
  P4  C5         seeded nucleation, b = 29
  P5  N1         homogeneous start, b = 29, Omega = 1e4
  P6  low Omega  {30, 10}, both b
Stream keys: gate 2; (b_idx, om_idx, stage, replicate). Usage: python src/run_mg2.py [phase ...]
"""
import os, sys, json, pickle, time
import numpy as np
from multiprocessing import Pool
sys.path.insert(0, os.path.dirname(__file__))
from runner import one_run, seed_for, reference
from stats_mg2 import bias_bounds, decide, a1_diagnostics
from measure import smooth_periodic, detect
from confirm import confirm_field

GATE = 2
B = [29.0, 29.16]
LADDER = [10, 30, 100, 300, 1000, 3000, 10000, 100000, 1000000, 100000000]
TAU = {0.2: 6.666666666666667e-4, 0.1: 1.6666666666666666e-4}
BUDGET = 12 * 3600.0   # PI decision 2026-10-28: raised from 8 h; low-Omega cells (p6) deferred as D3
OUT = "results/MG2"
CELLS = f"{OUT}/cells"
os.makedirs(CELLS, exist_ok=True)


def log(s):
    print(s, flush=True)
    with open(f"{OUT}/run_log.txt", "a") as f:
        f.write(s + "\n")


def used_core_seconds():
    tot = 0.0
    for fn in os.listdir(CELLS):
        c = pickle.load(open(f"{CELLS}/{fn}", "rb"))
        tot += sum(r.get("wall", 0.0) for r in c["runs"])
    return tot


def strip(r):
    keep = ("rng", "replicate", "eligible", "events", "blocks", "n_per_block", "farvar", "n_far_snap",
            "positions", "wall", "clip_frac", "n_elig_snap")
    return {k: r[k] for k in keep if k in r}


def mg2_cell(pool, name, b, Om, h=0.2, n_elig=12, cap=24, stage=0, start_rep=0, prior=None):
    """Attempt runs in stream order until n_elig eligible or cap attempted (cumulative with prior)."""
    path = f"{CELLS}/{name}.pkl"
    if os.path.exists(path) and prior is None:
        log(f"[skip] {name} (done)")
        return pickle.load(open(path, "rb"))
    if used_core_seconds() > BUDGET:
        log(f"[BUDGET] cap reached before {name}; stopping")
        raise SystemExit(2)
    b_idx = B.index(b) if b in B else 9
    om_idx = LADDER.index(Om) if Om in LADDER else 99
    runs = list(prior["runs"]) if prior else []
    elig = [r for r in runs if r.get("eligible")]
    target = len(elig) + n_elig if prior else n_elig
    rep = start_rep if not prior else max(r["replicate"] for r in runs) + 1
    attempts = 0
    while len(elig) < target and attempts < cap:
        nb = min(2, target - len(elig), cap - attempts)
        jobs = [dict(b=b, Om=float(Om), h=h, tau=TAU[h], T=100.0, start="soliton", kind="mg2",
                     rng=seed_for(GATE, b_idx, om_idx, stage + (10 if h == 0.1 else 0), rep + i)) for i in range(nb)]
        res = pool.map(one_run, jobs)
        for i, r in enumerate(res):
            r["replicate"] = rep + i
            runs.append(strip(r))
            if r.get("eligible"):
                elig.append(r)
        rep += nb; attempts += nb
    cell = dict(name=name, b=b, Om=Om, h=h, runs=runs, n_eligible=len(elig), n_attempted=len(runs))
    pickle.dump(cell, open(path, "wb"))
    ev = {}
    for r in runs:
        for e in (r.get("events") or []):
            ev[e[0]] = ev.get(e[0], 0) + 1
    log(f"[cell] {name}: b={b} Omega={Om:g} h={h}: {len(elig)} eligible of {len(runs)} attempted; events {ev}; "
        f"core-s {sum(r['wall'] for r in runs):.0f}; budget used {used_core_seconds()/3600:.2f} h")
    return cell


def cell_decision(cell, eps=0.01):
    ref = reference(cell["b"], cell["h"])
    win = np.abs(ref["x"]) < 10.0
    phis = ref["phi"][:, win]
    el = [r for r in cell["runs"] if r.get("eligible")]
    if len(el) < 12:
        return "INCONCLUSIVE", dict(reason="insufficient eligible runs", n=len(el))
    Bk = np.array([r["blocks"] for r in el])
    out = {}
    for s, nm in enumerate("GXY"):
        dg = a1_diagnostics(Bk[:, :, s], phis[s])
        bb = bias_bounds(Bk[:, :, s].mean(1), phis[s])
        out[nm] = dict(a1_ok=dg["ok"], decision=decide(bb, eps) if dg["ok"] else "INCONCLUSIVE",
                       Eb_point=bb["Eb_point"], E_total=bb["E_total"], lower=bb["lower"], upper=bb["upper"])
    ds = [out[nm]["decision"] for nm in "GXY"]
    v = "FAIL" if "FAIL" in ds else ("PASS" if all(d == "PASS" for d in ds) else "INCONCLUSIVE")
    return v, out


def p1_p2(pool):
    for (name, b, h) in (("C1_b29.0_1e8", 29.0, 0.2), ("C1_b29.16_1e8", 29.16, 0.2), ("C6_h0.1_b29.0_1e8", 29.0, 0.1)):
        cell = mg2_cell(pool, name, b, 1e8, h=h)
        v, det = cell_decision(cell)
        log(f"[decide] {name}: {v} " + json.dumps({k: (dict(dec=d['decision'], upper=round(100*d['upper'], 4)) if isinstance(d, dict) else d) for k, d in det.items()}))
        if v == "INCONCLUSIVE" and not cell.get("extended"):
            log(f"[extend] {name}: one extension of 12 eligible runs")
            cell = mg2_cell(pool, name, b, 1e8, h=h, prior=cell)
            cell["extended"] = True
            pickle.dump(cell, open(f"{CELLS}/{name}.pkl", "wb"))
            v, det = cell_decision(cell)
            log(f"[decide] {name} (after extension): {v}")


def p3(pool, oms=(10000, 3000, 1000, 300, 100, 100000, 1000000)):
    for Om in oms:
        for b in B:
            mg2_cell(pool, f"L_b{b}_{Om:g}", b, Om)


def p6(pool):
    p3(pool, oms=(30, 10))


def _c5_job(args):
    r = one_run(dict(args, kind="homog", keep_fields=True))
    ref = reference(29.0, 0.2)
    phi, times = r["phi"].astype(float), r["times"]
    ps = smooth_periodic(phi[:, 1, :], 0.2)
    first, pos = None, None
    for s in range(len(times)):
        d = detect(ps[s], ref["x"], 0.2, ref["thr"], ref["lam"])
        if d.size:
            first, pos = float(times[s]), float(d[0]); break
    out = dict(rng=args["rng"], seed_scale=args["seed_scale"], Om=args["Om"], first_detection=first, wall=r["wall"])
    if first is None:
        out.update(cls="none", nucleated=False)
        return out
    tc = min(first + 20.0, 60.0)
    s = int(np.argmin(np.abs(times - tc)))
    c = confirm_field(phi[s], 29.0, 0.2, ref["A"], ref["lam"], det_pos=pos)
    out.update(t_confirm=float(times[s]), cls=c["cls"], n_persistent=c["n_persistent"], near=c["near_detection"],
               nucleated=(c["cls"] == "persistent"), field=phi[s].astype(np.float32))
    return out


def p4(pool):
    cells = [(1.5, 300), (1.5, 1000), (1.5, 3000), (0.5, 3000), (1.0, 300), (1.0, 1000), (1.0, 3000)]
    for k, (sc, Om) in enumerate(cells):
        name = f"C5_s{sc}_{Om}"
        path = f"{CELLS}/{name}.pkl"
        if os.path.exists(path):
            log(f"[skip] {name}"); continue
        if used_core_seconds() > BUDGET:
            log("[BUDGET] cap reached in P4"); raise SystemExit(2)
        jobs = [dict(b=29.0, Om=float(Om), h=0.2, tau=TAU[0.2], T=60.0, start="homog", seed_scale=sc,
                     rng=seed_for(GATE, 0, LADDER.index(Om), 5, 100 * k + rep)) for rep in range(20)]
        res = pool.map(_c5_job, jobs)
        pickle.dump(dict(name=name, seed_scale=sc, Om=Om, runs=res), open(path, "wb"))
        n = sum(r["nucleated"] for r in res)
        log(f"[cell] {name}: nucleated {n}/20; classes {dict((c, sum(r['cls']==c for r in res)) for c in ('persistent','excursion','unresolved','none'))}; "
            f"budget used {used_core_seconds()/3600:.2f} h")


def p5(pool):
    name = "N1_b29.0_1e4"; path = f"{CELLS}/{name}.pkl"
    if os.path.exists(path):
        log(f"[skip] {name}"); return
    jobs = [dict(b=29.0, Om=1e4, h=0.2, tau=TAU[0.2], T=100.0, start="homog", kind="homog",
                 rng=seed_for(GATE, 0, LADDER.index(10000), 6, rep)) for rep in range(4)]
    res = pool.map(one_run, jobs)
    runs = [dict(rng=r["rng"], births=r["births"], max_dets=int(max(r["n_dets"])), wall=r["wall"]) for r in res]
    pickle.dump(dict(name=name, runs=runs), open(path, "wb"))
    log(f"[cell] {name}: births {[r['births'] for r in runs]}, max detections per snapshot {[r['max_dets'] for r in runs]}")


if __name__ == "__main__":
    phases = sys.argv[1:] or ["p1", "p3", "p4", "p5", "p6"]
    with Pool(int(os.environ.get("MODELG_WORKERS", "2"))) as pool:
        for ph in phases:
            log(f"=== phase {ph} start {time.strftime('%Y-%m-%d %H:%M:%S')}")
            dict(p1=p1_p2, p3=p3, p4=p4, p5=p5, p6=p6)[ph](pool)
            log(f"=== phase {ph} done")

"""MG7 Part 0 (proposal §7; GATES.md MG7 pre-run specification), on UNREGISTERED inputs except check 1, which only
reads and measures the registered start fields (no dynamics).

    $env:MODELG_WORKERS=8; python src/run_mg7.py part0

1. Endpoint inputs: SHA-256 of the 8 MG4c endpoint files against results/MG7/endpoint_manifest.csv; loaded array equals
   the file; each start detects single at theta = -4.45758; max nodal phi_G <= 0.14; P0_start(1.5) and M- recorded
   (kappa = 2, no shift, MG6 grid). Any failure: STOP.
   Then, per endpoint b, on the MG6 Part 0 graph (seed_for(6, 9, 0, 0, 0)) settled frozen at b from the soliton
   (Delta = 1, dt = 0.05, t = 300), with candidate stream seed_for(7, 9, 3, 1, ib), ib = 1 (27.5), 2 (28.5):
2. Step check (forced history): CP2, Delta 0.1, 100 units at dt 0.05; replay at dt 0.025 with the accepted edits
   imposed. Equal counts at every check, depth <= 2%, |D0bar(dt) - D0bar(dt/2)| <= 2 nodes (t = 50-100).
3. Tick check (free acceptance): CP2, Delta 0.1 against Delta 0.05, dt 0.05, same stream. Equal counts, depth <= 3%,
   |D0bar(0.1) - D0bar(0.05)| <= 6 nodes. Centre recorded, not scored (MG021).
4. Bound margin: max nodal phi_G <= 0.14 at every boundary of every Part 0 coupled run.
5. Timing and completion: every Part 0 run COMPLETE at its full horizon; projected cost of the 16 runs restated.
Outputs (committed): results/MG7/part0/decisions.json and check records. Files are never overwritten: an existing
output is a STOP (MG021 lesson), and attempt names carry b and the configuration.
"""
import os, sys, json, time, csv
import numpy as np
from concurrent.futures import ProcessPoolExecutor
sys.path.insert(0, os.path.dirname(__file__))
import graphs as gr
import graph_ode as go
import evolve as ev
import mg6
import run_mg5 as R5
import run_mg6 as R6
import run_mg6_part0 as P6
import run_mg7 as R
from lattice import P17
from runner import seed_for

P0 = os.path.join(R.OUT, "part0")
IB = {27.5: 1, 28.5: 2}
T_SETTLE, T_CHECK = 300.0, 100.0
WIN = (50.0, 100.0)
DEPTH_STEP, D0_STEP, DEPTH_TICK, D0_TICK = 0.02, 2.0, 0.03, 6.0
PART0_CAP = 30.0                                   # core-h, validation priced separately (proposal §6)
MG6_LEDGER = os.path.join(R.ROOT, "MG6", "budget_ledger.csv")


class Exists(RuntimeError):
    pass


def _fresh(path):
    if os.path.exists(path):
        raise Exists(f"{path} exists; Part 0 never overwrites an output (move it aside and record why)")
    return path


# ------------------------------------------------------------------ 1. endpoint inputs
def check_inputs():
    ok_files, rows = R.verify_endpoint_files()
    out = dict(files=rows, files_ok=ok_files, starts={})
    if not ok_files:
        out["ok"] = False
        return out
    good = True
    for b in R.B_END:
        for slot, pl in R.DRUNS:
            g, h, F, sha = R.start_state7(slot, pl, b)
            mesh = ev.Mesh(g["pos"], g["Lx"], g["Ly"])
            same_mesh = bool(np.array_equal(mesh.g["I"], g["I"]) and np.array_equal(mesh.g["J"], g["J"]))
            s = ev.smoothed_X(mesh.g, F)
            parts = ev.detect_full(mesh.pos, s, R.THETA, g["Lx"], g["Ly"])
            mx = float(F[0].max())
            row = dict(sha256=sha, h=h, n_parts=len(parts), max_phiG=mx, mesh_ok=same_mesh)
            if len(parts) == 1:
                grid = mg6.Grid(g["Lx"], g["Ly"], R.NX, R.NY)
                m = mg6.measure(mesh, F, parts[0]["centre"], R.KAPPA, None, grid, R.N0 / grid.A)
                row.update(centre=np.asarray(parts[0]["centre"]).tolist(), depth=float(parts[0]["depth"]),
                           P0_start_15=float(m["P0"][1]), M_minus=float(m["M_minus"]), D0_15=float(m["D0"][1]))
            row["ok"] = bool(len(parts) == 1 and mx <= R.MARGIN and same_mesh)
            good &= row["ok"]
            out["starts"][R.run_key("start", b, slot, pl)[6:]] = row
            print(f"  b {b} D{slot}_p{pl}: parts {len(parts)}, max phi_G {mx:.4f}, P0_start "
                  f"{row.get('P0_start_15', float('nan')):.2f}", flush=True)
    out["ok"] = bool(good)
    return out


# ------------------------------------------------------------------ settle and coupled runs (unregistered)
def settle(args):
    """The MG6 Part 0 graph settled frozen at endpoint b from the soliton (MG5 run_ticks, Delta 1, dt 0.05, t = 300)."""
    b = args["b"]
    path = _fresh(os.path.join(P0, f"settled_{R.bkey(b)}.npz"))
    g, h = P6.graph()
    box = {}

    def snap(t, pos, ids, F):
        box.update(t=t, pos=pos.copy(), F=F.copy())
    F_ic = go.soliton_ic(g, (g["Lx"] / 2, g["Ly"] / 2)).reshape(3, -1)
    res, checks, _ = R5.run_ticks(g["pos"], F_ic, h, R.params_for(b), R.THETA, 0.0, None, 1.0, 0.05, T_SETTLE,
                                  g["Lx"], g["Ly"], snap=snap)
    n_end = len(checks[-1]["parts"]) if checks else None
    os.makedirs(P0, exist_ok=True)
    if box:
        np.savez(path, **box)
    return dict(b=b, status=res["status"], wall_s=res["wall_s"], t_reached=res.get("t_reached"),
                n_parts_end=n_end, h=h, Lx=g["Lx"], Ly=g["Ly"], N=g["N"],
                sha256=R5.sha256(path) if box else None)


def coupled_run(a):
    """One Part 0 CP2 run (kappa 2, no shift) at endpoint b from that endpoint's settled state."""
    b = a["b"]
    tag = a["tag"]
    rec_path = _fresh(os.path.join(P0, f"check_{tag}.npz"))
    acc_path = _fresh(os.path.join(P0, f"accepted_{tag}.npz"))
    with np.load(os.path.join(P0, f"settled_{R.bkey(b)}.npz")) as z:
        pos, F = z["pos"], z["F"]
    Lx, Ly = a["Lx"], a["Ly"]
    cand = mg6.gen_candidates(R.N0, R.NU, a["t_end"], Lx, Ly, a["seed"])
    forced = None
    if a.get("forced_tag"):
        with np.load(os.path.join(P0, f"accepted_{a['forced_tag']}.npz")) as z:
            forced = np.unpackbits(z["bits"])[:int(z["n"])].astype(bool)
    co = mg6.Coupled(cand, R.KAPPA, None, forced=forced)
    grid = mg6.Grid(Lx, Ly, R.NX, R.NY)
    rho0 = R.N0 / grid.A
    res, checks, ticks = R6.run_ticks6(pos, F, a["h"], R.params_for(b), R.THETA, co, a["delta"], a["dt"], a["t_end"],
                                       Lx, Ly, R.KAPPA, None, grid, rho0, seed_audit=1)
    R6.save_checks6(rec_path, checks, ticks, co.accepted)
    np.savez_compressed(acc_path, bits=np.packbits(co.accepted), n=len(co.accepted))
    return dict(tag=tag, b=b, delta=a["delta"], dt=a["dt"], forced=a.get("forced_tag"), res=res,
                lists=[[dict(centre=np.asarray(p["centre"]).tolist(), depth=p["depth"]) for p in c["parts"]]
                       for c in checks],
                t=[c["t"] for c in checks],
                D0_15=[float(c["m"]["D0"][1]) if "m" in c else None for c in checks],
                P0_15=[float(c["m"]["P0"][1]) if "m" in c else None for c in checks],
                cand_sha256=mg6.candidates_hash(cand), record_sha256=R5.sha256(rec_path))


def window_mean(run, lo=WIN[0], hi=WIN[1]):
    v = [x for t, x in zip(run["t"], run["D0_15"]) if lo - 1e-9 <= t <= hi + 1e-9]
    return float(np.mean(v)) if v and None not in v else None


def compare(a, b, Lx, Ly, dtol, d0tol):
    """Equal particle counts at every check, relative depth difference <= dtol, |D0bar difference| <= d0tol nodes. The
    centre difference is recorded, not scored."""
    out = dict(status=[a["res"]["status"], b["res"]["status"]])
    if a["res"]["status"] != R.COMPLETE or b["res"]["status"] != R.COMPLETE:
        out.update(ok=False, why="a run is not COMPLETE")
        return out
    bad, cd, dd = [], [], []
    for i, (la, lb) in enumerate(zip(a["lists"], b["lists"])):
        if len(la) != len(lb):
            bad.append(dict(check=i, why="count", na=len(la), nb=len(lb)))
        elif len(la) >= 1:
            for pa, pb in zip(la, lb):
                cd.append(float(np.hypot(*gr.min_image(np.subtract(pa["centre"], pb["centre"]), Lx, Ly))))
                rel = abs(pa["depth"] - pb["depth"]) / abs(pa["depth"])
                dd.append(rel)
                if rel > dtol:
                    bad.append(dict(check=i, why="depth", depth_rel=rel))
    Da, Db = window_mean(a), window_mean(b)
    dD = None if Da is None or Db is None else abs(Da - Db)
    if dD is None:
        bad.append(dict(why="D0bar not available (no particle in a window check)"))
    elif dD > d0tol:
        bad.append(dict(why="D0bar", diff=dD))
    out.update(ok=not bad and len(a["lists"]) == len(b["lists"]), failures=bad[:20], D0_bar=[Da, Db], D0_diff=dD,
               max_depth_rel=max(dd, default=None), max_centre_diff_recorded=max(cd, default=None))
    return out


def mg6_production_rate():
    rows = list(csv.DictReader(open(MG6_LEDGER, encoding="utf-8")))
    r = {a: [float(x["wall_s"]) / R6.T_END for x in rows if x["key"].startswith(a + "_")] for a in ("CP2", "SH2")}
    return {a: float(np.mean(v)) for a, v in r.items()}


def authorise(dec, cutoff=R.CUTOFF):
    """(proceed, why) from the recorded checks. Pure function of the decision record."""
    if not dec.get("inputs", {}).get("ok"):
        return False, "check 1 (endpoint inputs) failed"
    for b in R.B_END:
        k = R.bkey(b)
        s = dec.get("settle", {}).get(k, {})
        if s.get("status") != R.COMPLETE or s.get("n_parts_end") != 1:
            return False, f"settle at {k} not COMPLETE and single at t = 300"
    for k, r in dec.get("runs", {}).items():
        if r["status"] != R.COMPLETE:
            return False, f"Part 0 run {k} not COMPLETE ({r['status']})"
    for b in R.B_END:
        k = R.bkey(b)
        if not dec.get("step_check", {}).get(k, {}).get("ok"):
            return False, f"check 2 (step) failed at {k}: PI decides"
        if not dec.get("tick_check", {}).get(k, {}).get("ok"):
            return False, f"check 3 (tick) failed at {k}: PI decides"
    if not dec.get("bound_margin", {}).get("ok"):
        return False, "check 4 (bound margin) failed"
    proj = dec.get("timing", {}).get("projected_core_h_16_runs")
    if proj is None or not np.isfinite(proj):
        return False, "projected cost not available"
    if proj > cutoff:
        return False, f"projected cost {proj:.1f} core-h exceeds the {cutoff:g} core-h cutoff: PI decision required"
    return True, None


def finish(dec):
    dec["part0_wall_s"] = time.time() - dec.pop("_t0", time.time())
    path = os.path.join(P0, "decisions.json")
    _fresh(path)
    json.dump(dec, open(path, "w"), indent=1, default=float)
    print("MG7 Part 0:", "PROCEED" if dec.get("proceed") else f"STOP ({dec.get('stop')})", flush=True)
    return dec


def main():
    os.makedirs(P0, exist_ok=True)
    _fresh(os.path.join(P0, "decisions.json"))
    go.continuum2d()
    dec = dict(design="GATES.md MG7 Part 0; unregistered graph seed_for(6, 9, 0, 0, 0); streams seed_for(7, 9, 3, 1, ib)",
               _t0=time.time(), platform=R5.platform_info(), code_revision=R5.code_revision(), workers=R.WORKERS)
    print("check 1: endpoint inputs", flush=True)
    dec["inputs"] = check_inputs()
    if not dec["inputs"]["ok"]:
        dec.update(proceed=False, stop="check 1 (endpoint inputs) failed")
        return finish(dec)

    with ProcessPoolExecutor(min(2, max(1, R.WORKERS))) as ex:
        st = list(ex.map(settle, [dict(b=b) for b in R.B_END]))
    dec["settle"] = {R.bkey(s["b"]): s for s in st}
    for s in st:
        print(f"settled at b {s['b']}: {s['status']}, {s['n_parts_end']} particle(s), {s['wall_s']:.0f} s", flush=True)
    if any(s["status"] != R.COMPLETE or s["n_parts_end"] != 1 for s in st):
        dec.update(proceed=False, stop="settle not COMPLETE and single at t = 300")
        return finish(dec)
    geo = dict(h=st[0]["h"], Lx=st[0]["Lx"], Ly=st[0]["Ly"], t_end=T_CHECK)

    def base(b):
        return dict(geo, b=b, seed=seed_for(R.GATE, 9, 3, 1, IB[b]))

    def tag(b, delta, dt, forced=False):
        return f"{R.bkey(b)}_d{delta:g}_dt{dt:g}" + ("_forced" if forced else "")

    runs = {}
    with ProcessPoolExecutor(min(4, max(1, R.WORKERS))) as ex:
        futs = {}
        for b in R.B_END:
            futs[tag(b, 0.1, 0.05)] = ex.submit(coupled_run, dict(base(b), tag=tag(b, 0.1, 0.05), delta=0.1, dt=0.05))
            futs[tag(b, 0.05, 0.05)] = ex.submit(coupled_run, dict(base(b), tag=tag(b, 0.05, 0.05), delta=0.05,
                                                                   dt=0.05))
        for k, f in futs.items():
            runs[k] = f.result()
            print(f"{k}: {runs[k]['res']['status']}, {runs[k]['res']['wall_s']:.0f} s", flush=True)
        wall1 = sum(r["res"]["wall_s"] for r in runs.values()) + sum(s["wall_s"] for s in st)
        est2 = 2 * 2.0 * float(np.mean([runs[tag(b, 0.1, 0.05)]["res"]["wall_s"] for b in R.B_END]))
        dec["part0_cost"] = dict(core_h_phase1=wall1 / 3600, est_core_h_phase2=est2 / 3600, cap=PART0_CAP)
        if (wall1 + est2) / 3600 > PART0_CAP:
            dec["runs"] = {k: dict(status=r["res"]["status"], wall_s=r["res"]["wall_s"]) for k, r in runs.items()}
            dec.update(proceed=False, stop=f"Part 0 projected {(wall1 + est2) / 3600:.1f} core-h > cap {PART0_CAP}")
            return finish(dec)
        futs = {tag(b, 0.1, 0.025, True): ex.submit(coupled_run, dict(base(b), tag=tag(b, 0.1, 0.025, True), delta=0.1,
                                                                      dt=0.025, forced_tag=tag(b, 0.1, 0.05)))
                for b in R.B_END}
        for k, f in futs.items():
            runs[k] = f.result()
            print(f"{k}: {runs[k]['res']['status']}, {runs[k]['res']['wall_s']:.0f} s", flush=True)

    Lx, Ly = geo["Lx"], geo["Ly"]
    dec["runs"] = {k: dict(status=r["res"]["status"], note=r["res"]["note"], wall_s=r["res"]["wall_s"],
                           t_reached=r["res"]["t_reached"], max_phiG=r["res"]["max_phiG"], b=r["b"],
                           delta=r["delta"], dt=r["dt"], forced=r["forced"], counts=[len(x) for x in r["lists"]],
                           t=r["t"], D0_15=r["D0_15"], P0_15=r["P0_15"], cand_sha256=r["cand_sha256"],
                           record_sha256=r["record_sha256"]) for k, r in runs.items()}
    dec["step_check"], dec["tick_check"] = {}, {}
    for b in R.B_END:
        A, B, C = runs[tag(b, 0.1, 0.05)], runs[tag(b, 0.1, 0.025, True)], runs[tag(b, 0.05, 0.05)]
        dec["step_check"][R.bkey(b)] = compare(A, B, Lx, Ly, DEPTH_STEP, D0_STEP)
        dec["tick_check"][R.bkey(b)] = compare(A, C, Lx, Ly, DEPTH_TICK, D0_TICK)
        print(f"b {b}: step {'ok' if dec['step_check'][R.bkey(b)]['ok'] else 'FAIL'} "
              f"(dD0 {dec['step_check'][R.bkey(b)].get('D0_diff')}), tick "
              f"{'ok' if dec['tick_check'][R.bkey(b)]['ok'] else 'FAIL'} "
              f"(dD0 {dec['tick_check'][R.bkey(b)].get('D0_diff')})", flush=True)
    mx = max(r["res"]["max_phiG"] for r in runs.values())
    dec["bound_margin"] = dict(max_phiG=mx, limit=R.MARGIN, ok=bool(mx <= R.MARGIN))
    rate6 = mg6_production_rate()
    spu_p0 = float(np.mean([runs[tag(b, 0.1, 0.05)]["res"]["wall_s"] / T_CHECK for b in R.B_END]))
    spu = max(np.mean(list(rate6.values())), spu_p0)
    proj = 16 * R.T_END * spu / 3600
    dec["timing"] = dict(mg6_production_s_per_unit=rate6, part0_s_per_unit_delta0p1=spu_p0, s_per_unit_used=spu,
                         projected_core_h_16_runs=proj, n_concurrent_part0=min(4, max(1, R.WORKERS)))
    dec["part0_cost"]["core_h_total"] = (sum(r["res"]["wall_s"] for r in runs.values())
                                         + sum(s["wall_s"] for s in st)) / 3600
    ok, why = authorise(dec)
    dec.update(proceed=ok, stop=why, delta=R.DELTA, dt=R.DT, t_end=R.T_END, window=[100.0, 200.0], nx=R.NX, ny=R.NY,
               cutoff_core_h=R.CUTOFF, run_cap=R.RUN_CAP)
    print(f"projected registered cost: {proj:.1f} core-h", flush=True)
    return finish(dec)

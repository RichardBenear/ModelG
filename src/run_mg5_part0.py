"""MG5 Part 0 (proposal §6; GATES.md MG5 pre-run specification): checks on an UNREGISTERED Delaunay graph before any
registered run, plus the timing and storage measurements used to restate the budget.

    $env:MODELG_WORKERS=9; python src/run_mg5.py part0

1. Solver check (frozen): ROS2 (dt) against MG4's BDF from the continuum initial condition at b = 28 for 50 time units;
   max |dphi| <= 1e-3.
2. Step check: from the settled state (frozen ROS2 to t = 300), nu = 0.2 for 100 time units, the same event list and
   ticks (Delta_fast) at dt and dt/2.
3. Tick check: the same event list at Delta_fast and Delta_fast/2 (both at dt).
Comparison (2, 3): particle counts equal at every check; where both have one particle, centre distance <= 0.1 and
|depth difference| <= 2% of the first run's depth.
Policy: a failure of 1 or 2 halves dt (all runs) and repeats the failed check once; a second failure stops for the PI.
A failure of 3 halves Delta_fast and repeats it once (Delta_fast against Delta_fast/2); a second failure stops.
Output: results/MG5/part0/decisions.json (committed) and part0 check records.
"""
import os, sys, json, time
import numpy as np
from concurrent.futures import ProcessPoolExecutor
sys.path.insert(0, os.path.dirname(__file__))
import graphs as gr
import graph_ode as go
import evolve as ev
import run_mg5 as R
from lattice import P17
from runner import seed_for

P0 = os.path.join(R.OUT, "part0")
B_IDX = 9                                   # stream index used by no registered study
NU_FAST, T_CHECK, T_SETTLE, T_SOLVER = 0.2, 100.0, 300.0, 50.0
SOLVER_TOL, CENTRE_TOL, DEPTH_TOL = 1e-3, 0.1, 0.02


def graph():
    tor = gr.torus(0.2, 122, 142)
    rng = np.random.default_rng(seed_for(ev.GATE, B_IDX, 0, 0, 0))
    pos = gr.poisson_points(R.N0, tor["Lx"], tor["Ly"], rng)
    g = ev.Mesh(pos, tor["Lx"], tor["Ly"]).g
    h = float(np.sqrt(gr.corrector(g)["h2"]))
    return g, h


def _final_capture():
    box = {}

    def snap(t, pos, ids, F):
        box.update(t=t, pos=pos.copy(), ids=ids.copy(), F=F.copy())
    return box, snap


def frozen(g, h, F0, t_end, dt):
    box, snap = _final_capture()
    res, checks, ticks = R.run_ticks(g["pos"], F0, h, {**P17, "b": 28.0}, R.THETA_28, 0.0, None, 1.0, dt, t_end,
                                     g["Lx"], g["Ly"], snap=snap)
    return res, box


def evolving(args):
    pos, F0, h, events, delta, dt, t_end, Lx, Ly, tag = args
    res, checks, ticks = R.run_ticks(pos, F0, h, {**P17, "b": 28.0}, R.THETA_28, NU_FAST, events, delta, dt, t_end,
                                     Lx, Ly, seed_audit=1)
    os.makedirs(P0, exist_ok=True)
    R.save_checks(os.path.join(P0, f"check_{tag}.npz"), checks, ticks)
    lists = [[dict(centre=np.asarray(p["centre"]).tolist(), depth=p["depth"]) for p in c["parts"]] for c in checks]
    return dict(tag=tag, res=res, lists=lists, n_ticks=len(ticks["t"]),
                bytes=os.path.getsize(os.path.join(P0, f"check_{tag}.npz")))


def compare(a, b, Lx, Ly):
    """Proposal §6 comparison rule for the paired replay checks."""
    bad = []
    if a["res"]["status"] != R.COMPLETE or b["res"]["status"] != R.COMPLETE:
        return dict(ok=False, why=f"status {a['res']['status']} / {b['res']['status']}")
    for i, (la, lb) in enumerate(zip(a["lists"], b["lists"])):
        if len(la) != len(lb):
            bad.append(dict(check=i, why="count", na=len(la), nb=len(lb)))
        elif len(la) == 1:
            d = float(np.hypot(*gr.min_image(np.subtract(la[0]["centre"], lb[0]["centre"]), Lx, Ly)))
            dd = abs(la[0]["depth"] - lb[0]["depth"]) / abs(la[0]["depth"])
            if d > CENTRE_TOL or dd > DEPTH_TOL:
                bad.append(dict(check=i, why="centre/depth", d=d, depth_rel=dd))
    md = [float(np.hypot(*gr.min_image(np.subtract(la[0]["centre"], lb[0]["centre"]), Lx, Ly)))
          for la, lb in zip(a["lists"], b["lists"]) if len(la) == 1 and len(lb) == 1]
    return dict(ok=not bad and len(a["lists"]) == len(b["lists"]), failures=bad[:20],
                max_centre_diff=max(md, default=None))


def main():
    os.makedirs(P0, exist_ok=True)
    go.continuum2d()
    t0 = time.time()
    g, h = graph()
    Lx, Ly = g["Lx"], g["Ly"]
    dec = dict(design="GATES.md MG5 Part 0; unregistered graph seed_for(5, 9, 0, 0, 0)", h=h, N=g["N"],
               platform=R.platform_info(), code_revision=R.code_revision())
    print(f"Part 0 graph: N = {g['N']}, h = {h:.6f}", flush=True)
    F_ic = go.soliton_ic(g, (Lx / 2, Ly / 2)).reshape(3, -1)
    P = {**P17, "b": 28.0}
    dt = 0.05

    # 1. solver check
    for attempt in range(2):
        tb = time.time()
        sol, _ = go.integrate(g, h, F_ic.ravel(), T_SOLVER, P=P)
        t_bdf = time.time() - tb
        res, box = frozen(g, h, F_ic, T_SOLVER, dt)
        diff = float(np.abs(box["F"].ravel() - sol.y[:, -1]).max())
        dec.setdefault("solver_check", []).append(dict(dt=dt, max_diff=diff, ok=diff <= SOLVER_TOL, bdf_wall_s=t_bdf,
                                                       ros2_wall_s=res["wall_s"]))
        print(f"solver check dt = {dt}: max|dphi| = {diff:.3e}", flush=True)
        if diff <= SOLVER_TOL:
            break
        if attempt == 0:
            dt /= 2
        else:
            dec.update(proceed=False, stop="solver check failed twice")
            return finish(dec)

    # settle (frozen ROS2) to t = 300
    res, box = frozen(g, h, F_ic, T_SETTLE, dt)
    dec["settle"] = dict(wall_s=res["wall_s"], status=res["status"], s_per_unit_frozen=res["wall_s"] / T_SETTLE)
    F_set, pos_set = box["F"], box["pos"]
    print(f"settled: {res['status']}, {res['wall_s']:.0f} s", flush=True)

    events = ev.gen_events(np.arange(len(pos_set)), R.N0, NU_FAST, T_CHECK, Lx, Ly,
                           seed_for(ev.GATE, B_IDX, 3, 1, 0))
    dec["check_events"] = dict(n=int(len(events)), sha256=ev.events_hash(events))
    delta = 0.2

    def pair(spec_a, spec_b):
        with ProcessPoolExecutor(min(2, R.WORKERS)) as ex:
            ra, rb = ex.map(evolving, [spec_a, spec_b])
        return ra, rb

    # 2. step check
    step_ok = False
    for attempt in range(2):
        ra, rb = pair((pos_set, F_set, h, events, delta, dt, T_CHECK, Lx, Ly, f"step_dt{dt:g}"),
                      (pos_set, F_set, h, events, delta, dt / 2, T_CHECK, Lx, Ly, f"step_dt{dt / 2:g}"))
        c = compare(ra, rb, Lx, Ly)
        c.update(dt=dt, delta=delta, wall_s=[ra["res"]["wall_s"], rb["res"]["wall_s"]], counts=[len(x) for x in ra["lists"]])
        dec.setdefault("step_check", []).append(c)
        print(f"step check dt {dt} vs {dt / 2}: {'ok' if c['ok'] else 'FAIL'}", flush=True)
        fast = ra
        if c["ok"]:
            step_ok = True
            break
        if attempt == 0:
            dt /= 2
    if not step_ok:
        dec.update(proceed=False, stop="step check failed twice")
        return finish(dec)

    # 3. tick check
    tick_ok = False
    for attempt in range(2):
        ra, rb = pair((pos_set, F_set, h, events, delta, dt, T_CHECK, Lx, Ly, f"tick_D{delta:g}"),
                      (pos_set, F_set, h, events, delta / 2, dt, T_CHECK, Lx, Ly, f"tick_D{delta / 2:g}"))
        c = compare(ra, rb, Lx, Ly)
        c.update(dt=dt, delta=delta, wall_s=[ra["res"]["wall_s"], rb["res"]["wall_s"]])
        dec.setdefault("tick_check", []).append(c)
        print(f"tick check Delta {delta} vs {delta / 2}: {'ok' if c['ok'] else 'FAIL'}", flush=True)
        if c["ok"]:
            tick_ok = True
            fast = ra
            break
        if attempt == 0:
            delta /= 2
    if not tick_ok:
        dec.update(proceed=False, stop="tick check failed twice")
        return finish(dec)

    # timing and storage, from the complete implementation (the accepted fast run of the tick check)
    s_tick_fast = fast["res"]["wall_s"] / (T_CHECK / delta)
    s_unit_frozen = dec["settle"]["s_per_unit_frozen"]
    s_unit_slow = s_unit_frozen + (fast["res"]["wall_s"] / T_CHECK - s_unit_frozen) * delta   # Delta = 1 at slow rates
    proj = {"A frozen (4 x 500)": 4 * 500 * s_unit_frozen, "A nu 0.002, 0.02 (8 x 500)": 8 * 500 * s_unit_slow,
            "A nu 0.2 (4 x 500)": 4 * 500 * fast["res"]["wall_s"] / T_CHECK,
            "B frozen (4 x 300)": 4 * 300 * s_unit_frozen, "B nu 0.2 (4 x 300)": 4 * 300 * fast["res"]["wall_s"] / T_CHECK}
    dec["timing"] = dict(s_per_tick_fast=s_tick_fast, s_per_unit_fast=fast["res"]["wall_s"] / T_CHECK,
                         s_per_unit_slow_est=s_unit_slow, s_per_unit_frozen=s_unit_frozen,
                         projected_core_h={k: v / 3600 for k, v in proj.items()},
                         projected_total_core_h=sum(proj.values()) / 3600)
    units = 16 * 500 + 8 * 300                   # simulated time units of all 24 registered runs
    dec["storage"] = dict(check_file_bytes_100_units_fast=fast["bytes"],
                          upper_estimate_committed_MB=fast["bytes"] / T_CHECK * units / 1e6,
                          note="fast-rate bytes per unit times all simulated units: an upper bound (frozen and slow "
                               "runs have fewer ticks)")
    dec.update(proceed=True, dt=dt, delta_fast=delta, part0_wall_s=time.time() - t0)
    return finish(dec)


def finish(dec):
    json.dump(dec, open(os.path.join(P0, "decisions.json"), "w"), indent=1, default=float)
    print("Part 0:", "PROCEED" if dec.get("proceed") else f"STOP ({dec.get('stop')})", flush=True)
    if dec.get("timing"):
        print(f"projected registered cost: {dec['timing']['projected_total_core_h']:.1f} core-h", flush=True)
    return dec


if __name__ == "__main__":
    main()

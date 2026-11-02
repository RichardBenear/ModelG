"""MG6 Part 0 (proposal §7; GATES.md MG6 pre-run specification): checks on an UNREGISTERED Delaunay graph before any
registered run, plus the full-pipeline timing with 8 concurrent runs (MG020a) used to restate the budget.

    $env:MODELG_WORKERS=8; python src/run_mg6.py part0

Checks (formulas in proposal §7):
1. acceptance rate and determinism (static synthetic field, kappa = 4, 20 units);
2. density equilibration (same field, 50 units, immigration-death mean and covariance);
3. step check: dt against dt/2 with the accepted edits of the dt run imposed (kappa = 4, 100 units);
4. tick check: Delta against Delta/2, fully coupled (kappa = 4, 100 units), parts (a) and (b);
5. bound margin: max nodal phi_G <= 0.7 * phi_lim in every coupled run;
6. transient: blocks B1 = {60, 70}, B2 = {90, 100} of D0(1.5) in the accepted configuration's run;
7. quadrature checks (i)-(iv);
8. timing: 8 concurrent coupled runs of 20 units (the production worker count, MG020a).
Output: results/MG6/part0/decisions.json and check records (committed).
"""
import os, sys, json, time
import numpy as np
from concurrent.futures import ProcessPoolExecutor
sys.path.insert(0, os.path.dirname(__file__))
import graphs as gr
import graph_ode as go
import evolve as ev
import mg6
import run_mg5 as R5
import run_mg6 as R
from lattice import P17
from runner import seed_for

P0 = os.path.join(R.OUT, "part0")
B_IDX = 9
KAP = 4.0
T_ACC, T_EQ, T_CHECK, T_SETTLE, T_TIME = 20.0, 50.0, 100.0, 300.0, 20.0
N_TIMING = 8                                # concurrent timing runs = production workers (MG020a)
CENTRE_STEP, DEPTH_STEP, CENTRE_TICK, DEPTH_TICK = 0.1, 0.02, 0.15, 0.03
D0_TOL, P0_MIN, MARGIN, DRIFT_TOL = 0.25, 10.0, 0.7, 0.25
EQ_EDGES = (0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0)


def graph():
    tor = gr.torus(0.2, 122, 142)
    rng = np.random.default_rng(seed_for(mg6.GATE, B_IDX, 0, 0, 0))
    pos = gr.poisson_points(R.N0, tor["Lx"], tor["Ly"], rng)
    g = ev.Mesh(pos, tor["Lx"], tor["Ly"]).g
    h = float(np.sqrt(gr.corrector(g)["h2"]))
    return g, h


def static_field(Lx, Ly):
    c = np.array([Lx / 2, Ly / 2])

    def G(pos):
        d = gr.min_image(np.asarray(pos) - c, Lx, Ly)
        return -0.25 * np.exp(-(d ** 2).sum(1) / 2.0)
    return dict(G=G, centre=c)


def settle(g, h):
    box = {}

    def snap(t, pos, ids, F):
        box.update(t=t, pos=pos.copy(), F=F.copy())
    F_ic = go.soliton_ic(g, (g["Lx"] / 2, g["Ly"] / 2)).reshape(3, -1)
    res, _, _ = R5.run_ticks(g["pos"], F_ic, h, {**P17, "b": 28.0}, R.THETA, 0.0, None, 1.0, 0.05, T_SETTLE,
                             g["Lx"], g["Ly"], snap=snap)
    return res, box


def coupled_run(args):
    """One Part 0 coupled run from the settled state. args: dict(tag, pos, F, h, delta, dt, t_end, seed, forced_tag,
    cross, keep_final)."""
    a = args
    Lx, Ly = a["Lx"], a["Ly"]
    cand = mg6.gen_candidates(R.N0, R.NU, a["t_end"], Lx, Ly, a["seed"])
    forced = None
    if a.get("forced_tag"):
        with np.load(os.path.join(P0, f"accepted_{a['forced_tag']}.npz")) as z:
            forced = np.unpackbits(z["bits"])[:int(z["n"])].astype(bool)
    co = mg6.Coupled(cand, KAP, None, forced=forced)
    grid = mg6.Grid(Lx, Ly)
    rho0 = R.N0 / grid.A
    box = {}

    def snap(t, pos, ids, F):
        if abs(t - a["t_end"]) < 1e-9:
            box.update(pos=pos.copy(), F=F.copy())
    res, checks, ticks = R.run_ticks6(a["pos"], a["F"], a["h"], {**P17, "b": 28.0}, R.THETA, co, a["delta"], a["dt"],
                                      a["t_end"], Lx, Ly, KAP, None, grid, rho0, seed_audit=1, snap=snap,
                                      cross=a.get("cross", True))
    os.makedirs(P0, exist_ok=True)
    R.save_checks6(os.path.join(P0, f"check_{a['tag']}.npz"), checks, ticks, co.accepted)
    np.savez_compressed(os.path.join(P0, f"accepted_{a['tag']}.npz"), bits=np.packbits(co.accepted),
                        n=len(co.accepted))
    out = dict(tag=a["tag"], res=res, n_ticks=len(ticks["t"]),
               lists=[[dict(centre=np.asarray(p["centre"]).tolist(), depth=p["depth"]) for p in c["parts"]]
                      for c in checks],
               t=[c["t"] for c in checks],
               D0_15=[float(c["m"]["D0"][1]) if "m" in c else None for c in checks],
               P0_15=[float(c["m"]["P0"][1]) if "m" in c else None for c in checks],
               cand_sha256=mg6.candidates_hash(cand))
    if a.get("keep_final") and box:
        np.savez(os.path.join(P0, f"final_{a['tag']}.npz"), **box)
    return out


def compare_tracks(a, b, Lx, Ly, ctol, dtol):
    if a["res"]["status"] != R.COMPLETE or b["res"]["status"] != R.COMPLETE:
        return dict(ok=False, why=f"status {a['res']['status']} / {b['res']['status']}")
    bad, md, mdd = [], [], []
    for i, (la, lb) in enumerate(zip(a["lists"], b["lists"])):
        if len(la) != len(lb):
            bad.append(dict(check=i, why="count", na=len(la), nb=len(lb)))
        elif len(la) == 1:
            d = float(np.hypot(*gr.min_image(np.subtract(la[0]["centre"], lb[0]["centre"]), Lx, Ly)))
            dd = abs(la[0]["depth"] - lb[0]["depth"]) / abs(la[0]["depth"])
            md.append(d); mdd.append(dd)
            if d > ctol or dd > dtol:
                bad.append(dict(check=i, why="centre/depth", d=d, depth_rel=dd))
    return dict(ok=not bad and len(a["lists"]) == len(b["lists"]), failures=bad[:20],
                max_centre_diff=max(md, default=None), max_depth_rel=max(mdd, default=None))


def window_mean(run, key, lo=50.0, hi=100.0):
    v = [x for t, x in zip(run["t"], run[key]) if lo - 1e-9 <= t <= hi + 1e-9]
    return float(np.mean(v)) if v and None not in v else None


def check_acceptance(g, Lx, Ly):
    """Check 1: |A - sum p| <= 4 sqrt(sum p(1-p)) over evaluated (non-ephemeral) candidates; determinism."""
    ff = static_field(Lx, Ly)
    out = []
    for rep in range(2):
        cand = mg6.gen_candidates(R.N0, R.NU, T_ACC, Lx, Ly, seed_for(mg6.GATE, B_IDX, 3, 1, 0))
        co = mg6.Coupled(cand, KAP, None, keep_p=True)
        grid = mg6.Grid(Lx, Ly)
        res, _, _ = R.run_ticks6(g["pos"], None, 0.0, None, R.THETA, co, 0.2, 0.05, T_ACC, Lx, Ly, KAP, None, grid,
                                 R.N0 / grid.A, measure_every=1e9, cross=False, fixed_field=ff)
        J = co.evaluated
        p = co.p[J]
        A = int(co.accepted.sum())
        mean, var = float(p.sum()), float((p * (1 - p)).sum())
        out.append(dict(status=res["status"], n_candidates=int(len(J)), n_evaluated=int(J.sum()),
                        n_ephemeral_frac=float(1 - J.mean()), accepted=A, expected=mean, sd=float(np.sqrt(var)),
                        z=float((A - mean) / np.sqrt(var)), bits_sha=__import__("hashlib").sha256(
                            np.packbits(co.accepted).tobytes()).hexdigest()))
    ok = (all(o["status"] == R.COMPLETE and abs(o["z"]) <= 4 for o in out) and out[0]["bits_sha"] == out[1]["bits_sha"])
    return dict(ok=ok, runs=out)


def check_equilibration(g, Lx, Ly):
    """Check 2: annular counts at t = 25..50 against the immigration-death mean and covariance."""
    ff = static_field(Lx, Ly)
    cand = mg6.gen_candidates(R.N0, R.NU, T_EQ, Lx, Ly, seed_for(mg6.GATE, B_IDX, 3, 1, 2))
    co = mg6.Coupled(cand, KAP, None)
    grid = mg6.Grid(Lx, Ly)
    rho0 = R.N0 / grid.A
    res, checks, _ = R.run_ticks6(g["pos"], None, 0.0, None, R.THETA, co, 0.2, 0.05, T_EQ, Lx, Ly, KAP, None, grid,
                                  rho0, measure_every=1.0, cross=False, fixed_field=ff)
    ts = [c["t"] for c in checks if 25 - 1e-9 <= c["t"] <= 50 + 1e-9]
    cnt = np.array([c["m"]["ann_count"][:len(EQ_EDGES) - 1] for c in checks if 25 - 1e-9 <= c["t"] <= 50 + 1e-9])
    gq = ff["G"](grid.pts)
    f = np.exp(KAP * gq)
    W = {R_: grid.disc_weights(ff["centre"], R_) for R_ in EQ_EDGES[1:]}
    W[0.0] = np.zeros(len(gq))
    out = []
    ts = np.array(ts)
    nu = R.NU
    for j, (lo, hi) in enumerate(zip(EQ_EDGES[:-1], EQ_EDGES[1:])):
        w = W[hi] - W[lo]
        a_f, a_1 = float((rho0 * f * w).sum()), float((rho0 * w).sum())
        m = a_f + (a_1 - a_f) * np.exp(-nu * ts)
        K = len(ts)
        S, T = np.meshgrid(ts, ts, indexing="ij")
        mm = a_f + (a_1 - a_f) * np.exp(-nu * np.minimum(S, T))
        var = float((mm * np.exp(-nu * np.abs(S - T))).sum() / K ** 2)
        obs, exp_ = float(cnt[:, j].mean()), float(m.mean())
        out.append(dict(annulus=[lo, hi], observed=obs, expected=exp_, sd=float(np.sqrt(var)),
                        z=float((obs - exp_) / np.sqrt(var))))
    return dict(ok=res["status"] == R.COMPLETE and all(abs(o["z"]) <= 4 for o in out), samples=len(ts),
                annuli=out, status=res["status"])


def check_quadrature(Lx, Ly, final):
    """Check 7: (i) constant integrates to A; (ii) disc areas; (iii) DEC zero; (iv) refinement on the coupled field."""
    grid = mg6.Grid(Lx, Ly)
    rng = np.random.default_rng(0)
    i_ok = abs(grid.dA * grid.nx * grid.ny - grid.A) <= 1e-12 * grid.A
    worst = {}
    for _ in range(20):
        c = rng.random(2) * (Lx, Ly)
        for Rr in (0.5,) + mg6.RADII:
            e = abs(grid.disc_weights(c, Rr).sum() - np.pi * Rr * Rr) / (np.pi * Rr * Rr)
            worst[Rr] = max(worst.get(Rr, 0.0), float(e))
    ii_ok = all(v <= (0.005 if r_ == 0.5 else 0.001) for r_, v in worst.items())
    mesh = ev.Mesh(final["pos"], Lx, Ly)
    F = final["F"]
    s = ev.smoothed_X(mesh.g, F)
    parts = ev.detect_full(mesh.pos, s, R.THETA, Lx, Ly)
    c = parts[0]["centre"] if parts else np.array([Lx / 2, Ly / 2])
    rho0 = R.N0 / grid.A
    md = mg6.measure(mesh, F, c, 0.0, None, grid, rho0)
    iii_ok = bool(np.all(md["P0"] == 0) and np.all(md["Pc"] == 0))
    m1 = mg6.measure(mesh, F, c, KAP, None, grid, rho0)
    g2 = mg6.Grid(Lx, Ly, 2 * mg6.NX, 2 * mg6.NY, 16)
    m2 = mg6.measure(mesh, F, c, KAP, None, g2, rho0)
    rel = dict(P0_15=abs(m2["P0"][1] - m1["P0"][1]) / abs(m1["P0"][1]),
               M_minus=abs(m2["M_minus"] - m1["M_minus"]) / abs(m1["M_minus"]))
    iv_ok = all(v <= 0.01 for v in rel.values())
    return dict(ok=bool(i_ok and ii_ok and iii_ok), refine_ok=bool(iv_ok), i=bool(i_ok),
                ii_worst_rel={str(k): v for k, v in worst.items()}, iii=iii_ok, iv_rel=rel,
                P0_15=[m1["P0"][1], m2["P0"][1]], M_minus=[m1["M_minus"], m2["M_minus"]])


def finish(dec):
    os.makedirs(P0, exist_ok=True)
    dec["part0_wall_s"] = time.time() - dec.pop("_t0", time.time())
    json.dump(dec, open(os.path.join(P0, "decisions.json"), "w"), indent=1, default=float)
    print("Part 0:", "PROCEED" if dec.get("proceed") else f"STOP ({dec.get('stop')})", flush=True)
    return dec


def main():
    os.makedirs(P0, exist_ok=True)
    go.continuum2d()
    dec = dict(design="GATES.md MG6 Part 0; unregistered graph seed_for(6, 9, 0, 0, 0)", _t0=time.time(),
               platform=R5.platform_info(), code_revision=R5.code_revision(), workers=R.WORKERS)
    g, h = graph()
    Lx, Ly = g["Lx"], g["Ly"]
    dec.update(h=h, N=g["N"], Lx=Lx, Ly=Ly)
    print(f"Part 0 graph: N = {g['N']}, h = {h:.6f}", flush=True)

    dec["acceptance"] = check_acceptance(g, Lx, Ly)
    print("check 1 acceptance:", "ok" if dec["acceptance"]["ok"] else "FAIL", flush=True)
    if not dec["acceptance"]["ok"]:
        dec.update(proceed=False, stop="check 1 (acceptance) failed")
        return finish(dec)
    dec["equilibration"] = check_equilibration(g, Lx, Ly)
    print("check 2 equilibration:", "ok" if dec["equilibration"]["ok"] else "FAIL", flush=True)
    if not dec["equilibration"]["ok"]:
        dec.update(proceed=False, stop="check 2 (density equilibration) failed")
        return finish(dec)

    res, box = settle(g, h)
    dec["settle"] = dict(status=res["status"], wall_s=res["wall_s"])
    print(f"settled: {res['status']}, {res['wall_s']:.0f} s", flush=True)
    base = dict(pos=box["pos"], F=box["F"], h=h, Lx=Lx, Ly=Ly, t_end=T_CHECK,
                seed=seed_for(mg6.GATE, B_IDX, 3, 1, 1))
    dt, delta = 0.05, 0.2
    runs, phis = {}, []
    with ProcessPoolExecutor(min(3, max(1, R.WORKERS))) as ex:
        for attempt in range(2):
            # 3 and 4 share the (delta, dt) run; the dt/2 replay needs its accepted edits, the delta/2 run does not
            tags = dict(A=f"d{delta:g}_dt{dt:g}", C=f"d{delta / 2:g}_dt{dt:g}", B=f"d{delta:g}_dt{dt / 2:g}_forced")
            futA = ex.submit(coupled_run, dict(base, tag=tags["A"], delta=delta, dt=dt, keep_final=True))
            futC = ex.submit(coupled_run, dict(base, tag=tags["C"], delta=delta / 2, dt=dt))
            runs["A"] = futA.result()
            runs["B"] = ex.submit(coupled_run, dict(base, tag=tags["B"], delta=delta, dt=dt / 2,
                                                    forced_tag=tags["A"])).result()
            runs["C"] = futC.result()
            phis += [r["res"]["max_phiG"] for r in runs.values()]
            step = compare_tracks(runs["A"], runs["B"], Lx, Ly, CENTRE_STEP, DEPTH_STEP)
            dec.setdefault("step_check", []).append(dict(dt=[dt, dt / 2], delta=delta, **step))
            tick_a = compare_tracks(runs["A"], runs["C"], Lx, Ly, CENTRE_TICK, DEPTH_TICK)
            D_a, D_c = window_mean(runs["A"], "D0_15"), window_mean(runs["C"], "D0_15")
            Pbar = window_mean(runs["C"], "P0_15")
            if Pbar is None or abs(Pbar) < P0_MIN:
                dec["tick_check"] = [dict(part_a=tick_a, D0_bar=[D_a, D_c], P0_bar=Pbar, ok=False)]
                dec.update(proceed=False, stop=f"tick check not scorable: |P0_bar| = {Pbar}")
                return finish(dec)
            b_ok = D_a is not None and D_c is not None and abs(D_a - D_c) <= D0_TOL * abs(Pbar)
            dec.setdefault("tick_check", []).append(dict(delta=[delta, delta / 2], dt=dt, part_a=tick_a,
                                                         D0_bar=[D_a, D_c], P0_bar=Pbar, part_b_ok=bool(b_ok),
                                                         ok=bool(tick_a["ok"] and b_ok)))
            print(f"step check dt {dt} vs {dt / 2}: {'ok' if step['ok'] else 'FAIL'}; tick check Delta {delta} vs "
                  f"{delta / 2}: {'ok' if dec['tick_check'][-1]['ok'] else 'FAIL'}", flush=True)
            if step["ok"] and dec["tick_check"][-1]["ok"]:
                break
            if attempt == 1:
                dec.update(proceed=False, stop="step or tick check failed twice")
                return finish(dec)
            if not step["ok"]:
                dt /= 2
            if not dec["tick_check"][-1]["ok"]:
                delta /= 2
        acc_run = runs["A"]                                   # the accepted configuration's run (delta, dt)
    dec["bound_margin"] = dict(max_phiG=max(phis), limit=MARGIN * mg6.phi_lim(KAP),
                               ok=max(phis) <= MARGIN * mg6.phi_lim(KAP))
    if not dec["bound_margin"]["ok"]:
        dec.update(proceed=False, stop="bound margin exceeded")
        return finish(dec)
    ser = dict(zip(acc_run["t"], acc_run["D0_15"]))
    b1, b2 = np.mean([ser[60.0], ser[70.0]]), np.mean([ser[90.0], ser[100.0]])
    Pbar = dec["tick_check"][-1]["P0_bar"]
    s = np.array([v for t, v in zip(acc_run["t"], acc_run["D0_15"]) if t > 0])
    lag1 = float(np.corrcoef(s[:-1], s[1:])[0, 1]) if len(s) > 2 else None
    drift = abs(b2 - b1) > DRIFT_TOL * abs(Pbar)
    dec["transient"] = dict(series=acc_run["D0_15"], t=acc_run["t"], B1=float(b1), B2=float(b2), lag1=lag1,
                            drift=bool(drift), window=[150.0, 250.0] if drift else [100.0, 200.0])
    with np.load(os.path.join(P0, f"final_{acc_run['tag']}.npz")) as z:
        final = {k: z[k] for k in z.files}
    dec["quadrature"] = check_quadrature(Lx, Ly, final)
    if not dec["quadrature"]["ok"]:
        dec.update(proceed=False, stop="quadrature check (i)-(iii) failed: implementation error")
        return finish(dec)
    nx, ny = (mg6.NX, mg6.NY) if dec["quadrature"]["refine_ok"] else (2 * mg6.NX, 2 * mg6.NY)

    # 8. timing: N_TIMING (= 8, the production worker count) concurrent coupled runs of 20 units at the accepted (delta, dt)
    tbase = dict(base, t_end=T_TIME, delta=delta, dt=dt, cross=True)
    with ProcessPoolExecutor(N_TIMING) as ex:
        tr = list(ex.map(coupled_run, [dict(tbase, tag=f"timing{i}", seed=seed_for(mg6.GATE, B_IDX, 3, 1, 10 + i))
                                       for i in range(N_TIMING)]))
    spu = [r["res"]["wall_s"] / T_TIME for r in tr]
    t_end = dec["transient"]["window"][1]
    proj = 20 * t_end * float(np.mean(spu)) / 3600
    dec["timing"] = dict(s_per_unit_concurrent=spu, n_concurrent=N_TIMING, projected_core_h_20_runs=proj, t_end=t_end)
    dec.update(proceed=True, dt=dt, delta=delta, t_end=t_end, window=dec["transient"]["window"], nx=nx, ny=ny,
               cutoff_core_h=R.CUTOFF)
    print(f"projected registered cost: {proj:.1f} core-h", flush=True)
    return finish(dec)


# ---------------------------------------------------------------------------------------------------------------------
# MG021 amendment (PI decision 2026-11-02 21:18 UTC, spec/MG6_amendment_1.md): Part 0 stopped on the tick check
# (decisions.json, kept unchanged). The PI accepted Delta = 0.1, dt = 0.05 with part (b) of the tick check governing.
# This completes checks 5-8 for that configuration and writes results/MG6/part0/decisions_mg021.json, which the
# driver reads in place of decisions.json.
AMEND_DELTA, AMEND_DT, AMEND_CUTOFF = 0.1, 0.05, 80.0
AMEND_RUNS = ("d0.2_dt0.05", "d0.2_dt0.025_forced", "d0.1_dt0.05", "d0.05_dt0.05", "d0.1_dt0.025_forced")
AMEND_ACC = "d0.1_dt0.05"                              # the accepted configuration's run (Delta 0.1, dt 0.05)
DEC_AMEND = os.path.join(P0, "decisions_mg021.json")


def _sha256(path):
    import hashlib
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for blk in iter(lambda: f.read(1 << 20), b""):
            h.update(blk)
    return h.hexdigest()


def amended_checks(old, final_path, manifest_sha=None):
    """Checks 5-7 from the committed Part 0 records (no new dynamics). Pure apart from file reads; tested."""
    tk = old["tick_check"]
    assert old.get("stop") == "step or tick check failed twice" and len(tk) == 2, "unexpected Part 0 record"
    assert tk[1]["delta"] == [AMEND_DELTA, AMEND_DELTA / 2] and tk[1]["part_b_ok"], "part (b) at Delta 0.1 not passed"
    assert all(s["ok"] for s in old["step_check"]) and old["step_check"][1]["delta"] == AMEND_DELTA
    out = {}
    phis = {}
    for tag in AMEND_RUNS:
        with np.load(os.path.join(P0, f"check_{tag}.npz")) as z:
            phis[tag] = float(np.max(z["tick_max_phiG"]))
    lim = MARGIN * mg6.phi_lim(KAP)
    out["bound_margin"] = dict(max_phiG_by_run=phis, max_phiG=max(phis.values()), limit=lim,
                               ok=bool(max(phis.values()) <= lim))
    with np.load(os.path.join(P0, f"check_{AMEND_ACC}.npz")) as z:
        t, D0 = z["t"].tolist(), z["m_D0"][:, 1].tolist()
    ser = dict(zip(t, D0))
    b1, b2 = np.mean([ser[60.0], ser[70.0]]), np.mean([ser[90.0], ser[100.0]])
    Pbar = tk[1]["P0_bar"]                                   # P0_bar of the passing comparison (Delta 0.05 run)
    s = np.array([v for tt, v in zip(t, D0) if tt > 0])
    drift = abs(b2 - b1) > DRIFT_TOL * abs(Pbar)
    out["transient"] = dict(run=AMEND_ACC, series=D0, t=t, B1=float(b1), B2=float(b2), P0_bar=Pbar,
                            lag1=float(np.corrcoef(s[:-1], s[1:])[0, 1]), drift=bool(drift),
                            window=[150.0, 250.0] if drift else [100.0, 200.0])
    sha = _sha256(final_path)
    out["final_state"] = dict(path=os.path.basename(final_path), sha256=sha,
                              matches_manifest=None if manifest_sha is None else sha == manifest_sha)
    return out


def main_amended():
    os.makedirs(P0, exist_ok=True)
    go.continuum2d()
    old = json.load(open(os.path.join(P0, "decisions.json")))
    dec = dict(design="spec/MG6_amendment_1.md (MG021): Part 0 checks 5-8 at Delta 0.1, dt 0.05", _t0=time.time(),
               platform=R5.platform_info(), code_revision=R5.code_revision(), workers=R.WORKERS,
               supersedes_outcome_of="decisions.json (STOP, kept unchanged)")
    final_path = os.path.join(P0, f"final_{AMEND_ACC}.npz")
    man = {}
    mpath = os.path.join(P0, "local_manifest.csv")
    if os.path.exists(mpath):
        import csv
        man = {os.path.basename(r["path"]): r["sha256"] for r in csv.DictReader(open(mpath))}
    dec.update(amended_checks(old, final_path, man.get(os.path.basename(final_path))))
    if dec["final_state"]["matches_manifest"] is False:
        dec.update(proceed=False, stop="final state differs from local_manifest.csv")
        return finish_amended(dec)
    print("check 5 bound margin:", "ok" if dec["bound_margin"]["ok"] else "FAIL", flush=True)
    if not dec["bound_margin"]["ok"]:
        dec.update(proceed=False, stop="bound margin exceeded")
        return finish_amended(dec)
    print(f"check 6 transient: window {dec['transient']['window']}", flush=True)
    g, h = graph()
    Lx, Ly = g["Lx"], g["Ly"]
    assert abs(h - old["h"]) < 1e-12 and g["N"] == old["N"], "Part 0 graph differs from the recorded one"
    with np.load(final_path) as z:
        final = {k: z[k] for k in z.files}
    dec["quadrature"] = check_quadrature(Lx, Ly, final)
    print("check 7 quadrature:", "ok" if dec["quadrature"]["ok"] else "FAIL", flush=True)
    if not dec["quadrature"]["ok"]:
        dec.update(proceed=False, stop="quadrature check (i)-(iii) failed: implementation error")
        return finish_amended(dec)
    nx, ny = (mg6.NX, mg6.NY) if dec["quadrature"]["refine_ok"] else (2 * mg6.NX, 2 * mg6.NY)
    res, box = settle(g, h)                                  # deterministic; the same settled state as the first Part 0
    dec["settle"] = dict(status=res["status"], wall_s=res["wall_s"])
    print(f"settled: {res['status']}, {res['wall_s']:.0f} s", flush=True)
    base = dict(pos=box["pos"], F=box["F"], h=h, Lx=Lx, Ly=Ly)
    tbase = dict(base, t_end=T_TIME, delta=AMEND_DELTA, dt=AMEND_DT, cross=True)
    with ProcessPoolExecutor(N_TIMING) as ex:
        tr = list(ex.map(coupled_run, [dict(tbase, tag=f"timing_mg021_{i}",
                                            seed=seed_for(mg6.GATE, B_IDX, 3, 1, 10 + i)) for i in range(N_TIMING)]))
    spu = [r["res"]["wall_s"] / T_TIME for r in tr]
    t_end = dec["transient"]["window"][1]
    proj = 20 * t_end * float(np.mean(spu)) / 3600
    dec["timing"] = dict(s_per_unit_concurrent=spu, n_concurrent=N_TIMING, projected_core_h_20_runs=proj, t_end=t_end,
                         max_phiG=[r["res"]["max_phiG"] for r in tr])
    dec.update(proceed=True, dt=AMEND_DT, delta=AMEND_DELTA, t_end=t_end, window=dec["transient"]["window"], nx=nx,
               ny=ny, cutoff_core_h=AMEND_CUTOFF)
    print(f"check 8 timing: projected registered cost {proj:.1f} core-h", flush=True)
    return finish_amended(dec)


def finish_amended(dec):
    dec["part0_wall_s"] = time.time() - dec.pop("_t0", time.time())
    json.dump(dec, open(DEC_AMEND, "w"), indent=1, default=float)
    print("Part 0 (MG021):", "PROCEED" if dec.get("proceed") else f"STOP ({dec.get('stop')})", flush=True)
    return dec

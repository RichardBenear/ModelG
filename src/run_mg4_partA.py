"""MG4 Part A driver (proposal §5; GATES.md MG4). Deterministic graph-ODE runs on T, the D and R panels, N3, and the
h = 0.3 check. Resumable at run level: one JSON summary + one compressed field file per run under results/MG4/partA/.
Workers: MODELG_WORKERS (default 2). Refuses to start without results/MG4/part0/decisions.json, or if Part 0 set
stop_for_PI, or if the box check moved Part A to the 30 x 30 torus (that needs the panel rebuilt; PI decision first).
    python src/run_mg4_partA.py            # everything, in the registered order
Order: T (8 runs; fixes the threshold theta), then D and R survival runs, D and R seeded runs, N3, the h = 0.3 check.
"""
import os, sys, json, time, platform, csv
import numpy as np
from concurrent.futures import ProcessPoolExecutor, as_completed
sys.path.insert(0, os.path.dirname(__file__))
import graphs as gr
import graph_ode as go
import measure_graph as mg
import run_mg4_part0 as p0
from runner import seed_for

ROOT = os.path.join(os.path.dirname(__file__), "..", "results", "MG4")
OUT = os.environ.get("MG4A_OUT") or os.path.join(ROOT, "partA")
DEC = os.environ.get("MG4_DEC") or os.path.join(ROOT, "part0", "decisions.json")
WORKERS = int(os.environ.get("MODELG_WORKERS", "2"))
CHECKS = np.arange(200.0, 300.0 + 1e-9, 10.0)
STATIONARY = 1e-8


# ------------------------------------------------------------------ graphs and placements
def get_graph(spec, dec):
    """spec: dict(fam, slot, torus) -> graph and its h. Rebuilt deterministically from the Part 0 records."""
    fam = spec["fam"]
    if spec.get("torus") == "h03":
        tor = gr.torus(p0.H_T3, p0.NX3, p0.NY3)
        if fam == "T":
            g = gr.triangular(tor)
        else:
            n = int(round(dec["c_F"]["D"] * tor["Lx"] * tor["Ly"] / p0.H_T3**2))
            g = p0.build("D", n, tor["Lx"], tor["Ly"], 200)          # h = 0.3 check graph: stream rep 200
        return g, float(np.sqrt(gr.corrector(g)["h2"]))
    tor = gr.torus(p0.H_T, p0.NX, p0.NY)
    if fam == "T":
        return gr.triangular(tor), float(dec["T"]["h"])
    if fam == "N3":
        return gr.random_regular(dec["N3"]["N"], 6, seed_for(p0.GATE, p0.B_IDX, p0.FAM["N3"], 0, spec["slot"])), \
            float(dec["N3"]["h"])
    pan = dec["panel"][f"{fam}{spec['slot']}"]
    return p0.build(fam, pan["n_points"], tor["Lx"], tor["Ly"], pan["rep"]), float(pan["h"])


def t_placements(g):
    """T: the node nearest the torus centre, the midpoint of its +x edge, the centroid of its upper-right triangle,
    one uniform random point (stage 1, replicate 0)."""
    Lx, Ly = g["Lx"], g["Ly"]
    i0 = int(np.argmin(np.hypot(*gr.min_image(g["pos"] - (Lx / 2, Ly / 2), Lx, Ly).T)))
    a = float(np.linalg.norm(g["dx"][0]))
    p = g["pos"][i0]
    up = np.array([0.5 * a, a * np.sqrt(3) / 2])
    rnd = np.random.default_rng(seed_for(p0.GATE, p0.B_IDX, p0.FAM["T"], 1, 0)).uniform((0, 0), (Lx, Ly))
    return [p, p + (a / 2, 0.0), p + (np.array([a, 0.0]) + up) / 3, rnd]


def placement(fam, slot, i, Lx, Ly):
    rng = np.random.default_rng(seed_for(p0.GATE, p0.B_IDX, p0.FAM[fam], 1, slot * 100 + i))
    return rng.uniform((0, 0), (Lx, Ly))


# ------------------------------------------------------------------ one run
def run_one(job):
    """job: dict(key, fam, slot, torus, centre, ic ('soliton'|'seed'), theta (None for T: classified later))."""
    path = os.path.join(OUT, "runs", job["key"] + ".json")
    if os.path.exists(path):
        return json.load(open(path))
    dec = json.load(open(DEC))
    g, h = get_graph(job, dec)
    N = g["N"]
    c = np.asarray(job["centre"], float)
    cont = go.continuum2d()
    t0 = time.time()
    if job["fam"] == "N3":                                   # no embedding: profile by hop distance x h (descriptive)
        from scipy.sparse.csgraph import shortest_path
        d = shortest_path(gr.laplacian(g) != 0, unweighted=True, indices=[int(job["slot_node"])])[0] * h
        y0 = np.concatenate([go.interp_radial(d, cont["r"], cont["phi"][s]) for s in range(3)])
        sol, rhs = go.integrate(g, h, y0, 300.0, t_eval=CHECKS)
    else:
        if job["ic"] == "soliton":
            sol, rhs = go.integrate(g, h, go.soliton_ic(g, c, cont), 300.0, t_eval=CHECKS)
        else:
            sol, rhs = go.integrate(g, h, np.zeros(3 * N), 300.0, t_eval=CHECKS, seed_centre=c, seed_scale=1.0)
    wall = time.time() - t0
    y = sol.y[:, -1]
    dydt = rhs(300.0, y)
    bal = go.balances(y, dydt, N)
    rec = dict(key=job["key"], fam=job["fam"], slot=job.get("slot"), torus=job.get("torus", "main"), ic=job["ic"],
               centre=c.tolist(), N=int(N), h=h, wall_s=wall, nfev=int(sol.nfev), max_dphidt=float(np.abs(dydt).max()),
               stationary=bool(np.abs(dydt).max() < STATIONARY),
               balances=dict(ok=bal["ok"], r1=bal["r1"], M1=bal["M1"], r2=bal["r2"], M2=bal["M2"],
                             I=[float(v) for v in bal["I"]], dI=[float(v) for v in bal["dI"]]),
               platform=dict(machine=platform.node(), os=platform.platform(), python=platform.python_version(),
                             numpy=np.__version__))
    if job["fam"] != "N3":
        S = mg.smoother(g)
        sX = np.array([S @ sol.y[N:2 * N, k] for k in range(len(CHECKS))])
        rec["min_smoothed_final"] = float(sX[-1].min())
        np.savez_compressed(os.path.join(OUT, "fields", job["key"] + ".npz"),
                            final=y.reshape(3, N).astype(np.float32), checks_sX=sX.astype(np.float32),
                            checks_X=sol.y[N:2 * N, :].T.astype(np.float32))
    else:
        np.savez_compressed(os.path.join(OUT, "fields", job["key"] + ".npz"), final=y.reshape(3, N).astype(np.float32))
        rec["min_phiX_final"] = float(y[N:2 * N].min())
    json.dump(rec, open(path, "w"), indent=1)
    ledger(job["key"], wall)
    return rec


def ledger(key, wall):
    path = os.path.join(OUT, "budget_ledger.csv")
    new = not os.path.exists(path)
    with open(path, "a", newline="") as f:
        w = csv.writer(f)
        if new:
            w.writerow(["item", "key", "wall_s", "workers", "machine", "os", "python", "numpy", "utc_time"])
        w.writerow(["partA", key, f"{wall:.2f}", WORKERS, platform.node(), platform.platform(),
                    platform.python_version(), np.__version__, time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())])


# ------------------------------------------------------------------ classification (after the threshold is known)
def analyse(rec, theta, dec):
    """Outcome class (§5.2), drift, radial profiles, depths, lambda. Uses the saved fields; no integration."""
    if rec["fam"] == "N3":
        return dict(cls="n/a")
    g, h = get_graph(dict(fam=rec["fam"], slot=rec["slot"], torus=rec["torus"]), dec)
    N = g["N"]
    f = np.load(os.path.join(OUT, "fields", rec["key"] + ".npz"))
    tree = mg._tree(g)
    checks = [mg.particles(g, f["checks_sX"][k].astype(float), theta, tree) for k in range(len(CHECKS))]
    cls, track = mg.classify(g, checks)
    out = dict(cls=cls, n_particles=[len(c) for c in checks])
    if cls == "single":
        cen = track[-1]
        out["centre_final"] = cen.tolist()
        out["drift"] = float(np.hypot(*gr.min_image(cen - np.asarray(rec["centre"]), g["Lx"], g["Ly"])))
        fin = f["final"].astype(float)
        S = mg.smoother(g)
        prof = {}
        for s, name in enumerate("GXY"):
            rc, p, cnt = mg.radial_profile(g, fin[s], cen)
            prof[name] = p
        rcs, ps, _ = mg.radial_profile(g, S @ fin[1], cen)
        out.update(rc=rc.tolist(), prof={k: v.tolist() for k, v in prof.items()}, counts=cnt.tolist(),
                   empty_bins_r6=int(((cnt == 0) & (rc < 6)).sum()),
                   depth_u=float(prof["X"][0]), depth_s=float(ps[0]), lam=mg.wavelength(rc, prof["X"]))
    return out


def c3b(fam, slot, centre, dec, cont):
    """C3(b) for one placement: e_S for S = G, X, Y, smoothed (criterion) and unsmoothed (reported)."""
    g, h = get_graph(dict(fam=fam, slot=slot), dec)
    L = gr.laplacian(g) / h**2
    S = mg.smoother(g)
    rn, _ = go.node_r(g, centre)
    m = rn < 6.0
    out = {}
    for s, name in enumerate("GXY"):
        fv = go.interp_radial(rn, cont["r"], cont["phi"][s])
        lc = go.interp_radial(rn, cont["r"], cont["lap"][s])
        Lf = L @ fv
        out[name] = float(np.linalg.norm((S @ Lf - S @ lc)[m]) / np.linalg.norm((S @ lc)[m]))
        out[name + "_unsmoothed"] = float(np.linalg.norm((Lf - lc)[m]) / np.linalg.norm(lc[m]))
    return out


# ------------------------------------------------------------------ main
def jobs_for(dec):
    tor = gr.torus(p0.H_T, p0.NX, p0.NY)
    Lx, Ly = tor["Lx"], tor["Ly"]
    gT = gr.triangular(tor)
    J = {"T": [], "main": [], "seed": [], "N3": [], "h03": []}
    for i, c in enumerate(t_placements(gT)):
        J["T"].append(dict(key=f"T_p{i}_sol", fam="T", slot=0, centre=list(c), ic="soliton"))
        J["seed"].append(dict(key=f"T_p{i}_seed", fam="T", slot=0, centre=list(c), ic="seed"))
    for fam in ("D", "R"):
        for slot in range(4):
            for i in range(8):
                c = placement(fam, slot, i, Lx, Ly)
                J["main"].append(dict(key=f"{fam}{slot}_p{i}_sol", fam=fam, slot=slot, centre=list(c), ic="soliton"))
                J["seed"].append(dict(key=f"{fam}{slot}_p{i}_seed", fam=fam, slot=slot, centre=list(c), ic="seed"))
    for slot in (0, 1):
        for i in range(2):
            J["N3"].append(dict(key=f"N3{slot}_n{i}", fam="N3", slot=slot, slot_node=i, centre=[0, 0], ic="soliton"))
    t3 = gr.torus(p0.H_T3, p0.NX3, p0.NY3)
    for i, c in enumerate(t_placements(gr.triangular(t3))):
        J["h03"].append(dict(key=f"h03_T_p{i}_sol", fam="T", slot=0, torus="h03", centre=list(c), ic="soliton"))
    for i in range(8):
        c = np.random.default_rng(seed_for(p0.GATE, p0.B_IDX, p0.FAM["D"], 1, 200 * 100 + i)).uniform(
            (0, 0), (t3["Lx"], t3["Ly"]))
        J["h03"].append(dict(key=f"h03_D_p{i}_sol", fam="D", slot=0, torus="h03", centre=list(c), ic="soliton"))
    return J


# Core-hour ceilings (GATES.md MG4 budget; stop triggers checked before every job is started). The ledger's summed
# per-run wall time is the measure. Re-stated values must be recorded in GATES.md before they are changed here.
CEILING = {"A": 80.0, "h03": 3.0}     # A re-stated 25 -> 80 by the PI, 2026-10-31 (GATES.md MG4)
# Which run groups are run, in order (the T runs always come first; they set the threshold). Fixed in GATES.md.
PARTS = ("main", "h03")      # PI decision 2026-10-31: seeded runs deferred, N3 dynamics dropped


def spent(group):
    path = os.path.join(OUT, "budget_ledger.csv")
    if not os.path.exists(path):
        return 0.0
    rows = list(csv.DictReader(open(path)))
    return sum(float(r["wall_s"]) for r in rows if r["key"].startswith("h03_") == (group == "h03")) / 3600


class BudgetStop(Exception):
    pass


def run_batch(jobs, group="A"):
    """Runs the jobs not yet recorded, keeping at most WORKERS in flight. Before each job is started the ledger is
    checked against the group's ceiling; once it is reached no further job starts (budget stop: INCONCLUSIVE)."""
    out, todo = {}, list(jobs)
    with ProcessPoolExecutor(WORKERS) as ex:
        live = {}
        while todo or live:
            while todo and len(live) < WORKERS:
                used = spent(group)
                if used >= CEILING[group]:
                    todo.clear()
                    print(f"BUDGET STOP: {group} ledger {used:.2f} core-h >= ceiling {CEILING[group]}", flush=True)
                    break
                j = todo.pop(0)
                live[ex.submit(run_one, j)] = j["key"]
            if not live:
                break
            done = next(as_completed(live))
            live.pop(done)
            r = done.result()
            out[r["key"]] = r
            print(f"{r['key']}: wall {r['wall_s']:.0f}s, max|dphi/dt| {r['max_dphidt']:.2e}", flush=True)
    if len(out) < len(jobs) and spent(group) >= CEILING[group]:
        raise BudgetStop(group)
    return out


def main():
    if not os.path.exists(DEC):
        sys.exit("decisions.json missing: run Part 0 first")
    dec = json.load(open(DEC))
    if dec.get("stop_for_PI"):
        sys.exit("Part 0 set stop_for_PI: MG4 stops for a PI decision (GATES.md MG4, proposal §3)")
    if dec.get("torus_for_partA") != "20x20":
        sys.exit("box check failed: Part A would need the 30 x 30 torus and a rebuilt panel (PI decision first)")
    for d in ("runs", "fields"):
        os.makedirs(os.path.join(OUT, d), exist_ok=True)
    go.continuum2d()
    J = jobs_for(dec)
    rT = run_batch(J["T"])
    theta = -0.5 * abs(float(np.mean([rT[j["key"]]["min_smoothed_final"] for j in J["T"]])))
    json.dump(dict(theta=theta, A_s_T=-2 * abs(theta), from_runs=[j["key"] for j in J["T"]]),
              open(os.path.join(OUT, "threshold.json"), "w"), indent=1)
    print("theta =", theta, flush=True)
    J["seed_TD"] = [j for j in J["seed"] if j["fam"] in ("T", "D")]
    J["seed_R"] = [j for j in J["seed"] if j["fam"] == "R"]
    for part in PARTS:
        run_batch(J[part], "h03" if part == "h03" else "A")
    print("Part A runs complete. Evaluate with: python src/evaluate_mg4.py", flush=True)


if __name__ == "__main__":
    try:
        main()
    except BudgetStop as e:
        sys.exit(f"Budget stop ({e}): MG4 is INCONCLUSIVE under proposal §7 item 2; all completed data are reported.")

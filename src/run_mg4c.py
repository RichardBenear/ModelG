"""MG4c driver (spec/MG4c_proposal.md draft 3, frozen at a946550; GATES.md MG4c and its pre-run implementation
specification). Single-particle window in b on unweighted Delaunay panel graphs D0 and D1, plus the cotangent-operator
diagnostic. Resumable at run level: one JSON summary and one field file per run under results/MG4c/.
    $env:MODELG_WORKERS=9; python src/run_mg4c.py          (PowerShell, PI's PC)
Order: Part 1 (40 runs) together with Part 3 (2 runs), then Part 2 rounds (each waits for its runs before the next b is
chosen), then the summary. Refuses to start a job once the ledger reaches the 20 core-h launch cutoff or the run cap.
"""
import os, sys, json, time, platform, csv
import numpy as np
from concurrent.futures import ProcessPoolExecutor, as_completed
sys.path.insert(0, os.path.dirname(__file__))
import graphs as gr
import graph_ode as go
import measure_graph as mg
import run_mg4_partA as pa
from lattice import P17

ROOT = os.path.join(os.path.dirname(__file__), "..", "results")
OUT = os.environ.get("MG4C_OUT") or os.path.join(ROOT, "MG4c")
DEC = os.path.join(ROOT, "MG4", "part0", "decisions.json")
MG4_EVAL = os.path.join(ROOT, "MG4", "evaluation.json")
TSTATE = os.path.join(ROOT, "MG4", "part0", "steps", "Tstate_122x142.npy")
WORKERS = int(os.environ.get("MODELG_WORKERS", "2"))
CHECKS = pa.CHECKS
GRID = [28.75, 28.5, 28.25, 28.0, 27.75, 27.5, 27.25, 27.0]          # Part 1; b = 29 comes from MG4 Part A
THETA_29, THETA_ELIG, ELIG_DEPTH = -5.117354421620525, -1.0, -2.0
STEPS, RUN_CAP, CUTOFF = 6, 102, 20.0
DRUNS = [(0, 0), (0, 1), (1, 0), (1, 1)]                          # (slot, placement): D0 p0, D0 p1, D1 p0, D1 p1
HP_DIGITS, HP_TOL = 40, 1e-25                                     # balance-flag recheck (pre-run spec)


def bkey(b):
    return f"b{b:.6f}"


def job(b, kind, slot=0, pl=0, op=None):
    key = f"{'w_' if op else ''}{bkey(b)}_{kind}{'' if kind == 'T' else slot}_p{pl}"
    return dict(key=key, b=float(b), kind=kind, slot=slot, pl=pl, op=op)


def jobs_at(b):
    return [job(b, "T")] + [job(b, "D", s, p) for s, p in DRUNS]


# ------------------------------------------------------------------ one run
def graph_for(j, dec):
    tor = gr.torus(pa.p0.H_T, pa.p0.NX, pa.p0.NY)
    if j["kind"] == "T":
        g, h = gr.triangular(tor), float(dec["T"]["h"])
        centre = pa.t_placements(g)[0]
    else:
        g, h = pa.get_graph(dict(fam="D", slot=j["slot"]), dec)
        centre = pa.placement("D", j["slot"], j["pl"], tor["Lx"], tor["Ly"])
    return g, h, np.asarray(centre, float)


def run_one(j):
    path = os.path.join(OUT, "runs", j["key"] + ".json")
    if os.path.exists(path):
        return json.load(open(path))
    dec = json.load(open(DEC))
    g, h, c = graph_for(j, dec)
    N = g["N"]
    P = {**P17, "b": j["b"]}
    op, w = None, None
    if j["op"] == "cotan":
        cot = gr.cotan(g)
        op, w = cot["L"], cot["A"]
    t0 = time.time()
    sol, rhs = go.integrate(g, h, go.soliton_ic(g, c), 300.0, t_eval=CHECKS, P=P, op=op)
    wall = time.time() - t0
    y = sol.y[:, -1]
    dydt = rhs(300.0, y)
    bal = go.balances(y, dydt, N, P, w=w)
    S = mg.smoother(g)
    sX = np.array([S @ sol.y[N:2 * N, k] for k in range(len(CHECKS))])
    np.savez_compressed(os.path.join(OUT, "fields", j["key"] + ".npz"), final=y.reshape(3, N),
                        checks_sX=sX.astype(np.float32))
    rec = dict(j, N=int(N), h=h, centre=c.tolist(), wall_s=wall, nfev=int(sol.nfev),
               max_dphidt=float(np.abs(dydt).max()), stationary=bool(np.abs(dydt).max() < pa.STATIONARY),
               min_smoothed_final=float(sX[-1].min()),
               balances=dict(ok=bal["ok"], r1=bal["r1"], M1=bal["M1"], r2=bal["r2"], M2=bal["M2"]),
               platform=dict(machine=platform.node(), os=platform.platform(), python=platform.python_version(),
                             numpy=np.__version__))
    json.dump(rec, open(path, "w"), indent=1)
    ledger(j["key"], wall)
    return rec


def ledger(key, wall):
    path = os.path.join(OUT, "budget_ledger.csv")
    new = not os.path.exists(path)
    with open(path, "a", newline="", encoding="utf-8") as f:
        wr = csv.writer(f)
        if new:
            wr.writerow(["key", "wall_s", "workers", "machine", "os", "python", "numpy", "utc_time"])
        wr.writerow([key, f"{wall:.2f}", WORKERS, platform.node(), platform.platform(), platform.python_version(),
                     np.__version__, time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())])


def spent(keys=False):
    path = os.path.join(OUT, "budget_ledger.csv")
    rows = list(csv.DictReader(open(path, encoding="utf-8"))) if os.path.exists(path) else []
    ch = sum(float(r["wall_s"]) for r in rows) / 3600
    return (ch, {r["key"] for r in rows}) if keys else (ch, len(rows))


def _worker(j):
    return run_one(j)                          # module-level name lookup (tests replace run_one)


class Stop(Exception):
    pass


def run_batch(jobs):
    """Runs the jobs not yet recorded, at most WORKERS in flight. Before each start: launch cutoff and run cap."""
    out, todo = {}, []
    for j in jobs:
        p = os.path.join(OUT, "runs", j["key"] + ".json")
        if os.path.exists(p):
            out[j["key"]] = json.load(open(p))
        else:
            todo.append(j)
    stopped = None
    with ProcessPoolExecutor(WORKERS) as ex:
        live = {}
        while todo or live:
            while todo and len(live) < WORKERS:
                ch, done_keys = spent(keys=True)
                n = len(done_keys | set(live.values()))    # a finished run not yet collected is counted once
                if ch >= CUTOFF or n >= RUN_CAP:
                    stopped = f"ledger {ch:.2f} core-h, {n} runs started"
                    todo.clear()
                    break
                j = todo.pop(0)
                live[ex.submit(_worker, j)] = j["key"]
            if not live:
                break
            done = next(as_completed(live))
            live.pop(done)
            r = done.result()
            out[r["key"]] = r
            print(f"{r['key']}: wall {r['wall_s']:.0f}s, max|dphi/dt| {r['max_dphidt']:.2e}, "
                  f"balance {'ok' if r['balances']['ok'] else 'FLAG'}", flush=True)
    if stopped:
        print("STOP (launch cutoff or run cap):", stopped, flush=True)
        raise Stop(stopped)
    return out


# ------------------------------------------------------------------ Part 3 construction checks (proposal §4)
def cotan_checks(g, fd_size=60):
    """The registered construction checks of the cotangent operator on graph g. Returns dict(name: bool, ok)."""
    Lx, Ly, N = g["Lx"], g["Ly"], g["N"]
    c = gr.cotan(g)
    ge = set(map(tuple, np.column_stack([g["I"], g["J"]])))
    fe, cnt = np.unique(c["edges"], axis=0, return_counts=True)
    K = c["K"]
    f = np.random.default_rng(0).normal(size=N)
    Lf = c["L"] @ f
    out = dict(faces_2N=len(c["faces"]) == 2 * N,
               face_area=bool(abs(c["face_area"].sum() / (Lx * Ly) - 1) < 1e-12),
               edges_in_two_faces=bool(set(map(tuple, fe)) == ge and np.all(cnt == 2)),
               weights_nonneg=bool(c["W"].data.min() >= -1e-12),
               areas=bool(c["A"].min() > 0 and abs(c["A"].sum() / (Lx * Ly) - 1) < 1e-12),
               K_symmetric=bool(abs(K - K.T).max() < 1e-12),
               annihilates_constants=bool(np.abs(c["L"] @ np.ones(N)).max() < 1e-10),
               weighted_sum_zero=bool(abs((c["A"] * Lf).sum()) <= 1e-12 * (c["A"] * np.abs(Lf)).sum()))
    # Jacobian against finite differences, on a small periodic Delaunay graph of the same construction
    gs = gr.delaunay(gr.poisson_points(fd_size, 3.0, 3.0, np.random.default_rng(1)), 3.0, 3.0)
    rhs, jac, _ = go.make_system(gs, None, op=gr.cotan(gs)["L"])
    y = np.random.default_rng(2).normal(0, 0.5, 3 * gs["N"])
    eps = 1e-6
    Jn = np.column_stack([(rhs(0, y + eps * e) - rhs(0, y - eps * e)) / (2 * eps) for e in np.eye(len(y))])
    J = jac(0.0, y).toarray()
    out["jacobian_fd"] = bool(np.abs(J - Jn).max() < 1e-5 * max(1.0, np.abs(J).max()))
    out["ok"] = all(out.values())
    return out


def part3_checks(dec):
    """Checks on D0 and on T (where L_w must equal MG4's L_G with h^2 = 1.5 a^2 entrywise to 1e-12)."""
    gD = graph_for(job(29.0, "D", 0, 0, "cotan"), dec)[0]
    gT = graph_for(job(29.0, "T", 0, 0, "cotan"), dec)[0]
    res = dict(D0=cotan_checks(gD), T=cotan_checks(gT))
    a = float(np.linalg.norm(gT["dx"][0]))
    LT = gr.laplacian(gT) / (1.5 * a * a)
    res["T"]["equals_LG"] = bool(abs(gr.cotan(gT)["L"] - LT).max() < 1e-12 * abs(LT).max())
    res["ok"] = res["D0"]["ok"] and res["T"]["ok"] and res["T"]["equals_LG"]
    return res


# ------------------------------------------------------------------ balance flags (high-precision recheck)
def hp_recheck(rec):
    """Balance residuals of the saved float64 final state in mpmath at HP_DIGITS (fallback: exact fractions). The
    operator is rebuilt from its off-diagonal weights, so its rows sum to zero exactly. Returns dict(r1, r2, M1, M2)."""
    dec = json.load(open(DEC))
    g, h, _ = graph_for(rec, dec)
    N = g["N"]
    y = np.load(os.path.join(OUT, "fields", rec["key"] + ".npz"))["final"].astype(float).ravel()
    try:
        import mpmath
        mpmath.mp.dps = HP_DIGITS
        F, Z = mpmath.mpf, mpmath.mpf(0)
    except ImportError:
        from fractions import Fraction as F
        Z = F(0)
    if rec["op"] == "cotan":
        cot = gr.cotan(g)
        Wc = cot["W"].tocoo()
        I, J, wv = Wc.row, Wc.col, [F(float(v)) for v in Wc.data]
        area = [F(float(a)) for a in cot["A"]]
        scale = [Z + 1 for _ in range(N)]
        node_w = area
    else:
        I = np.concatenate([g["I"], g["J"]]); J = np.concatenate([g["J"], g["I"]])
        wv = [F(1)] * len(I)
        scale = [F(1) / (F(float(h)) * F(float(h)))] * N
        area, node_w = None, [F(1)] * N
    P = {**P17, "b": rec["b"]}
    b, gg, p, q, dxc, dyc = (F(float(P[k])) for k in ("b", "g", "p", "q", "dx", "dy"))
    G0, X0, Y0 = (F(float(v)) for v in go.steady_state(P))
    f = [[F(float(v)) for v in y[k * N:(k + 1) * N]] for k in range(3)]

    def lap(u):                       # sum_i node_w_i (L u)_i, with (L u)_i = scale_i * (1/A_i) * sum_j w_ij (u_j - u_i)
        tot = Z
        for a, bb, ww in zip(I, J, wv):
            tot += (node_w[a] / (area[a] if area else 1)) * scale[a] * ww * (u[bb] - u[a])
        return tot

    sG, sX, sY = (sum((node_w[i] * f[k][i] for i in range(N)), Z) for k in range(3))
    r4 = sum((node_w[i] * ((X0 + f[1][i]) ** 2 * (Y0 + f[2][i]) - X0 ** 2 * Y0) for i in range(N)), Z)
    dG = lap(f[0]) - q * sG + gg * sX
    dX = dxc * lap(f[1]) + p * sG - (1 + b) * sX + r4
    dY = dyc * lap(f[2]) + b * sX - r4
    r1 = dG - (-q * sG + gg * sX)
    r2 = (dG + dX + dY) - (-(q - p) * sG - (1 - gg) * sX)
    M1 = abs(dG) + q * abs(sG) + gg * abs(sX)
    M2 = abs(dG) + abs(dX) + abs(dY) + abs(q - p) * abs(sG) + (1 - gg) * abs(sX)
    return dict(r1=float(r1), r2=float(r2), M1=float(M1), M2=float(M2))


def resolve_flags(recs):
    """Records the recheck of every flagged run in resolution.json before any class is read. Returns unresolved keys."""
    path = os.path.join(OUT, "resolution.json")
    res = json.load(open(path)) if os.path.exists(path) else {}
    for k, r in recs.items():
        if not r["balances"]["ok"] and k not in res:
            hp = hp_recheck(r)
            hp["resolved"] = bool(abs(hp["r1"]) <= HP_TOL * (1 + hp["M1"]) and abs(hp["r2"]) <= HP_TOL * (1 + hp["M2"]))
            res[k] = hp
            json.dump(res, open(path, "w"), indent=1)
            print(f"balance flag {k}: high-precision r1 {hp['r1']:.2e} r2 {hp['r2']:.2e} -> "
                  f"{'resolved' if hp['resolved'] else 'NOT RESOLVED'}", flush=True)
    return [k for k, r in recs.items() if not r["balances"]["ok"] and not res[k]["resolved"]]


# ------------------------------------------------------------------ classification (pure, given saved fields)
_graphs = {}


def classify_run(rec, theta):
    dec = json.load(open(DEC))
    gk = (rec["kind"], rec["slot"])
    if gk not in _graphs:
        g = graph_for(rec, dec)[0]
        _graphs[gk] = (g, mg._tree(g))
    g, tree = _graphs[gk]
    sX = np.load(os.path.join(OUT, "fields", rec["key"] + ".npz"))["checks_sX"].astype(float)
    return mg.classify(g, [mg.particles(g, sX[k], theta, tree) for k in range(len(CHECKS))])[0]


def threshold_at(Trec):
    """Proposal §3: fixed eligibility, then theta_b; otherwise theta_29 and an invalid T reference."""
    if Trec is None:
        return dict(theta=THETA_29, eligible=False, reason="no T run")
    c_elig = classify_run(Trec, THETA_ELIG)
    m = Trec["min_smoothed_final"]
    if c_elig == "single" and m <= ELIG_DEPTH:
        th = 0.5 * m
        cT = classify_run(Trec, th)
        if cT == "single":
            return dict(theta=th, eligible=True, T_class=cT, T_min=m)
        return dict(theta=THETA_29, eligible=False, reason=f"T {cT} under theta_b", T_min=m)
    return dict(theta=THETA_29, eligible=False, reason=f"T {c_elig} under theta_elig, min {m:.3f}", T_min=m,
                T_class_theta29=classify_run(Trec, THETA_29))


def table_at(b, recs):
    """Classes of the four D runs at b (MG4 records at b = 29)."""
    if abs(b - 29.0) < 1e-12:
        cl = json.load(open(MG4_EVAL))["classes"]
        return dict(theta=THETA_29, eligible=True, source="MG4 Part A",
                    D=[cl[f"D{s}_p{p}_sol"] for s, p in DRUNS])
    th = threshold_at(recs.get(f"{bkey(b)}_T_p0"))
    D = [classify_run(recs[f"{bkey(b)}_D{s}_p{p}"], th["theta"]) for s, p in DRUNS]
    D29 = [classify_run(recs[f"{bkey(b)}_D{s}_p{p}"], THETA_29) for s, p in DRUNS]
    mins = [recs[f"{bkey(b)}_D{s}_p{p}"]["min_smoothed_final"] for s, p in DRUNS]
    return dict(th, D=D, D_theta29=D29, D_min=mins)


all_single = lambda t: all(c == "single" for c in t["D"])
any_mult = lambda t: any(c == "multiplied" for c in t["D"])


# ------------------------------------------------------------------ Part 2 rules (pure functions of the grid table)
def plan(grid_tab):
    """grid_tab: {b: table} for b in [29] + GRID. Returns ('A', [(inside, outside), ...] or 'unbracketed' strings) or
    ('B', (lo, hi)) or ('B', None)."""
    bs = [29.0] + GRID                                              # descending
    single = [b for b in bs if all_single(grid_tab[b])]
    if single:
        top = single[0]
        i = bs.index(top)
        j = i
        while j + 1 < len(bs) and all_single(grid_tab[bs[j + 1]]):
            j += 1
        edges = []
        edges.append((top, bs[i - 1]) if i > 0 else f"upper edge > {top}")
        edges.append((bs[j], bs[j + 1]) if j + 1 < len(bs) else f"lower edge < {bs[j]}")
        return "A", edges
    nm = [b for b in bs if not any_mult(grid_tab[b])]
    if not nm:
        return "B", None
    hb = nm[0]
    return "B", (hb, bs[bs.index(hb) - 1])


def bisect_edge(inside, outside, pred, run_b, steps=STEPS, stop_on=None):
    """Bisection between inside (pred True) and outside (pred False). run_b(b) -> table. Returns bracket and visits."""
    visits = []
    for _ in range(steps):
        mid = round(0.5 * (inside + outside), 9)
        t = run_b(mid)
        visits.append((mid, t))
        if stop_on is not None and stop_on(t):
            return dict(bracket=None, stopped_at=mid, visits=visits)
        if pred(t):
            inside = mid
        else:
            outside = mid
    return dict(bracket=sorted([inside, outside]), visits=visits)


# ------------------------------------------------------------------ main
def main():
    for d in ("runs", "fields"):
        os.makedirs(os.path.join(OUT, d), exist_ok=True)
    go.continuum2d()
    log = dict(registered="GATES.md MG4c; spec/MG4c_proposal.md draft 3 (a946550)")
    status = "complete"
    recs = {}
    try:
        chk = part3_checks(json.load(open(DEC)))
        json.dump(chk, open(os.path.join(OUT, "part3_checks.json"), "w"), indent=1)
        log["part3_checks_ok"] = chk["ok"]
        p3 = [job(29.0, "D", 0, 0, "cotan"), job(29.0, "T", 0, 0, "cotan")] if chk["ok"] else []
        if not chk["ok"]:
            print("Part 3 construction checks FAILED: Part 3 not run (O-invalid). See part3_checks.json", flush=True)
        recs.update(run_batch([j for b in GRID for j in jobs_at(b)] + p3))
        if resolve_flags(recs):
            log["status"] = "NOT SCORED (unresolved balance flag)"
            return finish(log, recs)
        grid_tab = {b: table_at(b, recs) for b in [29.0] + GRID}
        log["grid"] = {bkey(b): t for b, t in grid_tab.items()}
        case, br = plan(grid_tab)
        log["case"] = case

        def run_b(b):
            recs.update(run_batch(jobs_at(b)))
            if resolve_flags(recs):
                raise Stop("unresolved balance flag")
            t = table_at(b, recs)
            log.setdefault("refined", {})[bkey(b)] = t
            return t

        if case == "A":
            log["edges"] = []
            for e in br:
                if isinstance(e, str):
                    log["edges"].append(dict(unbracketed=e))
                else:
                    log["edges"].append(bisect_edge(e[0], e[1], all_single, run_b))
        elif br is not None:
            log["case_B"] = bisect_edge(br[0], br[1], lambda t: not any_mult(t), run_b, stop_on=all_single)
        else:
            log["case_B"] = "multiplied at every tested b"
    except Stop as e:
        status = f"stopped: {e}"
    log["status"] = status
    return finish(log, recs)


def all_records():
    d = os.path.join(OUT, "runs")
    return {f[:-5]: json.load(open(os.path.join(d, f))) for f in sorted(os.listdir(d)) if f.endswith(".json")}


def finish(log, recs):
    recs = all_records()                       # includes runs finished in a batch that a stop interrupted
    ch, n = spent()
    log.update(core_h=ch, runs=n, platform=platform.platform())
    for k in ("w_b29.000000_D0_p0", "w_b29.000000_T_p0"):
        if k in recs:
            log.setdefault("part3", {})[k] = classify_run(recs[k], THETA_29)
    if "w_b29.000000_T_p0" in recs:
        yT = np.load(os.path.join(OUT, "fields", "w_b29.000000_T_p0.npz"))["final"]
        ref = np.load(TSTATE).reshape(yT.shape)               # MG4 Part 0 T reference: b = 29, centre node, same IC
        log["part3"]["T_control_maxdiff_X"] = float(np.abs(yT[1] - ref[1]).max())     # registered criterion: <= 1e-6
        log["part3"]["T_control_maxdiff_all"] = float(np.abs(yT - ref).max())
    json.dump(log, open(os.path.join(OUT, "driver_log.json"), "w"), indent=1, default=str)
    print("MG4c driver finished:", log.get("status"), f"{ch:.2f} core-h, {n} runs", flush=True)
    print("Evaluate with: python src/evaluate_mg4c.py", flush=True)
    return log


if __name__ == "__main__":
    main()

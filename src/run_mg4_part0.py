"""MG4 Part 0 driver (proposal §2–§4, §5.3 C3(a); implementation choices in GATES.md MG4, MG010).

Builds the pilot and panel graphs, calibrates h, computes spectra and the C3(a) ratios, runs T's positive control
(depth, existence window, box check) and writes results/MG4/part0/decisions.json. Resumable at step level (each step
writes its own JSON under results/MG4/part0/steps/). Workers: MODELG_WORKERS (default 2).
    python src/run_mg4_part0.py
"""
import os, sys, json, time, platform, csv
import numpy as np
from concurrent.futures import ProcessPoolExecutor
sys.path.insert(0, os.path.dirname(__file__))
import graphs as gr
import graph_ode as go
import measure_graph as mg
from runner import seed_for

OUT = os.environ.get("MG4_OUT") or os.path.join(os.path.dirname(__file__), "..", "results", "MG4", "part0")
STEPS = os.path.join(OUT, "steps")
GATE, B_IDX = 5, 0
FAM = dict(T=0, D=1, R=2, N3=3)
H_T, NX, NY = 0.2, 122, 142
if os.environ.get("MG4_SMOKE"):                     # plumbing test only: tiny torus, never a registered run
    NX, NY = 62, 70
H_T3, NX3, NY3 = 0.3, 82, 94
BOX_NX, BOX_NY = (184, 212) if not os.environ.get("MG4_SMOKE") else (74, 84)
KSTAR = 2 * np.pi / 3.08
CUTOFF = 1.5 * KSTAR**2
TOL_H = 0.03
C_GUESS = dict(D=2.39, R=12**2 / (8 * np.pi))
WORKERS = int(os.environ.get("MODELG_WORKERS", "2"))


def rng_for(fam, stage, rep):
    return np.random.default_rng(seed_for(GATE, B_IDX, FAM[fam], stage, rep))


def build(fam, n_pts, Lx, Ly, rep):
    rng = rng_for(fam, 0, rep)
    pos = gr.poisson_points(n_pts, Lx, Ly, rng)
    if fam == "D":
        return gr.delaunay(pos, Lx, Ly)
    rho = n_pts / (Lx * Ly)
    return gr.rgg(pos, Lx, Ly, np.sqrt(12 / (np.pi * rho)))


def ledger(item, key, wall):
    path = os.path.join(OUT, "budget_ledger.csv")
    new = not os.path.exists(path)
    with open(path, "a", newline="") as f:
        w = csv.writer(f)
        if new:
            w.writerow(["item", "key", "wall_s", "workers", "machine", "os", "python", "numpy", "scipy", "utc_time"])
        import scipy
        w.writerow([item, key, f"{wall:.2f}", WORKERS, platform.node(), platform.platform(), platform.python_version(),
                    np.__version__, scipy.__version__, time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())])


def step(name, fn, *a):
    path = os.path.join(STEPS, name + ".json")
    if os.path.exists(path):
        return json.load(open(path))
    t0 = time.time()
    out = fn(*a)
    out["wall_s"] = time.time() - t0
    json.dump(out, open(path, "w"), indent=1)
    ledger("part0", name, out["wall_s"])
    return out


# ------------------------------------------------------------------ graph steps
def graph_summary(g):
    c = gr.corrector(g)
    deg = gr.degrees(g)
    return dict(N=int(g["N"]), E=int(len(g["I"])), deg_min=int(deg.min()), deg_max=int(deg.max()),
                deg_mean=float(deg.mean()), D=c["D"].tolist(), eig=c["eig"].tolist(), angle=c["angle"],
                anisotropy=float((c["eig"][0] - c["eig"][1]) / (c["eig"][0] + c["eig"][1])),
                h=float(np.sqrt(c["h2"])), h_naive=float(np.sqrt(c["h2_naive"])),
                kept_fraction=float(g.get("kept_fraction", 1.0)))


def pilot(fam, rep, Lx, Ly, h_t):
    n = int(round(C_GUESS[fam] * Lx * Ly / h_t**2))
    g = build(fam, n, Lx, Ly, rep)
    s = graph_summary(g)
    s.update(n_points=n, c_F=n / (Lx * Ly) * s["h"] ** 2, rep=rep)
    return s


def panel_graph(fam, slot, c_F, Lx, Ly, h_t, spare):
    """Panel slot: rep = slot, then the family's next unused replacement reps (50, 51, ...) in order, until
    |h/h_t - 1| <= TOL_H. `spare` is the family's list of unused replacement reps (consumed in place)."""
    n = int(round(c_F * Lx * Ly / h_t**2))
    tried = []
    reps = [slot]
    while True:
        rep = reps.pop(0) if reps else spare.pop(0)
        g = build(fam, n, Lx, Ly, rep)
        s = graph_summary(g)
        s.update(n_points=n, rep=rep)
        tried.append(dict(rep=rep, h=s["h"]))
        if abs(s["h"] / h_t - 1) <= TOL_H:
            s["tried"] = tried
            return s


def spectral_step(fam, rep, n_pts, Lx, Ly, h):
    if fam == "T":
        g = gr.triangular(gr.torus(H_T, NX, NY))
    elif fam == "N3":
        g = gr.random_regular(n_pts, 6, seed_for(GATE, B_IDX, FAM["N3"], 0, rep))
    else:
        g = build(fam, n_pts, Lx, Ly, rep)
    nu = gr.spectrum(g, h, CUTOFF, k0=300)
    ks, r = gr.count_ratios(nu, Lx, Ly, KSTAR)
    return dict(fam=fam, rep=rep, h=h, n_eig=int(len(nu)), nu_head=nu[:10].tolist(), ks=ks.tolist(),
                ratios=r.tolist(), mean_ratio=float(r.mean()), pass_a=bool(abs(r.mean() - 1) <= 0.05))


# ------------------------------------------------------------------ T positive control
def t_reference(nx, ny):
    tor = gr.torus(H_T, nx, ny)
    g = gr.triangular(tor)
    h = np.sqrt(gr.corrector(g)["h2"])
    centre = np.array([tor["Lx"] / 2, tor["Ly"] / 2])
    i0 = int(np.argmin(np.hypot(*gr.min_image(g["pos"] - centre, tor["Lx"], tor["Ly"]).T)))
    c = g["pos"][i0]
    sol, rhs = go.integrate(g, h, go.soliton_ic(g, c), 300.0)
    y = sol.y[:, -1]
    N = g["N"]
    S = mg.smoother(g)
    sX = S @ y[N:2 * N]
    rc, prof, cnt = mg.radial_profile(g, y[N:2 * N], c)
    path = os.path.join(STEPS, f"Tstate_{nx}x{ny}.npy")
    np.save(path, y)
    return dict(nx=nx, ny=ny, N=N, h=h, centre=c.tolist(), depth=float(y[N:2 * N].min()),
                depth_s=float(sX.min()), max_dphidt=float(np.abs(rhs(300.0, y)).max()),
                balances_ok=go.balances(y, rhs(300.0, y), N)["ok"], rc=rc.tolist(),
                prof=np.where(np.isfinite(prof), prof, np.nan).tolist(), nfev=int(sol.nfev))


def window_edge(which, thr):
    """Bisection for b_s (exists: >= 1 particle at t = 300) or b_r (replicating: >= 2 particles at t = 300)."""
    tor = gr.torus(H_T, NX, NY)
    g = gr.triangular(tor)
    h = np.sqrt(gr.corrector(g)["h2"])
    y0 = np.load(os.path.join(STEPS, f"Tstate_{NX}x{NY}.npy"))
    S = mg.smoother(g)
    N = g["N"]
    from lattice import P17
    lo, hi, n = (26.0, 29.0, 10) if which == "b_s" else (29.0, 30.04, 9)
    if os.environ.get("MG4_SMOKE"):
        n = 2
    hist = []
    for _ in range(n):
        m = 0.5 * (lo + hi)
        sol, _ = go.integrate(g, h, y0, 300.0, P={**P17, "b": m})
        k = len(mg.particles(g, S @ sol.y[N:2 * N, -1], thr))
        hist.append(dict(b=m, n_particles=k))
        if which == "b_s":
            lo, hi = (lo, m) if k >= 1 else (m, hi)
        else:
            lo, hi = (m, hi) if k < 2 else (lo, m)
    return dict(which=which, interval=[lo, hi], history=hist)


def main():
    os.makedirs(STEPS, exist_ok=True)
    go.continuum2d(path=os.path.join(OUT, "continuum2d.npz"))
    tor = gr.torus(H_T, NX, NY)
    Lx, Ly = tor["Lx"], tor["Ly"]
    T = step("T_graph", lambda: graph_summary(gr.triangular(tor)))
    with ProcessPoolExecutor(WORKERS) as ex:
        futs = {f"pilot_{f}_{r}": ex.submit(step, f"pilot_{f}_{r}", pilot, f, r, Lx, Ly, H_T)
                for f in ("D", "R") for r in range(90, 94)}
        pilots = {k: v.result() for k, v in futs.items()}
    c_F = {f: float(np.mean([pilots[f"pilot_{f}_{r}"]["c_F"] for r in range(90, 94)])) for f in ("D", "R")}
    panel = {}
    for f in ("D", "R"):
        spare = list(range(50, 80))
        for slot in range(4):
            panel[f"{f}{slot}"] = step(f"panel_{f}_{slot}", panel_graph, f, slot, c_F[f], Lx, Ly, H_T, spare)
            done = {t["rep"] for t in panel[f"{f}{slot}"]["tried"]}
            spare[:] = [r for r in spare if r not in done]          # also correct when the step was resumed
    hD = float(np.mean([panel[f"D{s}"]["h"] for s in range(4)]))
    nD = int(round(np.mean([panel[f"D{s}"]["N"] for s in range(4)])))
    nD += nD % 2                                                    # k N must be even for a 6-regular graph
    jobs = [("T", 0, 0, T["h"])] + [(f, panel[f"{f}{s}"]["rep"], panel[f"{f}{s}"]["n_points"], panel[f"{f}{s}"]["h"])
                                     for f in ("D", "R") for s in range(4)] + [("N3", r, nD, hD) for r in (0, 1)]
    with ProcessPoolExecutor(WORKERS) as ex:
        futs = {f"spec_{f}_{r}": ex.submit(step, f"spec_{f}_{r}", spectral_step, f, r, n, Lx, Ly, h)
                for f, r, n, h in jobs}
        spec = {k: v.result() for k, v in futs.items()}
    ref = step(f"Tref_{NX}x{NY}", t_reference, NX, NY)
    thr = -0.5 * abs(ref["depth_s"])
    depth_ok = abs(ref["depth"] / -13.6 - 1) <= 0.20
    with ProcessPoolExecutor(WORKERS) as ex:
        fb = ex.submit(step, f"Tbox_{BOX_NX}x{BOX_NY}", t_reference, BOX_NX, BOX_NY)
        fs = ex.submit(step, "window_b_s", window_edge, "b_s", thr)
        fr = ex.submit(step, "window_b_r", window_edge, "b_r", thr)
        box, ws, wr = fb.result(), fs.result(), fr.result()
    rc = np.array(ref["rc"]); p1 = np.array(ref["prof"], float); p2 = np.array(box["prof"], float)
    box_l2 = mg.profile_l2(rc, p1, p2)
    in_window = ws["interval"][1] <= 29.0 - 0.1 and wr["interval"][0] >= 29.0 + 0.1
    dec = dict(
        design="spec/MG4_proposal.md draft 4 (773255d); GATES.md MG4",
        torus=dict(h_t=H_T, nx=NX, ny=NY, Lx=Lx, Ly=Ly, box=[BOX_NX, BOX_NY], h03=[NX3, NY3]),
        c_F=c_F, pilots=pilots, T=T, panel=panel, N3=dict(N=nD, h=hD, reps=[0, 1]),
        spectra={k: {kk: v[kk] for kk in ("fam", "rep", "h", "n_eig", "mean_ratio", "pass_a", "ratios")}
                 for k, v in spec.items()},
        positive_control=dict(depth=ref["depth"], depth_s=ref["depth_s"], depth_ok=depth_ok, threshold_node_ref=thr,
                              b_s=ws["interval"], b_r=wr["interval"], in_window=in_window,
                              box_l2=box_l2, box_ok=bool(box_l2 < 0.005),
                              T_max_dphidt=ref["max_dphidt"], T_balances_ok=ref["balances_ok"]),
        stop_for_PI=bool(not (depth_ok and in_window)),
        torus_for_partA="30x30" if box_l2 >= 0.005 else "20x20",
        platform=dict(machine=platform.node(), os=platform.platform(), python=platform.python_version(),
                      numpy=np.__version__, workers=WORKERS),
        note="T's registered threshold uses the mean of 4 placements (Part A); threshold_node_ref is the node-placement "
             "value used for the window bisection.")
    json.dump(dec, open(os.path.join(OUT, "decisions.json"), "w"), indent=1)
    print(json.dumps(dec["positive_control"], indent=1))
    print("stop_for_PI:", dec["stop_for_PI"], " torus:", dec["torus_for_partA"])


if __name__ == "__main__":
    main()

"""MG4 cost benchmark (timing only; never used in any estimate or criterion; GATES.md MG4, MG010 "Timing").

Runs one full-size graph-ODE particle run (t = 0..300, the registered solver settings) on each of T, D panel graph 0
and R panel graph 0, concurrently (3 processes), from the continuum soliton at the torus centre. Writes
results/MG4/timing/timing_<machine>.json with wall times and the platform.
    python src/time_mg4.py
"""
import os, sys, json, time, platform
import numpy as np
from concurrent.futures import ProcessPoolExecutor
sys.path.insert(0, os.path.dirname(__file__))
import graphs as gr
import graph_ode as go
import run_mg4_part0 as p0

OUT = os.path.join(os.path.dirname(__file__), "..", "results", "MG4", "timing")
GRAPHS = dict(T=None, D=(0, 18901), R=(0, 38595))      # (rep, n_points) of panel graphs D0 and R0 (Part 0 steps)


def one(fam):
    tor = gr.torus(p0.H_T, p0.NX, p0.NY)
    if fam == "T":
        g = gr.triangular(tor)
    else:
        rep, n = GRAPHS[fam]
        g = p0.build(fam, n, tor["Lx"], tor["Ly"], rep)
    h = float(np.sqrt(gr.corrector(g)["h2"]))
    c = np.array([tor["Lx"] / 2 + 0.031, tor["Ly"] / 2 + 0.047])
    t0 = time.time()
    sol, rhs = go.integrate(g, h, go.soliton_ic(g, c), 300.0)
    wall = time.time() - t0
    N = g["N"]
    return dict(fam=fam, N=int(N), h=h, wall_s=wall, nfev=int(sol.nfev), nlu=int(sol.nlu))


def main():
    os.makedirs(OUT, exist_ok=True)
    go.continuum2d()
    with ProcessPoolExecutor(3) as ex:
        res = list(ex.map(one, ["T", "D", "R"]))
    out = dict(runs=res, machine=platform.node(), os=platform.platform(), python=platform.python_version(),
               numpy=np.__version__, processes=3, utc=time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
               note="timing only; not used in any estimate")
    path = os.path.join(OUT, f"timing_{platform.node()}.json")
    json.dump(out, open(path, "w"), indent=1)
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()

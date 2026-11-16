"""EX001 (exploratory): MG192 gravity drift recheck at MG193's operating point a = 14, b = 29.20526.

Reuses src/mg192_mobility.py unchanged, with its parameter b set to the MG193 window centre.
"""
import json
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import mg192_mobility as m  # noqa: E402

BSTAR = 29.205258852005002
OUT = os.path.join(m.ROOT, "results", "EX001_gravity_recheck")


def dynamic_run_depth(grid, U, eps, x0, T=400.0, dt=0.002, sample=1.0):
    """Same as mg192_mobility.dynamic_run, but also records the core depth (min X minus local background)."""
    bg, a, bres = m.modulated_background(grid, eps)
    h0 = m.hom()
    u0 = bg + np.array([grid.shift(U[s], x0) for s in range(3)]) - h0[:, None]
    st = m.ETDRK4(grid, dt)
    vh = np.fft.fft(u0, axis=1)
    n_steps, every = int(round(T / dt)), int(round(sample / dt))
    ts, xs, ds = [0.0], [m.xmin_position(grid, u0[1], x0)], [float(np.min(u0[1] - bg[1]))]
    for n in range(1, n_steps + 1):
        vh = st.step(vh, a)
        if n % every == 0:
            u = np.real(np.fft.ifft(vh, axis=1))
            if not np.all(np.isfinite(u)) or np.any(u <= 0):
                raise RuntimeError("lost positivity/finiteness")
            ts.append(n * dt)
            xs.append(m.xmin_position(grid, u[1], xs[-1]))
            ds.append(float(np.min(u[1] - bg[1])))
    ts, ds = np.array(ts), np.array(ds)
    xs = np.unwrap(np.array(xs) * 2 * np.pi / grid.L) * grid.L / (2 * np.pi)
    sel = (ts >= 100.0) & (ts <= T)
    v = np.polyfit(ts[sel], xs[sel], 1)[0]
    dadx = float(m.P["a"] * eps * 2 * np.pi / grid.L * np.cos(2 * np.pi * x0 / grid.L))
    return dict(eps=eps, x0=x0, dt=dt, T=T, v=float(v), bg_residual=bres, dadx=dadx, depth_start=float(ds[0]),
                depth_end=float(ds[-1]), depth_min=float(ds.min()), depth_max=float(ds.max())), ts, xs, ds


def main():
    t0 = time.time()
    os.makedirs(OUT, exist_ok=True)
    g = m.Grid(100.0, 1024)
    U = np.load(os.path.join(m.ROOT, "results", "MG192_gravity_mobility", "soliton_adjoint_L100_N1024.npz"))["U"]
    # continue from b = 29 to b* in steps of 0.01
    nb = 21
    for i in range(1, nb + 1):
        m.P["b"] = 29.0 + (BSTAR - 29.0) * i / nb
        U, res, _ = m.newton_even(g, U, tol=1e-10)
    print(f"soliton at b*: residual {res:.2e}, depth {U[1].min() - m.hom()[1]:.4f}", flush=True)
    lam0, others, psi, Up, J = m.spectrum_and_adjoint(g, U)
    mob = m.mobilities(g, U, psi)
    lead = float(np.max(np.real(others)))
    print(f"translation eig {lam0:.2e}; least-stable other {lead:.4f}; mobilities {mob}", flush=True)
    out = dict(b=BSTAR, a=14.0, residual=res, translation_eig=float(np.real(lam0)), least_stable_other=lead,
               mobilities=mob, runs=[], platform=dict(numpy=np.__version__))
    series = {}
    L = g.L
    for eps, x0 in [(0.001, 0.0), (0.002, 0.0), (0.001, L / 8), (0.002, L / 8), (0.001, -L / 8), (0.002, -L / 8),
                    (0.0, 0.0), (0.0, L / 8)]:
        t1 = time.time()
        r, ts, xs, ds = dynamic_run_depth(g, U, eps, x0)
        r["predicted_v"] = mob["mu_a"] * r["dadx"]
        r["rel_err"] = (r["v"] - r["predicted_v"]) / r["predicted_v"] if eps > 0 else None
        r["local_a_rel"] = eps * np.sin(2 * np.pi * x0 / L)
        r["drift_total"] = float(xs[-1] - xs[0])
        r["wall_s"] = time.time() - t1
        out["runs"].append(r)
        series[f"eps{eps}_x0{x0}"] = np.vstack([ts, xs, ds])
        print(f"eps={eps} x0={x0:+.2f}: v={r['v']:.5e} pred={r['predicted_v']:.5e} rel={r['rel_err']} "
              f"drift={r['drift_total']:.4e} depth {r['depth_start']:.4f}->{r['depth_end']:.4f} ({r['wall_s']:.0f}s)", flush=True)
        json.dump(out, open(os.path.join(OUT, "results.json"), "w"), indent=1, default=float)
    np.savez_compressed(os.path.join(OUT, "trajectories.npz"), **series)
    out["wall_seconds"] = time.time() - t0
    json.dump(out, open(os.path.join(OUT, "results.json"), "w"), indent=1, default=float)
    print("done", out["wall_seconds"], flush=True)


if __name__ == "__main__":
    main()

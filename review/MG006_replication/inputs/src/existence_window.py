"""MG003: lower edge b_s of the stationary-soliton existence window (deterministic continuation in b).
For each b: start from the b=29 soliton, run seed-free to t=300, soliton 'exists' if phi_X(0) < -1.
Writes results/MG003/existence_window.txt."""
import sys
import numpy as np
sys.path.insert(0, "src")
from lattice import lattice_ode_1d, P17
from mg_continuum import run

lines = []
def bisect(exists, lo=28.5, hi=29.0, it=12):
    for _ in range(it):
        m = 0.5 * (lo + hi)
        if exists(m): hi = m
        else: lo = m
    return lo, hi

DIMS = [int(a) for a in sys.argv[1].split(",")] if len(sys.argv) > 1 else [1, 2, 3]
LATTICE = len(sys.argv) <= 1
for h in ((0.2, 0.1) if LATTICE else ()):
    N = int(round(100 / h))
    x, pG, pX, pY, _ = lattice_ode_1d(P17, N, h, 400.0)
    y0 = np.concatenate([pG, pX, pY])
    ex = lambda b: lattice_ode_1d({**P17, "b": b}, N, h, 300.0, seed_scale=0.0, y0=y0)[2][N // 2] < -1.0
    lo, hi = bisect(ex)
    lines.append(f"lattice h={h}: b_s in [{lo:.5f}, {hi:.5f}]"); print(lines[-1], flush=True)

for n in DIMS:
    g, P, s, rhs = run(n, h=0.025, t_end=300.0)
    y0 = s.y[:, -1]
    def ex(b, n=n, y0=y0):
        g2, P2, s2, _ = run(n, h=0.025, t_end=300.0, seed_scale=0.0, P={**P17, "b": b}, y0=y0)
        return s2.y[g2.N, -1] < -1.0
    lo, hi = bisect(ex, lo=(28.5 if n == 1 else 24.0), hi=29.0, it=(12 if n == 1 else 16))
    lines.append(f"continuum {n}D (h=0.025): b_s in [{lo:.5f}, {hi:.5f}]"); print(lines[-1], flush=True)
lines.append("Turing threshold b_c: continuum 30.04408; lattice h=0.1 30.0476, h=0.2 30.0458.")
open("results/MG003/existence_window.txt", "a").write("\n".join(lines) + "\n")

"""MG003: upper edge b_r of the single-soliton window (onset of self-replication), lattice h=0.2 and 0.1.
Single soliton = exactly one deep minimum (phi_X < -4) at t=400 after the standard seed.
Appends to results/MG003/existence_window.txt."""
import sys
import numpy as np
sys.path.insert(0, "src")
from lattice import lattice_ode_1d, P17

def nmin(b, h):
    N = int(round(100 / h))
    x, pG, pX, pY, m = lattice_ode_1d({**P17, "b": b}, N, h, 400.0)
    return int(np.sum((pX[1:-1] < pX[:-2]) & (pX[1:-1] < pX[2:]) & (pX[1:-1] < -4)))
out = []
for h in (0.2, 0.1):
    lo, hi = 29.3, 29.5
    for _ in range(9):
        m = 0.5 * (lo + hi)
        if nmin(m, h) > 1: hi = m
        else: lo = m
    out.append(f"lattice h={h}: self-replication onset b_r in [{lo:.4f}, {hi:.4f}] (single soliton below)")
    print(out[-1], flush=True)
open("results/MG003/replication_edge.txt", "w").write("\n".join(out) + "\n")

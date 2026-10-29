"""MG003: references at the h=0.2 window centre b_mid=29.16, and the b=29 vs b_mid profile difference
(expected plateau of the MG2 mis-specified-reference control N2). Writes results/MG003/lattice_reference_bmid.txt."""
import sys
import numpy as np
sys.path.insert(0, "src")
from lattice import lattice_ode_1d, lna_homogeneous, P17

h, N = 0.2, 500; i0 = N // 2
prof = dict(np.load("results/MG003/ref_profiles.npz"))
out = []
for h_, N_ in ((0.2, 500), (0.1, 1000)):
    P = {**P17, "b": 29.16}
    x, pG, pX, pY, mdt = lattice_ode_1d(P, N_, h_, 400.0)
    prof[f"h{h_}_b29.16"] = np.vstack([x, pG, pX, pY])
    out.append(f"b=29.16 h={h_}: amp G {pG[N_//2]:.4f} X {pX[N_//2]:.4f} Y {pY[N_//2]:.4f}; max|dphi/dt| {mdt:.1e}")
lna = lna_homogeneous({**P17, "b": 29.16}, N, h)
out.append(f"b=29.16 h=0.2: LNA Omega*Var G {lna[0]:.3f} X {lna[1]:.3f} Y {lna[2]:.3f}")
lo, hi = 0.0, 2.0
for _ in range(14):
    m = 0.5 * (lo + hi)
    if lattice_ode_1d({**P17, "b": 29.16}, N, h, 150.0, seed_scale=m)[2][i0] < -1.0: hi = m
    else: lo = m
out.append(f"b=29.16 h=0.2: seed threshold {lo:.4f}-{hi:.4f}")
A, B = prof["h0.2_b29.0"], prof["h0.2_b29.16"]
m = np.abs(A[0]) < 10
E = [np.linalg.norm(B[i][m] - A[i][m]) / np.linalg.norm(A[i][m]) for i in (1, 2, 3)]
out.append(f"h=0.2 L2 difference b=29 vs b=29.16 (|x|<10): G {100*E[0]:.2f}% X {100*E[1]:.2f}% Y {100*E[2]:.2f}% -> N2 expected plateau {100*max(E):.2f}%")
np.savez("results/MG003/ref_profiles.npz", **prof)
open("results/MG003/lattice_reference_bmid.txt", "w").write("\n".join(out) + "\n"); print("\n".join(out))

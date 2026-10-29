"""MG003: lattice references at the interior point b=29.5 (h=0.2, 0.1). Appends profiles to
results/MG003/ref_profiles.npz; writes results/MG003/lattice_reference_b29.5.txt."""
import sys
import numpy as np
sys.path.insert(0, "src")
from lattice import lattice_ode_1d, lna_homogeneous, P17
from mg_continuum import zeros

P = {**P17, "b": 29.5}
prof = dict(np.load("results/MG003/ref_profiles.npz"))
lines = []
for h in (0.2, 0.1):
    N = int(round(100 / h)); i0 = N // 2
    x, pG, pX, pY, mdt = lattice_ode_1d(P, N, h, 400.0)
    prof[f"h{h}_b29.5"] = np.vstack([x, pG, pX, pY])
    z = zeros(x[i0:], pX[i0:]); lam = 2 * np.diff(z[z < 30])[1:].mean()
    lo, hi = 0.0, 2.0
    for _ in range(14):
        m = 0.5 * (lo + hi)
        if lattice_ode_1d(P, N, h, 150.0, seed_scale=m)[2][i0] < -1.0: hi = m
        else: lo = m
    lna = lna_homogeneous(P, N, h)
    lines.append(f"b=29.5 h={h}: amp G {pG[i0]:.4f} X {pX[i0]:.4f} Y {pY[i0]:.4f}; lam {lam:.4f}; "
                 f"max|dphi/dt| {mdt:.1e}; seed threshold {lo:.4f}-{hi:.4f}; LNA Omega*Var G {lna[0]:.3f} X {lna[1]:.3f} Y {lna[2]:.3f}")
    print(lines[-1], flush=True)
np.savez("results/MG003/ref_profiles.npz", **prof)
open("results/MG003/lattice_reference_b29.5.txt", "w").write("\n".join(lines) + "\n")

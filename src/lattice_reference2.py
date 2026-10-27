"""MG003: save h=0.2 / h=0.1 lattice soliton reference profiles (b=29) and the mis-specified-reference
control profile (b=28.5). Writes results/MG003/ref_profiles.npz and ref_control.txt."""
import sys
import numpy as np
sys.path.insert(0, "src")
from lattice import lattice_ode_1d, P17

L = 100.0; prof = {}
for h in (0.2, 0.1):
    N = int(round(L / h))
    x, pG, pX, pY, mdt = lattice_ode_1d(P17, N, h, 400.0)
    prof[f"h{h}_b29.0"] = np.vstack([x, pG, pX, pY])
    # b=28.5: the standard seed does not nucleate; obtain the soliton (if it exists) by continuation
    x, pG2, pX2, pY2, mdt2 = lattice_ode_1d({**P17, "b": 28.5}, N, h, 400.0, seed_scale=0.0,
                                         y0=np.concatenate([pG, pX, pY]))
    prof[f"h{h}_b28.5"] = np.vstack([x, pG2, pX2, pY2])
    print(f"h={h}: b=28.5 continuation max|dphi/dt|={mdt2:.1e}, phiX(0)={pX2[N//2]:.4f} (b=29: {pX[N//2]:.4f})")
np.savez("results/MG003/ref_profiles.npz", **prof)
lines = []
for h in (0.2, 0.1):
    A, B = prof[f"h{h}_b29.0"], prof[f"h{h}_b28.5"]
    m = np.abs(A[0]) < 10
    E = {k: np.linalg.norm(B[i][m] - A[i][m]) / np.linalg.norm(A[i][m]) for i, k in ((1, "G"), (2, "X"), (3, "Y"))}
    lines.append(f"h={h}: L2 difference (|x|<10) between b=28.5 and b=29 solitons: " +
                 ", ".join(f"{k} {100*v:.2f}%" for k, v in E.items()) + f"; max {100*max(E.values()):.2f}%")
print("\n".join(lines)); open("results/MG003/ref_control.txt", "w").write("\n".join(lines) + "\n")

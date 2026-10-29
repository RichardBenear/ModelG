"""MG003: deterministic lattice references for the Phase 2 spec (no stochastic outcomes).
Writes results/MG003/lattice_reference.json and .txt."""
import json, sys
import numpy as np
sys.path.insert(0, "src")
from lattice import lattice_ode_1d, lna_homogeneous, steady_state, P17
from mg_continuum import zeros

P = P17; L = 100.0
MG1_1D = dict(G=-0.16101, X=-8.3645, Y=0.92991)   # MG002 1D continuum, h=0.0125, t=300 (stationary)
out, lines = {}, []
for h in (0.05, 0.1, 0.2):
    N = int(round(L / h))
    x, pG, pX, pY, mdt = lattice_ode_1d(P, N, h, 400.0)
    i0 = N // 2
    z = zeros(x[i0:], pX[i0:]); zz = np.diff(z[(z < 30)])[1:]
    lam = 2 * zz.mean()
    amp = dict(G=pG[i0], X=pX[i0], Y=pY[i0])
    # deterministic seed threshold (judged at t=150)
    lo, hi = 0.0, 2.0
    for _ in range(14):
        m = 0.5 * (lo + hi)
        if lattice_ode_1d(P, N, h, 150.0, seed_scale=m)[2][i0] < -1.0: hi = m
        else: lo = m
    # lattice Turing threshold b_c (max over lattice modes)
    G0, X0, _ = steady_state(P)
    def maxre(bb):
        Y = bb / X0
        J = np.array([[-1, .1, 0], [1, -(1 + bb) + 2 * X0 * Y, X0**2], [0, bb - 2 * X0 * Y, -X0**2]])
        ks = np.pi * np.arange(N // 2 + 1) * 2 / (N * h)
        lm = -(4 / h**2) * np.sin(ks * h / 2) ** 2
        return max(np.linalg.eigvals(J + l * np.diag([1, 1, 12])).real.max() for l in lm)
    from scipy.optimize import brentq
    bc = brentq(maxre, 29, 31, xtol=1e-5)
    lna = lna_homogeneous(P, N, h)
    out[str(h)] = dict(N=N, amp=amp, lam=lam, max_dphidt=mdt, seed_threshold=[lo, hi], b_c=bc,
                       lna_Omega_var=dict(G=lna[0], X=lna[1], Y=lna[2]),
                       amp_vs_continuum_pct={k: 100 * (amp[k] - MG1_1D[k]) / abs(MG1_1D[k]) for k in amp})
    lines.append(f"h={h} N={N}: amp G {amp['G']:.4f} X {amp['X']:.4f} Y {amp['Y']:.4f} "
                 f"(vs continuum {', '.join(f'{k} {v:+.2f}%' for k, v in out[str(h)]['amp_vs_continuum_pct'].items())}); "
                 f"lam {lam:.4f}; max|dphi/dt| {mdt:.1e}; seed threshold {lo:.4f}-{hi:.4f}; b_c {bc:.4f}; "
                 f"LNA Omega*Var(c) G {lna[0]:.3f} X {lna[1]:.3f} Y {lna[2]:.3f}")
    print(lines[-1], flush=True)
json.dump(out, open("results/MG003/lattice_reference.json", "w"), indent=1)
open("results/MG003/lattice_reference.txt", "w").write("\n".join(lines) + "\n")

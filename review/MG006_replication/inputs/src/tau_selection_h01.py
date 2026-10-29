"""MG003a: bounded tau selection for tau_leap_euler_1d (h = 0.1). For tau_k = 0.4 h^2/d_Y / 2^k, k=1..3:
 (i)  Omega->inf mean map (forward Euler) from the lattice-ODE soliton keeps it: E < 1e-6 at t = 100;
 (ii) discrete-time LNA single-node variance within 2% of continuous-time LNA (G, X, Y);
 (iii) max per-token exit-probability sum on the deterministic soliton state <= 0.8.
Both b = 29 and 29.16. Writes results/MG003a/tau_selection_h0.1.txt."""
import sys
import numpy as np
sys.path.insert(0, "src")
from lattice import euler_mean_map, lna_homogeneous, lna_discrete_euler, steady_state, P17

prof = np.load("results/MG003/ref_profiles.npz")
h, N = 0.1, 1000
lines = []
for k in (1, 2, 3):
    tau = 0.4 * h * h / 12 / 2**k
    ok = True
    for b in (29.0, 29.16):
        P = {**P17, "b": b}
        G0, X0, Y0 = steady_state(P)
        ref = prof[f"h0.1_b{b}"]
        G, X, Y = euler_mean_map(P, h, tau, 100.0, G0 + ref[1], X0 + ref[2], Y0 + ref[3])
        m = np.abs(ref[0]) < 10
        E = max(np.linalg.norm(o[m] - r[m]) / np.linalg.norm(r[m])
                for o, r in zip((G - G0, X - X0, Y - Y0), ref[1:]))
        lc = lna_homogeneous(P, N, h); ld = lna_discrete_euler(P, N, h, tau)
        dev = 100 * (ld - lc) / lc
        Xs, Ys = X0 + ref[2], Y0 + ref[3]
        psum = max((1 + b) * tau + 2 * tau / h**2, (P["p"] + P["q"] - P["p"]) * tau + 2 * tau / h**2,
                   float(np.max(Xs**2 * tau)) + 2 * 12 * tau / h**2)
        c = (E < 1e-6) and np.all(np.abs(dev) <= 2) and psum <= 0.8
        ok &= c
        lines.append(f"tau={tau:.4e} b={b}: mean-map E={E:.1e}; LNA discrete vs continuous G {dev[0]:+.2f}% X {dev[1]:+.2f}% Y {dev[2]:+.2f}% "
                     f"(discrete {ld[0]:.3f}/{ld[1]:.3f}/{ld[2]:.3f}); max exit-prob sum {psum:.3f} -> {'OK' if c else 'NO'}")
        print(lines[-1], flush=True)
    lines.append(f"tau={tau:.4e}: {'SELECTED' if ok else 'rejected'}")
    print(lines[-1], flush=True)
    if ok:
        break
open("results/MG003a/tau_selection_h0.1.txt", "w").write("\n".join(lines) + "\n")

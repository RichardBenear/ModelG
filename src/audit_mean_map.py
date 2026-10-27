"""MG003a: deterministic (Omega -> infinity) mean map of the MG003 tau-leap scheme (split, exponential
losses, explicit gains), started from the stationary lattice-ODE soliton. Reproduces review item 1."""
import sys
import numpy as np
sys.path.insert(0, "src")
from lattice import steady_state, P17

prof = np.load("results/MG003/ref_profiles.npz")
def mean_map_split(P, h, tau, T, y):
    G0, X0, Y0 = steady_state(P); b, g, p, q = P["b"], P["g"], P["p"], P["q"]
    G, X, Y = G0 + y[1], X0 + y[2], Y0 + y[3]
    PG, PX = 1 - np.exp(-q * tau), 1 - np.exp(-(1 + b) * tau)
    hp = lambda c, d: c * (1 - 2 * d * tau / h**2) + d * tau / h**2 * (np.roll(c, 1) + np.roll(c, -1))
    for _ in range(int(round(T / tau))):
        xl = PX * X; r4 = tau * X**2 * Y
        G, X, Y = (G + P["a"] * tau - PG * G + g / (1 + b) * xl,
                   X + p / q * PG * G - xl + r4,
                   Y + b / (1 + b) * xl - r4)
        G, X, Y = hp(G, 1.0), hp(X, P["dx"]), hp(Y, P["dy"])
    return G - G0, X - X0, Y - Y0

def E(ref, out):
    m = np.abs(ref[0]) < 10
    return max(np.linalg.norm(o[m] - r[m]) / np.linalg.norm(r[m]) for o, r in zip(out, ref[1:]))

h = 0.2; tau0 = 0.4 * h * h / 12
lines = []
for b in (29.0, 29.16):
    ref = prof[f"h0.2_b{b}"]
    for f in (1, 0.5):
        out = mean_map_split({**P17, "b": b}, h, tau0 * f, 100.0, ref)
        lines.append(f"split/exp scheme b={b} tau={tau0*f:.3e}: phiX(0) {out[1][250]:.3f} (ODE {ref[2][250]:.3f}), E={100*E(ref, out):.1f}%")
        print(lines[-1], flush=True)
open("results/MG003a/audit_split_scheme.txt", "w").write("\n".join(lines) + "\n")

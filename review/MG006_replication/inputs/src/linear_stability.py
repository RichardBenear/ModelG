"""MG001: homogeneous steady state and linear (Turing) dispersion of dimensionless Model G.

Equations: Pulver & LaViolette (2013), eqs. (7a)-(7c); parameters eq. (17).
Usage: python src/linear_stability.py > results/MG001/linear_stability.txt
"""
import numpy as np
from scipy.optimize import brentq

P17 = dict(a=14.0, b=29.0, g=0.1, p=1.0, q=1.0, s=0.0, u=0.0, w=0.0, dx=1.0, dy=12.0)


def steady_state(a, b, g, p, q, s, u, w, **_):
    G0 = (a + g * w) / (q - g * p)
    X0 = (p * a + q * w) / (q - g * p)
    Y0 = X0 * (s * X0**2 + b) / (X0**2 + u)
    return G0, X0, Y0


def rhs(G, X, Y, a, b, g, p, q, s, u, w, **_):
    return (-q * G + g * X + a,
            p * G - (1 + b) * X + u * Y - s * X**3 + X**2 * Y + w,
            b * X - u * Y + s * X**3 - X**2 * Y)


def jacobian(P):
    G0, X0, Y0 = steady_state(**P)
    b, g, p, q, s, u = (P[k] for k in "bgpqsu")
    return np.array([[-q, g, 0.0],
                     [p, -(1 + b) - 3 * s * X0**2 + 2 * X0 * Y0, u + X0**2],
                     [0.0, b + 3 * s * X0**2 - 2 * X0 * Y0, -u - X0**2]])


def dispersion(P, ks):
    J, D = jacobian(P), np.diag([1.0, P["dx"], P["dy"]])
    return np.array([np.linalg.eigvals(J - k * k * D).real.max() for k in ks])


def most_unstable(P, ks=np.linspace(0, 6, 60001)):
    sig = dispersion(P, ks)
    i = int(np.argmax(sig))
    return sig[i], ks[i]


if __name__ == "__main__":
    P = dict(P17)
    G0, X0, Y0 = steady_state(**P)
    print("Parameters (eq. 17):", P)
    print(f"Steady state: G0={G0:.6f} X0={X0:.6f} Y0={Y0:.6f}")
    print("RHS residual at steady state:", rhs(G0, X0, Y0, **P))
    sig, k = most_unstable(P)
    print(f"sigma(k=0) = {dispersion(P, [0.0])[0]:.6f}")
    print(f"Least-stable mode: k*={k:.4f}, lambda*=2pi/k*={2*np.pi/k:.4f}, max Re sigma={sig:.6f}")
    print("Homogeneous state is", "LINEARLY STABLE (subcritical)" if sig < 0 else "UNSTABLE")
    bc = brentq(lambda bb: most_unstable({**P, "b": bb})[0], 29.0, 32.0, xtol=1e-6)
    sc, kc = most_unstable({**P, "b": bc})
    print(f"Turing threshold b_c = {bc:.5f} (k_c={kc:.4f}, lambda_c={2*np.pi/kc:.4f}); "
          f"published b=29 is {100*(bc-29)/bc:.2f}% below threshold")
    print("Reported nonlinear soliton wavelength (Table 1): 3.08")

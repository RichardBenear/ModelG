"""MG192 separate check: mobility by co-moving continuation (no adjoint, no analytic Jacobian).

Solves  F(U) + kappa * x * W(x) * e_S + v * dU/dx = 0,  <U'_0, U - U_0> = 0
for (U, v) by Newton with a finite-difference Jacobian, starting from the lattice reference, and reports
dv/dkappa by central differences at two kappa sizes. W is a smooth window equal to 1 for |x| < 30.
Residual code is written independently of src/mg192_mobility.py; it shares only the model equations,
the reference profile file and NumPy/SciPy.
"""
import json
import os
import time

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "results", "MG192_gravity_mobility")
L, N = 100.0, 1024
DX, DY, A, B, GG = 1.0, 12.0, 14.0, 29.0, 0.1

h = L / N
x = -L / 2 + h * np.arange(N)
k = 2 * np.pi * np.fft.fftfreq(N, d=h)
kd = k.copy()
kd[N // 2] = 0.0
W = 0.5 * (np.tanh((x + 38.0) / 2.0) - np.tanh((x - 38.0) / 2.0))


def lap(f):
    return np.fft.irfft(-(np.fft.rfftfreq(N, d=h) * 2 * np.pi) ** 2 * np.fft.rfft(f), n=N)


def d1(f):
    return np.real(np.fft.ifft(1j * kd * np.fft.fft(f)))


def F(z, kap, species, U0p):
    G, X, Y = z[:N], z[N:2 * N], z[2 * N:3 * N]
    v = z[-1]
    src = kap * x * W
    eG = lap(G) - G + GG * X + A + v * d1(G)
    eX = DX * lap(X) + G - (1 + B) * X + X * X * Y + v * d1(X)
    eY = DY * lap(Y) + B * X - X * X * Y + v * d1(Y)
    if species == "G":
        eG = eG + src
    elif species == "X":
        eX = eX + src
    elif species == "b":          # perturbation of b: -db*X in X eq, +db*X in Y eq
        eX = eX - src * X
        eY = eY + src * X
    if U0p is None:
        phase = np.array([0.0])
    else:
        phase = np.array([np.dot(U0p, z[:3 * N] - Zref)])
    return np.concatenate([eG, eX, eY, phase])


def newton(z, kap, species, U0p, tol=1e-11):
    for it in range(30):
        r = F(z, kap, species, U0p)
        if np.max(np.abs(r)) < tol:
            return z, float(np.max(np.abs(r))), it
        J = np.empty((3 * N + 1, 3 * N + 1))
        r0 = r
        for j in range(3 * N + 1):
            e = 1e-7 * max(1.0, abs(z[j]))
            zp = z.copy()
            zp[j] += e
            J[:, j] = (F(zp, kap, species, U0p) - r0) / e
        z = z - np.linalg.solve(J, r0)
    return z, float(np.max(np.abs(F(z, kap, species, U0p)))), it


if __name__ == "__main__":
    t0 = time.time()
    d = np.load(os.path.join(ROOT, "results", "MG003", "ref_profiles.npz"))["h0.2_b29.0"]
    xr, phi = d[0], d[1:4]
    c = xr[np.argmin(phi[1])]
    xs = (xr - c + 50) % 100 - 50
    o = np.argsort(xs)
    h0 = np.array([140 / 9, 140 / 9, B / (140 / 9)])
    U = np.concatenate([h0[s] + np.interp(x, np.concatenate([xs[o] - 100, xs[o], xs[o] + 100]),
                                          np.concatenate([phi[s][o]] * 3)) for s in range(3)])
    # kappa = 0 soliton: fix translation with phase condition against the guess's derivative, v unknown
    Zref = U.copy()
    U0p = np.concatenate([d1(U[s * N:(s + 1) * N]) for s in range(3)])
    z0, r0, it0 = newton(np.concatenate([U, [0.0]]), 0.0, "G", U0p)
    Zref = z0[:3 * N].copy()
    U0p = np.concatenate([d1(Zref[s * N:(s + 1) * N]) for s in range(3)])
    out = dict(kappa0=dict(residual=r0, iterations=it0, v=float(z0[-1])), mobility={})
    for sp in ("G", "X", "b"):
        res = {}
        for kap in (1e-4, 2e-4):
            zp, rp, _ = newton(z0.copy(), kap, sp, U0p)
            zm, rm, _ = newton(z0.copy(), -kap, sp, U0p)
            # moving-frame equation uses +v U' with u(x - s(t)) convention: lab velocity = -v_frame? check sign:
            # u(x,t) = U(x - c t) gives u_t = -c U'; steady frame equation F + P + c U' = 0, so c = v.
            res[str(kap)] = dict(dvdk=float((zp[-1] - zm[-1]) / (2 * kap)), residuals=[rp, rm])
        out["mobility"][sp] = res
        print(sp, res)
    out["wall_seconds"] = time.time() - t0
    with open(os.path.join(OUT, "check_continuation.json"), "w") as f:
        json.dump(out, f, indent=1)
    print(json.dumps(out, indent=1))

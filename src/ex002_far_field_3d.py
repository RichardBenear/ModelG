"""EX002 (exploratory): 3D spherically symmetric Model G particle, net G integral and far-field decay.

Conservative finite-volume radial discretization: cell i has centre r_i = (i + 1/2) dr and exact shell volume
V_i = 4 pi (r_{i+1/2}^3 - r_{i-1/2}^3) / 3, so sum_i V_i Lap(u)_i equals the outer boundary flux exactly.
Outer boundary: 'neumann' (zero flux) or 'dirichlet' (u = homogeneous value at r = R).
Optional broadened core (Pulver & LaViolette 2013, eqs. 24-25): Laplacian multiplied by delta(r) = 1 + exp(-4 r^2/9).
"""
import json
import os
import sys
import time

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "results", "EX002_far_field_3d")
A, G_, P_, Q_, DX, DY = 14.0, 0.1, 1.0, 1.0, 1.0, 12.0


def hom(a, b):
    X0 = a / (Q_ - G_ * P_)
    return np.array([X0, X0, b / X0])


class Radial:
    def __init__(self, R, dr, bc, b, delta=False):
        self.R, self.dr, self.bc, self.b = R, dr, bc, b
        n = int(round(R / dr))
        self.n = n
        self.r = (np.arange(n) + 0.5) * dr
        rf = np.arange(n + 1) * dr
        self.V = 4 * np.pi * (rf[1:] ** 3 - rf[:-1] ** 3) / 3
        Af = 4 * np.pi * rf ** 2
        dfac = (1 + np.exp(-4 * self.r ** 2 / 9)) if delta else np.ones(n)
        main = np.zeros(n)
        up = np.zeros(n - 1)
        lo = np.zeros(n - 1)
        for i in range(n):
            if i < n - 1:
                c = Af[i + 1] / dr
                main[i] -= c
                up[i] += c
            if i > 0:
                c = Af[i] / dr
                main[i] -= c
                lo[i - 1] += c
        self.bvec = np.zeros(n)
        if bc == "dirichlet":
            c = Af[n] / (dr / 2)
            main[n - 1] -= c
            self.bvec[n - 1] = c      # times the boundary value
        Lap = sp.diags([lo, main, up], [-1, 0, 1]) / self.V[:, None] if False else \
            sp.diags(1 / self.V) @ sp.diags([lo, main, up], [-1, 0, 1])
        self.Lap = (sp.diags(dfac) @ Lap).tocsr()
        self.bterm = dfac * self.bvec / self.V
        self.h0 = hom(A, b)

    def lap(self, f, s):
        return self.Lap @ f + self.bterm * self.h0[s]

    def F(self, u):
        G, X, Y = u
        b = self.b
        return np.array([self.lap(G, 0) - Q_ * G + G_ * X + A,
                         DX * self.lap(X, 1) + P_ * G - (1 + b) * X + X * X * Y,
                         DY * self.lap(Y, 2) + b * X - X * X * Y])

    def J(self, u):
        G, X, Y = u
        b, n = self.b, self.n
        I = sp.identity(n)
        return sp.bmat([[self.Lap - Q_ * I, G_ * I, None],
                        [P_ * I, DX * self.Lap + sp.diags(-(1 + b) + 2 * X * Y), sp.diags(X * X)],
                        [None, sp.diags(b - 2 * X * Y), DY * self.Lap - sp.diags(X * X)]]).tocsc()

    def newton(self, u, tol=1e-10, maxit=30):
        for it in range(maxit):
            Fv = self.F(u).ravel()
            res = np.max(np.abs(Fv))
            if res < tol:
                return u, res, it, True
            du = spla.spsolve(self.J(u), -Fv)
            u = u + du.reshape(3, self.n)
            if np.any(u <= 0) or not np.all(np.isfinite(u)):
                return u, np.inf, it, False
        return u, float(np.max(np.abs(self.F(u)))), maxit, bool(np.max(np.abs(self.F(u))) < tol)

    def settle(self, u, dt=0.2, T=400.0):
        """Linearly implicit Euler in time, to relax the seed toward a stationary particle."""
        I = sp.identity(3 * self.n, format="csc")
        for k in range(int(T / dt)):
            Fv = self.F(u).ravel()
            du = spla.spsolve(I - dt * self.J(u), dt * Fv)
            u = u + du.reshape(3, self.n)
        return u

    def integral(self, f):
        return float(np.sum(self.V * f))


def seed(rad):
    """1D MG192 soliton profile (b = 29), mapped radially as a seed."""
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    d = np.load(os.path.join(ROOT, "results", "MG192_gravity_mobility", "soliton_adjoint_L100_N1024.npz"))["U"]
    x = -50 + (100 / 1024) * np.arange(1024)
    h1 = hom(14.0, 29.0)
    u = np.empty((3, rad.n))
    # the raw 1D profile decays in 3D; amplitude x4 and radial stretch x1.5 nucleate a 3D particle (found by trial,
    # seeds x2 / x4 unstretched decayed); positivity clipped at 0.05
    for s in range(3):
        u[s] = rad.h0[s] + 4.0 * np.interp(rad.r / 1.5, x[512:], d[s, 512:] - h1[s], right=0.0)
    return np.maximum(u, 0.05)


def tail_fit(rad, f, lo, hi):
    """Fit r*phi = e^{-kap r}(A cos kr + B sin kr) + C e^{-kap2 r} on [lo, hi] (least squares, weighted)."""
    from scipy.optimize import least_squares
    m = (rad.r >= lo) & (rad.r <= hi)
    r, y = rad.r[m], rad.r[m] * f[m]
    w = 1 / (np.abs(y).max() * np.exp(-0.46 * (r - lo)) + 1e-300)

    def mod(p):
        return np.exp(-p[0] * r) * (p[2] * np.cos(p[1] * r) + p[3] * np.sin(p[1] * r)) + p[4] * np.exp(-p[5] * r)

    best = None
    for ph in range(4):
        p0 = [0.46, 2.05, np.cos(ph * np.pi / 2) * abs(y).max() * np.exp(0.46 * lo),
              np.sin(ph * np.pi / 2) * abs(y).max() * np.exp(0.46 * lo), 0.0, 0.97]
        rr = least_squares(lambda p: (mod(p) - y) * w, p0, x_scale="jac", max_nfev=4000)
        if best is None or rr.cost < best.cost:
            best = rr
    p = best.x
    return dict(kappa=float(p[0]), k=float(p[1]), C=float(p[4]), kappa2=float(p[5]),
                rel_rms=float(np.sqrt(np.mean(((mod(p) - y) * w) ** 2)) / np.sqrt(np.mean((y * w) ** 2))))


def linear_rates(a, b):
    X0 = a / 0.9
    Y0 = b / X0
    J = np.array([[-1, 0.1, 0], [1, -(1 + b) + 2 * X0 * Y0, X0 ** 2], [0, b - 2 * X0 * Y0, -X0 ** 2]])
    D = np.diag([1, 1, 12.0])
    ss = np.arange(4)
    c = np.polyfit(ss, [np.linalg.det(J - D * s) for s in ss], 3)
    out = []
    for s in np.roots(c):
        k = np.sqrt(complex(s))
        out.append(dict(decay=abs(k.imag), wavenumber=abs(k.real)))
    return out


def main():
    os.makedirs(OUT, exist_ok=True)
    t0 = time.time()
    b = 29.0
    rec = dict(a=A, b=b, linear_rates=linear_rates(A, b), cases=[])
    print("linear rates", rec["linear_rates"], flush=True)
    # base solution: R = 40, Neumann, dr = 0.02, settled then Newton
    base = Radial(40.0, 0.02, "neumann", b)
    u = base.settle(seed(base))
    u, res, it, ok = base.newton(u)
    phiX0 = u[1, 0] - base.h0[1]
    print(f"base: ok {ok} res {res:.2e} it {it} core phi_X {phiX0:.4f} phi_G {u[0, 0] - base.h0[0]:.4f}", flush=True)
    if not ok or phiX0 > -1:
        rec["base_failed"] = dict(ok=ok, res=res, phiX0=phiX0)
        json.dump(rec, open(os.path.join(OUT, "results.json"), "w"), indent=1, default=float)
        return
    np.savez_compressed(os.path.join(OUT, "profile_R40_neumann.npz"), r=base.r, u=u)
    r_base, u_base = base.r, u
    for delta in (False, True):
        for R in (20.0, 30.0, 40.0, 60.0):
            for bc in ("neumann", "dirichlet"):
                for dr in ((0.02, 0.01) if (R == 40.0 and not delta) else (0.02,)):
                    rad = Radial(R, dr, bc, b, delta)
                    ui = np.array([np.interp(rad.r, r_base, u_base[s], right=rad.h0[s]) for s in range(3)])
                    if delta:
                        ui = rad.settle(ui, T=200.0)
                    ui, res, it, ok = rad.newton(ui)
                    phi = ui - rad.h0[:, None]
                    case = dict(delta=delta, R=R, bc=bc, dr=dr, ok=ok, residual=res,
                                core_phiX=float(phi[1, 0]), core_phiG=float(phi[0, 0]),
                                int_G=rad.integral(phi[0]), int_X=rad.integral(phi[1]), int_Y=rad.integral(phi[2]),
                                int_negG_core=rad.integral(np.minimum(phi[0], 0) * (rad.r < 3)),
                                rms_radius_G=float(np.sqrt(rad.integral(rad.r ** 2 * phi[0] ** 2) /
                                                           rad.integral(phi[0] ** 2))))
                    if R == 60.0 and ok:
                        case["tail_G"] = tail_fit(rad, phi[0], 8.0, 40.0)
                        case["tail_X"] = tail_fit(rad, phi[1], 8.0, 40.0)
                        np.savez_compressed(os.path.join(OUT, f"profile_R60_{bc}_{'delta' if delta else 'std'}.npz"),
                                            r=rad.r, u=ui)
                    rec["cases"].append(case)
                    print({k: (round(v, 6) if isinstance(v, float) else v) for k, v in case.items()
                           if not k.startswith("tail")}, flush=True)
                    if "tail_G" in case:
                        print("   tail G", case["tail_G"], "\n   tail X", case["tail_X"], flush=True)
                    json.dump(rec, open(os.path.join(OUT, "results.json"), "w"), indent=1, default=float)
    rec["wall_seconds"] = time.time() - t0
    json.dump(rec, open(os.path.join(OUT, "results.json"), "w"), indent=1, default=float)
    print("done", rec["wall_seconds"], flush=True)


if __name__ == "__main__":
    main()

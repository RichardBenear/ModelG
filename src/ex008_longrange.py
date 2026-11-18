"""EX008 (exploratory): designed test model for a long-range (1/r) G field. Our construction, not LaViolette's.

  X_t = Lap X + A - (1+b) X + X^2 Y + cx k G W
  Y_t = 12 Lap Y + b X - X^2 Y + cy k G W
  W_t = DW Lap W - W + c Y W^2 - W^3          (core-only species: W = 0 is the ether state)
  G_t = Lap G - k G W,  G(R) = G_inf          (no source or loss in the ether; supplied from far away)
Stationary particle: 3D spherically symmetric, conservative finite volumes (as EX002). Drift: l = 1 (dipole) linear
response to an imposed far gradient G -> G_inf + E z, v = -l^T P / l^T r1 with l, r1 the left/right null vectors of the
l = 1 Jacobian (r1 ~ U'(r), the translation mode).
"""
import json
import os
import sys
import time

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "results", "EX008_longrange")
A, B, DY, DW, C, GINF = 140 / 9, 29.0, 12.0, 0.1, 1.0, 1.0
DIFF = np.array([1.0, DY, DW, 1.0])          # X, Y, W, G


class Model:
    def __init__(self, R, dr, k, cx, cy):
        self.R, self.dr, self.k, self.cx, self.cy = R, dr, k, cx, cy
        n = int(round(R / dr))
        self.n = n
        self.r = (np.arange(n) + 0.5) * dr
        rf = np.arange(n + 1) * dr
        self.V = 4 * np.pi * (rf[1:] ** 3 - rf[:-1] ** 3) / 3
        Af = 4 * np.pi * rf ** 2
        main, up = np.zeros(n), np.zeros(n - 1)
        for i in range(n - 1):
            c = Af[i + 1] / dr
            main[i] -= c
            main[i + 1] -= c
            up[i] = c
        Kn = sp.diags([up, main, up], [-1, 0, 1])
        self.LapN = (sp.diags(1 / self.V) @ Kn).tocsr()
        cD = Af[n] / (dr / 2)
        mainD = main.copy()
        mainD[-1] -= cD
        self.LapD = (sp.diags(1 / self.V) @ sp.diags([up, mainD, up], [-1, 0, 1])).tocsr()
        self.bD = np.zeros(n)
        self.bD[-1] = cD / self.V[-1]
        self.h0 = np.array([A, B / A, 0.0, GINF])

    def F(self, u):
        X, Y, W, G = u
        s = self.k * G * W
        return np.array([self.LapN @ X + A - (1 + B) * X + X * X * Y + self.cx * s,
                         DY * (self.LapN @ Y) + B * X - X * X * Y + self.cy * s,
                         DW * (self.LapN @ W) - W + C * Y * W * W - W ** 3,
                         self.LapD @ G + self.bD * GINF - s])

    def J(self, u, l=0):
        X, Y, W, G = u
        k, n = self.k, self.n
        d = sp.diags
        if l == 0:
            LX, LY, LW, LG = self.LapN, DY * self.LapN, DW * self.LapN, self.LapD
        else:   # l = 1: Dirichlet-zero outer boundary for all, plus -2/r^2
            cen = d(-2.0 / self.r ** 2)
            LX, LY, LW, LG = self.LapD + cen, DY * (self.LapD + cen), DW * (self.LapD + cen), self.LapD + cen
        return sp.bmat([
            [LX + d(-(1 + B) + 2 * X * Y), d(X * X), d(self.cx * k * G), d(self.cx * k * W)],
            [d(B - 2 * X * Y), LY - d(X * X), d(self.cy * k * G), d(self.cy * k * W)],
            [None, d(C * W * W), LW + d(-1 + 2 * C * Y * W - 3 * W * W), None],
            [None, None, d(-k * G), LG + d(-k * W)]]).tocsc()

    def newton(self, u, tol=1e-10, maxit=40):
        for it in range(maxit):
            Fv = self.F(u).ravel()
            res = float(np.max(np.abs(Fv)))
            if res < tol:
                return u, res, True
            u = u + spla.spsolve(self.J(u), -Fv).reshape(4, self.n)
            if not np.all(np.isfinite(u)) or np.any(u[[0, 1, 3]] <= 0):
                return u, np.inf, False
        res = float(np.max(np.abs(self.F(u))))
        return u, res, res < tol

    def settle(self, u, dt=0.1, T=100.0):
        I = sp.identity(4 * self.n, format="csc")
        for _ in range(int(T / dt)):
            u = u + spla.spsolve(I - dt * self.J(u), dt * self.F(u).ravel()).reshape(4, self.n)
        return u

    def integral(self, f):
        return float(np.sum(self.V * f))

    def drift(self, U):
        """Mobility mu = v / E for an imposed far gradient G ~ G_inf + E z (v > 0: toward higher G)."""
        J1 = self.J(U, l=1)
        Up = np.array([np.gradient(U[s], self.r) for s in range(4)]).ravel()
        lu = spla.splu(J1)
        luT = spla.splu(J1.T.tocsc())
        rng = np.random.default_rng(0)
        rv, lv = Up.copy(), rng.standard_normal(J1.shape[0])
        lam_r = []
        for _ in range(8):               # inverse iteration for the eigenvalue nearest 0
            rv = lu.solve(rv)
            rv /= np.linalg.norm(rv)
            lv = luT.solve(lv)
            lv /= np.linalg.norm(lv)
        lam = float(rv @ (J1 @ rv))
        scale = (Up @ Up) / (rv @ Up)    # make rv ~ U' in amplitude
        rv = rv * scale
        X, Y, W, G = U
        # forcing from the harmonic part E r cos(theta) of the G field (per unit E): -k W r in G, +c k W r in X/Y
        P = np.zeros((4, self.n))
        P[0] = self.cx * self.k * W * self.r
        P[1] = self.cy * self.k * W * self.r
        P[3] = -self.k * W * self.r
        P = P.ravel()
        # J1 du + v rv + P = 0, with du the response; projecting on the left null vector:
        v = -float(lv @ P) / float(lv @ rv)
        cos_rU = float(rv @ Up / np.linalg.norm(rv) / np.linalg.norm(Up))
        return dict(mu=v, lambda_l1=lam, cos_translation=cos_rU)


def seed(model):
    d = np.load(os.path.join(ROOT, "results", "EX002_far_field_3d", "profile_R60_neumann_std.npz"))
    r0, u0 = d["r"], d["u"]
    h = np.array([140 / 9, 140 / 9, 29 / (140 / 9)])
    u = np.empty((4, model.n))
    u[0] = A + np.interp(model.r, r0, u0[1] - h[1], right=0.0)
    u[1] = B / A + np.interp(model.r, r0, u0[2] - h[2], right=0.0)
    u[2] = 3.0 * np.exp(-(model.r / 1.5) ** 2)
    u[3] = GINF
    return u


def leading_eigs(model, U, l=0, k=4):
    J = model.J(U, l=l)
    lu = spla.splu((J - 0.3 * sp.identity(J.shape[0])).tocsc())
    op = spla.LinearOperator(J.shape, matvec=lu.solve)
    mu = spla.eigs(op, k=k, which="LM", return_eigenvectors=False)
    return sorted([[float((0.3 + 1 / z).real), float((0.3 + 1 / z).imag)] for z in mu], reverse=True)


def particle(R, dr, k, cx, cy, prev=None):
    M = Model(R, dr, k, cx, cy)
    if prev is None:
        u = M.settle(seed(M), dt=0.05, T=20.0)
        ks = [k]
        u0 = u
        steps = np.linspace(0, k, 6)[1:] if k > 0 else []
        for kk in steps:
            Mk = Model(R, dr, kk, cx, cy)
            u0, res, ok = Mk.newton(u0)
            assert ok and u0[2, 0] > 1, (kk, res)
        u = u0
    else:
        rp, up = prev
        u = np.array([np.interp(M.r, rp, up[s]) for s in range(4)])
        u[3] = np.interp(M.r, rp, up[3] - GINF, right=0.0) + GINF
    u, res, ok = M.newton(u)
    return M, u, res, ok


def analyse(M, u):
    X, Y, W, G = u
    cons = M.integral(M.k * G * W)
    r = M.r
    m = (r > 15) & (r < M.R - 5)
    pred = cons / (4 * np.pi) * (1 / r[m] - 1 / M.R)
    dG = GINF - G[m]
    far = dict(rel_dev_from_point_sink=float(np.max(np.abs(dG - pred)) / np.max(np.abs(pred))),
               r_times_deficit_at_20=float(np.interp(20, r, r * (GINF - G))),
               M_over_4pi=cons / (4 * np.pi))
    return dict(core_X=float(X[0] - A), core_Y=float(Y[0]), core_W=float(W[0]), core_G=float(G[0]),
                consumption_M=cons, int_X=M.integral(X - A), int_Y=M.integral(Y - B / A),
                W_radius_rms=float(np.sqrt(M.integral(r ** 2 * W) / M.integral(W))), far=far)


def main():
    os.makedirs(OUT, exist_ok=True)
    t0 = time.time()
    rec = dict(note="EXPLORATORY", params=dict(A=A, b=B, DY=DY, DW=DW, c=C, G_inf=GINF), runs=[])
    for cx, cy, name in ((1, 0, "G feeds X"), (0, 1, "G feeds Y")):
        for k in (0.05, 0.1, 0.2, 0.5):
            for R in (60.0, 120.0):
                try:
                    M, u, res, ok = particle(R, 0.05, k, cx, cy)
                except AssertionError as ex:
                    rec["runs"].append(dict(variant=name, k=k, R=R, ok=False, why=str(ex)))
                    print(name, k, R, "FAILED", ex, flush=True)
                    continue
                run = dict(variant=name, cx=cx, cy=cy, k=k, R=R, ok=bool(ok), residual=res)
                if ok and u[2, 0] > 1:
                    run.update(analyse(M, u))
                    run["eigs_l0"] = leading_eigs(M, u, 0)
                    run["eigs_l1"] = leading_eigs(M, u, 1)
                    run["drift"] = M.drift(u)
                    run["pair_force_coeff"] = run["drift"]["mu"] * run["consumption_M"] / (4 * np.pi)
                    np.savez_compressed(os.path.join(OUT, f"particle_{'GX' if cx else 'GY'}_k{k:g}_R{R:g}.npz"),
                                        r=M.r, u=u)
                print({kk: v for kk, v in run.items()}, flush=True)
                rec["runs"].append(run)
                json.dump(rec, open(os.path.join(OUT, "results.json"), "w"), indent=1, default=float)
    rec["wall_seconds"] = time.time() - t0
    json.dump(rec, open(os.path.join(OUT, "results.json"), "w"), indent=1, default=float)
    print("done", rec["wall_seconds"])


if __name__ == "__main__":
    main()

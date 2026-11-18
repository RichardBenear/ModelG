"""EX009 part A: direct two-particle check of the EX008 1/d^2 attraction, axisymmetric (rho, z).

Two identical particles at z = +-d/2 on the axis; half domain z >= 0 with mirror symmetry at z = 0.
Co-moving steady state of the right particle: F(u) + v dz u = 0, with a phase condition fixing its position;
v (dz/dt of the right particle) comes out of a bordered Newton solve. v < 0 means approach.
Grid: uniform h near the particle, geometric stretching to far boundaries (G fixed there, others zero flux).
"""
import json
import os
import sys
import time

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ex008_longrange as e  # noqa: E402

OUT = os.path.join(e.ROOT, "results", "EX009_followups")
K, CX, CY = 0.5, 1, 0


def faces(core, h, far, grow=1.1):
    f = list(np.arange(0.0, core + 1e-9, h))
    w = h
    while f[-1] < far:
        w *= grow
        f.append(f[-1] + w)
    return np.array(f)


class Axi:
    def __init__(self, d, h=0.1, rho_core=8.0, pad=8.0, far=150.0, single=False):
        self.d = d
        zc = 0.0 if single else d / 2
        self.zp = zc
        rf = faces(rho_core, h, far)
        zf = faces(zc + pad, h, zc + far)
        self.rf, self.zf = rf, zf
        rc, zcen = 0.5 * (rf[1:] + rf[:-1]), 0.5 * (zf[1:] + zf[:-1])
        self.nr, self.nz = len(rc), len(zcen)
        self.R, self.Z = np.meshgrid(rc, zcen, indexing="ij")      # index [i_rho, j_z]
        dr, dz = np.diff(rf), np.diff(zf)
        ring = np.pi * (rf[1:] ** 2 - rf[:-1] ** 2)
        self.V = (ring[:, None] * dz[None, :]).ravel()
        n = self.nr * self.nz
        self.n = n
        idx = np.arange(n).reshape(self.nr, self.nz)
        rows, cols, vals = [], [], []
        diagN = np.zeros(n)
        diagD = np.zeros(n)
        bD = np.zeros(n)
        # rho faces (interior)
        for i in range(self.nr - 1):
            A = 2 * np.pi * rf[i + 1] * dz
            c = A / (rc[i + 1] - rc[i])
            a, b = idx[i], idx[i + 1]
            rows += [a, b]; cols += [b, a]; vals += [c, c]
            diagN[a] -= c; diagN[b] -= c
        # z faces (interior)
        for j in range(self.nz - 1):
            c = ring / (zcen[j + 1] - zcen[j])
            a, b = idx[:, j], idx[:, j + 1]
            rows += [a, b]; cols += [b, a]; vals += [c, c]
            diagN[a] -= c; diagN[b] -= c
        diagD[:] = diagN
        # outer boundaries: Dirichlet for G only
        cr = 2 * np.pi * rf[-1] * dz / (dr[-1] / 2)
        diagD[idx[-1]] -= cr
        bD[idx[-1]] += cr
        cz = ring / (dz[-1] / 2)
        diagD[idx[:, -1]] -= cz
        bD[idx[:, -1]] += cz
        rows, cols, vals = np.concatenate(rows), np.concatenate(cols), np.concatenate(vals)
        K0 = sp.csr_matrix((vals, (rows, cols)), shape=(n, n))
        iV = sp.diags(1 / self.V)
        self.LN = (iV @ (K0 + sp.diags(diagN))).tocsr()
        self.LD = (iV @ (K0 + sp.diags(diagD))).tocsr()
        self.bD = bD / self.V
        # d/dz, central on the nonuniform grid; mirror at z = 0 (value at -z1 equals z1), one-sided at the far end
        rr, cc, vv = [], [], []
        for j in range(self.nz):
            a = idx[:, j]
            if j == 0:
                continue                      # derivative of an even function at the first cell ~ 0
            if j == self.nz - 1:
                h_ = zcen[j] - zcen[j - 1]
                rr += [a, a]; cc += [a, idx[:, j - 1]]; vv += [np.full(self.nr, 1 / h_), np.full(self.nr, -1 / h_)]
                continue
            h_ = zcen[j + 1] - zcen[j - 1]
            rr += [a, a]; cc += [idx[:, j + 1], idx[:, j - 1]]; vv += [np.full(self.nr, 1 / h_), np.full(self.nr, -1 / h_)]
        self.Dz = sp.csr_matrix((np.concatenate(vv), (np.concatenate(rr), np.concatenate(cc))), shape=(n, n))

    def F(self, u):
        X, Y, W, G = u
        s = K * G * W
        return np.array([self.LN @ X + e.A - (1 + e.B) * X + X * X * Y + CX * s,
                         e.DY * (self.LN @ Y) + e.B * X - X * X * Y + CY * s,
                         e.DW * (self.LN @ W) - W + e.C * Y * W * W - W ** 3,
                         self.LD @ G + self.bD * e.GINF - s])

    def J(self, u):
        X, Y, W, G = u
        d = sp.diags
        return sp.bmat([
            [self.LN + d(-(1 + e.B) + 2 * X * Y), d(X * X), d(CX * K * G), d(CX * K * W)],
            [d(e.B - 2 * X * Y), e.DY * self.LN - d(X * X), d(CY * K * G), d(CY * K * W)],
            [None, d(e.C * W * W), e.DW * self.LN + d(-1 + 2 * e.C * Y * W - 3 * W * W), None],
            [None, None, d(-K * G), self.LD + d(-K * W)]]).tocsc()

    def guess(self, prof_r, prof_u, M):
        """Superpose radial profiles of the two particles (mirror image at -zp); G deficit continued as M/(4 pi r)."""
        out = np.empty((4, self.n))
        r1 = np.sqrt(self.R ** 2 + (self.Z - self.zp) ** 2).ravel()
        r2 = np.sqrt(self.R ** 2 + (self.Z + self.zp) ** 2).ravel()
        h0 = np.array([e.A, e.B / e.A, 0.0, e.GINF])
        def dev(s, r):
            if s < 3:
                return np.interp(r, prof_r, prof_u[s] - h0[s], right=0.0)
            near = np.interp(r, prof_r, prof_u[3] - h0[3])
            return np.where(r < 15, near, -M / (4 * np.pi * np.maximum(r, 1e-9)))
        for s in range(4):
            out[s] = h0[s] + dev(s, r1) + (dev(s, r2) if self.zp > 0 else 0.0)
        return out

    def solve(self, u, phase_ref=None, tol=1e-9, maxit=25, moving=True):
        n = self.n
        Dz4 = sp.block_diag([self.Dz] * 4).tocsr()
        ref = u.copy() if phase_ref is None else phase_ref
        phi = (Dz4 @ ref.ravel()) * np.repeat([1.0, 0, 0, 0], n)     # pin the X pattern
        v = 0.0
        for it in range(maxit):
            uf = u.ravel()
            Fv = self.F(u).ravel() + (v * (Dz4 @ uf) if moving else 0)
            g = float(phi @ (uf - ref.ravel()))
            res = float(np.max(np.abs(Fv)))
            if res < tol and abs(g) < 1e-9:
                return u, v, res, True, it
            Jm = self.J(u) + (v * Dz4 if moving else 0)
            if moving:
                col = (Dz4 @ uf)[:, None]
                Jb = sp.bmat([[Jm, sp.csr_matrix(col)], [sp.csr_matrix(phi[None, :]), None]]).tocsc()
                rhs = -np.concatenate([Fv, [g]])
                sol = spla.spsolve(Jb, rhs)
                u = u + sol[:-1].reshape(4, n)
                v = v + sol[-1]
            else:
                u = u + spla.spsolve(Jm.tocsc(), -Fv).reshape(4, n)
        res = float(np.max(np.abs(self.F(u).ravel() + (v * (Dz4 @ u.ravel()) if moving else 0))))
        return u, v, res, res < tol, maxit

    def integral(self, f):
        return float(np.sum(self.V * f))


def main():
    t0 = time.time()
    d = np.load(os.path.join(e.OUT, "particle_GX_k0.5_R120.npz"))
    pr, pu = d["r"], d["u"]
    rec = json.load(open(os.path.join(e.OUT, "results.json")))
    run = [x for x in rec["runs"] if x["variant"] == "G feeds X" and x["k"] == 0.5 and x["R"] == 120.0][0]
    Mr, mu = run["consumption_M"], run["drift"]["mu"]
    out = dict(note="EXPLORATORY", K=K, M_radial=Mr, mu_radial=mu, single=None, pairs=[])
    # single particle on the 2D grid (centred at z = 0): resolution check of M against the radial value
    a = Axi(0.0, single=True)
    u, v, res, ok, it = a.solve(a.guess(pr, pu, Mr), moving=False)
    M2 = a.integral(K * u[3] * u[2])
    out["single"] = dict(ok=bool(ok), residual=res, M=M2, M_rel_to_radial=M2 / Mr - 1, n=a.n)
    print("single", out["single"], round(time.time() - t0), flush=True)
    for dd in (8.0, 10.0, 12.0, 15.0, 20.0, 25.0, 30.0):
        a = Axi(dd)
        u0 = a.guess(pr, pu, Mr)
        u, v, res, ok, it = a.solve(u0)
        pred = mu * Mr / (4 * np.pi * dd ** 2)
        # local-G correction: partner lowers the ambient G by M/(4 pi d); first order, M and mu scale ~ with ambient G
        r_ = dict(d=dd, ok=bool(ok), residual=res, iters=it, v=v, v_pred=pred, ratio=v / pred, n=a.n,
                  M_each=a.integral(K * u[3] * u[2]) / 2, wall=round(time.time() - t0))
        print(r_, flush=True)
        out["pairs"].append(r_)
        json.dump(out, open(os.path.join(OUT, "partA_pair2d.json"), "w"), indent=1, default=float)
    # grid check at d = 12
    a = Axi(12.0, h=0.075)
    u, v, res, ok, it = a.solve(a.guess(pr, pu, Mr))
    out["grid_check"] = dict(d=12.0, h=0.075, ok=bool(ok), v=v, v_pred=mu * Mr / (4 * np.pi * 144))
    print("grid", out["grid_check"], flush=True)
    a = Axi(12.0, far=300.0)
    u, v, res, ok, it = a.solve(a.guess(pr, pu, Mr))
    out["far_check"] = dict(d=12.0, far=300.0, ok=bool(ok), v=v)
    print("far", out["far_check"], flush=True)
    out["wall_seconds"] = time.time() - t0
    json.dump(out, open(os.path.join(OUT, "partA_pair2d.json"), "w"), indent=1, default=float)


if __name__ == "__main__":
    main()

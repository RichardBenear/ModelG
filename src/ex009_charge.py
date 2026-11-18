"""EX009 part C: electric analogue design scan (radial 3D, one particle at a time).

Same X/Y/W particle as EX008, with a conserved signed quantity C (no ether source or loss, C = C_inf at R).
Core exchange Z + W <-> C + W with rate R = k W (a Z - b C), Z = X or Y:
  "+" particle: (a, b) = (1, eps)  -> net C source;   "-" particle: (a, b) = (eps, 1) -> net C sink.
  C_t = Lap C + R;  Z_t = ... - R.
Like charges repel and unlike attract if sources drift down the C gradient (mu_C < 0) and sinks up it (mu_C > 0).
"""
import json
import os
import sys

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ex008_longrange as e  # noqa: E402

OUT = os.path.join(e.ROOT, "results", "EX009_followups")
CINF = 1.0


class Charge(e.Model):
    def __init__(self, R, dr, k, a, b, z):
        super().__init__(R, dr, k, 0, 0)
        self.a, self.b, self.z = a, b, z
        self.h0 = np.array([e.A, e.B / e.A, 0.0, CINF])

    def F(self, u):
        X, Y, W, C = u
        Rr = self.k * W * (self.a * u[self.z] - self.b * C)
        f = np.array([self.LapN @ X + e.A - (1 + e.B) * X + X * X * Y,
                      e.DY * (self.LapN @ Y) + e.B * X - X * X * Y,
                      e.DW * (self.LapN @ W) - W + e.C * Y * W * W - W ** 3,
                      self.LapD @ C + self.bD * CINF + Rr])
        f[self.z] -= Rr
        return f

    def J(self, u, l=0):
        X, Y, W, C = u
        k, d, n = self.k, sp.diags, self.n
        if l == 0:
            LX, LY, LW, LC = self.LapN, e.DY * self.LapN, e.DW * self.LapN, self.LapD
        else:
            cen = d(-2.0 / self.r ** 2)
            LX, LY, LW, LC = self.LapD + cen, e.DY * (self.LapD + cen), e.DW * (self.LapD + cen), self.LapD + cen
        blocks = [[LX + d(-(1 + e.B) + 2 * X * Y), d(X * X), sp.csr_matrix((n, n)), sp.csr_matrix((n, n))],
                  [d(e.B - 2 * X * Y), LY - d(X * X), sp.csr_matrix((n, n)), sp.csr_matrix((n, n))],
                  [sp.csr_matrix((n, n)), d(e.C * W * W), LW + d(-1 + 2 * e.C * Y * W - 3 * W * W), sp.csr_matrix((n, n))],
                  [sp.csr_matrix((n, n)), sp.csr_matrix((n, n)), sp.csr_matrix((n, n)), LC]]
        dR = [None, None, None, None]                     # dR/d(species)
        dR[self.z] = d(k * W * self.a)
        dR[2] = d(k * (self.a * u[self.z] - self.b * C))
        dR[3] = d(-k * W * self.b)
        for j in range(4):
            if dR[j] is not None:
                blocks[3][j] = blocks[3][j] + dR[j]
                blocks[self.z][j] = blocks[self.z][j] - dR[j]
        return sp.bmat(blocks).tocsc()

    def drift(self, U):
        J1 = self.J(U, l=1)
        Up = np.array([np.gradient(U[s], self.r) for s in range(4)]).ravel()
        luT = spla.splu(J1.T.tocsc())
        lv = np.random.default_rng(0).standard_normal(J1.shape[0])
        for _ in range(8):
            lv = luT.solve(lv)
            lv /= np.linalg.norm(lv)
        W = U[2]
        P = np.zeros((4, self.n))
        P[3] = -self.k * self.b * W * self.r
        P[self.z] += self.k * self.b * W * self.r
        return -float(lv @ P.ravel()) / float(lv @ Up)


def run(k, eps, z, sign, R=60.0):
    a, b = (1.0, eps) if sign > 0 else (eps, 1.0)
    M = Charge(R, 0.05, k, a, b, z)
    d = np.load(os.path.join(e.OUT, "particle_GX_k0.05_R60.npz"))
    u = np.array([np.interp(M.r, d["r"], d["u"][s]) for s in range(4)])
    u[3] = CINF
    u, res, ok = M.newton(u)
    if not (ok and u[2, 0] > 1):
        return dict(k=k, eps=eps, z="XY"[z], sign=sign, ok=False, residual=res)
    Rr = k * u[2] * (a * u[z] - b * u[3])
    Q = M.integral(Rr)
    r = M.r
    m = (r > 15) & (r < R - 5)
    pred = Q / (4 * np.pi) * (1 / r[m] - 1 / R)
    dev = float(np.max(np.abs((u[3][m] - CINF) - pred)) / np.max(np.abs(pred)))
    mu = M.drift(u)
    eig0 = e.leading_eigs(M, u, 0)[0]
    return dict(k=k, eps=eps, z="XY"[z], sign=sign, ok=True, residual=res, Q=Q, far_dev=dev, mu_C=mu,
                eig_l0=eig0, core_C=float(u[3, 0]), wanted=("mu<0" if sign > 0 else "mu>0"),
                as_wanted=bool(mu < 0 if sign > 0 else mu > 0))


def main():
    out = []
    for z in (0, 1):
        for sign in (+1, -1):
            for k in (0.05, 0.2):
                r_ = run(k, 0.1, z, sign)
                print(r_, flush=True)
                out.append(r_)
    json.dump(out, open(os.path.join(OUT, "partC_charge.json"), "w"), indent=1, default=float)


if __name__ == "__main__":
    main()

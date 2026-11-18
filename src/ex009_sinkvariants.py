"""EX009 part C2: which C-consuming core reaction makes a sink drift UP the C gradient?

Sink reaction C + W -> W + P with rate k C W, product P in {X, Y, W (C feeds the core species), none}.
"""
import json
import os
import sys

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ex008_longrange as e  # noqa: E402
import ex009_charge as ch  # noqa: E402


class Sink(e.Model):
    def __init__(self, R, dr, k, p):
        super().__init__(R, dr, k, 0, 0)
        self.p = p          # 0 X, 1 Y, 2 W, None
        self.h0 = np.array([e.A, e.B / e.A, 0.0, ch.CINF])

    def F(self, u):
        X, Y, W, C = u
        Rr = self.k * C * W
        f = np.array([self.LapN @ X + e.A - (1 + e.B) * X + X * X * Y,
                      e.DY * (self.LapN @ Y) + e.B * X - X * X * Y,
                      e.DW * (self.LapN @ W) - W + e.C * Y * W * W - W ** 3,
                      self.LapD @ C + self.bD * ch.CINF - Rr])
        if self.p is not None:
            f[self.p] += Rr
        return f

    def J(self, u, l=0):
        X, Y, W, C = u
        k, d, n = self.k, sp.diags, self.n
        Z = sp.csr_matrix((n, n))
        if l == 0:
            LX, LY, LW, LC = self.LapN, e.DY * self.LapN, e.DW * self.LapN, self.LapD
        else:
            cen = d(-2.0 / self.r ** 2)
            LX, LY, LW, LC = self.LapD + cen, e.DY * (self.LapD + cen), e.DW * (self.LapD + cen), self.LapD + cen
        B = [[LX + d(-(1 + e.B) + 2 * X * Y), d(X * X), Z, Z],
             [d(e.B - 2 * X * Y), LY - d(X * X), Z, Z],
             [Z, d(e.C * W * W), LW + d(-1 + 2 * e.C * Y * W - 3 * W * W), Z],
             [Z, Z, d(-k * C), LC + d(-k * W)]]
        if self.p is not None:
            B[self.p][2] = B[self.p][2] + d(k * C)
            B[self.p][3] = B[self.p][3] + d(k * W)
        return sp.bmat(B).tocsc()

    def drift(self, U):
        J1 = self.J(U, l=1)
        Up = np.array([np.gradient(U[s], self.r) for s in range(4)]).ravel()
        luT = spla.splu(J1.T.tocsc())
        lv = np.random.default_rng(0).standard_normal(J1.shape[0])
        for _ in range(8):
            lv = luT.solve(lv)
            lv /= np.linalg.norm(lv)
        P = np.zeros((4, self.n))
        P[3] = -self.k * U[2] * self.r
        if self.p is not None:
            P[self.p] += self.k * U[2] * self.r
        return -float(lv @ P.ravel()) / float(lv @ Up)




def main():
    out = []
    for p, name in ((0, "X"), (1, "Y"), (2, "W"), (None, "none")):
        for k in (0.05, 0.2):
            M = Sink(60.0, 0.05, k, p)
            d = np.load(os.path.join(e.OUT, "particle_GX_k0.05_R60.npz"))
            u = np.array([np.interp(M.r, d["r"], d["u"][s]) for s in range(4)])
            u[3] = ch.CINF
            u, res, ok = M.newton(u)
            r_ = dict(product=name, k=k, ok=bool(ok and u[2, 0] > 1))
            if r_["ok"]:
                r_.update(Q=-M.integral(k * u[3] * u[2]), mu_C=M.drift(u), eig_l0=e.leading_eigs(M, u, 0)[0])
            print(r_, flush=True)
            out.append(r_)
    json.dump(out, open(os.path.join(ch.OUT, "partC_sink_variants.json"), "w"), indent=1, default=float)


if __name__ == "__main__":
    main()

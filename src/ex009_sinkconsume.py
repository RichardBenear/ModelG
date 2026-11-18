"""EX009 part C3: sink reaction C + Z + W -> W (removes C and Z = X or Y together), rate k C W Z."""
import json
import os
import sys

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ex008_longrange as e  # noqa: E402
import ex009_charge as ch  # noqa: E402
import ex009_sinkvariants as sv  # noqa: E402


class Sink2(sv.Sink):
    def __init__(self, R, dr, k, z):
        super().__init__(R, dr, k, None)
        self.z = z

    def F(self, u):
        X, Y, W, C = u
        Rr = self.k * C * W * u[self.z]
        f = np.array([self.LapN @ X + e.A - (1 + e.B) * X + X * X * Y,
                      e.DY * (self.LapN @ Y) + e.B * X - X * X * Y,
                      e.DW * (self.LapN @ W) - W + e.C * Y * W * W - W ** 3,
                      self.LapD @ C + self.bD * ch.CINF - Rr])
        f[self.z] -= Rr
        return f



def build_J(model, u, l):
    X, Y, W, C = u
    k, d, z, n = model.k, sp.diags, model.z, model.n
    Zm = sp.csr_matrix((n, n))
    if l == 0:
        LX, LY, LW, LC = model.LapN, e.DY * model.LapN, e.DW * model.LapN, model.LapD
    else:
        cen = d(-2.0 / model.r ** 2)
        LX, LY, LW, LC = model.LapD + cen, e.DY * (model.LapD + cen), e.DW * (model.LapD + cen), model.LapD + cen
    B = [[LX + d(-(1 + e.B) + 2 * X * Y), d(X * X), Zm, Zm],
         [d(e.B - 2 * X * Y), LY - d(X * X), Zm, Zm],
         [Zm, d(e.C * W * W), LW + d(-1 + 2 * e.C * Y * W - 3 * W * W), Zm],
         [Zm, Zm, Zm, LC]]
    Zv = u[z]
    dR = {z: k * C * W, 2: k * C * Zv, 3: k * W * Zv}
    for col, val in dR.items():
        B[3][col] = B[3][col] - d(val)
        B[z][col] = B[z][col] - d(val)
    return sp.bmat(B).tocsc()


Sink2.J = lambda self, u, l=0: build_J(self, u, l)


def drift(model, U):
    J1 = model.J(U, l=1)
    Up = np.array([np.gradient(U[s], model.r) for s in range(4)]).ravel()
    luT = spla.splu(J1.T.tocsc())
    lv = np.random.default_rng(0).standard_normal(J1.shape[0])
    for _ in range(8):
        lv = luT.solve(lv)
        lv /= np.linalg.norm(lv)
    P = np.zeros((4, model.n))
    v = model.k * U[2] * U[model.z] * model.r
    P[3] -= v
    P[model.z] -= v
    return -float(lv @ P.ravel()) / float(lv @ Up)


out = []
for z, name in ((0, "X"), (1, "Y")):
    for k in (0.005, 0.02):
        M = Sink2(60.0, 0.05, k, z)
        d = np.load(os.path.join(e.OUT, "particle_GX_k0.05_R60.npz"))
        u = np.array([np.interp(M.r, d["r"], d["u"][s]) for s in range(4)])
        u[3] = ch.CINF
        u, res, ok = M.newton(u)
        r_ = dict(consumes=name, k=k, ok=bool(ok and u[2, 0] > 1), residual=res)
        if r_["ok"]:
            r_.update(Q=-M.integral(k * u[3] * u[2] * u[z]), mu_C=drift(M, u), eig_l0=e.leading_eigs(M, u, 0)[0])
        print(r_, flush=True)
        out.append(r_)
json.dump(out, open(os.path.join(ch.OUT, "partC_sink_consume.json"), "w"), indent=1, default=float)

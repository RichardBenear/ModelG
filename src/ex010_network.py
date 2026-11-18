"""EX010 (exploratory): generic mass-action reaction network on the 3D radial grid, with l = 1 drift response.

Species: X, Y, W (the particle's core species), C (conserved charge carrier), G (conserved gravity carrier).
C and G have no ether reactions and are fixed at the outer boundary (C_inf = G_inf = 1); the others have zero flux.
Reactions: (k, {reactant: order}, {product: count}); rate k * prod u^order.
Drift: v = -l^T P / l^T U' for a far gradient E z in field f, P = (d reaction terms / d u_f) * r.
"""
import json
import os
import sys
import time

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
from scipy.optimize import fsolve

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ex008_longrange as e  # noqa: E402

OUT = os.path.join(e.ROOT, "results", "EX010_charge_gravity")
SP = ["X", "Y", "W", "C", "G"]
IX = {s: i for i, s in enumerate(SP)}
DIFF = np.array([1.0, 12.0, 0.1, 1.0, 1.0])
DIRICHLET = {3: 1.0, 4: 1.0}
A, B, CW = 140 / 9, 29.0, 1.0


def base():
    return [(A, {}, {"X": 1}), (B, {"X": 1}, {"Y": 1}), (1.0, {"X": 2, "Y": 1}, {"X": 3}), (1.0, {"X": 1}, {}),
            (1.0, {"W": 1}, {}), (CW, {"Y": 1, "W": 2}, {"Y": 1, "W": 3}), (1.0, {"W": 3}, {"W": 2})]


class Net(e.Model):
    def __init__(self, reactions, R=60.0, dr=0.05):
        super().__init__(R, dr, 0.0, 0, 0)
        self.rx = [(k, {IX[s]: o for s, o in re.items()}, {IX[s]: c for s, c in pr.items()}) for k, re, pr in reactions]
        self.ns = len(SP)

    def rates(self, u):
        out = []
        for k, re, pr in self.rx:
            r = np.full(self.n, k)
            for i, o in re.items():
                r = r * u[i] ** o
            out.append(r)
        return out

    def reac(self, u):
        f = np.zeros((self.ns, self.n))
        for (k, re, pr), r in zip(self.rx, self.rates(u)):
            for i, o in re.items():
                f[i] -= o * r
            for i, c in pr.items():
                f[i] += c * r
        return f

    def dreac(self, u):
        """D[i][j] = d f_i / d u_j (arrays)."""
        D = [[np.zeros(self.n) for _ in range(self.ns)] for _ in range(self.ns)]
        for k, re, pr in self.rx:
            for j, oj in re.items():
                dr = np.full(self.n, k * oj)
                for i, o in re.items():
                    dr = dr * (u[i] ** (o - 1 if i == j else o))
                for i, o in re.items():
                    D[i][j] -= o * dr
                for i, c in pr.items():
                    D[i][j] += c * dr
        return D

    def lap(self, s, f):
        if s in DIRICHLET:
            return self.LapD @ f + self.bD * DIRICHLET[s]
        return self.LapN @ f

    def F(self, u):
        return np.array([DIFF[s] * self.lap(s, u[s]) for s in range(self.ns)]) + self.reac(u)

    def J(self, u, l=0):
        D = self.dreac(u)
        cen = sp.diags(-2.0 / self.r ** 2)
        blocks = []
        for i in range(self.ns):
            row = []
            for j in range(self.ns):
                Mij = sp.diags(D[i][j])
                if i == j:
                    L = (self.LapD if (s_d := i in DIRICHLET) or l == 1 else self.LapN)
                    if l == 1:
                        L = self.LapD + cen
                    Mij = Mij + DIFF[i] * L
                row.append(Mij)
            blocks.append(row)
        return sp.bmat(blocks).tocsc()

    def newton(self, u, tol=1e-10, maxit=40):
        for it in range(maxit):
            Fv = self.F(u).ravel()
            res = float(np.max(np.abs(Fv)))
            if res < tol:
                return u, res, True
            u = u + spla.spsolve(self.J(u), -Fv).reshape(self.ns, self.n)
            if not np.all(np.isfinite(u)) or np.any(u[[0, 1]] <= 0):
                return u, np.inf, False
        res = float(np.max(np.abs(self.F(u))))
        return u, res, res < tol

    def drifts(self, U, fields=(3, 4)):
        J1 = self.J(U, l=1)
        Up = np.array([np.gradient(U[s], self.r) for s in range(self.ns)]).ravel()
        luT = spla.splu(J1.T.tocsc())
        lv = np.random.default_rng(0).standard_normal(J1.shape[0])
        for _ in range(8):
            lv = luT.solve(lv)
            lv /= np.linalg.norm(lv)
        D = self.dreac(U)
        out = {}
        for f in fields:
            P = np.array([D[i][f] * self.r for i in range(self.ns)]).ravel()
            out[SP[f]] = -float(lv @ P) / float(lv @ Up)
        return out

    def net_production(self, u, s):
        """Net reaction production of species s (= minus the boundary inflow at steady state)."""
        return self.integral(self.reac(u)[s])


def start(net):
    d = np.load(os.path.join(e.OUT, "particle_GX_k0.05_R60.npz"))
    u = np.ones((net.ns, net.n))
    for s in range(3):
        u[s] = np.interp(net.r, d["r"], d["u"][s])
    return u


def solve(reactions, u0=None):
    net = Net(reactions)
    u = start(net) if u0 is None else u0
    u, res, ok = net.newton(u)
    ok = bool(ok and u[2, 0] > 1)
    return net, u, res, ok


def plus(kf, eps):
    return [(kf, {"Y": 1, "W": 1}, {"W": 1, "C": 1}), (kf * eps, {"C": 1, "W": 1}, {"W": 1, "Y": 1})]


def minus(km):
    return [(km, {"C": 1, "X": 1, "W": 1}, {"W": 1})]


def grav(kg):
    return [(kg, {"G": 1, "W": 1}, {"W": 1, "X": 1})]

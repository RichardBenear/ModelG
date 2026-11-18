"""EX007 (exploratory): hybrid Model G-2 + G -> X (our construction).

G_t = G'' + beta - sigma G X - gamma G
X_t = X'' + alpha + gamma G - sigma G X + X^2 Y - X
Y_t = 12 Y'' + sigma G X - X^2 Y
Held fixed: X0 = 140/9, b_eff = sigma G0.  f = gamma G0 / X0 is the share of the X supply that comes from G:
alpha = X0 (1 - f), gamma = f X0 / G0, beta = G0 (sigma X0 + gamma).  f = 0 is EX006's G-2.
"""
import json
import os
import sys
import time
from multiprocessing import Pool

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ex006_g2_1d as w  # noqa: E402
import mg192_mobility as m  # noqa: E402

ROOT = w.ROOT
OUT = os.path.join(ROOT, "results", "EX007_g2_gx")
X0 = 140 / 9
SIGMA = 0.001
D = np.diag([1.0, 1.0, 12.0])


class G2X(w.G2):
    def __init__(self, sigma, beff, f):
        self.s, self.beff, self.f = sigma, beff, f
        self.G0 = beff / sigma
        self.gamma = f * X0 / self.G0
        self.alpha = X0 * (1 - f)
        self.beta = self.G0 * (sigma * X0 + self.gamma)
        self.h0 = np.array([self.G0, X0, beff / X0])

    def reac(self, u, beta=None):
        G, X, Y = u
        beta = self.beta if beta is None else beta
        sGX = self.s * G * X
        return np.array([beta - sGX - self.gamma * G, self.alpha + self.gamma * G - sGX + X * X * Y - X, sGX - X * X * Y])

    def jac(self, grid, u, D2):
        J = super().jac(grid, u, D2)
        N = grid.N
        J[:N, :N] -= self.gamma * np.eye(N)
        J[N:2 * N, :N] += self.gamma * np.eye(N)
        return J

    def dbeta_dG0(self):
        # beta(G0) with alpha, gamma, sigma fixed: X0 = alpha + gamma G0, beta = G0 (sigma X0 + gamma)
        return self.s * X0 + self.gamma + self.s * self.G0 * self.gamma

    def mobility(self, grid, psi):
        mu_b = -float(np.dot(psi[0], grid.x))
        return dict(mu_beta=mu_b, mu_G=mu_b * self.dbeta_dG0())


def jac_h(beff, f, sigma=SIGMA):
    G0 = beff / sigma
    g = f * X0 / G0
    return np.array([[-sigma * X0 - g, -sigma * G0, 0.0],
                     [g - sigma * X0, sigma * G0 - 1, X0 ** 2],
                     [sigma * X0, -sigma * G0, -X0 ** 2]])


def b_c(f, sigma=SIGMA):
    ks = np.linspace(0.0, 4, 801)
    lo, hi = 10.0, 60.0
    for _ in range(50):
        mid = 0.5 * (lo + hi)
        gr = max(np.max(np.linalg.eigvals(jac_h(mid, f, sigma) - D * k * k).real) for k in ks)
        lo, hi = (lo, mid) if gr > 0 else (mid, hi)
    return lo


def rates(beff, f, sigma=SIGMA):
    J = jac_h(beff, f, sigma)
    ss = np.arange(4.0)
    c = np.polyfit(ss, [np.linalg.det(J + D * s) for s in ss], 3)
    out = []
    for s in np.roots(c):
        k = np.sqrt(complex(s))
        out.append(dict(decay=float(abs(k.real)), wavenumber=float(abs(k.imag))))
    return sorted(out, key=lambda d: d["decay"])


FS = (0.0, 0.2, 0.4, 0.6, 0.8, 0.9, 0.95, 1.0)


def window(f, g, D2, u_seed, b_seed):
    """Existence window in b_eff at fixed f, by continuation up and down from (b_seed, u_seed)."""
    ends = []
    for sgn in (+1, -1):
        b, u, db = b_seed, u_seed, 0.02
        while db > 1e-4:
            M0, M = G2X(SIGMA, b, f), G2X(SIGMA, b + sgn * db, f)
            ut, res, ok = M.newton_even(g, u - M0.h0[:, None] + M.h0[:, None], D2=D2, maxit=25)
            if ok and ut[1].min() - M.h0[1] < -1:
                b, u = b + sgn * db, ut
            else:
                db /= 2
        ends.append(b)
    return ends[1], ends[0]


def embed(M, g1, u, g):
    phi = u - M.h0[:, None]
    big = np.repeat(M.h0[:, None], g.N, axis=1).astype(float)
    for s in range(3):
        big[s] += np.where(np.abs(g.x) < 49.9, np.interp(np.abs(g.x), g1.x, phi[s]), 0.0)
    U, res, ok = M.newton_even(g, big)
    assert ok, res
    return U, res


def scan():
    g1 = m.Grid(100.0, 1024)
    D2 = g1.D2()
    g = m.Grid(200.0, 2048)
    b0 = 29.221565
    M = G2X(SIGMA, b0, 0.0)
    u, res, ok = M.newton_even(g1, w.seed_from_model_g(g1, M), D2=D2)
    out = []
    fprev = 0.0
    for f in FS:
        # continue in f at b0
        for ff in np.linspace(fprev, f, 4)[1:] if f > fprev else []:
            Mp, Mn = G2X(SIGMA, b0, fprev), G2X(SIGMA, b0, ff)
            u, res, ok = Mn.newton_even(g1, u - Mp.h0[:, None] + Mn.h0[:, None], D2=D2)
            assert ok and u[1].min() - Mn.h0[1] < -1, (f, ff)
            fprev = ff
        fprev = f
        lo, hi = window(f, g1, D2, u, b0)
        bm = 0.5 * (lo + hi)
        M0, M = G2X(SIGMA, b0, f), G2X(SIGMA, bm, f)
        um, res, ok = M.newton_even(g1, u - M0.h0[:, None] + M.h0[:, None], D2=D2)
        assert ok
        U, res = embed(M, g1, um, g)
        psi, Up, lam, nr = M.adjoint(g, U)
        mob = M.mobility(g, psi)
        rt = rates(bm, f)
        km = rt[0]["decay"]
        x, phiG = g.x, U[0] - M.h0[0]
        msk = (x > 30) & (x < 70)
        basis = np.exp(-km * x[msk]) + np.exp(-km * (200 - x[msk]))
        A = float(np.sum(basis * phiG[msk]) / np.sum(basis * basis))
        r = dict(f=f, gamma=M.gamma, alpha=M.alpha, window=[lo, hi], beff=bm, b_c=b_c(f), residual=res,
                 core_X=float(U[1].min() - X0), core_G=float(phiG[g.N // 2]), int_G=float(phiG.sum() * g.h),
                 int_X=float((U[1] - X0).sum() * g.h), leading_eigs=[[float(z.real), float(z.imag)] for z in lam[:4]],
                 rates=rt, tail_A=A, **mob,
                 pred_far_force_amp=float(2 * mob["mu_G"] * A * (-km)))
        np.savez_compressed(os.path.join(OUT, f"particle_f{f:g}.npz"), U=U, x=g.x, h0=M.h0, beff=bm, f=f)
        print({k: v for k, v in r.items() if k != "rates"}, flush=True)
        out.append(r)
    return out


def pair_run(args):
    f, d, dt = args
    p = np.load(os.path.join(OUT, f"particle_f{f:g}.npz"))
    M = G2X(SIGMA, float(p["beff"]), f)
    g = m.Grid(200.0, 2048)
    phi = p["U"] - M.h0[:, None]
    u0 = M.h0[:, None] + np.array([g.shift(phi[s], d / 2) + g.shift(phi[s], -d / 2) for s in range(3)])
    st = w.ETD(M, g, dt)
    vh = np.fft.fft(u0, axis=1)
    ts, sep = [], []
    every = int(round(1.0 / dt))
    for n in range(int(round(100.0 / dt)) + 1):
        if n % every == 0:
            u = np.real(np.fft.ifft(vh, axis=1))
            xl, xr, core = w.two_min(g, u[1], X0)
            ts.append(n * dt)
            sep.append(xr - xl)
        vh = st.step(vh)
    ts, sep = np.array(ts), np.array(sep)
    sel = ts >= 20
    return dict(f=f, d=d, dt=dt, v_rel=float(np.polyfit(ts[sel], sep[sel], 1)[0]), core_end=core)


def main():
    os.makedirs(OUT, exist_ok=True)
    t0 = time.time()
    rec = dict(note="EXPLORATORY", sigma=SIGMA, X0=X0, scan=scan(), pairs=[])
    json.dump(rec, open(os.path.join(OUT, "results.json"), "w"), indent=1, default=float)
    jobs = [(r["f"], d, 0.005) for r in rec["scan"] for d in (40.0, 50.0, 60.0)]
    with Pool(int(os.environ.get("MODELG_WORKERS", "4"))) as pool:
        for r in pool.imap_unordered(pair_run, jobs):
            print("pair", r, flush=True)
            rec["pairs"].append(r)
    rec["wall_seconds"] = time.time() - t0
    json.dump(rec, open(os.path.join(OUT, "results.json"), "w"), indent=1, default=float)
    print("done", rec["wall_seconds"])


if __name__ == "__main__":
    main()

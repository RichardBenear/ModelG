"""EX006 (exploratory): Model G-2 in 1D, Fourier spectral (MG192/MG194 methods adapted). kappa = 0.

G_t = G'' + beta - sigma G X,  X_t = X'' + alpha - sigma G X + X^2 Y - X,  Y_t = 12 Y'' + sigma G X - X^2 Y.
Background: X0 = alpha, G0 = beta / (sigma alpha), Y0 = sigma G0 / alpha;  beff = sigma G0.
"""
import os
import sys

import numpy as np
import scipy.linalg as sla

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import mg192_mobility as m  # noqa: E402  (Grid, even_maps)

ROOT = m.ROOT
OUT = os.path.join(ROOT, "results", "EX006_g2_drift_pair")
ALPHA = 140 / 9
DIFF = np.array([1.0, 1.0, 12.0])


class G2:
    def __init__(self, sigma, beff):
        self.s, self.beff = sigma, beff
        self.G0 = beff / sigma
        self.beta = self.G0 * sigma * ALPHA
        self.h0 = np.array([self.G0, ALPHA, beff / ALPHA])

    def reac(self, u, beta=None):
        G, X, Y = u
        beta = self.beta if beta is None else beta
        sGX = self.s * G * X
        return np.array([beta - sGX, ALPHA - sGX + X * X * Y - X, sGX - X * X * Y])

    def rhs(self, grid, u, beta=None):
        return np.array([DIFF[i] * grid.lap(u[i]) for i in range(3)]) + self.reac(u, beta)

    def jac(self, grid, u, D2):
        G, X, Y = u
        s, N = self.s, grid.N
        I = np.eye(N)
        J = np.zeros((3 * N, 3 * N))
        J[:N, :N] = D2 + np.diag(-s * X)
        J[:N, N:2 * N] = np.diag(-s * G)
        J[N:2 * N, :N] = np.diag(-s * X)
        J[N:2 * N, N:2 * N] = D2 + np.diag(-s * G + 2 * X * Y - 1)
        J[N:2 * N, 2 * N:] = np.diag(X * X)
        J[2 * N:, :N] = np.diag(s * X)
        J[2 * N:, N:2 * N] = np.diag(s * G - 2 * X * Y)
        J[2 * N:, 2 * N:] = 12.0 * D2 - np.diag(X * X)
        del I
        return J

    def newton_even(self, grid, u0, tol=1e-9, maxit=40, D2=None):
        N = grid.N
        D2 = grid.D2() if D2 is None else D2
        E, rows = m.even_maps(N)
        u = u0.copy()
        scale = np.array([self.G0, ALPHA, 1.0])[:, None]
        for it in range(maxit):
            F = self.rhs(grid, u)
            res = float(np.max(np.abs(F / np.maximum(scale, 1.0))))
            if res < tol:
                return u, res, True
            J = self.jac(grid, u, D2)
            dr = np.linalg.solve(J[rows] @ E, -F.ravel()[rows])
            u = u + (E @ dr).reshape(3, N)
            if np.any(u <= 0) or not np.all(np.isfinite(u)):
                return u, np.inf, False
        F = self.rhs(grid, u)
        res = float(np.max(np.abs(F / np.maximum(scale, 1.0))))
        return u, res, res < tol

    def adjoint(self, grid, U, nev=6):
        """Translation adjoint psi (J^T psi = 0, <psi, U'> = 1) and the leading eigenvalues of J."""
        D2 = grid.D2()
        J = self.jac(grid, U, D2)
        Up = np.array([grid.ddx(U[s]) for s in range(3)]).ravel()
        lu = sla.lu_factor(J.T - 1e-10 * np.eye(J.shape[0]))
        psi = np.random.default_rng(1).standard_normal(J.shape[0])
        for _ in range(6):
            psi = sla.lu_solve(lu, psi)
            psi /= np.linalg.norm(psi)
        psi = psi / np.dot(psi, Up)
        # leading eigenvalues by shift-invert Arnoldi around +0.5
        import scipy.sparse.linalg as spla
        lus = sla.lu_factor(J - 0.5 * np.eye(J.shape[0]))
        op = spla.LinearOperator(J.shape, matvec=lambda v: sla.lu_solve(lus, v))
        mu = spla.eigs(op, k=nev, which="LM", return_eigenvectors=False)
        lam = np.sort_complex(0.5 + 1 / mu)[::-1]
        return psi.reshape(3, grid.N), Up.reshape(3, grid.N), lam, float(np.max(np.abs(J @ Up.ravel())))

    def mobility(self, grid, psi):
        """v = mu_beta * dbeta/dx (uniform source gradient); per background-G gradient: mu_G = mu_beta * sigma * alpha."""
        mu_b = -float(np.dot(psi[0], grid.x))
        return dict(mu_beta=mu_b, mu_G=mu_b * self.s * ALPHA)


class ETD:
    """ETDRK4 (Kassam & Trefethen), as MG192, with G-2 reactions and an optional beta(x)."""

    def __init__(self, model, grid, dt, M=64):
        self.m, self.grid, self.dt = model, grid, dt
        Lop = -DIFF[:, None] * grid.k2[None, :]
        self.E = np.exp(dt * Lop)
        self.E2 = np.exp(dt * Lop / 2)
        r = np.exp(1j * np.pi * (np.arange(1, M + 1) - 0.5) / M)
        LR = dt * Lop[:, :, None] + r[None, None, :]
        self.Q = dt * np.real(np.mean((np.exp(LR / 2) - 1) / LR, axis=2))
        self.f1 = dt * np.real(np.mean((-4 - LR + np.exp(LR) * (4 - 3 * LR + LR ** 2)) / LR ** 3, axis=2))
        self.f2 = dt * np.real(np.mean((2 + LR + np.exp(LR) * (-2 + LR)) / LR ** 3, axis=2))
        self.f3 = dt * np.real(np.mean((-4 - 3 * LR - LR ** 2 + np.exp(LR) * (4 - LR)) / LR ** 3, axis=2))

    def nl(self, vh, beta):
        return np.fft.fft(self.m.reac(np.real(np.fft.ifft(vh, axis=1)), beta), axis=1)

    def step(self, vh, beta=None):
        Nv = self.nl(vh, beta)
        av = self.E2 * vh + self.Q * Nv
        Na = self.nl(av, beta)
        bv = self.E2 * vh + self.Q * Na
        Nb = self.nl(bv, beta)
        cv = self.E2 * av + self.Q * (2 * Nb - Nv)
        Nc = self.nl(cv, beta)
        return self.E * vh + Nv * self.f1 + 2 * (Na + Nb) * self.f2 + Nc * self.f3


def seed_from_model_g(grid, model):
    """Model G 1D soliton (MG194, b* = 29.205) X and Y deviations on the G-2 background; G deviation 0."""
    U = np.load(os.path.join(ROOT, "results", "MG194_two_particle", "soliton_bstar.npz"))["U"]
    g0 = m.Grid(100.0, 1024)
    h = np.array([140 / 9, 140 / 9, 29.205258852005002 / (140 / 9)])
    phi = U - h[:, None]
    u = np.repeat(model.h0[:, None], grid.N, axis=1).astype(float)
    for s in (1, 2):
        u[s] += np.interp(grid.x, g0.x, phi[s], left=0.0, right=0.0)
    return u


# ---------------- driver ----------------
import json  # noqa: E402
import time  # noqa: E402
from multiprocessing import Pool  # noqa: E402

CASES = {0.0003: dict(window=(29.03625, 29.47453), L=300.0, N=3072),
         0.001: dict(window=(29.04813, 29.39500), L=200.0, N=2048),
         0.003: dict(window=(28.69734, 29.17484), L=200.0, N=2048)}
DT, T_PAIR, T_FIT = 0.005, 100.0, 20.0
D_LIST = np.concatenate([np.arange(8.0, 20.0, 0.5), np.arange(20.0, 70.01, 1.0)])


def beff_of(sigma):
    lo, hi = CASES[sigma]["window"]
    return 0.5 * (lo + hi)


def particle(sigma, L, N):
    model = w_model = G2(sigma, beff_of(sigma))
    g1 = m.Grid(100.0, 1024)
    # continue from the Model G seed at 29.2 to the window midpoint
    u, res, ok = G2(sigma, 29.2).newton_even(g1, seed_from_model_g(g1, G2(sigma, 29.2)))
    bs = np.linspace(29.2, model.beff, 6)[1:]
    prev = G2(sigma, 29.2)
    for b in bs:
        cur = G2(sigma, b)
        u, res, ok = cur.newton_even(g1, u - prev.h0[:, None] + cur.h0[:, None])
        prev = cur
    assert ok and u[1].min() - model.h0[1] < -1, (sigma, res)
    g = m.Grid(L, N)
    phi = u - model.h0[:, None]
    big = np.repeat(model.h0[:, None], N, axis=1).astype(float)
    for s in range(3):
        big[s] += np.where(np.abs(g.x) < 49.9, np.interp(np.abs(g.x), g1.x, phi[s]), 0.0)
    U, res, ok = w_model.newton_even(g, big)
    assert ok, res
    return model, g, U, res


def xmin(g, X, guess):
    return m.xmin_position(g, X, guess)


def two_min(g, X, X0):
    left, right = g.x < 0, g.x > 0
    jl = np.where(left)[0][np.argmin(X[left])]
    jr = np.where(right)[0][np.argmin(X[right])]
    return xmin(g, X, g.x[jl]), xmin(g, X, g.x[jr]), float(min(X[jl], X[jr]) - X0)


def pair_run(args):
    sigma, d, L, N, dt = args
    model = G2(sigma, beff_of(sigma))
    g = m.Grid(L, N)
    U = np.load(os.path.join(OUT, f"particle_s{sigma:g}.npz"))["U"]
    if L != CASES[sigma]["L"]:
        g0 = m.Grid(CASES[sigma]["L"], CASES[sigma]["N"])
        phi0 = U - model.h0[:, None]
        lim = g0.L / 2 - g0.h
        U = model.h0[:, None] + np.array([np.where(np.abs(g.x) < lim, np.interp(np.abs(g.x), g0.x, phi0[s]), 0.0) for s in range(3)])
        U, _, _ = model.newton_even(g, U)
    phi = U - model.h0[:, None]
    if d > 0:
        u0 = model.h0[:, None] + np.array([g.shift(phi[s], d / 2) + g.shift(phi[s], -d / 2) for s in range(3)])
    else:
        u0 = U.copy()
    st = ETD(model, g, dt)
    vh = np.fft.fft(u0, axis=1)
    X0 = model.h0[1]
    every = int(round(1.0 / dt))
    ts, sep, dep = [], [], []
    xl, xr = (-d / 2, d / 2) if d > 0 else (0.0, 0.0)
    t0 = time.time()
    for n in range(int(round(T_PAIR / dt)) + 1):
        if n % every == 0:
            u = np.real(np.fft.ifft(vh, axis=1))
            if not np.all(np.isfinite(u)) or np.any(u <= 0):
                return dict(sigma=sigma, d=d, L=L, error="lost positivity", t=n * dt)
            if d > 0:
                xl, xr, core = two_min(g, u[1], X0)
            else:
                xr = xmin(g, u[1], xr)
                core = float(u[1].min() - X0)
            ts.append(n * dt)
            sep.append(xr - xl)
            dep.append(core)
        vh = st.step(vh)
    ts, sep = np.array(ts), np.array(sep)
    sel = ts >= T_FIT
    v = float(np.polyfit(ts[sel], sep[sel], 1)[0])
    return dict(sigma=sigma, d=d, L=L, N=N, dt=dt, v_rel=v, delta_d=float(sep[-1] - sep[0]),
                merged=bool(max(dep) > -1.0 or (d > 0 and sep.min() < 2.0)), depth_end=float(dep[-1]),
                wall_s=time.time() - t0)


def drift_check(sigma, eps=0.002, L=100.0, N=1024, T=200.0, x0=12.5):
    """Direct run with beta(x) = beta (1 + eps sin(2 pi x / L)); background by full Newton."""
    model = G2(sigma, beff_of(sigma))
    g = m.Grid(L, N)
    model, g, U, _ = particle(sigma, L, N)
    psi, Up, lam, null_res = model.adjoint(g, U)
    mob = model.mobility(g, psi)
    beta = model.beta * (1 + eps * np.sin(2 * np.pi * g.x / L))
    D2 = g.D2()
    bg = np.repeat(model.h0[:, None], N, axis=1).astype(float)
    bg[0] = beta / (sigma * ALPHA)
    for _ in range(30):
        F = model.rhs(g, bg, beta).ravel()
        if np.max(np.abs(F)) < 1e-9 * model.G0:
            break
        bg = bg + np.linalg.solve(model.jac(g, bg, D2), -F).reshape(3, N)
    u0 = bg + np.array([g.shift(U[s] - model.h0[s], x0) for s in range(3)])
    st = ETD(model, g, DT)
    vh = np.fft.fft(u0, axis=1)
    ts, xs = [], []
    x = x0
    for n in range(int(round(T / DT)) + 1):
        if n % 200 == 0:
            u = np.real(np.fft.ifft(vh, axis=1))
            x = xmin(g, u[1], x)
            ts.append(n * DT)
            xs.append(x)
        vh = st.step(vh, beta)
    ts, xs = np.array(ts), np.array(xs)
    sel = ts >= 50
    v = float(np.polyfit(ts[sel], xs[sel], 1)[0])
    dbdx = model.beta * eps * 2 * np.pi / L * np.cos(2 * np.pi * x0 / L)
    pred = mob["mu_beta"] * dbdx
    return dict(sigma=sigma, eps=eps, x0=x0, L=L, v=v, v_pred=pred, rel_err=(v - pred) / pred,
                dG0dx=dbdx / (sigma * ALPHA), mobility_L100=mob, moved=float(xs[-1] - xs[0]))


def main():
    os.makedirs(OUT, exist_ok=True)
    t0 = time.time()
    rec = dict(note="EXPLORATORY", cases={str(k): v for k, v in CASES.items()}, partA=[], partB=[])
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    import ex005_model_g2 as g5
    for sigma, c in CASES.items():
        model, g, U, res = particle(sigma, c["L"], c["N"])
        np.savez_compressed(os.path.join(OUT, f"particle_s{sigma:g}.npz"), U=U, x=g.x, h0=model.h0)
        psi, Up, lam, null_res = model.adjoint(g, U)
        mob = model.mobility(g, psi)
        rt = g5.rates(model.beff, sigma, kappa=0.0)
        phi = U - model.h0[:, None]
        a = dict(sigma=sigma, beff=model.beff, G0=model.G0, residual=res, core_X=float(phi[1].min()),
                 core_G=float(phi[0][g.N // 2]), int_G=float(np.sum(phi[0]) * g.h), int_X=float(np.sum(phi[1]) * g.h),
                 leading_eigs=[[float(z.real), float(z.imag)] for z in lam], null_residual=null_res,
                 rates=rt, **mob)
        print("A", a, flush=True)
        rec["partA"].append(a)
        json.dump(rec, open(os.path.join(OUT, "results.json"), "w"), indent=1, default=float)
    for sigma in (0.001, 0.003):
        dc = drift_check(sigma)
        print("A-check", dc, flush=True)
        rec.setdefault("drift_check", []).append(dc)
        json.dump(rec, open(os.path.join(OUT, "results.json"), "w"), indent=1, default=float)
    jobs = []
    for sigma, c in CASES.items():
        jobs += [(sigma, 0.0, c["L"], c["N"], DT)] + [(sigma, float(d), c["L"], c["N"], DT) for d in D_LIST]
        jobs += [(sigma, 40.0, 2 * c["L"], 2 * c["N"], DT), (sigma, 40.0, c["L"], c["N"], DT / 2)]
    with Pool(int(os.environ.get("MODELG_WORKERS", "4"))) as pool:
        for r in pool.imap_unordered(pair_run, jobs):
            print("B", {k: (round(v, 12) if isinstance(v, float) else v) for k, v in r.items()}, flush=True)
            rec["partB"].append(r)
            json.dump(rec, open(os.path.join(OUT, "results.json"), "w"), indent=1, default=float)
    rec["wall_seconds"] = time.time() - t0
    json.dump(rec, open(os.path.join(OUT, "results.json"), "w"), indent=1, default=float)
    print("done", rec["wall_seconds"], flush=True)


if __name__ == "__main__":
    main()

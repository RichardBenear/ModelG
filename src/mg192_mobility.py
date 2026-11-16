"""MG192: SQK gravity mobility by linear response (continuum 1D, Fourier spectral).

Frozen specification: spec/MG192_gravity_mobility.md (registered a2ba0fa).
Run:  python src/mg192_mobility.py            (all stages)
Outputs go to results/MG192_gravity_mobility/.
"""
import json
import os
import platform
import sys
import time

import numpy as np
import scipy
import scipy.linalg as sla

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
OUT = os.path.join(ROOT, "results", "MG192_gravity_mobility")
REF = os.path.join(ROOT, "results", "MG003", "ref_profiles.npz")

# Published eq. 17 parameters (u = s = w = 0, p = q = 1)
P = dict(dx=1.0, dy=12.0, a=14.0, b=29.0, g=0.1, p=1.0, q=1.0)
PUBLISHED_CORE = dict(G=-0.161, X=-8.36, Y=0.930)
PUBLISHED_MU = -0.0599


def hom(a=P["a"], b=P["b"], w=0.0):
    G0 = (a + P["g"] * w) / (P["q"] - P["g"] * P["p"])
    X0 = (P["p"] * a + P["q"] * w) / (P["q"] - P["g"] * P["p"])
    Y0 = b / X0
    return np.array([G0, X0, Y0])


class Grid:
    def __init__(self, L, N):
        self.L, self.N = L, N
        self.h = L / N
        self.x = -L / 2 + self.h * np.arange(N)          # x[N//2] == 0
        self.k = 2 * np.pi * np.fft.fftfreq(N, d=self.h)
        self.k2 = self.k ** 2

    def lap(self, f):
        return np.real(np.fft.ifft(-self.k2 * np.fft.fft(f)))

    def ddx(self, f):
        kk = self.k.copy()
        kk[self.N // 2] = 0.0                              # Nyquist derivative set to zero
        return np.real(np.fft.ifft(1j * kk * np.fft.fft(f)))

    def D2(self):
        I = np.eye(self.N)
        return np.real(np.fft.ifft(-self.k2[:, None] * np.fft.fft(I, axis=0), axis=0))

    def shift(self, f, s):
        """f(x - s) by Fourier shift."""
        return np.real(np.fft.ifft(np.fft.fft(f) * np.exp(-1j * self.k * s)))


def rhs(grid, u, a=None, b=None, w=None):
    G, X, Y = u
    a = P["a"] if a is None else a
    b = P["b"] if b is None else b
    w = 0.0 if w is None else w
    fG = grid.lap(G) - P["q"] * G + P["g"] * X + a
    fX = P["dx"] * grid.lap(X) + P["p"] * G - (1 + b) * X + X * X * Y + w
    fY = P["dy"] * grid.lap(Y) + b * X - X * X * Y
    return np.array([fG, fX, fY])


def jac(grid, u, D2=None, b=None):
    G, X, Y = u
    N = grid.N
    b = P["b"] if b is None else b
    D2 = grid.D2() if D2 is None else D2
    I = np.eye(N)
    J = np.zeros((3 * N, 3 * N))
    J[:N, :N] = D2 - P["q"] * I
    J[:N, N:2 * N] = P["g"] * I
    J[N:2 * N, :N] = P["p"] * I
    J[N:2 * N, N:2 * N] = P["dx"] * D2 + np.diag(-(1 + b) + 2 * X * Y)
    J[N:2 * N, 2 * N:] = np.diag(X * X)
    J[2 * N:, N:2 * N] = np.diag(b - 2 * X * Y)
    J[2 * N:, 2 * N:] = P["dy"] * D2 - np.diag(X * X)
    return J


def even_maps(N):
    """Expansion E (full <- reduced) and row selection for even functions about x = 0 (index N//2)."""
    M = N // 2 + 1
    m = np.empty(N, dtype=int)
    for j in range(N):
        m[j] = j - N // 2 if j >= N // 2 else N // 2 - j
    E1 = np.zeros((N, M))
    E1[np.arange(N), m] = 1.0
    rows1 = np.array([N // 2 + i for i in range(N // 2)] + [0])  # reduced index i -> a full row with m == i
    E = sla.block_diag(E1, E1, E1)
    rows = np.concatenate([rows1, rows1 + N, rows1 + 2 * N])
    return E, rows


def initial_guess(grid):
    d = np.load(REF)["h0.2_b29.0"]
    xr, phi = d[0], d[1:4]
    c = xr[np.argmin(phi[1])]
    xs = (xr - c + 50.0) % 100.0 - 50.0
    o = np.argsort(xs)
    xs, phi = xs[o], phi[:, o]
    u = np.empty((3, grid.N))
    h0 = hom()
    for s in range(3):
        xp = np.concatenate([xs - 100, xs, xs + 100])
        fp = np.concatenate([phi[s]] * 3)
        xx = np.abs(grid.x) if grid.L <= 100 else np.minimum(np.abs(grid.x), 49.9)
        v = np.interp(xx, xp, fp) if grid.L <= 100 else np.where(np.abs(grid.x) < 49.9, np.interp(xx, xp, fp), 0.0)
        u[s] = h0[s] + v
    return u


def newton_even(grid, u0, tol=1e-10, maxit=40, log=None):
    N = grid.N
    D2 = grid.D2()
    E, rows = even_maps(N)
    u = u0.copy()
    hist = []
    for it in range(maxit):
        F = rhs(grid, u).ravel()
        res = np.max(np.abs(F))
        hist.append(float(res))
        if log:
            log(f"    Newton it {it}: max|F| = {res:.3e}")
        if res <= tol * 1e-2:
            break
        J = jac(grid, u, D2)
        Jr = J[rows] @ E
        dr = np.linalg.solve(Jr, -F[rows])
        u = u + (E @ dr).reshape(3, N)
        if np.any(u <= 0):
            raise RuntimeError("positivity lost in Newton")
    F = rhs(grid, u).ravel()
    return u, float(np.max(np.abs(F))), hist


def spectrum_and_adjoint(grid, U):
    N = grid.N
    J = jac(grid, U)
    lam, vl, vr = sla.eig(J, left=True, right=True)
    i0 = int(np.argmin(np.abs(lam)))
    lam0 = lam[i0]
    others = np.delete(lam, i0)
    psi = np.real(vl[:, i0])
    # refine adjoint by inverse iteration on J^T (tiny shift)
    A = J.T - (np.real(lam0) + 1e-12) * np.eye(3 * N)
    lu = sla.lu_factor(A)
    for _ in range(3):
        psi = sla.lu_solve(lu, psi)
        psi /= np.linalg.norm(psi)
    Up = np.array([grid.ddx(U[s]) for s in range(3)]).ravel()
    psi = psi / np.dot(psi, Up)
    return lam0, others, psi.reshape(3, N), Up.reshape(3, N), J


def mobilities(grid, U, psi):
    x = grid.x
    mu_a = -np.dot(psi[0], x)
    mu_w = -np.dot(psi[1], x)
    X = U[1]
    mu_b = -np.dot(psi[1] * (-X) + psi[2] * X, x)
    return dict(mu_a=float(mu_a), mu_b=float(mu_b), mu_w=float(mu_w))


# ---------------- dynamics (ETDRK4, Kassam & Trefethen 2005) ----------------

class ETDRK4:
    def __init__(self, grid, dt, M=64):
        self.grid, self.dt = grid, dt
        D = np.array([1.0, P["dx"], P["dy"]])
        Lop = -D[:, None] * grid.k2[None, :]
        self.E = np.exp(dt * Lop)
        self.E2 = np.exp(dt * Lop / 2)
        r = np.exp(1j * np.pi * (np.arange(1, M + 1) - 0.5) / M)
        LR = dt * Lop[:, :, None] + r[None, None, :]
        self.Q = dt * np.real(np.mean((np.exp(LR / 2) - 1) / LR, axis=2))
        self.f1 = dt * np.real(np.mean((-4 - LR + np.exp(LR) * (4 - 3 * LR + LR ** 2)) / LR ** 3, axis=2))
        self.f2 = dt * np.real(np.mean((2 + LR + np.exp(LR) * (-2 + LR)) / LR ** 3, axis=2))
        self.f3 = dt * np.real(np.mean((-4 - 3 * LR - LR ** 2 + np.exp(LR) * (4 - LR)) / LR ** 3, axis=2))

    def nl(self, vh, a):
        G, X, Y = np.real(np.fft.ifft(vh, axis=1))
        b = P["b"]
        n = np.array([-P["q"] * G + P["g"] * X + a,
                      P["p"] * G - (1 + b) * X + X * X * Y,
                      b * X - X * X * Y])
        return np.fft.fft(n, axis=1)

    def step(self, vh, a):
        Nv = self.nl(vh, a)
        av = self.E2 * vh + self.Q * Nv
        Na = self.nl(av, a)
        bv = self.E2 * vh + self.Q * Na
        Nb = self.nl(bv, a)
        cv = self.E2 * av + self.Q * (2 * Nb - Nv)
        Nc = self.nl(cv, a)
        return self.E * vh + Nv * self.f1 + 2 * (Na + Nb) * self.f2 + Nc * self.f3


def xmin_position(grid, X, guess):
    """Minimum of the Fourier interpolant of X near guess (Newton on X')."""
    Xh = np.fft.fft(X) / grid.N
    k = grid.k.copy()
    k[grid.N // 2] = 0.0
    j = int(np.argmin(X))
    s = grid.x[j]
    if guess is not None and abs(((s - guess + grid.L / 2) % grid.L) - grid.L / 2) > 1.0:
        s = guess
    for _ in range(30):
        e = np.exp(1j * k * (s + grid.L / 2))
        d1 = np.real(np.sum(1j * k * Xh * e))
        d2 = np.real(np.sum(-(k ** 2) * Xh * e))
        ds = -d1 / d2
        s += ds
        if abs(ds) < 1e-14:
            break
    return s


def modulated_background(grid, eps):
    """Particle-free steady state with a(x) = a(1 + eps sin(2 pi x / L)), full-space Newton."""
    a = P["a"] * (1 + eps * np.sin(2 * np.pi * grid.x / grid.L))
    u = np.array([hom(a=ai) for ai in a]).T
    D2 = grid.D2()
    for it in range(30):
        F = rhs(grid, u, a=a).ravel()
        if np.max(np.abs(F)) < 1e-12:
            break
        J = jac(grid, u, D2)
        u = u + np.linalg.solve(J, -F).reshape(3, grid.N)
    return u, a, float(np.max(np.abs(rhs(grid, u, a=a))))


def dynamic_run(grid, U, eps, x0, T=400.0, dt=0.002, sample=1.0):
    bg, a, bres = modulated_background(grid, eps)
    h0 = hom()
    u0 = bg + np.array([grid.shift(U[s], x0) for s in range(3)]) - h0[:, None]
    stepper = ETDRK4(grid, dt)
    vh = np.fft.fft(u0, axis=1)
    nsteps = int(round(T / dt))
    every = int(round(sample / dt))
    ts, xs = [0.0], [xmin_position(grid, u0[1], x0)]
    for n in range(1, nsteps + 1):
        vh = stepper.step(vh, a)
        if n % every == 0:
            u = np.real(np.fft.ifft(vh, axis=1))
            if not np.all(np.isfinite(u)) or np.any(u <= 0):
                raise RuntimeError("dynamic run lost positivity/finiteness")
            ts.append(n * dt)
            xs.append(xmin_position(grid, u[1], xs[-1]))
    ts, xs = np.array(ts), np.unwrap(np.array(xs) * 2 * np.pi / grid.L) * grid.L / (2 * np.pi)
    m = (ts >= 100.0) & (ts <= T)
    v = np.polyfit(ts[m], xs[m], 1)[0]
    return dict(eps=eps, x0=x0, dt=dt, T=T, v=float(v), bg_residual=bres,
                dadx=float(P["a"] * eps * 2 * np.pi / grid.L * np.cos(2 * np.pi * x0 / grid.L))), ts, xs


# ---------------- main ----------------

def main():
    os.makedirs(OUT, exist_ok=True)
    t_start = time.time()
    logf = open(os.path.join(OUT, "run_log.txt"), "w")

    def log(s):
        print(s)
        logf.write(s + "\n")
        logf.flush()

    log(f"MG192 run started {time.strftime('%Y-%m-%d %H:%M:%S UTC', time.gmtime())}")
    log(f"platform {platform.platform()} python {sys.version.split()[0]} numpy {np.__version__} scipy {scipy.__version__}")
    rec = dict(spec="spec/MG192_gravity_mobility.md", registration_commit="a2ba0fa", params=P, cases={})

    sols = {}
    for (L, N) in [(100.0, 1024), (100.0, 2048), (200.0, 2048)]:
        tag = f"L{int(L)}_N{N}"
        log(f"[{tag}] stationary soliton")
        g = Grid(L, N)
        u0 = initial_guess(g)
        U, res, hist = newton_even(g, u0, log=log)
        log(f"[{tag}] residual {res:.3e}; eigen-decomposition (3N = {3 * N})")
        t0 = time.time()
        lam0, others, psi, Up, _ = spectrum_and_adjoint(g, U)
        log(f"[{tag}] eig took {time.time() - t0:.1f} s")
        mu = mobilities(g, U, psi)
        h0 = hom()
        core = dict(G=float(U[0, N // 2] - h0[0]), X=float(U[1, N // 2] - h0[1]), Y=float(U[2, N // 2] - h0[2]))
        top = others[np.argsort(-np.real(others))[:6]]
        edge = float(np.max(np.abs(psi[:, np.abs(g.x) > 0.4 * L])) / np.max(np.abs(psi)))
        case = dict(L=L, N=N, newton_history=hist, residual=res, lambda0=[float(np.real(lam0)), float(np.imag(lam0))],
                    max_re_other=float(np.max(np.real(others))),
                    leading_other=[[float(np.real(z)), float(np.imag(z))] for z in top],
                    core_amplitudes=core, mobilities=mu, adjoint_edge_ratio=edge)
        rec["cases"][tag] = case
        sols[tag] = (g, U, psi)
        log(f"[{tag}] lambda0 = {lam0:.3e}; max Re(other) = {case['max_re_other']:.6f}")
        log(f"[{tag}] core {core}; mobilities {mu}; adjoint edge ratio {edge:.2e}")

    # ---- P0 ----
    c1, c2, c3 = rec["cases"]["L100_N1024"], rec["cases"]["L100_N2048"], rec["cases"]["L200_N2048"]
    g1, U1, _ = sols["L100_N1024"]
    g2, U2, _ = sols["L100_N2048"]
    prof_diff = float(np.max(np.abs(U2[:, ::2] - U1)))
    p0a = c1["residual"] <= 1e-10 and c2["residual"] <= 1e-10 and prof_diff <= 1e-8
    p0b = all(abs(complex(*c["lambda0"])) <= 1e-7 and c["max_re_other"] <= -1e-4 for c in (c1, c2))
    p0c = all(abs(c1["core_amplitudes"][s] - PUBLISHED_CORE[s]) / abs(PUBLISHED_CORE[s]) <= 0.02 for s in "GXY")
    reldN = {m: abs(c2["mobilities"][m] - c1["mobilities"][m]) / abs(c1["mobilities"][m]) for m in c1["mobilities"]}
    reldL = {m: abs(c3["mobilities"][m] - c1["mobilities"][m]) / abs(c1["mobilities"][m]) for m in c1["mobilities"]}
    p0d = all(v <= 1e-4 for v in reldN.values()) and all(v <= 0.01 for v in reldL.values())
    P0 = dict(a=bool(p0a), b=bool(p0b), c=bool(p0c), d=bool(p0d), profile_diff_N=prof_diff,
              mobility_rel_diff_N=reldN, mobility_rel_diff_L=reldL,
              core_rel_err={s: abs(c1["core_amplitudes"][s] - PUBLISHED_CORE[s]) / abs(PUBLISHED_CORE[s]) for s in "GXY"})
    P0["pass"] = bool(p0a and p0b and p0c and p0d)
    rec["P0"] = P0
    log(f"P0: {P0}")

    mu_a = c1["mobilities"]["mu_a"]
    unc = max(abs(c2["mobilities"]["mu_a"] - mu_a), abs(c3["mobilities"]["mu_a"] - mu_a))
    rec["P1"] = dict(mu_a=mu_a, numerical_uncertainty=unc, ratio=abs(mu_a) / unc if unc > 0 else float("inf"),
                     label=("SUPPORTED" if mu_a < 0 and abs(mu_a) > 100 * unc else
                            "CONTRADICTED" if mu_a > 0 and abs(mu_a) > 100 * unc else "NOT RESOLVED"))
    rec["P3"] = dict(mu_a=mu_a, published=PUBLISHED_MU, rel_diff=abs(mu_a - PUBLISHED_MU) / abs(PUBLISHED_MU),
                     label="CONSISTENT" if abs(mu_a - PUBLISHED_MU) / abs(PUBLISHED_MU) <= 0.15 else "DIFFERENT")
    mb, mw = c1["mobilities"]["mu_b"], c1["mobilities"]["mu_w"]
    # background responses per unit slope of each parameter (slow modulation): dG0/da=1/0.9, dX0/da=1/0.9;
    # dG0/dw=0.1/0.9, dX0/dw=1/0.9; b leaves G0, X0 unchanged (dY0/db = 1/X0).
    rec["P4"] = dict(mu_b=mb, mu_w=mw, mu_b_over_mu_a=mb / mu_a, mu_w_over_mu_a=mw / mu_a,
                     mobility_per_unit_background_G_slope_a_arm=-mu_a * 0.9,
                     note="v = mu * d(param)/dx; for the a arm v = -mu_a*0.9 * (-dG_bg/dx)")
    log(f"P1: {rec['P1']}")
    log(f"P3: {rec['P3']}")
    log(f"P4: {rec['P4']}")
    np.savez_compressed(os.path.join(OUT, "soliton_adjoint_L100_N1024.npz"), x=g1.x, U=U1, psi=sols["L100_N1024"][2])

    # ---- P2 dynamics ----
    if not P0["pass"]:
        log("P0 failed: STOPPED before dynamics, as registered.")
        rec["outcome"] = "STOPPED"
    else:
        runs = []
        series = {}
        L = g1.L
        for eps in (0.001, 0.002, 0.0):
            for x0 in (0.0, L / 8):
                t0 = time.time()
                r, ts, xs = dynamic_run(g1, U1, eps, x0)
                r["predicted_v"] = mu_a * r["dadx"]
                r["wall_s"] = time.time() - t0
                runs.append(r)
                series[f"eps{eps}_x0{x0}"] = np.vstack([ts, xs])
                log(f"dynamic eps={eps} x0={x0}: v={r['v']:.6e} predicted={r['predicted_v']:.6e} ({r['wall_s']:.0f} s)")
        t0 = time.time()
        rdt, ts, xs = dynamic_run(g1, U1, 0.001, 0.0, dt=0.001)
        rdt["wall_s"] = time.time() - t0
        series["dtcheck_eps0.001_x00.0"] = np.vstack([ts, xs])
        log(f"dt check (dt=0.001): v={rdt['v']:.6e}")
        np.savez_compressed(os.path.join(OUT, "trajectories.npz"), **series)
        act = [r for r in runs if r["eps"] > 0]
        nul = [r for r in runs if r["eps"] == 0]
        within = [abs(r["v"] - r["predicted_v"]) / abs(r["predicted_v"]) for r in act]
        smallest = min(abs(r["predicted_v"]) for r in act)
        null_ok = all(abs(r["v"]) < 0.01 * smallest for r in nul)
        lin = [next(r for r in act if r["eps"] == 0.002 and r["x0"] == x0)["v"] /
               next(r for r in act if r["eps"] == 0.001 and r["x0"] == x0)["v"] for x0 in (0.0, L / 8)]
        lin_ok = all(1.8 <= q <= 2.2 for q in lin)
        p2 = all(w <= 0.10 for w in within) and null_ok and lin_ok
        rec["P2"] = dict(runs=runs, dt_check=rdt, rel_err=within, null_ok=bool(null_ok), linearity_ratios=lin,
                         linearity_ok=bool(lin_ok), dt_rel_change=abs(rdt["v"] - runs[0]["v"]) / abs(runs[0]["v"]),
                         **{"pass": bool(p2)})
        log(f"P2: rel_err {within}, null_ok {null_ok}, linearity {lin}, pass {p2}")
        if rec["P1"]["label"] == "SUPPORTED" and p2:
            rec["outcome"] = "SUPPORTED"
        elif mu_a > 0 and p2:
            rec["outcome"] = "CONTRADICTED"
        else:
            rec["outcome"] = "UNRESOLVED"
    rec["wall_seconds"] = time.time() - t_start
    rec["platform"] = dict(platform=platform.platform(), python=sys.version.split()[0], numpy=np.__version__,
                           scipy=scipy.__version__, cores_used=1)
    with open(os.path.join(OUT, "results.json"), "w") as f:
        json.dump(rec, f, indent=1)
    log(f"OUTCOME: {rec['outcome']}  (wall {rec['wall_seconds']:.0f} s)")


if __name__ == "__main__":
    main()

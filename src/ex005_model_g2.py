"""EX005 (exploratory): LaViolette's Model G-2 (A -> X, B <-> G, G + X -> Y + Z, 2X + Y -> 3X, X -> Omega).

Dimensionless (D_G = D_X = 1, D_Y = 12, X decay rate 1, as Model G):
  G_t = Lap G + beta - kappa G - sigma G X
  X_t = Lap X + alpha - sigma G X + X^2 Y - X
  Y_t = 12 Lap Y + sigma G X - X^2 Y
Background: X0 = alpha, G0 = beta / (kappa + sigma X0), Y0 = sigma G0 / X0.  B_eff = sigma G0 plays Brusselator b.
Stationary localized budget: Int phi_X = 0 exactly and (kappa + sigma X0) Int phi_G = -sigma Int phi_G phi_X.
"""
import json
import os
import sys
import time

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ex002_far_field_3d as e  # noqa: E402

ROOT = e.ROOT
OUT = os.path.join(ROOT, "results", "EX005_model_g2")
ALPHA, KAPPA = 140 / 9, 1.0
D = np.diag([1.0, 1.0, 12.0])
FRAC = 29.0 / 30.04408          # same fraction of the Turing threshold as the Model G reference point


def hom(beff, sigma, alpha=ALPHA):
    G0 = beff / sigma
    return np.array([G0, alpha, beff / alpha])


def jac(beff, sigma, kappa=KAPPA, alpha=ALPHA):
    G0, X0, Y0 = hom(beff, sigma, alpha)
    return np.array([[-kappa - sigma * X0, -sigma * G0, 0.0],
                     [-sigma * X0, -(sigma * G0 + 1) + 2 * X0 * Y0, X0 ** 2],
                     [sigma * X0, sigma * G0 - 2 * X0 * Y0, -X0 ** 2]])


def growth(beff, sigma, ks=np.linspace(0.0, 4, 801), **kw):
    J = jac(beff, sigma, **kw)
    return np.array([np.max(np.linalg.eigvals(J - D * k * k).real) for k in ks]), ks


def b_c(sigma, **kw):
    lo, hi = 10.0, 60.0
    for _ in range(50):
        mid = 0.5 * (lo + hi)
        lo, hi = (lo, mid) if growth(mid, sigma, **kw)[0].max() > 0 else (mid, hi)
    return lo


def rates(beff, sigma, **kw):
    """Spatial decay rates of the linearized stationary far field (roots of det(J + D s) = 0, s = -k^2)."""
    J = jac(beff, sigma, **kw)
    ss = np.arange(4.0)
    c = np.polyfit(ss, [np.linalg.det(J + D * s) for s in ss], 3)
    out = []
    for s in np.roots(c):
        k = np.sqrt(complex(s))
        out.append(dict(decay=float(abs(k.real)), wavenumber=float(abs(k.imag))))
    return sorted(out, key=lambda d: d["decay"])


class Radial2(e.Radial):
    def __init__(self, R, dr, bc, beff, sigma, kappa=KAPPA, alpha=ALPHA):
        self.sigma, self.kappa, self.alpha, self.beff = sigma, kappa, alpha, beff
        super().__init__(R, dr, bc, beff)
        self.h0 = hom(beff, sigma, alpha)
        self.beta = self.h0[0] * (kappa + sigma * alpha)

    def F(self, u):
        G, X, Y = u
        s = self.sigma
        return np.array([self.lap(G, 0) + self.beta - self.kappa * G - s * G * X,
                         self.lap(X, 1) + self.alpha - s * G * X + X * X * Y - X,
                         12.0 * self.lap(Y, 2) + s * G * X - X * X * Y])

    def J(self, u):
        G, X, Y = u
        s, I = self.sigma, sp.identity(self.n)
        dg = sp.diags
        return sp.bmat([[self.Lap + dg(-self.kappa - s * X), dg(-s * G), None],
                        [dg(-s * X), self.Lap + dg(-s * G + 2 * X * Y - 1), dg(X * X)],
                        [dg(s * X), dg(s * G - 2 * X * Y), 12.0 * self.Lap - dg(X * X)]]).tocsc()


def stepper(rad, u, dt, T):
    I = sp.identity(3 * rad.n, format="csc")
    for _ in range(int(round(T / dt))):
        Fv = rad.F(u).ravel()
        u = u + spla.spsolve(I - dt * rad.J(u), dt * Fv).reshape(3, rad.n)
        if not np.all(np.isfinite(u)):
            break
    return u


def analyse(rad, u, kap_m, rt):
    phi = u - rad.h0[:, None]
    r = rad.r
    I = [rad.integral(phi[s]) for s in range(3)]
    s = rad.sigma
    lhs = (rad.kappa + s * rad.alpha) * I[0]
    rhs = -s * rad.integral(phi[0] * phi[1])
    w = rad.V
    def corr(a, b):
        return float(np.sum(w * a * b) / np.sqrt(np.sum(w * a * a) * np.sum(w * b * b)))
    out = dict(core=dict(G=float(phi[0, 0]), X=float(phi[1, 0]), Y=float(phi[2, 0])),
               int_G=I[0], int_X=I[1], int_Y=I[2], identity_lhs=lhs, identity_rhs=rhs,
               int_G_rel=I[0] / rad.integral(np.abs(phi[0])),
               corr_GY=corr(phi[0], phi[2]), corr_GX=corr(phi[0], phi[1]),
               G0=float(rad.h0[0]), rel_core_G=float(phi[0, 0] / rad.h0[0]))
    # G tail: linear least squares for r phi_G = C e^{-k_m r} + e^{-k_o r}(A cos w r + B sin w r), rates from linear theory
    osc = [x for x in rt if x["wavenumber"] > 1e-6][0]
    lo, hi = 8.0, rad.R - 8.0
    m = (r > lo) & (r < hi)
    y = r[m] * phi[0][m]
    basis = np.vstack([np.exp(-kap_m * r[m]), np.exp(-osc["decay"] * r[m]) * np.cos(osc["wavenumber"] * r[m]),
                       np.exp(-osc["decay"] * r[m]) * np.sin(osc["wavenumber"] * r[m])]).T
    wts = 1 / (np.abs(basis).max(axis=1) + 1e-300)
    coef = np.linalg.lstsq(basis * wts[:, None], y * wts, rcond=None)[0]
    fit = basis @ coef
    out["tail_G"] = dict(C_monotone=float(coef[0]), osc_amp=float(np.hypot(coef[1], coef[2])), kappa_monotone=kap_m,
                         kappa_osc=osc["decay"], k_osc=osc["wavenumber"], window=[lo, hi],
                         rel_rms=float(np.sqrt(np.mean(((fit - y) * wts) ** 2) / np.mean((y * wts) ** 2))),
                         rphiG_at_20=float(np.interp(20.0, r, r * phi[0])))
    return out


def main():
    os.makedirs(OUT, exist_ok=True)
    t0 = time.time()
    d = np.load(os.path.join(ROOT, "results", "EX002_far_field_3d", "profile_R60_neumann_std.npz"))
    hg = e.hom(14.0, 29.0)
    rec = dict(note="EXPLORATORY", alpha=ALPHA, frac=FRAC, runs=[])
    for kappa, sigmas in ((0.0, (0.001, 0.002, 0.003, 0.005, 0.01, 0.02, 0.03, 0.05, 0.1, 0.15, 0.2, 0.3, 0.5, 0.7, 1.0)), (1.0, (0.01, 0.1, 1.0))):
        prev = None
        for sigma in sigmas:
            bc = b_c(sigma, kappa=kappa)
            beff = FRAC * bc
            rt = rates(beff, sigma, kappa=kappa)
            mono = [x for x in rt if x["wavenumber"] < 1e-6]
            kap_m = mono[0]["decay"] if mono else rt[0]["decay"]
            R = float(min(150.0, max(60.0, 25 + 8 / kap_m)))
            rad = Radial2(R, 0.05, "neumann", beff, sigma, kappa=kappa)
            if prev is None:
                rp, pp = d["r"], d["u"] - hg[:, None]
                pp[0] = 0.0
            else:
                rp, pp = prev
            u = np.array([rad.h0[s] + np.interp(rad.r, rp, pp[s], right=0.0) for s in range(3)])
            u = np.maximum(u, 1e-3)
            u, res, it, ok = rad.newton(u, tol=1e-9 * max(1.0, rad.h0[0] / 100), maxit=40)
            how = "newton"
            if not ok or u[1, 0] - rad.h0[1] > -1:
                u = np.array([rad.h0[s] + np.interp(rad.r, rp, pp[s], right=0.0) for s in range(3)])
                u = stepper(rad, np.maximum(u, 1e-3), 0.1, 200.0)
                u, res, it, ok = rad.newton(u, tol=1e-9 * max(1.0, rad.h0[0] / 100), maxit=40)
                how = "settle+newton"
            alive = bool(ok and u[1, 0] - rad.h0[1] < -1)
            run = dict(kappa=kappa, sigma=sigma, b_c=bc, beff=beff, beta=rad.beta, R=R, rates=rt, ok=bool(ok),
                       residual=float(res), alive=alive, how=how)
            if alive:
                run.update(analyse(rad, u, kap_m, rt))
                prev = (rad.r, u - rad.h0[:, None])
                np.savez_compressed(os.path.join(OUT, f"profile_k{kappa:g}_s{sigma:g}.npz"), r=rad.r, u=u, h0=rad.h0)
            print({k: v for k, v in run.items() if k != "rates"}, flush=True)
            rec["runs"].append(run)
            json.dump(rec, open(os.path.join(OUT, "results.json"), "w"), indent=1, default=float)
    rec["wall_seconds"] = time.time() - t0
    json.dump(rec, open(os.path.join(OUT, "results.json"), "w"), indent=1, default=float)
    print("done", rec["wall_seconds"])


if __name__ == "__main__":
    main()

"""MG6 in-run quadrature cross-check (proposal §8, audit scope): a second, separately written implementation of
P0(1.5), Pc(1.5), M- and M+. Differences from mg6: scipy LinearNDInterpolator on the 3 x 3 tiled node positions (its own
Qhull triangulation) instead of the driver's mesh point location; a 2x finer midpoint grid (796 x 804); point-mask discs
(a cell counts if its midpoint is inside). Shares no code with mg6 apart from numpy and scipy."""
import numpy as np
from scipy.interpolate import LinearNDInterpolator

NXC, NYC, R_CHECK = 796, 804, 1.5


def _wrap(d, L):
    return d - L * np.round(d / L)


def cross_check(pos, phiG, Lx, Ly, c, kappa, shift, rho0, R=R_CHECK):
    pos = np.mod(np.asarray(pos, float), (Lx, Ly))
    tiles = [(a * Lx, b * Ly) for a in (-1, 0, 1) for b in (-1, 0, 1)]
    P = np.concatenate([pos + t for t in tiles])
    V = np.tile(np.asarray(phiG, float), 9)
    interp = LinearNDInterpolator(P, V)
    dx, dy = Lx / NXC, Ly / NYC
    xs = (np.arange(NXC) + 0.5) * dx
    ys = (np.arange(NYC) + 0.5) * dy
    X, Y = np.meshgrid(xs, ys, indexing="ij")
    X, Y = X.ravel(), Y.ravel()
    dA = Lx * Ly / (NXC * NYC)
    g = interp(X, Y)
    inside = np.hypot(_wrap(X - c[0], Lx), _wrap(Y - c[1], Ly)) < R
    out = dict(M_minus=float(np.sum(np.maximum(-g, 0)) * dA), M_plus=float(np.sum(np.maximum(g, 0)) * dA))
    if kappa == 0:
        out.update(P0=0.0, Pc=0.0)
        return out
    if shift is None:
        ge = g
    else:
        ge = interp(np.mod(X + shift[0], Lx), np.mod(Y + shift[1], Ly))
    f = np.exp(kappa * ge)
    out["P0"] = float(rho0 * np.sum(f[inside] - 1.0) * dA)
    out["Pc"] = float(rho0 * np.sum(f[inside]) * dA - inside.sum() * dA * rho0 * np.sum(f) * dA / (Lx * Ly))
    return out

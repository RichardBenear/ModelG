"""MG4 (Phase 3): the deterministic graph ODE (Omega -> infinity mean field of rules v1.1 on a graph) and the
continuum 2D soliton used as the initial condition and in C3(b).

Equations (potentials phi_S = S - S0), with L_G = (A - D_deg) / h^2:
  dphi_G/dt = L_G phi_G - q phi_G + g phi_X
  dphi_X/dt = dx L_G phi_X + p phi_G - (1+b) phi_X + (X^2 Y - X0^2 Y0) + chi
  dphi_Y/dt = dy L_G phi_Y + b phi_X - (X^2 Y - X0^2 Y0)
chi is the published seed (eq. 18) at a centre c, used only by the descriptive seeded-formation runs.
"""
import os, sys
import numpy as np
import scipy.sparse as sp
from scipy.integrate import solve_ivp
sys.path.insert(0, os.path.dirname(__file__))
from lattice import steady_state, P17
import graphs as gr

CONT_FILE = os.path.join(os.path.dirname(__file__), "..", "results", "MG4", "part0", "continuum2d.npz")


def continuum2d(R=50.0, h=0.025, t_end=300.0, path=CONT_FILE):
    """MG1 continuum 2D stationary soliton at b = 29 (src/mg_continuum.run, seed scale 1), cached.
    Returns dict(r, phi (3, n), lap (3, n)) with lap = the solver's own conservative radial Laplacian of phi."""
    if os.path.exists(path):
        d = np.load(path)
        return {k: d[k] for k in d.files}
    from mg_continuum import run
    grid, P, sol, rhs = run(2, R=R, h=h, t_end=t_end)
    N = grid.N
    y = sol.y[:, -1]
    phi = np.vstack([y[:N], y[N:2 * N], y[2 * N:]])
    lap = np.vstack([grid.lap(f) for f in phi])
    out = dict(r=grid.ri.copy(), phi=phi, lap=lap, max_dphidt=np.array(np.abs(rhs(t_end, y)).max()))
    os.makedirs(os.path.dirname(path), exist_ok=True)
    np.savez(path, **out)
    return out


def node_r(g, c):
    """Periodic Euclidean distance of every node from centre c, and the minimal-image offset vectors."""
    d = gr.min_image(g["pos"] - np.asarray(c, float), g["Lx"], g["Ly"])
    return np.hypot(d[:, 0], d[:, 1]), d


def interp_radial(r_nodes, r, f):
    """Cubic interpolation of a radial profile onto node distances (beyond the grid: 0)."""
    from scipy.interpolate import CubicSpline
    cs = CubicSpline(r, f, bc_type=((1, 0.0), "not-a-knot"))
    out = cs(np.clip(r_nodes, 0, r[-1]))
    out[r_nodes > r[-1]] = 0.0
    return out


def soliton_ic(g, c, cont=None):
    """Initial condition: continuum 2D profiles of G, X, Y at distance |x - c| (proposal §5)."""
    cont = cont or continuum2d()
    rn, _ = node_r(g, c)
    return np.concatenate([interp_radial(rn, cont["r"], cont["phi"][s]) for s in range(3)])


def make_system(g, h, P=None, seed_centre=None, seed_scale=0.0):
    P = dict(P17 if P is None else P)
    G0, X0, Y0 = steady_state(P)
    b, gg, p, q = P["b"], P["g"], P["p"], P["q"]
    N = g["N"]
    L = (gr.laplacian(g) / h**2).tocsr()
    rs = node_r(g, seed_centre)[0] if seed_scale > 0 else None

    def rhs(t, y):
        pG, pX, pY = y[:N], y[N:2 * N], y[2 * N:]
        X, Y = X0 + pX, Y0 + pY
        r4 = X**2 * Y - X0**2 * Y0
        dX = P["dx"] * (L @ pX) + p * pG - (1 + b) * pX + r4
        if seed_scale > 0:
            dX = dX - seed_scale * np.exp(-rs**2 / 2) * np.exp(-(t - 10.0) ** 2 / 18.0)
        return np.concatenate([L @ pG - q * pG + gg * pX, dX, P["dy"] * (L @ pY) + b * pX - r4])

    I = sp.identity(N, format="csr")

    def jac(t, y):
        X, Y = X0 + y[N:2 * N], Y0 + y[2 * N:]
        return sp.bmat([[L - q * I, gg * I, None],
                        [p * I, P["dx"] * L + sp.diags(-(1 + b) + 2 * X * Y), sp.diags(X**2)],
                        [None, sp.diags(b - 2 * X * Y), P["dy"] * L - sp.diags(X**2)]], format="csc")

    return rhs, jac, P


def integrate(g, h, y0, t_end=300.0, t_eval=None, P=None, seed_centre=None, seed_scale=0.0,
              rtol=1e-8, atol=1e-10):
    """BDF with sparse Jacobian, as lattice_ode_1d (with its first-step retries). Returns (sol, rhs)."""
    rhs, jac, P = make_system(g, h, P, seed_centre, seed_scale)
    sol = None
    for fs in (None, 1e-3, 1e-2):
        sol = solve_ivp(rhs, (0.0, t_end), y0, method="BDF", jac=jac, rtol=rtol, atol=atol,
                        max_step=0.5, first_step=fs, t_eval=t_eval)
        sol.first_step_used = fs
        if sol.success:
            break
    if not sol.success:
        raise RuntimeError(sol.message)
    return sol, rhs


def balances(y, dydt, N, P=None):
    """§5.4 node-sum balances. Returns the two residuals and their term-magnitude sums M:
       r1 = dI_G/dt - (-q I_G + g I_X);  r2 = dI_tot/dt - (-(q-p) I_G - (1-g) I_X)."""
    P = dict(P17 if P is None else P)
    q, gg, p = P["q"], P["g"], P["p"]
    I = [y[k * N:(k + 1) * N].sum() for k in range(3)]
    dI = [dydt[k * N:(k + 1) * N].sum() for k in range(3)]
    r1 = dI[0] - (-q * I[0] + gg * I[1])
    M1 = abs(dI[0]) + q * abs(I[0]) + gg * abs(I[1])
    r2 = sum(dI) - (-(q - p) * I[0] - (1 - gg) * I[1])
    M2 = sum(abs(v) for v in dI) + abs(q - p) * abs(I[0]) + (1 - gg) * abs(I[1])
    return dict(I=I, dI=dI, r1=r1, M1=M1, r2=r2, M2=M2,
                ok=bool(abs(r1) <= 1e-10 + 1e-9 * M1 and abs(r2) <= 1e-10 + 1e-9 * M2))

"""MG003: lattice tools for Phase 2 (rules v1.1 on a fixed periodic 1D lattice).

- `tau_leap_1d`: binomial/Poisson tau-leaping of rules R1, R1r, R2, R2r, R3, R4, R5 + hopping D_S,
  plus the optional seed-forcing rule R_seed ({v,X} -> 0, zero-order propensity Omega*|chi(v,t)|).
  Lie splitting: reactions, then hopping. Hopping: each token hops left/right with probability
  p_S = d_S*tau/h^2 per direction (Binomial(n, 2p_S), then split 50/50). The mean field is
  exactly the lattice Laplacian. Requires 2 p_S <= 1 (spec uses <= 0.8).
- `lattice_ode_*`: the deterministic mean-field lattice ODE (same lattice, same Laplacian).
- `lna_homogeneous`: linear-noise prediction of the stationary single-node variance about the
  homogeneous state (Lyapunov equation per lattice Fourier mode).
Concentrations are n/Omega; Omega = tokens per node at unit concentration.
"""
import numpy as np
import numba as nb
import scipy.sparse as sp
from scipy.integrate import solve_ivp
from scipy.linalg import solve_continuous_lyapunov

P17 = dict(a=14.0, b=29.0, g=0.1, p=1.0, q=1.0, s=0.0, u=0.0, w=0.0, dx=1.0, dy=12.0)


def steady_state(P):
    G0 = (P["a"] + P["g"] * P["w"]) / (P["q"] - P["g"] * P["p"])
    X0 = (P["p"] * P["a"] + P["q"] * P["w"]) / (P["q"] - P["g"] * P["p"])
    Y0 = X0 * (P["s"] * X0**2 + P["b"]) / (X0**2 + P["u"])
    return G0, X0, Y0


# ------------------------------------------------------------------ stochastic kernel
@nb.njit(cache=True)
def _seed(x, t, scale):
    return scale * np.exp(-x * x / 2.0) * np.exp(-(t - 10.0) ** 2 / 18.0)   # |chi|


@nb.njit(cache=True)
def tau_leap_1d(nG, nX, nY, Om, a, b, g, p, q, dx, dy, h, tau, nsteps, t0,
                seed_scale, xs, rec_every, rng_seed):
    """Advance counts in place; returns recorded X snapshots (every rec_every steps) and times."""
    np.random.seed(rng_seed)
    N = nG.size
    pg, px, py = tau / h**2, dx * tau / h**2, dy * tau / h**2
    nrec = nsteps // rec_every
    recX = np.zeros((nrec, N), dtype=np.int64)
    recG = np.zeros((nrec, N), dtype=np.int64)
    recY = np.zeros((nrec, N), dtype=np.int64)
    rect = np.zeros(nrec)
    PG2X = 1.0 - np.exp(-q * tau)                 # G loss: p (to X) + (q-p) (to reservoir)
    fracGX = p / q
    PXL = 1.0 - np.exp(-(1.0 + b) * tau)         # X loss: g (to G) + b (to Y) + (1-g) (sink)
    fG, fY = g / (1.0 + b), b / (1.0 + b)
    inv_Om2 = 1.0 / (Om * Om)
    lamA = Om * a * tau
    dG = np.zeros(N, dtype=np.int64); dXa = np.zeros(N, dtype=np.int64); dYa = np.zeros(N, dtype=np.int64)
    t = t0
    k = 0
    for step in range(nsteps):
        # ---- reactions
        for i in range(N):
            gi, xi, yi = nG[i], nX[i], nY[i]
            src = np.random.poisson(lamA)                                   # R1
            gl = np.random.binomial(gi, PG2X) if gi > 0 else 0              # R2 + R1r
            g2x = np.random.binomial(gl, fracGX) if gl > 0 else 0
            xl = np.random.binomial(xi, PXL) if xi > 0 else 0               # R2r, R3, R5
            x2g = np.random.binomial(xl, fG) if xl > 0 else 0
            rest = xl - x2g
            x2y = np.random.binomial(rest, fY / (1.0 - fG)) if rest > 0 else 0
            lam4 = xi * (xi - 1) * yi * inv_Om2 * tau                       # R4: Y -> X
            y2x = np.random.poisson(lam4) if lam4 > 0 else 0
            if y2x > yi:
                y2x = yi
            sd = 0
            if seed_scale > 0.0:                                            # R_seed
                lam_s = Om * _seed(xs[i], t, seed_scale) * tau
                if lam_s > 1e-300:
                    sd = np.random.poisson(lam_s)
            nG[i] = gi + src - gl + x2g
            nX[i] = xi + g2x - xl + y2x
            nY[i] = yi + x2y - y2x
            if sd > nX[i]:
                sd = nX[i]
            nX[i] -= sd
        # ---- hopping (periodic)
        for arr_id in range(3):
            if arr_id == 0:
                arr = nG; pp = pg; dd = dG
            elif arr_id == 1:
                arr = nX; pp = px; dd = dXa
            else:
                arr = nY; pp = py; dd = dYa
            for i in range(N):
                dd[i] = 0
            for i in range(N):
                ni = arr[i]
                if ni > 0:
                    mv = np.random.binomial(ni, 2.0 * pp)
                    if mv > 0:
                        lf = np.random.binomial(mv, 0.5)
                        dd[i] -= mv
                        dd[(i - 1) % N] += lf
                        dd[(i + 1) % N] += mv - lf
            for i in range(N):
                arr[i] += dd[i]
        t += tau
        if (step + 1) % rec_every == 0 and k < nrec:
            recX[k, :] = nX; recG[k, :] = nG; recY[k, :] = nY; rect[k] = t
            k += 1
    return recG, recX, recY, rect


# ------------------------------------------------------------------ deterministic lattice ODE
def lap_periodic(N, h):
    e = np.ones(N)
    L = sp.diags([e[:-1], -2 * e, e[:-1]], [-1, 0, 1], format="lil")
    L[0, N - 1] = 1; L[N - 1, 0] = 1
    return L.tocsr() / h**2


def lattice_ode_1d(P, N, h, t_end, seed_scale=1.0, y0=None, rtol=1e-8, atol=1e-10):
    """Mean-field of tau_leap_1d (potentials, periodic, x centred at 0)."""
    G0, X0, Y0 = steady_state(P)
    x = (np.arange(N) - N // 2) * h
    L = lap_periodic(N, h)
    b, g, p, q = P["b"], P["g"], P["p"], P["q"]

    def rhs(t, y):
        pG, pX, pY = y[:N], y[N:2 * N], y[2 * N:]
        X, Y = X0 + pX, Y0 + pY
        chi = -seed_scale * np.exp(-x**2 / 2) * np.exp(-(t - 10) ** 2 / 18)
        return np.concatenate([L @ pG - q * pG + g * pX,
                               P["dx"] * (L @ pX) + p * pG - (1 + b) * pX + X**2 * Y - X0**2 * Y0 + chi,
                               P["dy"] * (L @ pY) + b * pX - (X**2 * Y - X0**2 * Y0)])

    def jac(t, y):
        X, Y = X0 + y[N:2 * N], Y0 + y[2 * N:]
        I = sp.identity(N, format="csr")
        return sp.bmat([[L - q * I, g * I, None],
                        [p * I, P["dx"] * L + sp.diags(-(1 + b) + 2 * X * Y), sp.diags(X**2)],
                        [None, sp.diags(b - 2 * X * Y), P["dy"] * L - sp.diags(X**2)]], format="csc")

    y0 = np.zeros(3 * N) if y0 is None else y0
    sol = None
    for fs in (None, 1e-3, 1e-2):
        sol = solve_ivp(rhs, (0, t_end), y0, method="BDF", jac=jac, rtol=rtol, atol=atol,
                        max_step=0.5, first_step=fs)
        if sol.success:
            break
    y = sol.y[:, -1]
    return x, y[:N], y[N:2 * N], y[2 * N:], float(np.abs(rhs(t_end, y)).max())


# ------------------------------------------------------------------ linear noise approximation
def lna_homogeneous(P, N, h):
    """Stationary single-node variance (concentration units, times Omega) about the homogeneous
    state on a periodic 1D lattice, including reaction and hopping noise."""
    G0, X0, Y0 = steady_state(P)
    a, b, g, p, q = P["a"], P["b"], P["g"], P["p"], P["q"]
    J = np.array([[-q, g, 0.0], [p, -(1 + b) + 2 * X0 * Y0, X0**2], [0.0, b - 2 * X0 * Y0, -X0**2]])
    # reaction channels: (stoichiometry, propensity / Omega)
    ch = [((1, 0, 0), a), ((-1, 0, 0), (q - p) * G0), ((-1, 1, 0), p * G0), ((1, -1, 0), g * X0),
          ((0, -1, 1), b * X0), ((0, -1, 0), (1 - g) * X0), ((0, 1, -1), X0**2 * Y0)]
    Br = sum(np.outer(v, v) * w for v, w in ch)
    D = np.diag([1.0, P["dx"], P["dy"]]); c = np.array([G0, X0, Y0])
    var = np.zeros(3)
    ks = 2 * np.pi * np.arange(N) / (N * h)
    lam = -(4 / h**2) * np.sin(ks * h / 2) ** 2          # lattice Laplacian eigenvalues
    for lm in lam:
        Jm = J + lm * D
        Bm = Br + np.diag(2 * np.diag(D) * c * (-lm))     # hopping noise
        C = solve_continuous_lyapunov(Jm, -Bm)
        var += np.diag(C)
    return var / N   # Omega * Var(concentration) per node, for G, X, Y


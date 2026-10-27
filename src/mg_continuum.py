"""MG002: continuum solver for dimensionless Model G (Pulver & LaViolette 2013, eqs. 13a-c).

Radially symmetric method of lines, conservative finite-volume Laplacian:
    lap(phi)_i = [A_{i+1/2}(phi_{i+1}-phi_i) - A_{i-1/2}(phi_i-phi_{i-1})] / (h V_i),
    A = r^(n-1), V_i = control-volume measure / omega_n.
Nodes r_i = i*h, i = 0..N, r_N = R (Dirichlet phi = 0). Node 0 control volume [0, h/2]
(symmetry: zero flux at r = 0). n = 1 is the symmetric half of the 1D domain [-R, R].

The discrete scheme conserves exactly, so the integrated balances
    dI_G/dt = F_G - q I_G + g I_X  and  d(sum I)/dt = sum F - (q-p) I_G - (1-g) I_X
hold to round-off at every instant, with F_S the discrete outer-face flux. What the
diagnostic measures is therefore the size of the dI/dt and F terms, not discretization error
in the balance itself (truncation error is estimated by grid refinement).

Integrals use geometric weights omega_n * V_i with omega_1 = 2 (full line), omega_2 = 2 pi,
omega_3 = 4 pi.
"""
import json
import numpy as np
import scipy.sparse as sp
from scipy.integrate import solve_ivp

P17 = dict(a=14.0, b=29.0, g=0.1, p=1.0, q=1.0, s=0.0, u=0.0, w=0.0, dx=1.0, dy=12.0)
OMEGA = {1: 2.0, 2: 2 * np.pi, 3: 4 * np.pi}


def steady_state(P):
    a, b, g, p, q, s, u, w = (P[k] for k in "abgpqsuw")
    G0 = (a + g * w) / (q - g * p)
    X0 = (p * a + q * w) / (q - g * p)
    Y0 = X0 * (s * X0**2 + b) / (X0**2 + u)
    return G0, X0, Y0


class Grid:
    def __init__(self, n, R, h):
        self.n, self.R, self.h = n, R, h
        self.N = int(round(R / h))
        self.r = np.arange(self.N + 1) * h          # includes boundary node N
        ri = self.r[:-1]                             # interior unknown nodes 0..N-1
        rf = (np.arange(self.N) + 0.5) * h           # faces i+1/2, i = 0..N-1
        self.Af = rf ** (n - 1)
        lo = np.concatenate([[0.0], rf[:-1]])
        self.V = (rf**n - lo**n) / n                 # exact measure / omega_n
        self.w = OMEGA[n] * self.V                   # quadrature weights for interior nodes
        self.ri = ri

    def lap(self, f):
        """f: values at nodes 0..N-1 (node N is 0 by Dirichlet)."""
        fe = np.append(f, 0.0)
        flux = self.Af * (fe[1:] - fe[:-1]) / self.h  # outward flux at faces
        out = flux.copy()
        out[1:] -= flux[:-1]
        return out / self.V

    def lap_matrix(self):
        N, h = self.N, self.h
        main = np.zeros(N); up = np.zeros(N - 1); dn = np.zeros(N - 1)
        main -= self.Af / h
        main[1:] -= self.Af[:-1] / h
        up += self.Af[:-1] / h
        dn += self.Af[:-1] / h
        L = sp.diags([dn, main, up], [-1, 0, 1], format="csr")
        return sp.diags(1.0 / self.V) @ L

    def outer_flux(self, f, d):
        """d * outward flux through r = R - h/2 (times omega_n): F = d * omega * A * dphi/dr."""
        return d * OMEGA[self.n] * self.Af[-1] * (0.0 - f[-1]) / self.h


def seed(r, t, scale=1.0):
    return -scale * np.exp(-r**2 / 2) * np.exp(-(t - 10.0) ** 2 / 18.0)


def make_rhs(grid, P, seed_scale):
    G0, X0, Y0 = steady_state(P)
    b, g, p, q, s, u, dx, dy = (P[k] for k in ("b", "g", "p", "q", "s", "u", "dx", "dy"))
    N = grid.N

    def reactions(pG, pX, pY):
        G, X, Y = G0 + pG, X0 + pX, Y0 + pY
        fG = -q * G + g * X + P["a"]
        fX = p * G - (1 + b) * X + u * Y - s * X**3 + X**2 * Y + P["w"]
        fY = b * X - u * Y + s * X**3 - X**2 * Y
        return fG, fX, fY

    def rhs(t, y):
        pG, pX, pY = y[:N], y[N:2 * N], y[2 * N:]
        fG, fX, fY = reactions(pG, pX, pY)
        dG = grid.lap(pG) + fG
        dX = dx * grid.lap(pX) + fX + seed(grid.ri, t, seed_scale)
        dY = dy * grid.lap(pY) + fY
        return np.concatenate([dG, dX, dY])

    Lm = grid.lap_matrix()

    def jac(t, y):
        pX, pY = y[N:2 * N], y[2 * N:]
        X, Y = X0 + pX, Y0 + pY
        I = sp.identity(N, format="csr")
        dfXdX = sp.diags(-(1 + b) - 3 * s * X**2 + 2 * X * Y)
        dfXdY = sp.diags(u + X**2)
        dfYdX = sp.diags(b + 3 * s * X**2 - 2 * X * Y)
        dfYdY = sp.diags(-u - X**2)
        return sp.bmat([[Lm - q * I, g * I, None],
                        [p * I, dx * Lm + dfXdX, dfXdY],
                        [None, dfYdX, dy * Lm + dfYdY]], format="csc")

    return rhs, jac, reactions


def run(n, R=50.0, h=0.05, t_end=100.0, seed_scale=1.0, rtol=1e-8, atol=1e-10,
        P=None, t_eval=None, y0=None, t_start=0.0):
    P = dict(P17 if P is None else P)
    grid = Grid(n, R, h)
    rhs, jac, _ = make_rhs(grid, P, seed_scale)
    if y0 is None:
        y0 = np.zeros(3 * grid.N)
    sol = None
    # scipy BDF step control can fail when the state sits at round-off level (zero-seed runs).
    # Retry with explicit first steps; the retry used is recorded on the solution object.
    for fs in (None, 1e-3, 1e-2, 1e-4):
        sol = solve_ivp(rhs, (t_start, t_end), y0, method="BDF", jac=jac, rtol=rtol, atol=atol,
                        max_step=0.5, t_eval=t_eval, first_step=fs)
        sol.first_step_used = fs
        if sol.success:
            break
    if not sol.success:
        raise RuntimeError(sol.message)
    return grid, P, sol, rhs


# ---------------------------------------------------------------- measurements
def zeros(r, f):
    """Linearly interpolated zero crossings of f(r)."""
    sgn = np.sign(f)
    idx = np.where(sgn[:-1] * sgn[1:] < 0)[0]
    return r[idx] - f[idx] * (r[idx + 1] - r[idx]) / (f[idx + 1] - f[idx])


def core_measures(grid, f):
    r = grid.ri
    z = zeros(r, f)
    r0 = z[0] if len(z) else np.nan
    m = r <= r0
    # partial last cell: integrate linear interpolation inside [r_k, r0] by trapezoid on a fine sub-grid
    rr = np.linspace(0, r0, 4001)
    ff = np.interp(rr, r, f)
    wgt = OMEGA[grid.n] * rr ** (grid.n - 1)
    core_int = np.trapezoid(ff * wgt, rr)
    rms = np.sqrt(np.trapezoid(ff * rr**2 * wgt, rr) / core_int)
    full_int = float(np.sum(grid.w * f))
    return dict(amp=float(f[0]), r0=float(r0), rms=float(rms), core_int=float(core_int),
                full_int=full_int, zeros=z)


def wavelength(grid, fields, r_max_frac=0.6):
    """Turing wavelength = 2 x mean spacing of successive zeros, using zeros beyond the core
    (excluding the first zero) and inside r < r_max_frac * R, pooled over the three fields."""
    sp_all, per = [], {}
    for k, f in fields.items():
        z = zeros(grid.ri, f)
        z = z[(z < r_max_frac * grid.R)]
        d = np.diff(z)[1:]           # drop the core-to-first-ring spacing
        per[k] = 2 * float(np.mean(d)) if len(d) else np.nan
        sp_all.extend(d)
    sp_all = np.array(sp_all)
    return 2 * float(np.mean(sp_all)), per, float(2 * np.std(sp_all))


def analyze(grid, P, sol, rhs):
    N = grid.N
    y = sol.y[:, -1]
    F = dict(G=y[:N], X=y[N:2 * N], Y=y[2 * N:])
    dydt = rhs(sol.t[-1], y)
    dF = dict(G=dydt[:N], X=dydt[N:2 * N], Y=dydt[2 * N:])
    out = dict(n=grid.n, R=grid.R, h=grid.h, t=float(sol.t[-1]))
    out["fields"] = {k: {kk: vv for kk, vv in core_measures(grid, f).items() if kk != "zeros"}
                     for k, f in F.items()}
    lam, per, lam_sd = wavelength(grid, F)
    out["lambda"] = lam; out["lambda_per_field"] = per; out["lambda_spread_2sd"] = lam_sd
    out["max_abs_dphidt"] = float(max(np.abs(v).max() for v in dF.values()))
    d = dict(G=1.0, X=P["dx"], Y=P["dy"])
    I = {k: float(np.sum(grid.w * F[k])) for k in F}
    dI = {k: float(np.sum(grid.w * dF[k])) for k in F}
    Fl = {k: float(grid.outer_flux(F[k], d[k])) for k in F}
    g, p, q = P["g"], P["p"], P["q"]
    out["integrals"] = dict(
        I=I, dIdt=dI, flux=Fl,
        balance_G_residual=dI["G"] - (Fl["G"] - q * I["G"] + g * I["X"]),
        balance_sum_residual=sum(dI.values()) - (sum(Fl.values()) - (q - p) * I["G"] - (1 - g) * I["X"]))
    G0, X0, Y0 = steady_state(P)
    out["min_concentration"] = dict(G=float(G0 + F["G"].min()), X=float(X0 + F["X"].min()),
                                    Y=float(Y0 + F["Y"].min()))
    return out, F


if __name__ == "__main__":
    import sys
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 1
    h = float(sys.argv[2]) if len(sys.argv) > 2 else 0.05
    g_, P_, s_, rhs_ = run(n, h=h)
    res, _ = analyze(g_, P_, s_, rhs_)
    print(json.dumps(res, indent=1))

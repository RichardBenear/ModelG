"""MG6 (Phase 5, first campaign): G-coupled spatial growth around a particle.

Registered design: spec/MG6_proposal.md draft 3 (frozen at 106b3e2), GATES.md MG6.

Contents
- gen_candidates / candidates_hash / regen_hash_fresh_process: the dominating candidate process with Exp(nu) lifetimes
  for the initial nodes and every candidate (proposal §3.1), and its canonical SHA-256.
- draw_shift: the fixed periodic offset of the shifted null (§3.3).
- Coupled: the per-run coupled birth-death state (alive ids, death times, acceptance bitmap) and tick_edits(), which
  turns one tick into MG5's event format for evolve.apply_tick (deaths first, accepted births in candidate order).
- Grid / disc_weights / field_on_grid: the coordinate-area quadrature of §4.3 (nx x ny midpoints that tile the torus;
  exact weights for whole cells, an 8 x 8 sub-lattice for cells a disc boundary cuts).
- measure: the per-check measurements of §4.2 (counts, both baselines, annuli, radial phi_G, M-, M+, M, P0, Pc).
"""
import hashlib, os, sys
import numpy as np
sys.path.insert(0, os.path.dirname(__file__))
import graphs as gr

GATE = 6
PHI_CAP = 0.10
KAPPA_MAX = 4.0
FBOUND = float(np.exp(KAPPA_MAX * PHI_CAP))         # e^0.4: one dominating bound for every arm
SHIFT_MIN = 6.0
RADII = (1.0, 1.5, 2.0, 3.0, 5.0, 8.0)
ANN = (0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0, 8.0)
FINE = tuple(np.round(np.arange(0.0, 3.0001, 0.25), 2))
NX, NY, SUB = 398, 402, 8


def phi_lim(kappa):
    """Largest nodal phi_G allowed at a boundary (proposal §3.2): exp(kappa phi) <= F  <=>  phi <= 0.4 / kappa."""
    return np.inf if kappa == 0 else KAPPA_MAX * PHI_CAP / kappa


# ------------------------------------------------------------------ candidate process (§3.1)
def gen_candidates(N0, nu, t_end, Lx, Ly, seed, Fb=FBOUND):
    """Initial lifetimes, then the candidates of a Poisson process at rate nu * rho0 * F per unit area and time on the
    torus x [0, t_end] (expected number nu * N0 * F * t_end). Candidate j has id N0 + j. Draw order is fixed:
    life0 (N0), n ~ Poisson, tb (sorted uniforms), x, y, u, ell."""
    rng = np.random.default_rng(seed)
    life0 = rng.exponential(1.0 / nu, int(N0))
    n = int(rng.poisson(nu * N0 * Fb * t_end))
    tb = np.sort(rng.random(n) * t_end)
    x = rng.random(n) * Lx
    y = rng.random(n) * Ly
    u = rng.random(n)
    ell = rng.exponential(1.0 / nu, n)
    return dict(N0=int(N0), life0=life0, tb=tb, x=x, y=y, u=u, ell=ell)


CAND_FIELDS = ("life0", "tb", "x", "y", "u", "ell")


def candidates_hash(c):
    """SHA-256 over N0 (little-endian int64) and the six float64 arrays, little-endian, in CAND_FIELDS order."""
    h = hashlib.sha256()
    h.update(np.array([c["N0"]], dtype="<i8").tobytes())
    for f in CAND_FIELDS:
        h.update(np.ascontiguousarray(c[f], dtype="<f8").tobytes())
    return h.hexdigest()


def regen_hash_fresh_process(args):
    """Regenerate the candidate list in a fresh Python process and return its hash. args = (N0, nu, t_end, Lx, Ly,
    seed)."""
    import subprocess
    N0, nu, t_end, Lx, Ly, seed = args
    cmd = [sys.executable, os.path.abspath(__file__), "--regen-hash", str(int(N0)), repr(float(nu)), repr(float(t_end)),
           repr(float(Lx)), repr(float(Ly)), str(int(seed))]
    return subprocess.check_output(cmd).decode().strip().splitlines()[-1]


def draw_shift(seed, Lx, Ly, smin=SHIFT_MIN):
    """Uniform offset on the torus conditioned on minimum-image length >= smin (first accepted draw)."""
    rng = np.random.default_rng(seed)
    while True:
        s = np.array([rng.random() * Lx, rng.random() * Ly])
        if np.hypot(*gr.min_image(s, Lx, Ly)) >= smin:
            return s


class Coupled:
    """State of the coupled birth-death process for one run. kappa: coupling strength (0 = DEC); shift: None (CP, DEC)
    or the offset s (SH). forced: optional boolean array over candidates imposing the acceptances (Part 0 step check)."""

    def __init__(self, cand, kappa, shift=None, forced=None, Fb=FBOUND, keep_p=False):
        self.c, self.kappa, self.Fb = cand, float(kappa), Fb
        self.shift = None if shift is None else np.asarray(shift, float)
        self.forced = forced
        N0, n = cand["N0"], len(cand["tb"])
        self.death = np.concatenate([cand["life0"], np.full(n, np.inf)])
        self.accepted = np.zeros(n, bool)
        self.evaluated = np.zeros(n, bool)
        self.keep_p = keep_p
        self.p = np.full(n, np.nan) if keep_p else None
        self.max_p = 0.0

    def tick_edits(self, mesh, F, ids, t0, t1):
        """Edits of the tick (t0, t1]: deaths of alive nodes, then the candidates born in the interval that are not
        ephemeral (tb + ell > t1), accepted iff u < f(phi_G*(x [+ s])) / F with phi_G* interpolated on the pre-edit
        mesh. Returns MG5's event array (deaths first, accepted births in candidate order) and a dict of counts."""
        import evolve as ev
        c, N0 = self.c, self.c["N0"]
        dt_ = self.death[ids]
        dead = ids[(dt_ > t0) & (dt_ <= t1)]
        a, b = np.searchsorted(c["tb"], [t0, t1], side="right")
        j = np.arange(a, b)
        eph = c["tb"][j] + c["ell"][j] <= t1
        j = j[~eph]
        if len(j):
            q = np.column_stack([c["x"][j], c["y"][j]])
            if self.shift is not None:
                q = q + self.shift
            if self.kappa == 0.0:
                p = np.full(len(j), 1.0 / self.Fb)
            else:
                phi = mesh.interpolate(F, np.mod(q, (mesh.Lx, mesh.Ly)))[0][0]
                p = np.exp(self.kappa * phi) / self.Fb
            self.max_p = max(self.max_p, float(p.max()))
            acc = (c["u"][j] < p) if self.forced is None else self.forced[j]
            self.evaluated[j] = True
            if self.keep_p:
                self.p[j] = p
            ja = j[acc]
        else:
            ja = j
        self.accepted[ja] = True
        self.death[N0 + ja] = c["tb"][ja] + c["ell"][ja]
        E = np.empty(len(dead) + len(ja), dtype=ev.EV_DTYPE)
        E["t"][:len(dead)] = self.death[dead]
        E["type"][:len(dead)] = 0
        E["id"][:len(dead)] = dead
        E["x"][:len(dead)] = np.nan
        E["y"][:len(dead)] = np.nan
        E["t"][len(dead):] = c["tb"][ja]
        E["type"][len(dead):] = 1
        E["id"][len(dead):] = N0 + ja
        E["x"][len(dead):] = c["x"][ja]
        E["y"][len(dead):] = c["y"][ja]
        return E, dict(n_cand=int(b - a), n_eph=int(eph.sum()), n_eval=int(len(j)), n_acc=int(len(ja)),
                       n_dead=int(len(dead)))


# ------------------------------------------------------------------ quadrature (§4.3)
class Grid:
    """nx x ny cells tiling the torus [0, Lx) x [0, Ly); sample points at the cell midpoints; dA = A / (nx ny)."""

    def __init__(self, Lx, Ly, nx=NX, ny=NY, sub=SUB):
        self.Lx, self.Ly, self.nx, self.ny, self.sub = float(Lx), float(Ly), int(nx), int(ny), int(sub)
        self.dx, self.dy = self.Lx / nx, self.Ly / ny
        self.A = self.Lx * self.Ly
        self.dA = self.A / (nx * ny)
        ix, iy = np.meshgrid(np.arange(nx), np.arange(ny), indexing="ij")
        self.pts = np.column_stack([(ix.ravel() + 0.5) * self.dx, (iy.ravel() + 0.5) * self.dy])
        o = (np.arange(sub) + 0.5) / sub - 0.5
        self.sub_off = np.array([(a * self.dx, b * self.dy) for a in o for b in o])

    def disc_weights(self, c, R):
        """Area weights (cells,) of the disc |x - c| < R (minimum image): dA for cells wholly inside, 0 wholly outside,
        dA x (fraction of the sub x sub sub-lattice inside) for cells the boundary cuts."""
        d = gr.min_image(self.pts - np.asarray(c, float), self.Lx, self.Ly)
        ax, ay = np.abs(d[:, 0]), np.abs(d[:, 1])
        hx, hy = 0.5 * self.dx, 0.5 * self.dy
        far = np.hypot(ax + hx, ay + hy)
        near = np.hypot(np.maximum(ax - hx, 0.0), np.maximum(ay - hy, 0.0))
        w = np.zeros(len(d))
        w[far <= R] = 1.0
        cut = np.flatnonzero((far > R) & (near < R))
        if len(cut):
            sx = d[cut, 0][:, None] + self.sub_off[None, :, 0]
            sy = d[cut, 1][:, None] + self.sub_off[None, :, 1]
            w[cut] = (np.hypot(sx, sy) < R).mean(1)
        return w * self.dA

    def field(self, mesh, F, row=0, shift=None):
        """phi (row) of F at the midpoints (+ shift), barycentric on mesh."""
        q = self.pts if shift is None else np.mod(self.pts + np.asarray(shift, float), (self.Lx, self.Ly))
        out = np.empty(len(q))
        for a in range(0, len(q), 200000):
            out[a:a + 200000] = mesh.interpolate(F[row:row + 1], q[a:a + 200000])[0][0]
        return out


def _dist(pos, c, Lx, Ly):
    return np.hypot(*gr.min_image(np.asarray(pos) - np.asarray(c, float), Lx, Ly).T)


def measure(mesh, F, c, kappa, shift, grid, rho0):
    """Per-check measurements around centre c (§4.2): node counts N(<R), D0, Dc; annular counts and densities
    (coarse and fine edges) against both baselines; radial phi_G by quadrature; M-, M+, M; quasistatic predictions
    P0, Pc for discs and annuli. kappa = 0 gives f = 1 exactly (DEC)."""
    Lx, Ly, A = mesh.Lx, mesh.Ly, grid.A
    N = mesh.n
    rbar = N / A
    r = _dist(mesh.pos, c, Lx, Ly)
    out = dict(N=N, rho_bar=rbar, rho0=rho0)
    NR = np.array([(r < R).sum() for R in RADII])
    out["N_R"] = NR
    out["D0"] = NR - rho0 * np.pi * np.square(RADII)
    out["Dc"] = NR - rbar * np.pi * np.square(RADII)
    g = grid.field(mesh, F, 0)
    out["M_minus"] = float(np.maximum(-g, 0.0).sum() * grid.dA)
    out["M_plus"] = float(np.maximum(g, 0.0).sum() * grid.dA)
    out["M"] = float(g.sum() * grid.dA)
    if kappa == 0:
        f = np.ones_like(g)
    else:
        ge = g if shift is None else grid.field(mesh, F, 0, shift)
        f = np.exp(kappa * ge)
    fmean = rho0 * f.sum() * grid.dA / A                     # (1/A) ∫ rho0 f dA
    W = {R: grid.disc_weights(c, R) for R in sorted(set(RADII) | set(ANN[1:]) | set(FINE[1:]))}
    W[0.0] = np.zeros(len(g))
    area = {R: W[R].sum() for R in W}
    P0 = {R: float(rho0 * ((f - 1.0) * W[R]).sum()) if kappa != 0 else 0.0 for R in W}
    Pc = {R: float(rho0 * (f * W[R]).sum() - area[R] * fmean) if kappa != 0 else 0.0 for R in W}
    gint = {R: float((g * W[R]).sum()) for R in W}
    out["P0"] = np.array([P0[R] for R in RADII])
    out["Pc"] = np.array([Pc[R] for R in RADII])
    out["quad_area"] = np.array([area[R] for R in RADII])
    for name, edges in (("ann", ANN), ("fine", FINE)):
        lo, hi = np.array(edges[:-1]), np.array(edges[1:])
        cnt = np.array([((r >= a) & (r < b)).sum() for a, b in zip(lo, hi)])
        ex_area = np.pi * (hi ** 2 - lo ** 2)
        q_area = np.array([area[b] - area[a] for a, b in zip(lo, hi)])
        out[name + "_count"] = cnt
        out[name + "_dn0"] = cnt / ex_area - rho0
        out[name + "_dnc"] = cnt / ex_area - rbar
        out[name + "_P0"] = np.array([P0[b] - P0[a] for a, b in zip(lo, hi)]) / ex_area
        out[name + "_Pc"] = np.array([Pc[b] - Pc[a] for a, b in zip(lo, hi)]) / ex_area
        out[name + "_phiG"] = np.array([gint[b] - gint[a] for a, b in zip(lo, hi)]) / q_area
    return out


if __name__ == "__main__" and len(sys.argv) > 1 and sys.argv[1] == "--regen-hash":
    N0, nu, t_end, Lx, Ly, seed = sys.argv[2:8]
    print(candidates_hash(gen_candidates(int(N0), float(nu), float(t_end), float(Lx), float(Ly), int(seed))))

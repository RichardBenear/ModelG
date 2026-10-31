"""MG4 (Phase 3): frozen periodic 2D graphs, the graph Laplacian, homogenized calibration of h, spectra.

Registered design: spec/MG4_proposal.md draft 4 (frozen at 773255d), GATES.md MG4.

A graph is a dict with
  pos  (N, 2) node positions on the torus [0, Lx) x [0, Ly)  (measurement only; dynamics never sees them)
  I, J (E,)   edge endpoints, I < J, each undirected edge once
  dx   (E, 2) minimal-image displacement x_J - x_I
  Lx, Ly, kind, plus family-specific metadata.
The unit-rate Laplacian is L1 = A - D_deg (negative semidefinite; rows sum to 0). The dynamics use L_G = L1 / h^2.
"""
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
from scipy.spatial import Delaunay, cKDTree
from scipy.sparse.csgraph import connected_components

SQ3 = np.sqrt(3.0)


# ------------------------------------------------------------------ torus and helpers
def torus(h_t, nx, ny):
    """Rectangular torus commensurate with the triangular lattice of target h_t (a = h_t / sqrt(1.5))."""
    assert ny % 2 == 0, "n_y must be even for a periodic triangular lattice"
    a = h_t / np.sqrt(1.5)
    return dict(a=a, nx=nx, ny=ny, Lx=nx * a, Ly=ny * a * SQ3 / 2)


def min_image(d, Lx, Ly):
    d = np.array(d, dtype=float)
    d[..., 0] -= Lx * np.round(d[..., 0] / Lx)
    d[..., 1] -= Ly * np.round(d[..., 1] / Ly)
    return d


def _finish(pos, I, J, Lx, Ly, kind, **meta):
    I, J = np.asarray(I, dtype=np.int64), np.asarray(J, dtype=np.int64)
    lo, hi = np.minimum(I, J), np.maximum(I, J)
    key = np.unique(lo * len(pos) + hi)
    I, J = key // len(pos), key % len(pos)
    assert np.all(I != J)
    dx = min_image(pos[J] - pos[I], Lx, Ly)
    return dict(pos=pos, I=I, J=J, dx=dx, Lx=Lx, Ly=Ly, kind=kind, N=len(pos), **meta)


def degrees(g):
    return np.bincount(np.concatenate([g["I"], g["J"]]), minlength=g["N"])


# ------------------------------------------------------------------ families
def triangular(tor):
    """T: equilateral triangular lattice on the commensurate torus (degree 6)."""
    a, nx, ny = tor["a"], tor["nx"], tor["ny"]
    j, i = np.divmod(np.arange(nx * ny), nx)
    pos = np.column_stack([(i + 0.5 * (j % 2)) * a, j * a * SQ3 / 2])
    idx = lambda ii, jj: (jj % ny) * nx + (ii % nx)
    n0 = idx(i, j)
    up_l = np.where(j % 2 == 0, idx(i - 1, j + 1), idx(i, j + 1))
    up_r = np.where(j % 2 == 0, idx(i, j + 1), idx(i + 1, j + 1))
    I = np.concatenate([n0, n0, n0])
    J = np.concatenate([idx(i + 1, j), up_l, up_r])
    return _finish(pos, I, J, tor["Lx"], tor["Ly"], "T", a=a)


def poisson_points(n, Lx, Ly, rng):
    return np.column_stack([rng.uniform(0, Lx, n), rng.uniform(0, Ly, n)])


def delaunay(pos, Lx, Ly):
    """D: periodic Delaunay triangulation (3x3 tiling; keep edges touching the central tile)."""
    n = len(pos)
    shifts = np.array([(sx, sy) for sx in (-1, 0, 1) for sy in (-1, 0, 1)], dtype=float) * (Lx, Ly)
    Q = np.concatenate([pos + s for s in shifts])
    centre = int(np.flatnonzero((shifts == 0).all(1))[0])
    tri = Delaunay(Q).simplices
    e = np.concatenate([tri[:, [0, 1]], tri[:, [1, 2]], tri[:, [0, 2]]])
    keep = ((e // n) == centre).any(1)
    e = e[keep] % n
    return _finish(pos, e[:, 0], e[:, 1], Lx, Ly, "D")


def rgg(pos, Lx, Ly, radius):
    """R: periodic random geometric graph, largest connected component only (nodes relabelled)."""
    tree = cKDTree(pos, boxsize=(Lx, Ly))
    pr = tree.query_pairs(radius, output_type="ndarray")
    n = len(pos)
    A = sp.coo_matrix((np.ones(len(pr)), (pr[:, 0], pr[:, 1])), shape=(n, n))
    nc, lab = connected_components(A, directed=False)
    big = np.argmax(np.bincount(lab))
    keep = np.flatnonzero(lab == big)
    new = -np.ones(n, dtype=np.int64)
    new[keep] = np.arange(len(keep))
    m = (lab[pr[:, 0]] == big) & (lab[pr[:, 1]] == big)
    return _finish(pos[keep], new[pr[m, 0]], new[pr[m, 1]], Lx, Ly, "R",
                   radius=radius, n_points=n, kept_fraction=len(keep) / n)


def random_regular(n, k, seed):
    """N3: random k-regular graph (no geometry). Positions are dummies (zeros); never used for measurement."""
    import networkx as nx
    G = nx.random_regular_graph(k, n, seed=int(seed))
    e = np.array(G.edges(), dtype=np.int64)
    I, J = np.minimum(e[:, 0], e[:, 1]), np.maximum(e[:, 0], e[:, 1])
    return dict(pos=np.zeros((n, 2)), I=I, J=J, dx=np.zeros((len(I), 2)), Lx=np.nan, Ly=np.nan,
                kind="N3", N=n, k=k)


# ------------------------------------------------------------------ operators
def laplacian(g):
    """Unit-rate graph Laplacian L1 = A - D_deg (CSR, symmetric)."""
    n = g["N"]
    A = sp.coo_matrix((np.ones(len(g["I"])), (g["I"], g["J"])), shape=(n, n))
    A = (A + A.T).tocsr()
    return (A - sp.diags(np.asarray(A.sum(1)).ravel())).tocsr()


def incidence(g):
    """Signed incidence B (E x N): (B c)_e = c_J - c_I."""
    E, n = len(g["I"]), g["N"]
    rows = np.concatenate([np.arange(E), np.arange(E)])
    cols = np.concatenate([g["J"], g["I"]])
    vals = np.concatenate([np.ones(E), -np.ones(E)])
    return sp.csr_matrix((vals, (rows, cols)), shape=(E, n))


def corrector(g):
    """Homogenized diffusivity tensor of the unit-rate Laplacian (proposal §4).

    For a macroscopic gradient e, D_eff(e, e) = (1/N) min_chi sum_edges (e.dx + chi_J - chi_I)^2; the minimiser solves
    (B^T B) chi = -B^T (dx e) with one node pinned (B^T B = -L1 is singular only along constants on a connected graph).
    Returns dict(D (2x2), eig (e1 >= e2), angle (radians, of e1's eigenvector), h2 = (e1+e2)/2, h2_naive)."""
    B = incidence(g)
    M = (B.T @ B).tocsc()[1:, 1:]
    lu = spla.splu(M)
    cols = []
    for a in range(2):
        rhs = -(B.T @ g["dx"][:, a])[1:]
        chi = np.concatenate([[0.0], lu.solve(rhs)])
        cols.append(g["dx"][:, a] + B @ chi)
    F = np.column_stack(cols)
    D = F.T @ F / g["N"]
    D = 0.5 * (D + D.T)
    w, v = np.linalg.eigh(D)
    naive = float((g["dx"] ** 2).sum() / g["N"] / 2.0)     # mean over nodes of sum_j |dx|^2 / 4 (each edge at 2 nodes)
    return dict(D=D, eig=w[::-1].copy(), angle=float(np.arctan2(v[1, 1], v[0, 1])), h2=float(w.sum() / 2),
                h2_naive=naive)


def cotan(g):
    """MG4c Part 3 diagnostic operator: cotangent (P1 finite-element, lumped-mass) Laplacian on the periodic Delaunay
    triangulation of g's positions (spec/MG4c_proposal.md §4). Faces come from the 3x3 tiling; a face is kept once, iff
    its centroid lies in the central tile (origin offset slightly so that lattice-symmetric centroids never sit on the
    tile border; a torus triangulation has exactly 2N faces, checked in the tests). Returns dict(K (CSR, symmetric, rows sum to 0), A (barycentric node areas),
    L = diag(1/A) K, faces (F, 3) node indices, face_area, edges (E, 2) sorted node pairs of the faces' edges)."""
    pos, Lx, Ly, n = g["pos"], g["Lx"], g["Ly"], g["N"]
    shifts = np.array([(sx, sy) for sx in (-1, 0, 1) for sy in (-1, 0, 1)], dtype=float) * (Lx, Ly)
    Q = np.concatenate([pos + s for s in shifts])
    tri = Delaunay(Q).simplices
    cen = Q[tri].mean(1) - (1.234567e-4 * Lx, 2.345678e-4 * Ly)   # tile origin offset: no centroid sits on its border
    keep = (cen[:, 0] >= 0) & (cen[:, 0] < Lx) & (cen[:, 1] >= 0) & (cen[:, 1] < Ly)
    tri = tri[keep]
    P = Q[tri]                                              # (F, 3, 2) true (unwrapped) geometry
    faces = tri % n
    rows, cols, vals = [], [], []
    A = np.zeros(n)
    e1, e2 = P[:, 1] - P[:, 0], P[:, 2] - P[:, 0]
    area = 0.5 * np.abs(e1[:, 0] * e2[:, 1] - e1[:, 1] * e2[:, 0])
    for k in range(3):                                      # angle at vertex k is opposite edge (k+1, k+2)
        i, j = (k + 1) % 3, (k + 2) % 3
        u, v = P[:, i] - P[:, k], P[:, j] - P[:, k]
        cot = (u * v).sum(1) / np.abs(u[:, 0] * v[:, 1] - u[:, 1] * v[:, 0])
        rows += [faces[:, i], faces[:, j]]
        cols += [faces[:, j], faces[:, i]]
        vals += [0.5 * cot, 0.5 * cot]
        np.add.at(A, faces[:, k], area / 3)
    W = sp.coo_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))), shape=(n, n)).tocsr()
    K = (W - sp.diags(np.asarray(W.sum(1)).ravel())).tocsr()
    ed = np.concatenate([faces[:, [0, 1]], faces[:, [1, 2]], faces[:, [0, 2]]])
    ed = np.sort(ed, 1)
    return dict(K=K, A=A, L=(sp.diags(1.0 / A) @ K).tocsr(), faces=faces, face_area=area, edges=ed, W=W)


def nd_order(g, leaf=64):
    """Geometric nested-dissection node order (fill-reducing order for the sparse LU in the ODE solver; MG010).
    Recursive median bisection of the positions, alternating axes; the nodes of the lower half that touch the upper
    half form the separator and are ordered after both halves. Changes only the order of elimination, not the system."""
    A = (laplacian(g) != 0).tocsr()
    pos, out = g["pos"], []

    def rec(nodes, depth):
        if len(nodes) <= leaf:
            out.append(nodes)
            return
        x = pos[nodes, depth % 2]
        m = np.median(x)
        lo, hi = nodes[x < m], nodes[x >= m]
        in_hi = np.zeros(g["N"], bool)
        in_hi[hi] = True
        sub = A[lo]
        rows = np.repeat(np.arange(len(lo)), np.diff(sub.indptr))
        sep = np.zeros(len(lo), bool)
        sep[rows[in_hi[sub.indices]]] = True
        rec(lo[~sep], depth + 1)
        rec(hi, depth + 1)
        out.append(lo[sep])

    rec(np.arange(g["N"]), 0)
    o = np.concatenate(out)
    assert len(o) == g["N"] and len(np.unique(o)) == g["N"]
    return o


def spectrum(g, h, cutoff, k0=None, vectors=False):
    """Nonnegative eigenvalues nu of -L_G = -(A - D_deg)/h^2 up to `cutoff` (inclusive), ascending.
    Shift-invert Lanczos; k grows until the largest returned eigenvalue exceeds the cutoff."""
    K = (-laplacian(g) / h**2).tocsc()
    k = k0 or 64
    while True:
        k = min(k, g["N"] - 2)
        w, v = spla.eigsh(K, k=k, sigma=-1e-3, which="LM")
        o = np.argsort(w)
        w, v = w[o], v[:, o]
        if w[-1] > cutoff or k >= g["N"] - 2:
            m = w <= cutoff
            return (w[m], v[:, m]) if vectors else w[m]
        k = int(k * 1.5) + 16


def torus_count(k, Lx, Ly):
    """Exact finite-torus count N_c(k) = #{q = (2 pi m/Lx, 2 pi n/Ly): |q| <= k} (includes q = 0)."""
    mx, my = int(np.ceil(k * Lx / (2 * np.pi))) + 1, int(np.ceil(k * Ly / (2 * np.pi))) + 1
    qx = 2 * np.pi * np.arange(-mx, mx + 1) / Lx
    qy = 2 * np.pi * np.arange(-my, my + 1) / Ly
    q2 = qx[:, None] ** 2 + qy[None, :] ** 2
    return int((q2 <= k * k * (1 + 1e-12)).sum())


def count_ratios(nu, Lx, Ly, kstar, n=21, band=0.1):
    """C3(a): N_graph(k^2) / N_c(k) at n equally spaced k in [(1-band) k*, (1+band) k*]; returns (ks, ratios)."""
    ks = np.linspace((1 - band) * kstar, (1 + band) * kstar, n)
    r = np.array([np.sum(nu <= k * k * (1 + 1e-12)) / torus_count(k, Lx, Ly) for k in ks])
    return ks, r

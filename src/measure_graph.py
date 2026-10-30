"""MG4 (Phase 3): measurement on embedded graphs (proposal §5.1–5.3). The embedding is used here only.

- smoother(g): sparse row-normalised Gaussian kernel, sigma = 0.4, cutoff 4 sigma, periodic distances.
- particles(...): candidates, single-linkage merging at 1.5, centroid refinement.
- classify(...): outcome classes with precedence multiplied > decayed > single > unresolved.
- radial_profile, wavelength: bins [0, 0.5) then 0.25 wide to 8; empty bins dropped; lambda from interpolated zeros.
"""
import numpy as np
import scipy.sparse as sp
from scipy.spatial import cKDTree
from scipy.sparse.csgraph import connected_components
import graphs as gr

SIGMA, CUT = 0.4, 1.6
CAND_R, MERGE_R, CENT_R, LINK_R = 1.0, 1.5, 0.6, 1.0
BIN_EDGES = np.concatenate([[0.0], np.arange(0.5, 8.0 + 1e-9, 0.25)])


def _tree(g):
    return cKDTree(np.mod(g["pos"], (g["Lx"], g["Ly"])), boxsize=(g["Lx"], g["Ly"]))


def smoother(g, sigma=SIGMA, cut=CUT):
    tree = _tree(g)
    M = tree.sparse_distance_matrix(tree, cut, output_type="coo_matrix")
    off = M.row != M.col                               # self pairs may or may not be stored; add them once below
    w = np.exp(-M.data[off]**2 / (2 * sigma**2))
    S = sp.coo_matrix((w, (M.row[off], M.col[off])), shape=(g["N"], g["N"])).tocsr()
    S = S + sp.identity(g["N"], format="csr")
    return sp.diags(1.0 / np.asarray(S.sum(1)).ravel()) @ S


def particles(g, s, thr, tree=None):
    """Distinct particles in a smoothed field s. Returns list of dict(node, centre, depth), deepest first."""
    tree = tree or _tree(g)
    below = np.flatnonzero(s < thr)
    cand = []
    for i in below:
        nb = tree.query_ball_point(g["pos"][i] % (g["Lx"], g["Ly"]), CAND_R)
        nb = np.asarray(nb)
        m = s[nb].min()
        if s[i] < m or (s[i] == m and i == nb[s[nb] == m].min()):
            cand.append(i)
    if not cand:
        return []
    cand = np.array(cand)
    ct = cKDTree(np.mod(g["pos"][cand], (g["Lx"], g["Ly"])), boxsize=(g["Lx"], g["Ly"]))
    pr = ct.query_pairs(MERGE_R, output_type="ndarray")
    A = sp.coo_matrix((np.ones(len(pr)), (pr[:, 0], pr[:, 1])), shape=(len(cand),) * 2)
    _, lab = connected_components(A, directed=False)
    out = []
    for c in np.unique(lab):
        mem = cand[lab == c]
        i0 = mem[np.argmin(s[mem])]
        nb = np.asarray(tree.query_ball_point(g["pos"][i0] % (g["Lx"], g["Ly"]), CENT_R))
        w = np.maximum(thr - s[nb], 0.0)
        d = gr.min_image(g["pos"][nb] - g["pos"][i0], g["Lx"], g["Ly"])
        cen = g["pos"][i0] + (w[:, None] * d).sum(0) / w.sum()
        out.append(dict(node=int(i0), centre=np.mod(cen, (g["Lx"], g["Ly"])), depth=float(s[i0])))
    return sorted(out, key=lambda p: p["depth"])


def classify(g, checks):
    """checks: list (in time order) of particle lists at t = 200, 210, ..., 300. Returns (class, track)."""
    if any(len(c) >= 2 for c in checks):
        return "multiplied", None
    if len(checks[-1]) == 0:
        return "decayed", None
    if all(len(c) == 1 for c in checks):
        track = [checks[0][0]["centre"]]
        for c in checks[1:]:
            d = gr.min_image(c[0]["centre"] - track[-1], g["Lx"], g["Ly"])
            if np.hypot(*d) > LINK_R:
                return "unresolved", None
            track.append(c[0]["centre"])
        return "single", np.array(track)
    return "unresolved", None


def radial_profile(g, f, centre, edges=BIN_EDGES):
    """Mean of f over nodes in each ring; returns (bin centres, means (nan where empty), counts)."""
    d = gr.min_image(g["pos"] - centre, g["Lx"], g["Ly"])
    r = np.hypot(d[:, 0], d[:, 1])
    k = np.digitize(r, edges) - 1
    ok = (k >= 0) & (k < len(edges) - 1)
    cnt = np.bincount(k[ok], minlength=len(edges) - 1)
    tot = np.bincount(k[ok], weights=f[ok], minlength=len(edges) - 1)
    mean = np.where(cnt > 0, tot / np.maximum(cnt, 1), np.nan)
    return 0.5 * (edges[1:] + edges[:-1]), mean, cnt


def wavelength(rc, prof, lo=1.0, hi=7.0):
    """2 x mean spacing of zero crossings in lo < r < hi (linear interpolation between non-empty bins); None if < 2."""
    m = np.isfinite(prof)
    r, f = rc[m], prof[m]
    z = []
    for i in range(len(f) - 1):
        if f[i] == 0.0:
            z.append(r[i])
        elif f[i] * f[i + 1] < 0:
            z.append(r[i] - f[i] * (r[i + 1] - r[i]) / (f[i + 1] - f[i]))
    z = np.array([v for v in z if lo < v < hi])
    return None if len(z) < 2 else float(2 * np.mean(np.diff(z)))


def profile_l2(rc, prof, ref, rmax=6.0):
    """C2(a): area-weighted (r dr) relative L2 difference over r < rmax, non-empty bins of both."""
    m = (rc < rmax) & np.isfinite(prof) & np.isfinite(ref)
    w = rc[m]
    return float(np.sqrt((w * (prof[m] - ref[m]) ** 2).sum() / (w * ref[m] ** 2).sum()))

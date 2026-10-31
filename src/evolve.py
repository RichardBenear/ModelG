"""MG5 (Phase 4): a particle on an evolving periodic Delaunay graph.

Registered design: spec/MG5_proposal.md draft 3 (frozen at e2a25e3), GATES.md MG5.

Contents
- Mesh: periodic Delaunay (3x3 tiling, as graphs.delaunay) that also keeps the triangulation for interpolation.
- Ros2: fixed-step, linearly implicit Rosenbrock method ROS2 (Verwer et al. 1999, gamma = 1 + 1/sqrt 2) for the graph ODE
  of graph_ode.make_system; one Jacobian + LU per tick (proposal §4).
- gen_events / events_hash: the continuous-time birth-death event list (proposal §3.1) and its canonical SHA-256.
- apply_tick: the tick-synchronous edit of proposal §3.2 (barycentric births on the pre-edit mesh, ephemeral nodes
  excluded, effective deaths removed, deterministic row order) with jump accounting (§7).
- audit_births / exact_births: the independent-interpolator audit and its exact-rational resolution (§7).
- detect_full / sparse_record / detect_sparse: detection on the full field, the per-check sparse record, and exact
  re-detection from the sparse record alone (§7, "Saved data").
- classify_run: per-run class and mobility (§7).
"""
import hashlib
import os, sys
from fractions import Fraction
import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import splu
from scipy.spatial import Delaunay, cKDTree
sys.path.insert(0, os.path.dirname(__file__))
import graphs as gr
import graph_ode as go
import measure_graph as mg

GATE = 5
GAMMA = 1.0 + 1.0 / np.sqrt(2.0)
SPARSE_CUT = -1.0            # sparse records keep every node with smoothed phi_X < -1.0
JUMP_TOL = 1e-10             # jump-accounting identity tolerance (relative to 1 + |I| + |R| + |B|)
AUDIT_TOL = 1e-9             # independent-interpolator audit flag level
EXACT_TOL = 1e-12            # exact-recheck agreement (relative to 1 + |value|)
AUDIT_N = 20                 # births audited per tick
AUDIT_WIN = 1.0              # local window (radius) of the independent interpolator
SPLIT_R, STEP_R, MOBILE_R = 4.0, 2.0, 1.5


# ------------------------------------------------------------------ mesh
class Mesh:
    """Periodic Delaunay triangulation of pos (N, 2) in [0, Lx) x [0, Ly). Builds the same graph as graphs.delaunay
    (3x3 tiling, edges touching the central tile) and keeps the tiled triangulation for point location."""

    def __init__(self, pos, Lx, Ly):
        self.pos = np.mod(np.asarray(pos, float), (Lx, Ly))
        self.Lx, self.Ly, self.n = Lx, Ly, len(pos)
        shifts = np.array([(sx, sy) for sx in (-1, 0, 1) for sy in (-1, 0, 1)], dtype=float) * (Lx, Ly)
        Q = np.concatenate([self.pos + s for s in shifts])
        self.tri = Delaunay(Q)
        n = self.n
        centre = int(np.flatnonzero((shifts == 0).all(1))[0])
        simp = self.tri.simplices
        e = np.concatenate([simp[:, [0, 1]], simp[:, [1, 2]], simp[:, [0, 2]]])
        keep = ((e // n) == centre).any(1)
        e = e[keep] % n
        self.g = gr._finish(self.pos, e[:, 0], e[:, 1], Lx, Ly, "D")
        self.Q = Q

    def locate(self, q):
        """Triangle (tiled vertex indices) and barycentric weights for points q in [0, Lx) x [0, Ly)."""
        q = np.mod(np.atleast_2d(np.asarray(q, float)), (self.Lx, self.Ly))
        s = self.tri.find_simplex(q)
        if np.any(s < 0):
            raise RuntimeError("point outside the tiled triangulation")
        T = self.tri.transform[s]
        b = np.einsum("ijk,ik->ij", T[:, :2], q - T[:, 2])
        w = np.column_stack([b, 1.0 - b.sum(1)])
        return self.tri.simplices[s], w

    def interpolate(self, F, q):
        """Barycentric interpolation of F (3, N) at q (k, 2). Returns (3, k), the vertex rows (k, 3) and weights."""
        v, w = self.locate(q)
        rows = v % self.n
        return np.einsum("kj,skj->sk", w, F[:, rows]), rows, w


# ------------------------------------------------------------------ ROS2
class Ros2:
    """ROS2 on one fixed graph. factor(y) evaluates the Jacobian at y and factorises (I - gamma dt J) once;
    advance(y, n) takes n steps of dt with that factorisation (a W-method: order 2 for any approximate Jacobian)."""

    def __init__(self, g, h, P, dt):
        self.rhs, self.jac, self.P = go.make_system(g, h, P)
        self.N, self.dt = g["N"], float(dt)
        o = gr.nd_order(g)
        N = self.N
        self.perm = np.column_stack([o, N + o, 2 * N + o]).ravel()
        self.inv = np.argsort(self.perm)
        self.lu = None
        self.n_factor = 0

    def factor(self, y):
        J = self.jac(0.0, y)
        M = (sp.identity(3 * self.N, format="csc") - GAMMA * self.dt * J).tocsc()
        M = M[self.perm][:, self.perm].tocsc()
        self.lu = splu(M, permc_spec="NATURAL")
        self.n_factor += 1

    def _solve(self, v):
        return self.lu.solve(v[self.perm])[self.inv]

    def advance(self, y, n):
        dt = self.dt
        for _ in range(int(n)):
            k1 = self._solve(self.rhs(0.0, y))
            k2 = self._solve(self.rhs(0.0, y + dt * k1) - 2.0 * k1)
            y = y + 1.5 * dt * k1 + 0.5 * dt * k2
        return y


def n_steps(span, dt):
    """Integer number of steps of dt in span; the registered Delta/dt combinations are exact multiples."""
    n = int(round(span / dt))
    if abs(n * dt - span) > 1e-9 * max(1.0, span):
        raise ValueError(f"span {span} is not a multiple of dt {dt}")
    return n


# ------------------------------------------------------------------ events
EV_DTYPE = [("t", "<f8"), ("type", "i1"), ("id", "<i8"), ("x", "<f8"), ("y", "<f8")]   # type 0 death, 1 birth


def gen_events(ids0, N0, nu, t_end, Lx, Ly, seed):
    """Exact continuous-time spatial birth-death list (proposal §3.1). Deaths at rate nu per alive node, births at rate
    nu * N0 (= nu rho A) at uniform positions. A death picks uniformly among the nodes alive at its time.
    New ids are N0 + (number of births so far). Returns a structured array in time order."""
    rng = np.random.default_rng(seed)
    ids0 = np.asarray(ids0, dtype=np.int64)
    cap = len(ids0) + int(nu * N0 * t_end * 1.5) + 1000
    alive = np.empty(cap, dtype=np.int64)
    alive[:len(ids0)] = ids0
    na = len(ids0)
    next_id = int(N0) if len(ids0) == 0 else max(int(N0), int(ids0.max()) + 1)
    out_t, out_type, out_id, out_x, out_y = [], [], [], [], []
    t = 0.0
    birth_rate = nu * N0
    while True:
        R = nu * na + birth_rate
        if R <= 0:
            break
        t += rng.exponential(1.0 / R)
        if t > t_end:
            break
        if rng.random() * R < nu * na:
            k = int(rng.integers(na))
            out_t.append(t); out_type.append(0); out_id.append(int(alive[k])); out_x.append(np.nan); out_y.append(np.nan)
            na -= 1
            alive[k] = alive[na]
        else:
            x, y = rng.random() * Lx, rng.random() * Ly
            if na >= len(alive):
                alive = np.concatenate([alive, np.empty(len(alive), dtype=np.int64)])
            alive[na] = next_id
            na += 1
            out_t.append(t); out_type.append(1); out_id.append(next_id); out_x.append(x); out_y.append(y)
            next_id += 1
    ev = np.empty(len(out_t), dtype=EV_DTYPE)
    ev["t"], ev["type"], ev["id"], ev["x"], ev["y"] = out_t, out_type, out_id, out_x, out_y
    return ev


def events_hash(ev):
    """SHA-256 of the canonical serialization: the five columns, little-endian, in the order t, type, id, x, y."""
    h = hashlib.sha256()
    for name in ("t", "type", "id", "x", "y"):
        h.update(np.ascontiguousarray(ev[name]).astype(np.dtype(ev.dtype[name]).newbyteorder("<")).tobytes())
    return h.hexdigest()


def regen_hash_fresh_process(args):
    """Regenerate the event list in a fresh Python process (a subprocess, so it also works inside a pool worker) and
    return its hash (proposal §3.1). args = (ids0, N0, nu, t_end, Lx, Ly, seed) with ids0 = arange(len(ids0))."""
    import subprocess
    ids0, N0, nu, t_end, Lx, Ly, seed = args
    ids0 = np.asarray(ids0)
    if not np.array_equal(ids0, np.arange(len(ids0))):
        raise ValueError("fresh-process regeneration supports ids0 = arange(n) only")
    cmd = [sys.executable, os.path.abspath(__file__), "--regen-hash", str(len(ids0)), str(int(N0)), repr(float(nu)),
           repr(float(t_end)), repr(float(Lx)), repr(float(Ly)), str(int(seed))]
    return subprocess.check_output(cmd).decode().strip().splitlines()[-1]


def tick_slices(ev_t, delta, t_end):
    """Index ranges of the events in (t_(k-1), t_k] for k = 1 .. K (K = t_end / delta)."""
    K = n_steps(t_end, delta)
    edges = np.arange(K + 1) * delta
    idx = np.searchsorted(ev_t, edges, side="right")
    return [(int(idx[k - 1]), int(idx[k])) for k in range(1, K + 1)]


# ------------------------------------------------------------------ one tick of edits
def apply_tick(mesh, F, ids, ev):
    """Proposal §3.2 steps 3-6. mesh: pre-edit Mesh (positions mesh.pos, rows = ids order); F (3, N) boundary field;
    ids (N,) persistent ids; ev: this tick's events (time order). Returns (pos, ids, F, acct, births) with births =
    dict(q, rows, w, val) for the effective births (used by the audit)."""
    N = len(ids)
    deaths = ev["id"][ev["type"] == 0]
    bmask = ev["type"] == 1
    born_ids = ev["id"][bmask]
    eph = np.isin(born_ids, deaths)                          # born and deleted within this tick
    eff_born = ~eph
    q = np.column_stack([ev["x"][bmask][eff_born], ev["y"][bmask][eff_born]])
    new_ids = born_ids[eff_born]
    eff_deaths = deaths[~np.isin(deaths, born_ids)]
    row_of = {int(i): r for r, i in enumerate(ids)} if len(eff_deaths) else {}
    dead_rows = np.array([row_of[int(i)] for i in eff_deaths], dtype=np.int64)
    if len(q):
        val, vrows, w = mesh.interpolate(F, q)
    else:
        val, vrows, w = np.zeros((3, 0)), np.zeros((0, 3), dtype=np.int64), np.zeros((0, 3))
    I_before = F.sum(1)
    R = F[:, dead_rows].sum(1) if len(dead_rows) else np.zeros(3)
    B = val.sum(1)
    keep = np.ones(N, bool)
    keep[dead_rows] = False
    F_new = np.concatenate([F[:, keep], val], axis=1)
    pos_new = np.concatenate([mesh.pos[keep], np.mod(q, (mesh.Lx, mesh.Ly))])
    ids_new = np.concatenate([ids[keep], new_ids])
    I_after = F_new.sum(1)
    resid = I_after - (I_before - R + B)
    scale = 1.0 + np.abs(I_before) + np.abs(R) + np.abs(B)
    acct = dict(n_prop_births=int(bmask.sum()), n_prop_deaths=int(len(deaths)), n_eff_births=int(len(new_ids)),
                n_eff_deaths=int(len(dead_rows)), n_ephemeral=int(eph.sum()),
                I_before=I_before.tolist(), R=R.tolist(), B=B.tolist(), I_after=I_after.tolist(),
                jump_resid=float(np.max(np.abs(resid) / scale)), jump_ok=bool(np.all(np.abs(resid) <= JUMP_TOL * scale)))
    births = dict(q=q, rows=vrows, w=w, val=val)
    return pos_new, ids_new, F_new, acct, births


# ------------------------------------------------------------------ interpolation audit and exact resolution
def audit_births(mesh, F, births, rng=None, n=AUDIT_N):
    """Independent interpolator: scipy LinearNDInterpolator on the periodically padded points within AUDIT_WIN of each
    audited birth. Returns (max |difference|, audited indices)."""
    from scipy.interpolate import LinearNDInterpolator
    k = len(births["q"])
    if k == 0:
        return 0.0, np.zeros(0, dtype=np.int64)
    sel = np.arange(k) if k <= n else np.sort((rng or np.random.default_rng(0)).choice(k, n, replace=False))
    Lx, Ly = mesh.Lx, mesh.Ly
    tree = cKDTree(mesh.pos, boxsize=(Lx, Ly))
    worst = 0.0
    for i in sel:
        q = births["q"][i]
        nb = np.asarray(tree.query_ball_point(np.mod(q, (Lx, Ly)), AUDIT_WIN))
        d = gr.min_image(mesh.pos[nb] - q, Lx, Ly)
        ip = LinearNDInterpolator(d, F[:, nb].T)
        v = ip(np.zeros((1, 2)))[0]
        if not np.all(np.isfinite(v)):
            return np.inf, sel
        worst = max(worst, float(np.max(np.abs(v - births["val"][:, i]))))
    return worst, sel


def _orient(ax, ay, bx, by, cx, cy):
    return (bx - ax) * (cy - ay) - (by - ay) * (cx - ax)


def exact_births(mesh, F, births):
    """Exact rational recheck of every effective birth of a tick (proposal §7): exhaustive search over all triangles of
    the pre-edit tiled mesh whose bounding box contains the point (a necessary condition), exact orientation tests and
    exact barycentric weights on the float64 positions. Returns dict(ok, max_rel, n_multi)."""
    simp, Q, n = mesh.tri.simplices, mesh.Q, mesh.n
    P = Q[simp]                                              # (T, 3, 2)
    lo, hi = P.min(1), P.max(1)
    max_rel, n_multi, ok = 0.0, 0, True
    for i, q in enumerate(births["q"]):
        qx, qy = float(q[0]), float(q[1])
        cand = np.flatnonzero((lo[:, 0] <= qx) & (hi[:, 0] >= qx) & (lo[:, 1] <= qy) & (hi[:, 1] >= qy))
        Qx, Qy = Fraction(qx), Fraction(qy)
        vals = []
        for c in cand:
            (ax, ay), (bx, by), (cx, cy) = [(Fraction(float(p[0])), Fraction(float(p[1]))) for p in P[c]]
            D = _orient(ax, ay, bx, by, cx, cy)
            w0 = _orient(Qx, Qy, bx, by, cx, cy) / D
            w1 = _orient(ax, ay, Qx, Qy, cx, cy) / D
            w2 = 1 - w0 - w1
            if w0 < 0 or w1 < 0 or w2 < 0:
                continue
            rows = simp[c] % n
            vals.append([w0 * Fraction(float(F[s, rows[0]])) + w1 * Fraction(float(F[s, rows[1]]))
                         + w2 * Fraction(float(F[s, rows[2]])) for s in range(3)])
        if not vals:
            return dict(ok=False, max_rel=np.inf, n_multi=n_multi, reason=f"birth {i}: no containing triangle")
        if len(vals) > 1:
            n_multi += 1
            for v in vals[1:]:
                for s in range(3):
                    if abs(float(v[s] - vals[0][s])) > EXACT_TOL * (1 + abs(float(vals[0][s]))):
                        return dict(ok=False, max_rel=np.inf, n_multi=n_multi, reason=f"birth {i}: containing triangles disagree")
        for s in range(3):
            ex = float(vals[0][s])
            rel = abs(births["val"][s, i] - ex) / (1 + abs(ex))
            max_rel = max(max_rel, rel)
    ok = max_rel <= EXACT_TOL
    return dict(ok=bool(ok), max_rel=float(max_rel), n_multi=n_multi)


# ------------------------------------------------------------------ detection and sparse records
def smoothed_X(g, F):
    return mg.smoother(g) @ F[1]


def sparse_record(pos, ids, s, Lx, Ly):
    """Per-check record (proposal §7): every node with smoothed phi_X < SPARSE_CUT, with its row index (increasing),
    position, id and smoothed value; the node count; the global minimum and its position."""
    rows = np.flatnonzero(s < SPARSE_CUT)
    j = int(np.argmin(s))
    return dict(rows=rows.astype(np.int64), pos=np.asarray(pos)[rows].astype(np.float64),
                ids=np.asarray(ids)[rows].astype(np.int64), s=np.asarray(s)[rows].astype(np.float64),
                N=int(len(s)), smin=float(s[j]), smin_pos=np.asarray(pos)[j].astype(np.float64), Lx=Lx, Ly=Ly)


def detect_full(pos, s, thr, Lx, Ly):
    g = dict(pos=np.asarray(pos), Lx=Lx, Ly=Ly, N=len(s))
    return mg.particles(g, np.asarray(s), thr)


def detect_sparse(rec, thr):
    """Exact re-detection from a sparse record at any thr <= SPARSE_CUT. Stored nodes are in increasing row order, so
    MG4's lowest-index tie break is preserved."""
    if thr > SPARSE_CUT:
        raise ValueError("sparse records support thresholds <= SPARSE_CUT only")
    if len(rec["s"]) == 0:
        return []
    g = dict(pos=rec["pos"], Lx=rec["Lx"], Ly=rec["Ly"], N=len(rec["s"]))
    out = mg.particles(g, rec["s"], thr)
    for p in out:
        p["row"] = int(rec["rows"][p["node"]])
    return out


# ------------------------------------------------------------------ per-run classification (proposal §7)
def classify_run(counts_centres, c0, Lx, Ly):
    """counts_centres: list over checks (t = 10, 20, ..., t_end) of particle lists (each with 'centre'), in time order;
    c0: starting centre (t = 0 detection). Returns dict(cls, sub, mobility, max_d, track)."""
    lists = counts_centres
    n = [len(c) for c in lists]
    pairs = list(zip(n[:-1], n[1:]))
    if any(a >= 2 and b >= 2 for a, b in pairs):
        k = next(i for i, (a, b) in enumerate(pairs) if a >= 2 and b >= 2)
        # split vs remote: the first check with >= 2 particles; the "second" particle is the one farther from the
        # original's last known single position
        first = next(i for i, m in enumerate(n) if m >= 2)
        ref = np.asarray(c0, float)
        for i in range(first - 1, -1, -1):
            if n[i] == 1:
                ref = lists[i][0]["centre"]
                break
        d = [np.hypot(*gr.min_image(np.asarray(p["centre"]) - ref, Lx, Ly)) for p in lists[first]]
        d.sort()
        sub = "split" if d[1] <= SPLIT_R else "remote"
        return dict(cls="multiplied", sub=sub, first_multi_check=first, first_qualifying=k)
    if any(a == 0 and b == 0 for a, b in pairs):
        return dict(cls="dissolved", sub=None)
    if all(m == 1 for m in n):
        track = [np.asarray(c0, float)]
        unwrapped = [np.asarray(c0, float)]
        for c in lists:
            step = gr.min_image(np.asarray(c[0]["centre"]) - track[-1], Lx, Ly)
            if np.hypot(*step) > STEP_R:
                return dict(cls="unresolved", sub="step")
            track.append(np.asarray(c[0]["centre"], float))
            unwrapped.append(unwrapped[-1] + step)
        U = np.array(unwrapped)
        d = np.hypot(*(U - U[0]).T)
        return dict(cls="single", sub=None, mobility="mobile" if d.max() >= MOBILE_R else "pinned",
                    max_d=float(d.max()), track=U.tolist())
    return dict(cls="unresolved", sub=None)


def msd_dp(track, dt_check=10.0, lags=range(1, 11)):
    """Descriptive D_p: time-averaged MSD over lags 10-100 fitted as MSD = 4 D_p tau (least squares through 0)."""
    U = np.asarray(track, float)
    tau, msd = [], []
    for L in lags:
        if L >= len(U):
            break
        d = U[L:] - U[:-L]
        tau.append(L * dt_check)
        msd.append(float(np.mean((d ** 2).sum(1))))
    if not tau:
        return None
    tau, msd = np.array(tau), np.array(msd)
    return float((tau * msd).sum() / (4 * (tau ** 2).sum()))


# ------------------------------------------------------------------ flag rechecks (proposal §7)
HP_DIGITS, HP_TOL = 40, 1e-25


def hp_balance(g, h, P, y):
    """MG4c-style recheck of the node-sum balances of state y on graph g (unweighted operator) in mpmath at HP_DIGITS
    (exact fractions if mpmath is missing). The operator is rebuilt from its edges, so rows sum to zero exactly.
    Returns dict(r1, r2, M1, M2, resolved)."""
    try:
        import mpmath
        mpmath.mp.dps = HP_DIGITS
        Fm, Z = mpmath.mpf, mpmath.mpf(0)
    except ImportError:
        Fm, Z = Fraction, Fraction(0)
    N = g["N"]
    I = np.concatenate([g["I"], g["J"]]); J = np.concatenate([g["J"], g["I"]])
    inv_h2 = 1 / (Fm(float(h)) * Fm(float(h)))
    b, gg, p, q, dxc, dyc = (Fm(float(P[k])) for k in ("b", "g", "p", "q", "dx", "dy"))
    G0, X0, Y0 = (Fm(float(v)) for v in go.steady_state(P))
    f = [[Fm(float(v)) for v in y[k * N:(k + 1) * N]] for k in range(3)]

    def lap(u):
        tot = Z
        for a, c in zip(I, J):
            tot += inv_h2 * (u[c] - u[a])
        return tot

    sG, sX = sum(f[0], Z), sum(f[1], Z)
    r4 = sum(((X0 + f[1][i]) ** 2 * (Y0 + f[2][i]) - X0 ** 2 * Y0 for i in range(N)), Z)
    dG = lap(f[0]) - q * sG + gg * sX
    dX = dxc * lap(f[1]) + p * sG - (1 + b) * sX + r4
    dY = dyc * lap(f[2]) + b * sX - r4
    r1 = dG - (-q * sG + gg * sX)
    r2 = (dG + dX + dY) - (-(q - p) * sG - (1 - gg) * sX)
    M1 = abs(dG) + q * abs(sG) + gg * abs(sX)
    M2 = abs(dG) + abs(dX) + abs(dY) + abs(q - p) * abs(sG) + (1 - gg) * abs(sX)
    out = dict(r1=float(r1), r2=float(r2), M1=float(M1), M2=float(M2))
    out["resolved"] = bool(abs(out["r1"]) <= HP_TOL * (1 + out["M1"]) and abs(out["r2"]) <= HP_TOL * (1 + out["M2"]))
    return out


def exact_jump(F_before, dead_rows, births_val, F_after):
    """Exact recheck of a jump-accounting flag: the four node sums per species with math.fsum (correctly rounded).
    Resolved iff the identity holds to 1e-13 relative."""
    import math
    out = []
    for s in range(3):
        Ib = math.fsum(F_before[s]); R = math.fsum(F_before[s, dead_rows]); B = math.fsum(births_val[s])
        Ia = math.fsum(F_after[s])
        out.append(abs(Ia - (Ib - R + B)) / (1 + abs(Ib) + abs(R) + abs(B)))
    return dict(rel=[float(v) for v in out], resolved=bool(max(out) <= 1e-13))


if __name__ == "__main__" and len(sys.argv) == 9 and sys.argv[1] == "--regen-hash":
    n, N0_, nu_, te_, Lx_, Ly_, seed_ = sys.argv[2:]
    print(events_hash(gen_events(np.arange(int(n)), int(N0_), float(nu_), float(te_), float(Lx_), float(Ly_),
                                 int(seed_))))

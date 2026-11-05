"""MG8 (Phase 6): scheduler sensitivity of the 1D lattice approximation.

Registered design: spec/MG8_proposal.md draft 4 (frozen at 3ab4ec8), GATES.md MG8.

Schedulers (proposal §3.1). Every node update takes ONE snapshot (g, x, y) of the node's counts, draws every exit of every
species from it in a fixed sequence (G: R2, R1r, hopL, hopR; X: R2r, R3, R5, hopL, hopR; Y: R4, hopL, hopR; then the R1
source), and applies all of the node's changes at once.
  J = 0  synchronous: snapshots at start of step, changes applied after the last node (== lattice.tau_leap_euler_1d).
  R = 1  rightward sweep from a seam s ~ U{0..N-1} drawn each step: s, s+1, ..., s-1 (mod N).
  L = 2  mirror of R: seam s drawn as for R, s' = N-1-s, order s', s'-1, ...; the two hop draws swap direction.
  P = 3  sequential in a uniform random permutation drawn each step.
Per step the random draws are consumed as: the order draw (seam or permutation; none for J), then node updates in
processing order.  Deterministic versions (Omega -> infinity) apply the same orders to expected concentrations
(exit = c * prob, same rescaling rule, R4 per-Y probability x^2 tau); in them only the order draw is random.
"""
import os, sys
import numpy as np
import numba as nb
sys.path.insert(0, os.path.dirname(__file__))
from lattice import _exits, P17, steady_state
from measure import smooth_periodic, fourier_shift

GATE = 8
TAU0 = 6.667e-4
H = 0.2
L_DOMAIN = 100.0
J, R, L, P = 0, 1, 2, 3
NAMES = {J: "J", R: "R", L: "L", P: "P"}


# ------------------------------------------------------------------ order draws
@nb.njit(cache=True)
def _order(sched, N, order):
    """Fill `order` with this step's node sequence, consuming the order draw. Returns the mirror flag."""
    if sched == R or sched == L:
        s = np.random.randint(0, N)
        if sched == R:
            for k in range(N):
                order[k] = (s + k) % N
            return 0
        sp_ = N - 1 - s
        for k in range(N):
            order[k] = (sp_ - k) % N
        return 1
    if sched == P:
        perm = np.random.permutation(N)
        for k in range(N):
            order[k] = perm[k]
        return 0
    for k in range(N):
        order[k] = k
    return 0


# ------------------------------------------------------------------ stochastic kernel
@nb.njit(cache=True)
def run_stoch(nG, nX, nY, sched, Om, a, b, g, p, q, dx, dy, h, tau, nsteps, rec_every, rng_seed):
    """Advance integer counts in place for nsteps under scheduler `sched`. Records counts every rec_every steps.
    rng_seed < 0 continues the generator state of the previous call (segmented runs).
    Returns recG, recX, recY, rec_step, clips[3] (rescaled node-updates per species)."""
    if rng_seed >= 0:
        np.random.seed(rng_seed)
    N = nG.size
    hg, hx, hy = tau / h**2, dx * tau / h**2, dy * tau / h**2
    nrec = nsteps // rec_every
    recG = np.zeros((nrec, N), dtype=np.int64)
    recX = np.zeros((nrec, N), dtype=np.int64)
    recY = np.zeros((nrec, N), dtype=np.int64)
    recs = np.zeros(nrec, dtype=np.int64)
    inv_Om2 = 1.0 / (Om * Om)
    lamA = Om * a * tau
    PG = np.array([p * tau, (q - p) * tau, hg, hg])
    PX = np.array([g * tau, b * tau, (1.0 - g) * tau, hx, hx])
    PY = np.zeros(3)
    oG = np.zeros(4, dtype=np.int64); oX = np.zeros(5, dtype=np.int64); oY = np.zeros(3, dtype=np.int64)
    dG = np.zeros(N, dtype=np.int64); dX = np.zeros(N, dtype=np.int64); dY = np.zeros(N, dtype=np.int64)
    order = np.zeros(N, dtype=np.int64)
    clips = np.zeros(3, dtype=np.int64)
    k = 0
    for step in range(nsteps):
        mirror = _order(sched, N, order)
        if sched == J:
            for i in range(N):
                dG[i] = 0; dX[i] = 0; dY[i] = 0
        for kk in range(N):
            i = order[kk]
            gi, xi, yi = nG[i], nX[i], nY[i]
            if mirror == 0:
                il, ir = (i - 1) % N, (i + 1) % N
            else:
                il, ir = (i + 1) % N, (i - 1) % N
            clips[0] += _exits(gi, PG, oG)
            clips[1] += _exits(xi, PX, oX)
            PY[0] = xi * (xi - 1) * inv_Om2 * tau
            PY[1] = hy; PY[2] = hy
            clips[2] += _exits(yi, PY, oY)
            src = np.random.poisson(lamA)
            cg = -(oG[0] + oG[1] + oG[2] + oG[3]) + oX[0] + src
            cx = oG[0] - (oX[0] + oX[1] + oX[2] + oX[3] + oX[4]) + oY[0]
            cy = oX[1] - (oY[0] + oY[1] + oY[2])
            if sched == J:
                dG[i] += cg; dX[i] += cx; dY[i] += cy
                dG[il] += oG[2]; dG[ir] += oG[3]
                dX[il] += oX[3]; dX[ir] += oX[4]
                dY[il] += oY[1]; dY[ir] += oY[2]
            else:
                nG[i] += cg; nX[i] += cx; nY[i] += cy
                nG[il] += oG[2]; nG[ir] += oG[3]
                nX[il] += oX[3]; nX[ir] += oX[4]
                nY[il] += oY[1]; nY[ir] += oY[2]
        if sched == J:
            for i in range(N):
                nG[i] += dG[i]; nX[i] += dX[i]; nY[i] += dY[i]
        if (step + 1) % rec_every == 0 and k < nrec:
            recG[k, :] = nG; recX[k, :] = nX; recY[k, :] = nY; recs[k] = step + 1
            k += 1
    return recG, recX, recY, recs, clips


# ------------------------------------------------------------------ deterministic kernel
@nb.njit(cache=True)
def _exits_mean(c, probs, out):
    s = 0.0
    for j in range(probs.size):
        s += probs[j]
    scale = 1.0
    clip = 0
    if s > 1.0:
        clip = 1
        scale = 1.0 / s
    for j in range(probs.size):
        out[j] = c * probs[j] * scale
    return clip


@nb.njit(cache=True)
def run_det(G, X, Y, sched, a, b, g, p, q, dx, dy, h, tau, nsteps, rng_seed):
    """Deterministic (expected-concentration) map under scheduler `sched`, nsteps steps, in place.
    Only the order draw is random. rng_seed < 0 continues the previous generator state. Returns clips[3]."""
    if rng_seed >= 0:
        np.random.seed(rng_seed)
    N = G.size
    hg, hx, hy = tau / h**2, dx * tau / h**2, dy * tau / h**2
    PG = np.array([p * tau, (q - p) * tau, hg, hg])
    PX = np.array([g * tau, b * tau, (1.0 - g) * tau, hx, hx])
    PY = np.zeros(3)
    oG = np.zeros(4); oX = np.zeros(5); oY = np.zeros(3)
    dG = np.zeros(N); dX = np.zeros(N); dY = np.zeros(N)
    order = np.zeros(N, dtype=np.int64)
    clips = np.zeros(3, dtype=np.int64)
    srcA = a * tau
    for step in range(nsteps):
        mirror = _order(sched, N, order)
        if sched == J:
            for i in range(N):
                dG[i] = 0.0; dX[i] = 0.0; dY[i] = 0.0
        for kk in range(N):
            i = order[kk]
            gi, xi, yi = G[i], X[i], Y[i]
            if mirror == 0:
                il, ir = (i - 1) % N, (i + 1) % N
            else:
                il, ir = (i + 1) % N, (i - 1) % N
            clips[0] += _exits_mean(gi, PG, oG)
            clips[1] += _exits_mean(xi, PX, oX)
            PY[0] = xi * xi * tau
            PY[1] = hy; PY[2] = hy
            clips[2] += _exits_mean(yi, PY, oY)
            cg = -(oG[0] + oG[1] + oG[2] + oG[3]) + oX[0] + srcA
            cx = oG[0] - (oX[0] + oX[1] + oX[2] + oX[3] + oX[4]) + oY[0]
            cy = oX[1] - (oY[0] + oY[1] + oY[2])
            if sched == J:
                dG[i] += cg; dX[i] += cx; dY[i] += cy
                dG[il] += oG[2]; dG[ir] += oG[3]
                dX[il] += oX[3]; dX[ir] += oX[4]
                dY[il] += oY[1]; dY[ir] += oY[2]
            else:
                G[i] += cg; X[i] += cx; Y[i] += cy
                G[il] += oG[2]; G[ir] += oG[3]
                X[il] += oX[3]; X[ir] += oX[4]
                Y[il] += oY[1]; Y[ir] += oY[2]
        if sched == J:
            for i in range(N):
                G[i] += dG[i]; X[i] += dX[i]; Y[i] += dY[i]
    return clips


# ------------------------------------------------------------------ diffusion-only conditional maps (Part 0)
@nb.njit(cache=True)
def diffusion_map(c, order, r, sync):
    """One step of the diffusion-only deterministic map for one species (h = 1, D = 1, r = tau) under a given order.
    sync = True gives J. Returns the new vector (input not modified)."""
    N = c.size
    u = c.copy()
    d = np.zeros(N)
    for kk in range(N):
        i = order[kk]
        ci = u[i]
        il, ir = (i - 1) % N, (i + 1) % N
        if sync:
            d[i] -= 2 * r * ci; d[il] += r * ci; d[ir] += r * ci
        else:
            u[i] -= 2 * r * ci; u[il] += r * ci; u[ir] += r * ci
    if sync:
        for i in range(N):
            u[i] += d[i]
    return u


def map_matrix(order, r, sync=False):
    """T_o built by applying the production-style diffusion map to basis vectors (column j = M_o e_j)."""
    N = len(order)
    T = np.zeros((N, N))
    o = np.asarray(order, dtype=np.int64)
    for j in range(N):
        e = np.zeros(N); e[j] = 1.0
        T[:, j] = diffusion_map(e, o, r, sync)
    return T


def chain_matrix(order, r):
    """Independent construction of T_o from per-token chain rules: a token at node i is processed at its node's
    position in the order; it stays (1-2r) or hops; a token landing on a node processed LATER is processed again."""
    N = len(order)
    pos = np.empty(N, dtype=int)
    for k, i in enumerate(order):
        pos[i] = k
    T = np.zeros((N, N))
    for j in range(N):
        # distribution over (node, processed-position) frontier
        stack = [(j, 1.0, -1)]          # node, mass, position of the update that delivered it (-1: start of step)
        while stack:
            i, m, came = stack.pop()
            if pos[i] > came:           # will still be processed this step
                T[i, j] += m * (1 - 2 * r)
                for nb_ in ((i - 1) % N, (i + 1) % N):
                    stack.append((nb_, m * r, pos[i]))
            else:
                T[i, j] += m
    return T


def displacement_dist(order, r, start, dmax):
    """Augmented benchmark: expected fraction of tokens starting at `start` ending with unwrapped displacement d,
    d in [-dmax, dmax], by cohort propagation of the chain rules. Returns array indexed d + dmax."""
    N = len(order)
    pos = np.empty(N, dtype=int)
    for k, i in enumerate(order):
        pos[i] = k
    out = np.zeros(2 * dmax + 1)
    stack = [(start, 0, 1.0, -1)]
    while stack:
        i, d, m, came = stack.pop()
        if m < 1e-300:
            continue
        if pos[i] > came:
            out[d + dmax] += m * (1 - 2 * r)
            if abs(d) < dmax:
                stack.append(((i - 1) % N, d - 1, m * r, pos[i]))
                stack.append(((i + 1) % N, d + 1, m * r, pos[i]))
            else:
                out[d + dmax] += m * 2 * r      # truncation (never reached for dmax >= N with r <= 0.2 to 1e-12)
        else:
            out[d + dmax] += m
    return out


@nb.njit(cache=True)
def bench_one_step(n0, N, r, sched, rng_seed, dmax):
    """Stochastic benchmark: n0 tokens at node 0, one step, cohorts labelled by unwrapped displacement.
    Returns (counts over d in [-dmax, dmax], the order used)."""
    np.random.seed(rng_seed)
    order = np.zeros(N, dtype=np.int64)
    _order(sched, N, order)
    W = 2 * dmax + 1
    cnt = np.zeros((N, W), dtype=np.int64)       # tokens at node i with displacement d
    cnt[0, dmax] = n0
    probs = np.array([r, r])
    o = np.zeros(2, dtype=np.int64)
    if sched == J:
        new = np.zeros((N, W), dtype=np.int64)
        for i in range(N):
            for w in range(W):
                n = cnt[i, w]
                if n == 0:
                    continue
                _exits(n, probs, o)
                new[i, w] += n - o[0] - o[1]
                if w > 0:
                    new[(i - 1) % N, w - 1] += o[0]
                if w < W - 1:
                    new[(i + 1) % N, w + 1] += o[1]
        cnt = new
    else:
        for kk in range(N):
            i = order[kk]
            for w in range(W):
                n = cnt[i, w]
                if n == 0:
                    continue
                _exits(n, probs, o)
                cnt[i, w] -= o[0] + o[1]
                if w > 0:
                    cnt[(i - 1) % N, w - 1] += o[0]
                if w < W - 1:
                    cnt[(i + 1) % N, w + 1] += o[1]
    out = np.zeros(W, dtype=np.int64)
    for i in range(N):
        for w in range(W):
            out[w] += cnt[i, w]
    return out, order


# ------------------------------------------------------------------ measurements (proposal §3.3)
def centre(phiX_s, x, h):
    """Parabolic sub-grid minimum of smoothed phi_X (wrapped to [-L/2, L/2))."""
    N = phiX_s.size
    Ld = N * h
    i = int(np.argmin(phiX_s))
    ym, y0, yp = phiX_s[(i - 1) % N], phiX_s[i], phiX_s[(i + 1) % N]
    den = ym - 2 * y0 + yp
    off = 0.5 * (ym - yp) / den if den != 0 else 0.0
    return ((x[i] + off * h + Ld / 2) % Ld) - Ld / 2


def unwrap(c, Ld):
    c = np.asarray(c, dtype=float)
    out = c.copy()
    for k in range(1, c.size):
        d = c[k] - c[k - 1]
        d -= Ld * np.round(d / Ld)
        out[k] = out[k - 1] + d
    return out


def slope(t, y):
    t = np.asarray(t, float); y = np.asarray(y, float)
    tm = t.mean()
    return float(np.sum((t - tm) * (y - y.mean())) / np.sum((t - tm) ** 2))


def nyquist_free(f):
    """Remove the Nyquist Fourier mode (even N) so that sub-grid Fourier shifts are exactly invertible
    (pre-run specification, spec/MG8_prerun.md)."""
    N = f.shape[-1]
    if N % 2:
        return np.asarray(f, float)
    F = np.fft.fft(f, axis=-1)
    F[..., N // 2] = 0.0
    return np.real(np.fft.ifft(F, axis=-1))


def _err(phiA, phiB, nB, h, s):
    """max over species of ||phiA(x+s) - phiB(x)|| / ||phiB||."""
    sh = fourier_shift(phiA, h, s)
    return float(np.max(np.sqrt(np.sum((sh - phiB) ** 2, axis=-1)) / nB))


def shape_error(phiA, phiB, h, tol=1e-10):
    """Common-shift alignment (proposal §3.3, draft 4): sigma* = argmin_sigma max_s e_s(sigma); returns
    (E, sigma*, per-species e_s(sigma*)). phiA, phiB: arrays (3, N). Integer-node search, then golden section on
    [sigma_int - h/2, sigma_int + h/2]; ties -> smallest |sigma|, then smaller sigma."""
    N = phiA.shape[-1]
    phiA, phiB = nyquist_free(phiA), nyquist_free(phiB)
    nB = np.sqrt(np.sum(phiB ** 2, axis=-1))
    best = None
    for m in range(N):
        sft = (m if m <= N // 2 else m - N) * h
        e = _err(phiA, phiB, nB, h, sft)
        key = (e, abs(sft), sft)
        if best is None or key < best[0]:
            best = (key, sft)
    s0 = best[1]
    lo, hi = s0 - h / 2, s0 + h / 2
    gr = (np.sqrt(5) - 1) / 2
    c1, c2 = hi - gr * (hi - lo), lo + gr * (hi - lo)
    f1, f2 = _err(phiA, phiB, nB, h, c1), _err(phiA, phiB, nB, h, c2)
    while hi - lo > tol:
        if f1 < f2 or (f1 == f2 and (abs(c1), c1) <= (abs(c2), c2)):
            hi, c2, f2 = c2, c1, f1
            c1 = hi - gr * (hi - lo); f1 = _err(phiA, phiB, nB, h, c1)
        else:
            lo, c1, f1 = c1, c2, f2
            c2 = lo + gr * (hi - lo); f2 = _err(phiA, phiB, nB, h, c2)
    sg = 0.5 * (lo + hi)
    cand = [(_err(phiA, phiB, nB, h, sg), abs(sg), sg), best[0]]
    e, _, sg = min(cand)
    sh = fourier_shift(phiA, h, sg)
    per = np.sqrt(np.sum((sh - phiB) ** 2, axis=-1)) / nB
    return float(e), float(sg), per


# ------------------------------------------------------------------ labels (proposal §4.2, §4.3)
def settled(E, E1, v, v1):
    return abs(E - E1) <= 0.01 * E + 1e-7 and abs(v - v1) <= 0.1 * abs(v) + 1e-7


def label_production(single, E, Efloor, v, vfloor, is_settled):
    """Per-replicate labels of §4.2. Returns dict(fate, shape, drift, summary)."""
    if not single:
        return dict(fate="CHANGED", shape=None, drift=None, summary="SENSITIVE")
    if 10 * Efloor >= 1e-3:
        shape = "INSUFFICIENT RESOLUTION"
    elif not is_settled:
        shape = "UNSETTLED"
    elif E >= 1e-2:
        shape = "MATERIAL"
    elif E >= 1e-3:
        shape = "SMALL"
    else:
        shape = "NEGLIGIBLE"
    if 10 * vfloor >= 1e-6:
        drift = "INSUFFICIENT RESOLUTION"
    elif not is_settled:
        drift = "UNSETTLED"
    elif abs(v) >= 1e-6:
        drift = "DRIFTING"
    else:
        drift = "NOT DRIFTING AT THIS HORIZON"
    if shape == "MATERIAL" or drift == "DRIFTING":
        summ = "SENSITIVE"
    elif shape == "NEGLIGIBLE" and drift == "NOT DRIFTING AT THIS HORIZON":
        summ = "PRACTICALLY INSENSITIVE"
    elif shape == "SMALL" and drift == "NOT DRIFTING AT THIS HORIZON":
        summ = "MINOR"
    elif "INSUFFICIENT RESOLUTION" in (shape, drift):
        summ = "INSUFFICIENT RESOLUTION"
    else:
        summ = "UNRESOLVED"
    return dict(fate="SAME", shape=shape, drift=drift, summary=summ)


def combine_replicates(labels):
    s = [l["summary"] for l in labels]
    return s[0] if all(x == s[0] for x in s) else "SEAM-DEPENDENT"


def label_refinement(levels, horizon):
    """§4.3. levels: dict k -> dict(Ek, Efloor, S, all_single, all_settled) for k = 0..6.
    horizon: dict(single, settled, E300, E600, Efloor). Returns (label, info)."""
    from scipy import stats
    hz_rel = horizon["single"] and horizon["settled"] and \
        abs(horizon["E600"] - horizon["E300"]) <= 0.05 * horizon["E300"] + 1e-7
    elig = [k for k in range(2, 7) if levels[k]["all_single"] and levels[k]["all_settled"]
            and levels[k]["S"] <= 0.25 and levels[k]["Ek"] >= 10 * levels[k]["Efloor"]]
    info = dict(eligible=elig)
    if len(elig) >= 4 and (max(elig) - min(elig)) >= 3:
        lt = np.log([TAU0 / 2 ** k for k in elig]); le = np.log([levels[k]["Ek"] for k in elig])
        res = stats.linregress(lt, le)
        n = len(elig)
        lo = res.slope - stats.t.ppf(0.975, n - 2) * res.stderr
        resid = le - (res.intercept + res.slope * lt)
        info.update(p=res.slope, p_lo=lo, resid=resid.tolist())
        ok = (res.slope >= 0.7 and lo >= 0.5 and np.all(np.abs(resid) <= 0.15)
              and levels[max(elig)]["Ek"] <= 0.1 * levels[min(elig)]["Ek"] and hz_rel)
        if ok:
            return "CONVERGES TO J", info
    below = [k for k in range(2, 7) if levels[k]["Ek"] < 10 * levels[k]["Efloor"]]
    if len(elig) < 4 and below:
        floors_ok = all(10 * levels[k]["Efloor"] < 1e-4 for k in range(2, 7))
        all_ok = all(levels[k]["all_single"] and levels[k]["all_settled"] for k in range(7))
        seam_ok = all(levels[k]["S"] <= 0.25 for k in range(7) if levels[k]["Ek"] >= 10 * levels[k]["Efloor"])
        mono = all(levels[k + 1]["Ek"] <= levels[k]["Ek"] for k in range(6))
        hz = horizon["single"] and horizon["settled"] and (
            hz_rel if horizon["E300"] >= 10 * horizon["Efloor"] else horizon["E600"] < 1e-4)
        if floors_ok and all_ok and seam_ok and mono and levels[6]["Ek"] < 1e-4 and hz:
            return "BELOW MEASUREMENT RESOLUTION", info
    return "NOT RESOLVED BY THIS LADDER", info


# ------------------------------------------------------------------ helpers
def params(b=29.0):
    return {**P17, "b": b}


def reference_start(b=29.0):
    from runner import reference
    ref = reference(b, H)
    return ref, ref["S0"][:, None] + ref["phi"]


def steps_per_unit(tau):
    return int(round(1.0 / tau))

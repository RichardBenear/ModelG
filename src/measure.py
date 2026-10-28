"""MG004: measurement pipeline of protocol v1.3 §5.1–5.2 (frozen).

Operates on one run's recorded snapshots (counts n_S[t, i] on a periodic ring) and produces:
  - per-snapshot particle detections (smoothed phi_X minima below -0.5|A*|, merged within lambda*/2),
  - tracks and events (loss, birth, replication, merge; persistence 20 time units),
  - eligibility of snapshots and of the run,
  - Fourier-aligned run-level and K-block profiles on |x| < 10,
  - far-field Omega*Var for G, X, Y (nodes at circular distance > 20 from every particle).
"""
import numpy as np

SIGMA_SMOOTH = 0.4      # physical units (protocol §5.1)
PERSIST = 20.0          # time units
LINK = 1.5              # max circular link distance between snapshots
WINDOW = 10.0           # |x| < WINDOW for profiles
FAR = 20.0              # far-field circular distance


def circ(a, b, L):
    d = np.abs(a - b) % L
    return np.minimum(d, L - d)


def smooth_periodic(f, h, sigma=SIGMA_SMOOTH):
    N = f.shape[-1]
    k = 2 * np.pi * np.fft.fftfreq(N, d=h)
    return np.real(np.fft.ifft(np.fft.fft(f, axis=-1) * np.exp(-0.5 * (k * sigma) ** 2), axis=-1))


def fourier_shift(f, h, s):
    """Return g(x) = f(x + s) on the periodic lattice (band-limited interpolation)."""
    N = f.shape[-1]
    k = 2 * np.pi * np.fft.fftfreq(N, d=h)
    return np.real(np.fft.ifft(np.fft.fft(f, axis=-1) * np.exp(1j * k * s), axis=-1))


def detect(phiX_s, x, h, thr, lam):
    """Local minima (periodic) of smoothed phi_X below thr; merge within lam/2 keeping the deeper;
    sub-node position by parabola. Returns sorted array of positions."""
    N = phiX_s.size
    L = N * h
    l, r = np.roll(phiX_s, 1), np.roll(phiX_s, -1)
    idx = np.where((phiX_s < l) & (phiX_s <= r) & (phiX_s < thr))[0]
    if idx.size == 0:
        return np.zeros(0)
    order = idx[np.argsort(phiX_s[idx])]
    kept = []
    for i in order:
        if all(circ(x[i], x[j], L) >= lam / 2 for j in kept):
            kept.append(i)
    pos = []
    for i in kept:
        ym, y0, yp = phiX_s[(i - 1) % N], phiX_s[i], phiX_s[(i + 1) % N]
        den = ym - 2 * y0 + yp
        off = 0.5 * (ym - yp) / den if den != 0 else 0.0
        p = x[i] + off * h
        pos.append(((p + L / 2) % L) - L / 2)
    return np.sort(np.array(pos))


def track(dets, times, L, x0=0.0):
    """Greedy nearest-neighbour tracking. Track 0 is the original particle (expected near x0 at t0;
    pass x0=None for runs without an initial particle). Returns dict with per-snapshot index of track 0's
    detection (or -1), the list of tracks (start, end, positions), and events."""
    tracks = []      # each: dict(id, t_start, t_last, pos_last, hits:list of (snap, det_idx))
    orig_hit = -np.ones(len(times), dtype=int)
    if x0 is not None:
        tracks.append(dict(id=0, t_start=times[0], t_last=times[0] - 1e-9, pos=x0, hits=[], alive=True))
    for s, (t, D) in enumerate(zip(times, dets)):
        used = set()
        pairs = []
        for tr in tracks:
            if not tr["alive"]:
                continue
            for j, p in enumerate(D):
                pairs.append((circ(p, tr["pos"], L), tr["id"], j))
        pairs.sort()
        matched_tracks = {}
        det_owner = {}
        for dist, tid, j in pairs:
            if dist > LINK:
                break
            if tid in matched_tracks:
                continue
            if j in det_owner:
                # second track claiming the same detection -> merge candidate
                tr = next(q for q in tracks if q["id"] == tid)
                tr.setdefault("merge_into", det_owner[j]); tr.setdefault("merge_t", t)
                continue
            matched_tracks[tid] = j
            det_owner[j] = tid
        for tr in tracks:
            if not tr["alive"]:
                continue
            if tr["id"] in matched_tracks:
                j = matched_tracks[tr["id"]]
                tr["pos"], tr["t_last"] = D[j], t
                tr["hits"].append((s, j))
                used.add(j)
                if tr["id"] == 0:
                    orig_hit[s] = j
            elif t - tr["t_last"] >= PERSIST:
                tr["alive"] = False
        for j, p in enumerate(D):
            if j not in used and j not in det_owner:
                tracks.append(dict(id=len(tracks), t_start=t, t_last=t, pos=p, hits=[(s, j)], alive=True))
    # events (persistence-qualified)
    ev = []
    t_end = times[-1]
    for tr in tracks:
        dur = tr["t_last"] - tr["t_start"]
        if tr["id"] == 0:
            if t_end - tr["t_last"] >= PERSIST:
                ev.append(("loss", tr["t_last"]))
            if "merge_t" in tr:
                ev.append(("merge", tr["merge_t"]))
        elif dur >= PERSIST:
            # birth; replication if born within 2 lambda of a living particle (decided by caller via lam)
            ev.append(("birth", tr["t_start"], tr["hits"][0]))
    return dict(orig_hit=orig_hit, tracks=tracks, events=ev)


def analyze_run(nG, nX, nY, times, Om, h, S0, Astar, lam,  # Astar: SMOOTHED reference depth (v1.3.1) t_discard=10.0, K=6, T=100.0,
                min_snap_per_block=20, x0=0.0):
    """Full §5.1–5.2 analysis of one MG2-type run. Returns dict with eligibility, block profiles
    (K, 3, n_window) aligned to the original particle, far-field Omega*Var per field, and events."""
    N = nX.shape[1]
    L = N * h
    x = (np.arange(N) - N // 2) * h
    phi = np.stack([nG / Om - S0[0], nX / Om - S0[1], nY / Om - S0[2]], axis=1)   # (T_snap, 3, N)
    phiX_s = smooth_periodic(phi[:, 1, :], h)
    thr = -0.5 * abs(Astar)
    dets = [detect(phiX_s[s], x, h, thr, lam) for s in range(len(times))]
    tr = track(dets, times, L, x0=x0)
    events = []
    for e in tr["events"]:
        if e[0] == "birth":
            s, j = e[2]
            # replication if a pre-existing tracked particle is within 2 lambda at birth time
            pos_b = dets[s][j]
            near = [p for k, p in enumerate(dets[s]) if k != j and circ(p, pos_b, L) < 2 * lam]
            events.append(("replication" if near else "birth", float(e[1])))
        else:
            events.append((e[0], float(e[1])))
    win = np.abs(x) < WINDOW
    # eligible snapshot: exactly one detection, and it is the original track
    elig = np.array([(len(dets[s]) == 1 and tr["orig_hit"][s] == 0) for s in range(len(times))])
    use = times >= t_discard - 1e-9
    edges = np.linspace(t_discard, T, K + 1)
    blocks = np.full((K, 3, int(win.sum())), np.nan)
    nb = np.zeros(K, dtype=int)
    fs1 = np.zeros(3); fs2 = np.zeros(3); fcount = 0
    for s in np.where(use & elig)[0]:
        p = dets[s][0]
        al = fourier_shift(phi[s], h, p)          # particle moved to x = 0
        b = min(int(np.searchsorted(edges, times[s], side="right") - 1), K - 1)
        if nb[b] == 0:
            blocks[b] = 0.0
        blocks[b] += al[:, win]
        nb[b] += 1
    for s in np.where(use)[0]:
        far = np.ones(N, dtype=bool)
        for p in dets[s]:
            far &= circ(x, p, L) > FAR
        if far.sum() > 0:   # v1.3.1: snapshots without any particle contribute all nodes (C3 data rule)
            v = phi[s][:, far]
            fs1 += v.sum(1); fs2 += (v**2).sum(1); fcount += v.shape[1]
    for b in range(K):
        if nb[b] > 0:
            blocks[b] /= nb[b]
    # protocol §5.2: a run with any loss, replication or merge before T is not eligible
    disq = [e for e in events if e[0] in ("loss", "replication", "merge") and e[1] <= T]
    eligible = (len(disq) == 0) and bool(np.all(nb >= min_snap_per_block))
    # far-field variance about the run's own far-field mean (pooled over far nodes and snapshots)
    farvar = Om * (fs2 / fcount - (fs1 / fcount) ** 2) if fcount else np.full(3, np.nan)
    return dict(eligible=eligible, blocks=blocks, n_per_block=nb.tolist(), events=events,
                farvar=farvar.tolist(), n_elig_snap=int((use & elig).sum()),
                n_dets=[len(d) for d in dets])


def analyze_homogeneous(nG, nX, nY, times, Om, h, S0, Astar, lam, t_start=0.0):
    """Detection-only analysis for runs started homogeneous (MG2-N1, MG2-C5, MG3): per-snapshot detections
    and persistent births (time of first threshold crossing of a track that persists >= 20)."""
    N = nX.shape[1]
    L = N * h
    x = (np.arange(N) - N // 2) * h
    phiX = nX / Om - S0[1]
    phiX_s = smooth_periodic(phiX, h)
    thr = -0.5 * abs(Astar)
    dets = [detect(phiX_s[s], x, h, thr, lam) for s in range(len(times))]
    tr = track(dets, times, L, x0=None)
    births = sorted(float(e[1]) for e in tr["events"] if e[0] == "birth")
    first = None
    for t_ in births:
        if t_ >= t_start:
            first = t_
            break
    return dict(births=births, first_birth=first, n_dets=[len(d) for d in dets], dets=dets)

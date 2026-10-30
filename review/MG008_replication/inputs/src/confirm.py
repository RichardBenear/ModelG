"""MG004g: deterministic confirmation (protocol v1.3.3 §5.4, Amendment 3 as finalized).

Field-level test. From a field phi (3, N) at confirmation time, run the seed-free lattice ODE to t = 300 and
check t = 200, 210, ..., 300 (11 checks):
  - deep minima: local minima of the (unsmoothed) phi_X below -0.5 |A*| (A* = unsmoothed reference depth);
  - IDENTITY: a structure is persistent if a chain of deep minima can be linked across all 11 checks with
    each step's circular displacement <= lambda* (greedy nearest link starting from every deep minimum at
    t = 200).
Outcome classes:
  'persistent'  : >= 1 persistent particle-depth structure over continuation times 200-300
                  ("deterministically persistent particle-depth structure"; interpreted as basin membership,
                  with the finite-horizon limitation stated in the protocol);
  'excursion'   : max |phi| < 1e-3 at t = 300 (relaxed to the homogeneous state);
  'unresolved'  : otherwise.
Also returned: number of persistent structures, number of deep minima at t = 300, and whether a persistent
structure lies within 2 lambda* of a given detection position (detection-level attribution)."""
import numpy as np
from lattice import lattice_ode_1d, P17

CHECKS = np.arange(200.0, 300.0 + 1e-9, 10.0)


def _deep_minima(f, thr):
    l, r = np.roll(f, 1), np.roll(f, -1)
    return np.where((f < l) & (f <= r) & (f < thr))[0]


def confirm_field(phi, b, h, A_unsmoothed, lam, det_pos=None, t_end=300.0):
    N = phi.shape[1]
    L = N * h
    x = (np.arange(N) - N // 2) * h
    y0 = np.concatenate([phi[0], phi[1], phi[2]])
    xx, ts, Y, mdt = lattice_ode_1d({**P17, "b": b}, N, h, t_end, seed_scale=0.0, y0=y0, t_eval=CHECKS)
    thr = -0.5 * abs(A_unsmoothed)
    mins = [x[_deep_minima(Y[N:2 * N, k], thr)] for k in range(Y.shape[1])]
    circ = lambda a, c: np.minimum(np.abs(a - c) % L, L - np.abs(a - c) % L)
    persistent_ends = []
    for start in mins[0]:
        pos, ok = start, True
        for k in range(1, len(mins)):
            if mins[k].size == 0:
                ok = False; break
            d = circ(mins[k], pos)
            j = int(np.argmin(d))
            if d[j] > lam:
                ok = False; break
            pos = mins[k][j]
        if ok:
            persistent_ends.append(float(pos))
    persistent_ends = sorted(set(np.round(persistent_ends, 6)))
    yend = Y[:, -1]
    if persistent_ends:
        cls = "persistent"
    elif np.max(np.abs(yend)) < 1e-3:
        cls = "excursion"
    else:
        cls = "unresolved"
    near = None
    if det_pos is not None and persistent_ends:
        near = bool(min(circ(np.array(persistent_ends), det_pos)) <= 2 * lam)
    return dict(cls=cls, n_persistent=len(persistent_ends), n_deep_end=int(mins[-1].size),
                min_phiX_end=float(yend[N:2 * N].min()), maxdt_end=mdt, near_detection=near)

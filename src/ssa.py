"""MG004: exact Gillespie direct-method SSA of rules v1.1 on a periodic ring (validation V1 only).
Channels per node: R1 (0->G, Om*a), R2 (G->X, p nG), R1r (G->0, (q-p) nG), R2r (X->G, g nX), R3 (X->Y, b nX),
R5 (X->0, (1-g) nX), R4 (Y->X, nX(nX-1)nY/Om^2), hops of G, X, Y left/right (d_S/h^2 per token per direction).
Records counts on a regular time grid."""
import numpy as np
import numba as nb


@nb.njit(cache=True)
def _node_props(g_, x_, y_, Om, a, b, g, p, q, dx, dy, h, out):
    ih2 = 1.0 / (h * h)
    out[0] = Om * a
    out[1] = p * g_
    out[2] = (q - p) * g_
    out[3] = g * x_
    out[4] = b * x_
    out[5] = (1.0 - g) * x_
    out[6] = x_ * (x_ - 1) * y_ / (Om * Om)
    out[7] = g_ * ih2; out[8] = g_ * ih2
    out[9] = dx * x_ * ih2; out[10] = dx * x_ * ih2
    out[11] = dy * y_ * ih2; out[12] = dy * y_ * ih2
    s = 0.0
    for k in range(13):
        s += out[k]
    return s


@nb.njit(cache=True)
def ssa_ring(nG, nX, nY, Om, a, b, g, p, q, dx, dy, h, T, dt_rec, seed):
    np.random.seed(seed)
    N = nG.size
    props = np.zeros((N, 13))
    tot = np.zeros(N)
    for i in range(N):
        tot[i] = _node_props(nG[i], nX[i], nY[i], Om, a, b, g, p, q, dx, dy, h, props[i])
    nrec = int(T / dt_rec)
    rec = np.zeros((nrec, 3, N), dtype=np.int64)
    t = 0.0
    k = 0
    next_rec = dt_rec
    nev = 0
    while k < nrec:
        A = 0.0
        for i in range(N):
            A += tot[i]
        t += -np.log(np.random.random()) / A
        while k < nrec and t >= next_rec:
            for i in range(N):
                rec[k, 0, i] = nG[i]; rec[k, 1, i] = nX[i]; rec[k, 2, i] = nY[i]
            k += 1
            next_rec += dt_rec
        if k >= nrec:
            break
        r = np.random.random() * A
        i = 0
        while i < N - 1 and r >= tot[i]:
            r -= tot[i]
            i += 1
        c = 0
        while c < 12 and r >= props[i, c]:
            r -= props[i, c]
            c += 1
        il, ir = (i - 1) % N, (i + 1) % N
        j = -1
        if c == 0:
            nG[i] += 1
        elif c == 1:
            nG[i] -= 1; nX[i] += 1
        elif c == 2:
            nG[i] -= 1
        elif c == 3:
            nX[i] -= 1; nG[i] += 1
        elif c == 4:
            nX[i] -= 1; nY[i] += 1
        elif c == 5:
            nX[i] -= 1
        elif c == 6:
            nY[i] -= 1; nX[i] += 1
        elif c == 7:
            nG[i] -= 1; nG[il] += 1; j = il
        elif c == 8:
            nG[i] -= 1; nG[ir] += 1; j = ir
        elif c == 9:
            nX[i] -= 1; nX[il] += 1; j = il
        elif c == 10:
            nX[i] -= 1; nX[ir] += 1; j = ir
        elif c == 11:
            nY[i] -= 1; nY[il] += 1; j = il
        else:
            nY[i] -= 1; nY[ir] += 1; j = ir
        tot[i] = _node_props(nG[i], nX[i], nY[i], Om, a, b, g, p, q, dx, dy, h, props[i])
        if j >= 0:
            tot[j] = _node_props(nG[j], nX[j], nY[j], Om, a, b, g, p, q, dx, dy, h, props[j])
        nev += 1
    return rec, nev

"""MG042 deterministic explicit reservoirs. P17 only; no finite stochastic claim.

Rows: G,X,Y,uA,uB,F1,F2,F2r,F3,F4,F5,F1r,HG,HX,HY.
uA/uB are depleted fractions; rA=1-uA, rB=1-uB.
F are cumulative reaction concentration fluxes. HS is cumulative total
outgoing hopping flux, so net influx is (HS_left+HS_right)/2-HS.
Z=F3, W=F5. All initial fluxes are zero. Stocks use the Omega->infinity
normalization A0=C*14*300, B0=C*29*(140/9)*300, without integer rounding.
"""
from pathlib import Path
import numpy as np
import numba as nb
import scipy.sparse as sp
from scipy.integrate import BDF

ROOT = Path(__file__).resolve().parents[1]
CAPACITIES = (1, 10, 100, 1000, 10000, 100000, 1000000)
STEPS = {"E0": 1500, "E1": 3000, "E2": 6000}
S0 = np.array([140/9, 140/9, 29/(140/9)])
H = .2
DIFFUSION = np.array([1., 1., 12.])
# Species order for integer stoichiometry: G,X,Y,A,B,Z,W.
STOICH = np.array([[1,-1,1,0,0,0,-1], [0,1,-1,-1,1,-1,0],
                   [0,0,0,1,-1,0,0], [-1,0,0,0,0,0,1],
                   [0,0,0,-1,0,0,0], [0,0,0,1,0,0,0],
                   [0,0,0,0,0,1,0]], dtype=np.int64)


def stocks(capacity):
    return capacity * 14 * 300., capacity * 29 * S0[1] * 300.


def initial(start, n=500):
    state = np.zeros((15, n))
    state[:3] = S0[:, None]
    state[3:5] = 0.
    if start == "S":
        with np.load(ROOT / "results/MG003/ref_profiles.npz") as z:
            phi = z["h0.2_b29.0"][1:]
        if n != phi.shape[1]:
            raise ValueError("Reference start requires N=500")
        state[:3] += phi
    elif start != "H":
        raise ValueError(start)
    return state


@nb.njit(cache=True)
def advance_euler(state, capacity, finite, steps, nsteps):
    """All increments read the same active-state snapshot. No clipping."""
    n = state.shape[1]
    nxt = np.empty((5, n))
    dt = 1. / steps
    a0 = capacity * 4200.
    b0 = capacity * 29. * (140./9.) * 300.
    for _ in range(nsteps):
        for i in range(n):
            left, right = (i-1) % n, (i+1) % n
            g, x, y, ra, rb = state[0,i], state[1,i], state[2,i], 1.-state[3,i], 1.-state[4,i]
            f1, f2, f2r = 14.*ra, g, .1*x
            f3, f4, f5 = 29.*rb*x, x*x*y, .9*x
            lg = (state[0,left]+state[0,right]-2*g) / (.2*.2)
            lx = (state[1,left]+state[1,right]-2*x) / (.2*.2)
            ly = 12.*(state[2,left]+state[2,right]-2*y) / (.2*.2)
            nxt[0,i] = g + dt*(lg + f1-f2+f2r)
            nxt[1,i] = x + dt*(lx + f2-f2r-f3+f4-f5)
            nxt[2,i] = y + dt*(ly + f3-f4)
            nxt[3,i] = state[3,i]+dt*f1/a0 if finite else state[3,i]
            nxt[4,i] = state[4,i]+dt*f3/b0 if finite else state[4,i]
            state[5,i] += dt*f1
            state[6,i] += dt*f2
            state[7,i] += dt*f2r
            state[8,i] += dt*f3
            state[9,i] += dt*f4
            state[10,i] += dt*f5
            # F1r stays zero because q=p at P17.
            state[12,i] += dt*2*g/(.2*.2)
            state[13,i] += dt*2*x/(.2*.2)
            state[14,i] += dt*24*y/(.2*.2)
        state[:5] = nxt
    return state


def ode_rhs(t, flat, capacity, finite):
    """Independently expressed continuous RHS; never calls Euler increments."""
    z = flat.reshape(15, -1)
    g, x, y, ua, ub = z[:5]
    ra, rb = 1-ua, 1-ub
    out = np.zeros_like(z)
    a0, b0 = stocks(capacity)
    lap = (np.roll(z[:3], 1, axis=1) + np.roll(z[:3], -1, axis=1) - 2*z[:3]) / H**2
    out[0] = lap[0] - g + .1*x + 14*ra
    out[1] = lap[1] + g - (1+29*rb)*x + x**2*y
    out[2] = 12*lap[2] + 29*rb*x - x**2*y
    if finite:
        out[3] = 14*ra/a0
        out[4] = 29*rb*x/b0
    out[5] = 14*ra
    out[6] = g
    out[7] = .1*x
    out[8] = 29*rb*x
    out[9] = x**2*y
    out[10] = .9*x
    out[12:15] = 2*DIFFUSION[:, None]*z[:3]/H**2
    return out.ravel()


def ode_jac(t, flat, capacity, finite):
    z = flat.reshape(15, -1)
    n = z.shape[1]
    x, y, rb = z[1], z[2], 1-z[4]
    a0, b0 = stocks(capacity)
    # All blocks are N-by-N; row-major field ordering is frozen in pre-run spec.
    rows, cols, vals = [], [], []
    nodes = np.arange(n)
    def diagonal(r, c, value):
        rows.extend(r*n+nodes); cols.extend(c*n+nodes)
        vals.extend(np.broadcast_to(value, (n,)))
    for r, d in enumerate(DIFFUSION):
        for offset, value in ((-1,d/H**2), (0,-2*d/H**2), (1,d/H**2)):
            rows.extend(r*n+nodes); cols.extend(r*n+(nodes+offset)%n)
            vals.extend(np.full(n, value))
    diagonal(0,0,-1); diagonal(0,1,.1); diagonal(0,3,-14)
    diagonal(1,0,1); diagonal(1,1,-1-29*rb+2*x*y)
    diagonal(1,2,x*x); diagonal(1,4,29*x)
    diagonal(2,1,29*rb-2*x*y); diagonal(2,2,-x*x); diagonal(2,4,-29*x)
    if finite:
        diagonal(3,3,-14/a0)
        diagonal(4,1,29*rb/b0); diagonal(4,4,-29*x/b0)
    diagonal(5,3,-14); diagonal(6,0,1); diagonal(7,1,.1)
    diagonal(8,1,29*rb); diagonal(8,4,-29*x)
    diagonal(9,1,2*x*y); diagonal(9,2,x*x); diagonal(10,1,.9)
    for r,d in enumerate(DIFFUSION):
        diagonal(12+r,r,2*d/H**2)
    return sp.coo_matrix((vals, (rows, cols)), shape=(15*n,15*n)).tocsc()


def valid_state(z):
    # Never clip a numerically negative population or availability.
    return bool(np.all(np.isfinite(z)) and np.all(z[:5] >= 0) and np.all(z[3:5] <= 1) and np.all(z[5:] >= 0))


def trajectory(job, checkpoint=None):
    state = initial(job["start"], job.get("n", 500))
    horizon = job.get("T", 300)
    finite = job["mode"] == "FINITE"
    cap = job["capacity"]
    saved = [state.copy()]
    times = [0.]
    if job["method"] in STEPS:
        steps = STEPS[job["method"]]
        for t in range(1, horizon+1):
            advance_euler(state, cap, finite, steps, steps)
            if not valid_state(state):
                raise FloatingPointError(f"Invalid Euler state at t={t}")
            saved.append(state.copy()); times.append(float(t))
            if checkpoint and (t % 25 == 0 or t == horizon):
                checkpoint(np.array(times), np.array(saved))
        return np.array(times), np.array(saved)
    if job["method"] != "O":
        raise ValueError(job["method"])
    n = state.shape[1]
    # Eliminate identically constant coordinates, rather than letting sparse LU
    # roundoff perturb zero-depletion chemostats or the inactive R1r ledger.
    rows = [i for i in range(15) if i != 11 and (finite or i not in (3,4))]
    indices = np.concatenate([np.arange(i*n,(i+1)*n) for i in rows])
    def unpack(packed):
        full = np.zeros(15*n)
        full[indices] = packed
        return full
    atol = np.repeat([1e-11]*3+[1e-12]*2+[1e-8]*10, n)[indices]
    solver = BDF(lambda t,z: ode_rhs(t,unpack(z),cap,finite)[indices], 0., state.ravel()[indices], float(horizon),
                 rtol=1e-9, atol=atol, max_step=.5,
                 jac=lambda t,z: ode_jac(t,unpack(z),cap,finite)[indices,:][:,indices])
    next_time = 1
    while solver.status == "running":
        solver.step()
        if solver.status == "failed":
            raise RuntimeError("BDF step failed")
        if next_time <= solver.t:
            dense = solver.dense_output()
            while next_time <= solver.t and next_time <= horizon:
                snapshot = unpack(dense(next_time)).reshape(15, n)
                if not valid_state(snapshot):
                    raise FloatingPointError(f"Invalid ODE state at t={next_time}")
                saved.append(snapshot.copy()); times.append(float(next_time))
                if checkpoint and (next_time % 25 == 0 or next_time == horizon):
                    checkpoint(np.array(times), np.array(saved))
                next_time += 1
    return np.array(times), np.array(saved)


def accounting(times, data, job):
    """Reconstruct changes in seven species without huge-stock cancellation."""
    a0, b0 = stocks(job["capacity"])
    f = data[:, 5:12]
    hops = data[:, 12:15]
    net = (np.roll(hops,1,axis=2)+np.roll(hops,-1,axis=2))*.5-hops
    expected = np.einsum("sr,trn->tsn", STOICH, f)
    expected[:,:3] += net
    throughput = np.einsum("sr,trn->tsn", np.abs(STOICH), f)
    throughput[:,:3] += hops + (np.roll(hops,1,axis=2)+np.roll(hops,-1,axis=2))*.5
    external = np.zeros_like(expected)
    if job["mode"] == "CHEMOSTAT":
        external[:,3] = f[:,0]-f[:,6]
        external[:,4] = f[:,3]
    expected += external
    throughput += np.abs(external)
    actual = np.zeros_like(expected)
    actual[:,:3] = data[:,:3]-data[0,:3]
    actual[:,3] = -(data[:,3]-data[0,3])*a0
    actual[:,4] = -(data[:,4]-data[0,4])*b0
    actual[:,5] = data[:,8]
    actual[:,6] = data[:,10]
    residual = actual-expected
    per_site = np.max(np.abs(residual)/(1+throughput))
    global_ratio = np.max(np.abs(residual.sum(axis=2))/(1+throughput.sum(axis=2)))
    if job["mode"] == "FINITE":
        if job["method"] in STEPS:
            steps = STEPS[job["method"]]
            oracle = np.exp(times*steps*np.log1p(-1/(steps*job["capacity"]*300)))
        else:
            oracle = np.exp(-times/(job["capacity"]*300))
    else:
        oracle = np.ones(len(times))
    a_error = np.max(np.abs((1-data[:,3])-oracle[:,None]))
    return {"per_site_scaled_residual": float(per_site), "global_scaled_residual": float(global_ratio),
            "availability_oracle_error": float(a_error),
            "external_A_total": float(external[-1,3].sum()),
            "external_B_total": float(external[-1,4].sum()),
            "pass": bool(per_site <= 1e-8 and global_ratio <= 1e-8 and a_error <= 1e-8)}

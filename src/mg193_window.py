"""MG193: continuum 1D existence/stability window of the Model G soliton in (a, b).

Frozen specification: spec/MG193_existence_window.md (registered 7586da8).
Reuses the MG192 discretization (src/mg192_mobility.py) unchanged.
Run: python src/mg193_window.py
"""
import json
import os
import platform
import sys
import time
from multiprocessing import Pool

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import scipy
import mg192_mobility as m

OUT = os.path.join(m.ROOT, "results", "MG193_existence_window")
SEED = 20261116


# ---------------- discretization helpers ----------------

def odd_maps(N):
    M = N // 2 - 1
    E1 = np.zeros((N, M))
    for j in range(N):
        if j > N // 2:
            E1[j, j - N // 2 - 1] = 1.0
        elif 0 < j < N // 2:
            E1[j, N // 2 - j - 1] = -1.0
    rows1 = np.arange(N // 2 + 1, N)
    import scipy.linalg as sla
    E = sla.block_diag(E1, E1, E1)
    rows = np.concatenate([rows1, rows1 + N, rows1 + 2 * N])
    return E, rows


class Ctx:
    def __init__(self, L=100.0, N=1024):
        self.g = m.Grid(L, N)
        self.D2 = self.g.D2()
        self.Ee, self.re = m.even_maps(N)
        self.Eo, self.ro = odd_maps(N)


def depth(ctx, u, a):
    return float(u[1, ctx.g.N // 2] - m.hom(a=a)[1])


def newton(ctx, u, a, b, maxit=20, tol=1e-10):
    g = ctx.g
    for it in range(maxit):
        F = m.rhs(g, u, a=a, b=b).ravel()
        r = float(np.max(np.abs(F)))
        if r <= tol:
            ok = bool(np.all(u > 0)) and depth(ctx, u, a) <= -1.0
            return u, r, ok, ("ok" if ok else "decay")
        J = m.jac(g, u, ctx.D2, b=b)
        Jr = J[ctx.re] @ ctx.Ee
        try:
            dr = np.linalg.solve(Jr, -F[ctx.re])
        except np.linalg.LinAlgError:
            return u, r, False, "singular"
        u = u + (ctx.Ee @ dr).reshape(3, g.N)
        if not np.all(np.isfinite(u)) or np.any(u <= 0):
            return u, np.inf, False, "positivity"
    return u, r, False, "noconv"


def eigs(ctx, u, b):
    J = m.jac(ctx.g, u, ctx.D2, b=b)
    le = np.linalg.eigvals(J[ctx.re] @ ctx.Ee)
    lo = np.linalg.eigvals(J[ctx.ro] @ ctx.Eo)
    i0 = int(np.argmin(np.abs(lo)))
    lam0 = lo[i0]
    lo = np.delete(lo, i0)
    me, mo = float(np.max(le.real)), float(np.max(lo.real))
    ie, io = int(np.argmax(le.real)), int(np.argmax(lo.real))
    return dict(max_re_even=me, im_even=float(le[ie].imag), max_re_odd=mo, im_odd=float(lo[io].imag),
                translation=float(abs(lam0)))


def turing_threshold(a):
    X0 = m.hom(a=a)[1]
    G0 = X0 * 0 + m.hom(a=a)[0]

    def maxre(b):
        Y0 = b / X0
        best = -np.inf
        for k in np.linspace(0.01, 6, 2000):
            A = np.array([[-k * k - 1, 0.1, 0],
                          [1, -k * k - (1 + b) + 2 * X0 * Y0, X0 * X0],
                          [0, b - 2 * X0 * Y0, -12 * k * k - X0 * X0]])
            best = max(best, np.max(np.linalg.eigvals(A).real))
        return best

    lo, hi = 1.0, 200.0
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        if maxre(mid) > 0:
            hi = mid
        else:
            lo = mid
    return 0.5 * (lo + hi)


# ---------------- continuation ----------------

def param(p, base, s):
    """p='b': absolute b = base + s; p='a': a = base * (1 + s)."""
    return base + s if p == "b" else base * (1 + s)


def continue_line(ctx, u0, a0, b0, p, direction, step0, smax, stop_value=None, log=print, tol=1e-10, use_eigs=True):
    """Natural continuation in parameter p from (a0, b0). Returns rows and the edge record."""
    base = b0 if p == "b" else a0
    s, step, u = 0.0, step0, u0.copy()
    rows = []

    def ab(sv):
        return (a0, param("b", b0, sv)) if p == "b" else (param("a", a0, sv), b0)

    a, b = ab(0.0)
    u, r, ok, why = newton(ctx, u, a, b, tol=tol)
    if not ok:
        return rows, dict(edge="START_FAILED", why=why)
    ev = eigs(ctx, u, b) if use_eigs else dict(max_re_even=-1.0, max_re_odd=-1.0, im_even=0.0, im_odd=0.0, translation=0.0)
    stable = ev["max_re_even"] < 0 and ev["max_re_odd"] < 0
    rows.append(dict(s=0.0, a=a, b=b, residual=r, depth=depth(ctx, u, a), stable=stable, **ev))
    prev_u, prev_s, prev_stable = u, 0.0, stable
    while True:
        s_try = prev_s + direction * step
        if abs(s_try) > smax:
            return rows, dict(edge="LIMIT", last_s=prev_s)
        a, b = ab(s_try)
        if stop_value is not None and p == "b" and direction > 0 and b >= stop_value:
            return rows, dict(edge="TURING", last_s=prev_s, b_c=stop_value)
        u, r, ok, why = newton(ctx, prev_u.copy(), a, b, tol=tol)
        if not ok:
            if step < 1e-5:
                lab = "DECAY" if why == "decay" else "FOLD"
                return rows, dict(edge=lab, why=why, last_ok_s=prev_s, bracket=[prev_s, s_try],
                                  last_ok_param=ab(prev_s), fail_param=[a, b])
            step /= 2
            continue
        ev = eigs(ctx, u, b) if use_eigs else dict(max_re_even=-1.0, max_re_odd=-1.0, im_even=0.0, im_odd=0.0, translation=0.0)
        stable = ev["max_re_even"] < 0 and ev["max_re_odd"] < 0
        if stable != prev_stable:
            # bisect the stability crossing
            lo_s, hi_s, ulo = prev_s, s_try, prev_u
            while abs(hi_s - lo_s) > 1e-5:
                mid = 0.5 * (lo_s + hi_s)
                am, bm = ab(mid)
                um, rm, okm, _ = newton(ctx, ulo.copy(), am, bm)
                if not okm:
                    hi_s = mid
                    continue
                evm = eigs(ctx, um, bm)
                stm = evm["max_re_even"] < 0 and evm["max_re_odd"] < 0
                if stm == prev_stable:
                    lo_s, ulo = mid, um
                else:
                    hi_s = mid
            rows.append(dict(s=s_try, a=a, b=b, residual=r, depth=depth(ctx, u, a), stable=stable, **ev))
            if prev_stable and not stable:
                which = "even" if ev["max_re_even"] >= 0 else "odd"
                return rows, dict(edge="STABILITY_LOSS", mode=which, bracket=[lo_s, hi_s],
                                  last_stable_param=ab(lo_s), unstable_param=ab(hi_s), last_ok_s=lo_s)
        rows.append(dict(s=s_try, a=a, b=b, residual=r, depth=depth(ctx, u, a), stable=stable, **ev))
        log(f"  {p} dir {direction:+d} s={s_try:.6f} a={a:.6f} b={b:.6f} depth={rows[-1]['depth']:.4f} "
            f"maxRe e/o {ev['max_re_even']:.4f}/{ev['max_re_odd']:.4f}")
        prev_u, prev_s, prev_stable = u, s_try, stable
        step = min(step * 1.5, step0 * 4)


def solution_at(ctx, s_target, rows_unused, u0, a0, b0, p):
    """Re-solve by stepping from (a0,b0) to s_target in small steps (helper for checks)."""
    n = max(1, int(abs(s_target) / (0.01 if p == "b" else 2e-3)))
    u = u0.copy()
    for i in range(1, n + 1):
        sv = s_target * i / n
        a, b = (a0, b0 + sv) if p == "b" else (a0 * (1 + sv), b0)
        u, r, ok, why = newton(ctx, u, a, b)
        if not ok:
            return None
    return u


# ---------------- time-dependent checks ----------------

class Stepper(m.ETDRK4):
    def __init__(self, grid, dt, b):
        super().__init__(grid, dt)
        self.b = b

    def nl(self, vh, a):
        G, X, Y = np.real(np.fft.ifft(vh, axis=1))
        b = self.b
        n = np.array([-G + 0.1 * X + a, G - (1 + b) * X + X * X * Y, b * X - X * X * Y])
        return np.fft.fft(n, axis=1)


def evolve(grid, u0, a, b, T, dt=0.002):
    st = Stepper(grid, dt, b)
    vh = np.fft.fft(u0, axis=1)
    out = []
    n = int(round(T / dt))
    for i in range(1, n + 1):
        vh = st.step(vh, a)
        if i % int(round(50 / dt)) == 0:
            u = np.real(np.fft.ifft(vh, axis=1))
            X0 = m.hom(a=a, b=b)[1]
            dX = u[1] - X0
            mins = int(np.sum((dX < -1.0) & (dX <= np.roll(dX, 1)) & (dX <= np.roll(dX, -1))))
            out.append(dict(t=i * dt, depth=float(u[1, grid.N // 2] - X0), min_depth=float(dX.min()),
                            deep_minima=mins))
    return out


# ---------------- tasks for the pool ----------------

def task_bline(args):
    a_target, bstar, U14 = args
    ctx = Ctx()
    logs = []
    log = lambda s: logs.append(s)
    res = dict(a=a_target)
    # move from a=14 to a_target at b = bstar
    u = U14.copy()
    if abs(a_target - 14.0) > 1e-12:
        n = int(np.ceil(abs(a_target - 14.0) / 0.02))
        for i in range(1, n + 1):
            ai = 14.0 + (a_target - 14.0) * i / n
            u, r, ok, why = newton(ctx, u, ai, bstar)
            if not ok:
                res["reach"] = dict(ok=False, failed_at_a=ai, why=why)
                return res, logs
    res["reach"] = dict(ok=True, b=bstar)
    bc = turing_threshold(a_target)
    res["b_c"] = bc
    up_rows, up_edge = continue_line(ctx, u, a_target, bstar, "b", +1, 0.02, 5.0, stop_value=bc, log=log)
    dn_rows, dn_edge = continue_line(ctx, u, a_target, bstar, "b", -1, 0.02, 5.0, log=log)
    res.update(up=dict(edge=up_edge, rows=up_rows), down=dict(edge=dn_edge, rows=dn_rows))
    return res, logs


def save(rec):
    with open(os.path.join(OUT, "window.json"), "w") as f:
        json.dump(rec, f, indent=1, default=float)


def main():
    os.makedirs(OUT, exist_ok=True)
    t_start = time.time()
    logf = open(os.path.join(OUT, "run_log.txt"), "w")

    def log(s):
        print(s, flush=True)
        logf.write(s + "\n")
        logf.flush()

    log(f"MG193 start {time.strftime('%Y-%m-%d %H:%M:%S UTC', time.gmtime())} numpy {np.__version__} scipy {scipy.__version__}")
    U14 = np.load(os.path.join(m.ROOT, "results", "MG192_gravity_mobility", "soliton_adjoint_L100_N1024.npz"))["U"]
    ctx = Ctx()
    rec = dict(spec="spec/MG193_existence_window.md", registration_commit="7586da8")

    # 1. a = 14 line in b, from b = 29
    bc14 = turing_threshold(14.0)
    log(f"b_c(14) = {bc14:.5f}")
    up_rows, up_edge = continue_line(ctx, U14, 14.0, 29.0, "b", +1, 0.02, 5.0, stop_value=bc14, log=log)
    dn_rows, dn_edge = continue_line(ctx, U14, 14.0, 29.0, "b", -1, 0.002, 5.0, log=log)
    log(f"a=14 up edge {up_edge}")
    log(f"a=14 down edge {dn_edge}")

    def edge_value(edge, base):
        s = edge.get("last_ok_s", edge.get("last_s"))
        if edge["edge"] == "STABILITY_LOSS":
            s = 0.5 * sum(edge["bracket"])
        elif "bracket" in edge:
            s = 0.5 * sum(edge["bracket"])
        return base + s

    b_lo14 = edge_value(dn_edge, 29.0)
    b_hi14 = edge_value(up_edge, 29.0)
    bstar = 0.5 * (b_lo14 + b_hi14)
    rec["a14"] = dict(b_c=bc14, up=dict(edge=up_edge, rows=up_rows), down=dict(edge=dn_edge, rows=dn_rows),
                      b_low=b_lo14, b_high=b_hi14, b_star=bstar)
    log(f"a=14 window [{b_lo14:.5f}, {b_hi14:.5f}], centre b* = {bstar:.5f}")

    # 2. a-line at b = 29 (and at b*)
    aline = {}
    for bb in (29.0, bstar):
        u = U14 if bb == 29.0 else solution_at(ctx, bstar - 29.0, None, U14, 14.0, 29.0, "b")
        aup, aup_e = continue_line(ctx, u, 14.0, bb, "a", +1, 2e-5 if bb == 29.0 else 2e-3, 0.5, log=log)
        adn, adn_e = continue_line(ctx, u, 14.0, bb, "a", -1, 2e-3, 0.5, log=log)
        aline[f"b{bb:.5f}"] = dict(b=bb, up=dict(edge=aup_e, rows=aup), down=dict(edge=adn_e, rows=adn),
                                   a_high=14.0 * (1 + 0.5 * sum(aup_e["bracket"])) if "bracket" in aup_e else None,
                                   a_low=14.0 * (1 + 0.5 * sum(adn_e["bracket"])) if "bracket" in adn_e else None)
        log(f"a-line at b={bb:.5f}: up {aup_e.get('edge')} {aline[f'b{bb:.5f}']['a_high']}, "
            f"down {adn_e.get('edge')} {aline[f'b{bb:.5f}']['a_low']}")
    rec["a_lines"] = aline
    save(rec)

    # 3. b-lines at other a (pool), starting at b*
    with Pool(4) as pool:
        outs = pool.map(task_bline, [(a, bstar, U14) for a in (12.0, 13.0, 15.0, 16.0)])
    rec["b_lines"] = {}
    for res, logs in outs:
        for s in logs:
            logf.write(s + "\n")
        rec["b_lines"][f"a{res['a']}"] = res
        if res["reach"]["ok"]:
            lo = edge_value(res["down"]["edge"], bstar)
            hi = edge_value(res["up"]["edge"], bstar)
            res["b_low"], res["b_high"] = lo, hi
            log(f"a={res['a']}: window [{lo:.5f}, {hi:.5f}] edges {res['down']['edge']['edge']}/{res['up']['edge']['edge']}")
        else:
            log(f"a={res['a']}: could not reach from a=14 at b*={bstar:.5f}: {res['reach']}")
    save(rec)

    # 4. qualification: a=14 edges at N=2048 and L=200
    qual = {}
    for (L, N) in [(100.0, 2048), (200.0, 2048)]:
        cq = Ctx(L, N)
        g1 = ctx.g
        # interpolate the b* solution onto the new grid
        ustar = solution_at(ctx, bstar - 29.0, None, U14, 14.0, 29.0, "b")
        xq = cq.g.x
        uq = np.empty((3, N))
        for s_ in range(3):
            Fh = np.fft.fft(ustar[s_]) / g1.N
            kk = g1.k.copy()
            kk[g1.N // 2] = 0.0
            inside = np.abs(xq) < 50.0
            vals = np.real(np.exp(1j * np.outer(xq[inside] + 50.0, kk)) @ Fh)
            uq[s_] = m.hom()[s_]
            uq[s_][inside] = vals
        b_c = bc14
        # Deviation D1 (recorded in the report): Newton-only continuation, residual tolerance 1e-9 at N=2048,
        # because the FFT round-off floor there is ~1e-10 and the 1e-10 test failed intermittently in attempt 1.
        upq, upq_e = continue_line(cq, uq, 14.0, bstar, "b", +1, 0.02, 5.0, stop_value=b_c, log=log, tol=1e-9, use_eigs=False)
        dnq, dnq_e = continue_line(cq, uq, 14.0, bstar, "b", -1, 0.02, 5.0, log=log, tol=1e-9, use_eigs=False)
        lo, hi = edge_value(dnq_e, bstar), edge_value(upq_e, bstar)
        qual[f"L{int(L)}_N{N}"] = dict(b_low=lo, b_high=hi, up_edge=upq_e, down_edge=dnq_e,
                                      diff_low=abs(lo - b_lo14), diff_high=abs(hi - b_hi14),
                                      ok=bool(abs(lo - b_lo14) <= 1e-3 and abs(hi - b_hi14) <= 1e-3))
        log(f"qualification {L},{N}: [{lo:.5f},{hi:.5f}] ok={qual[f'L{int(L)}_N{N}']['ok']}")
    rec["qualification"] = qual
    save(rec)

    # 5. time-dependent checks (separate check)
    g = ctx.g
    checks = {}
    jobs = []
    off_b = 0.01
    jobs.append(("b_low_inside", 14.0, b_lo14 + off_b))
    jobs.append(("b_low_outside", 14.0, b_lo14 - off_b))
    jobs.append(("b_high_inside", 14.0, b_hi14 - off_b))
    jobs.append(("b_high_outside", 14.0, b_hi14 + off_b))
    al = aline["b29.00000"]
    if al["a_high"]:
        jobs.append(("a_high_inside_b29", al["a_high"] * (1 - 5e-4), 29.0))
        jobs.append(("a_high_outside_b29", al["a_high"] * (1 + 5e-4), 29.0))
    if al["a_low"]:
        jobs.append(("a_low_inside_b29", al["a_low"] * (1 + 5e-4), 29.0))
        jobs.append(("a_low_outside_b29", al["a_low"] * (1 - 5e-4), 29.0))
    for name, a, b in jobs:
        # start from the nearest converged stationary profile (inside), placed on the target parameters
        if "a_" in name:
            s_in = (a / 14.0 - 1) if "inside" in name else ((a / 14.0 - 1) - np.sign(a / 14.0 - 1) * 1e-3)
            u0 = solution_at(ctx, s_in, None, U14, 14.0, 29.0, "a")
        else:
            b_in = b if "inside" in name else (b + off_b * 2 if "low" in name else b - off_b * 2)
            u0 = solution_at(ctx, b_in - 29.0, None, U14, 14.0, 29.0, "b")
        if u0 is None:
            checks[name] = dict(a=a, b=b, result="NO_START_PROFILE")
            continue
        t0 = time.time()
        hist = evolve(g, u0, a, b, 500.0)
        d0 = float(u0[1, g.N // 2] - m.hom(a=a, b=b)[1])
        last = hist[-1]
        if "inside" in name:
            verdict = "PERSISTS" if abs(last["depth"] - d0) <= 0.01 * abs(d0) and last["deep_minima"] == 1 else "CHANGED"
        else:
            if last["min_depth"] > -1.0:
                verdict = "DISSOLVED"
            elif last["deep_minima"] > 1 or abs(last["depth"] - d0) > 0.1 * abs(d0):
                verdict = "CHANGED"
            else:
                verdict = "INCONCLUSIVE"
        checks[name] = dict(a=a, b=b, start_depth=d0, history=hist, verdict=verdict, wall_s=time.time() - t0)
        log(f"check {name}: a={a:.6f} b={b:.5f} -> {verdict} (depth {d0:.3f} -> {last['depth']:.3f}, "
            f"minima {last['deep_minima']})")
    rec["time_checks"] = checks
    save(rec)

    # 6. recommended point and perturbed stability
    rec_pt = dict(a=14.0, b=bstar, margin_b_low=bstar - b_lo14, margin_b_high=b_hi14 - bstar,
                  a_window_at_bstar=[aline[f"b{bstar:.5f}"]["a_low"], aline[f"b{bstar:.5f}"]["a_high"]])
    ustar = solution_at(ctx, bstar - 29.0, None, U14, 14.0, 29.0, "b")
    rng = np.random.default_rng(SEED)
    pert = np.zeros_like(ustar)
    for s_ in range(3):
        coef = rng.normal(size=8) + 1j * rng.normal(size=8)
        f = np.zeros(g.N, dtype=complex)
        f[1:9] = coef
        r = np.real(np.fft.ifft(f))
        r /= np.max(np.abs(r))
        pert[s_] = 0.01 * np.max(np.abs(ustar[s_] - m.hom(b=bstar)[s_])) * r
    hist = evolve(g, ustar + pert, 14.0, bstar, 200.0)
    d0 = float(ustar[1, g.N // 2] - m.hom(b=bstar)[1])
    dmin = min(h_["min_depth"] for h_ in hist[-1:])
    rec_pt["perturbed_run"] = dict(start_depth=d0, final_min_depth=hist[-1]["min_depth"],
                                   deep_minima=hist[-1]["deep_minima"],
                                   ok=bool(abs(hist[-1]["min_depth"] - d0) <= 0.01 * abs(d0) and hist[-1]["deep_minima"] == 1))
    rec["recommended_point"] = rec_pt
    log(f"recommended point {rec_pt}")

    rec["wall_seconds"] = time.time() - t_start
    rec["platform"] = dict(platform=platform.platform(), python=sys.version.split()[0], numpy=np.__version__,
                           scipy=scipy.__version__)
    with open(os.path.join(OUT, "window.json"), "w") as f:
        json.dump(rec, f, indent=1, default=float)
    log(f"done, wall {rec['wall_seconds']:.0f} s")


if __name__ == "__main__":
    main()

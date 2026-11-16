"""MG194: gravity drift confirmation (Part A) and two-particle interaction (Part B), continuum 1D.

Specification: spec/MG194_two_particle_gravity.md. Numerics reuse src/mg192_mobility.py unchanged, with b set to the
MG193 window centre. Run: python src/mg194_two_particle.py   (MODELG_WORKERS processes, default 4)
"""
import json
import os
import platform
import sys
import time
from multiprocessing import Pool

import numpy as np
import scipy
import scipy.linalg as sla
from scipy.optimize import least_squares

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import mg192_mobility as m  # noqa: E402

BSTAR = 29.205258852005002
m.P["b"] = BSTAR
OUT = os.path.join(m.ROOT, "results", "MG194_two_particle")
KAPPA0, K0, KAPPA_REAL = 0.4608152900342353, 2.04914, 0.9657108197751885
HALF_WL = np.pi / K0
WORKERS = int(os.environ.get("MODELG_WORKERS", "4"))


def soliton(grid):
    """Stationary soliton at b* on grid (L=100 N=1024 from MG192's profile, continued in b)."""
    U = np.load(os.path.join(m.ROOT, "results", "MG192_gravity_mobility", "soliton_adjoint_L100_N1024.npz"))["U"]
    nb = 21
    for i in range(1, nb + 1):
        m.P["b"] = 29.0 + (BSTAR - 29.0) * i / nb
        U, res, _ = m.newton_even(grid, U, tol=1e-10)
    m.P["b"] = BSTAR
    if grid.N != 1024 or grid.L != 100.0:
        raise ValueError("soliton() builds the L=100, N=1024 profile; use embed() for other grids")
    return U, res


def embed(g_from, U, g_to):
    """Put the L=100 soliton on a larger periodic grid (homogeneous outside |x| < 50), then refine by Newton."""
    h0 = m.hom()
    Fh = np.fft.fft(U - h0[:, None], axis=1) / g_from.N
    kk = g_from.k.copy()
    kk[g_from.N // 2] = 0.0
    inside = np.abs(g_to.x) < g_from.L / 2
    u = np.repeat(h0[:, None], g_to.N, axis=1).astype(float)
    ph = np.exp(1j * np.outer(g_to.x[inside] + g_from.L / 2, kk))
    for s in range(3):
        u[s, inside] += np.real(ph @ Fh[s])
    u, res, _ = m.newton_even(g_to, u, tol=1e-10)
    return u, res


def pair_state(grid, U, d):
    h0 = m.hom()
    phi = U - h0[:, None]
    return h0[:, None] + np.array([grid.shift(phi[s], d / 2) + grid.shift(phi[s], -d / 2) for s in range(3)])


def refine_min(grid, X, s0):
    Xh = np.fft.fft(X) / grid.N
    k = grid.k.copy()
    k[grid.N // 2] = 0.0
    s = s0
    for _ in range(30):
        e = np.exp(1j * k * (s + grid.L / 2))
        d1 = np.real(np.sum(1j * k * Xh * e))
        d2 = np.real(np.sum(-(k ** 2) * Xh * e))
        ds = -d1 / d2
        s += ds
        if abs(ds) < 1e-14:
            break
    return s


def two_minima(grid, X, X0):
    """Positions of the deepest minimum in x<0 and in x>0, refined; and their depths."""
    left, right = grid.x < 0, grid.x > 0
    jl = np.where(left)[0][np.argmin(X[left])]
    jr = np.where(right)[0][np.argmin(X[right])]
    xl, xr = refine_min(grid, X, grid.x[jl]), refine_min(grid, X, grid.x[jr])
    return xl, xr, float(X[jl] - X0), float(X[jr] - X0)


# ---------------- Part A ----------------

def part_a_run(args):
    eps, x0 = args
    m.P["b"] = BSTAR
    g = m.Grid(100.0, 1024)
    U = np.load(os.path.join(OUT, "soliton_bstar.npz"))["U"]
    bg, a, bres = m.modulated_background(g, eps)
    u0 = bg + np.array([g.shift(U[s], x0) for s in range(3)]) - m.hom()[:, None]
    st = m.ETDRK4(g, 0.002)
    vh = np.fft.fft(u0, axis=1)
    ts, xs, ds = [0.0], [m.xmin_position(g, u0[1], x0)], [float(np.min(u0[1] - bg[1]))]
    t0 = time.time()
    for n in range(1, 200001):
        vh = st.step(vh, a)
        if n % 500 == 0:
            u = np.real(np.fft.ifft(vh, axis=1))
            if not np.all(np.isfinite(u)) or np.any(u <= 0):
                return dict(eps=eps, x0=x0, error="positivity/finiteness lost", t=n * 0.002)
            ts.append(n * 0.002)
            xs.append(m.xmin_position(g, u[1], xs[-1]))
            ds.append(float(np.min(u[1] - bg[1])))
    ts, ds = np.array(ts), np.array(ds)
    xs = np.unwrap(np.array(xs) * 2 * np.pi / g.L) * g.L / (2 * np.pi)
    sel = ts >= 100.0
    v = float(np.polyfit(ts[sel], xs[sel], 1)[0])
    dadx = float(14.0 * eps * 2 * np.pi / g.L * np.cos(2 * np.pi * x0 / g.L))
    return dict(eps=eps, x0=x0, v=v, dadx=dadx, bg_residual=bres, depth_start=float(ds[0]), depth_end=float(ds[-1]),
                depth_maxdev=float(np.max(np.abs(ds - ds[0])) / abs(ds[0])), wall_s=time.time() - t0,
                series=np.vstack([ts, xs, ds]).tolist())


# ---------------- Part B ----------------

def pair_run(args):
    d, L, N = args
    m.P["b"] = BSTAR
    g = m.Grid(L, N)
    U = np.load(os.path.join(OUT, f"soliton_bstar{'' if L == 100.0 else '_L' + str(int(L))}.npz"))["U"]
    u0 = pair_state(g, U, d) if d > 0 else U.copy()
    X0 = m.hom()[1]
    st = m.ETDRK4(g, 0.002)
    vh = np.fft.fft(u0, axis=1)
    a = np.full(g.N, 14.0)
    t0 = time.time()
    if d > 0:
        xl, xr, dl, dr = two_minima(g, u0[1], X0)
    else:
        xl, xr, dl, dr = 0.0, m.xmin_position(g, u0[1], 0.0), 0.0, float(u0[1].min() - X0)
    ts, sep, depth = [0.0], [xr - xl], [max(dl, dr) if d > 0 else dr]
    merged = False
    for n in range(1, 50001):
        vh = st.step(vh, a)
        if n % 500 == 0:
            u = np.real(np.fft.ifft(vh, axis=1))
            if not np.all(np.isfinite(u)) or np.any(u <= 0):
                return dict(d=d, L=L, N=N, error="positivity/finiteness lost", t=n * 0.002)
            if d > 0:
                xl, xr, dl, dr = two_minima(g, u[1], X0)
                core = max(dl, dr)
                if core > -1.0 or xr - xl < 2.0:
                    merged = True
            else:
                xr = m.xmin_position(g, u[1], xr)
                core = float(u[1].min() - X0)
            ts.append(n * 0.002)
            sep.append(xr - xl)
            depth.append(core)
    ts, sep = np.array(ts), np.array(sep)
    sel = ts >= 20.0
    v = float(np.polyfit(ts[sel], sep[sel], 1)[0])
    return dict(d=d, L=L, N=N, v_rel=v, delta_d=float(sep[-1] - sep[0]), moving=bool(abs(sep[-1] - sep[0]) > 0.25),
                merged=merged, depth_start=float(depth[0]), depth_end=float(depth[-1]), wall_s=time.time() - t0,
                series=np.vstack([ts, sep, depth]).tolist())


def lin_curve(g, U, psi, ds):
    h0 = m.hom()
    X, Y = U[1], U[2]
    X0, Y0 = h0[1], h0[2]
    phi = U - h0[:, None]
    out = []
    for d in ds:
        du = np.array([g.shift(phi[s], d) for s in range(3)])
        PX = (2 * X * Y - 2 * X0 * Y0) * du[1] + (X * X - X0 * X0) * du[2]
        v = -(np.dot(psi[1], PX) + np.dot(psi[2], -PX))
        out.append(-2.0 * v)
    return np.array(out)


def pair_equilibrium(args):
    d0 = args
    m.P["b"] = BSTAR
    g = m.Grid(100.0, 1024)
    U = np.load(os.path.join(OUT, "soliton_bstar.npz"))["U"]
    X0 = m.hom()[1]
    try:
        u, res, hist = m.newton_even(g, pair_state(g, U, d0), tol=1e-10, maxit=40)
    except RuntimeError as e:
        return dict(d0=d0, converged=False, why=str(e))
    xl, xr, dl, dr = two_minima(g, u[1], X0)
    ok = res <= 1e-10 and max(dl, dr) < -1.0 and xr - xl > 2.0
    rec = dict(d0=d0, converged=bool(ok), residual=res, d_star=float(xr - xl), depths=[dl, dr])
    if not ok:
        return rec
    J = m.jac(g, u)
    lam = sla.eigvals(J)
    # translation mode: eigenvalue nearest 0 whose eigenvector is odd -> identify via odd-subspace spectrum
    Up = np.array([g.ddx(u[s]) for s in range(3)]).ravel()
    tr = np.linalg.norm(J @ Up) / np.linalg.norm(Up)
    i0 = int(np.argmin(np.abs(lam)))
    rest = np.delete(lam, i0)
    lead = rest[np.argmax(np.real(rest))]
    rec.update(translation_eig=complex(lam[i0]).real, translation_residual=float(tr),
               least_stable=[float(np.real(lead)), float(np.imag(lead))], stable=bool(np.real(lead) < 0))
    np.savez_compressed(os.path.join(OUT, f"pair_eq_{xr - xl:.4f}.npz"), u=u)
    return rec


def zeros_of(ds, v, mask=None):
    z = []
    for i in range(len(ds) - 1):
        if mask is not None and not (mask[i] and mask[i + 1]):
            continue
        if v[i] == 0 or v[i] * v[i + 1] < 0:
            z.append(float(ds[i] - v[i] * (ds[i + 1] - ds[i]) / (v[i + 1] - v[i])))
    return z


def classify(ds, v, resolved):
    dr, vr = np.array(ds)[resolved], np.array(v)[resolved]
    if len(dr) < 10:
        return "UNRESOLVED", [], []
    z = zeros_of(np.array(ds), np.array(v), resolved)
    gaps = np.diff(z).tolist()
    if len(z) >= 3 and all(abs(gp - HALF_WL) <= 0.25 * HALF_WL for gp in gaps):
        return "OSCILLATORY", z, gaps
    if np.all(vr < 0):
        return "MONOTONIC ATTRACTION", z, gaps
    if np.all(vr > 0):
        return "MONOTONIC REPULSION", z, gaps
    return "UNRESOLVED", z, gaps


def fit_decay(ds, v, with_real=False):
    ds, v = np.asarray(ds), np.asarray(v)
    w = np.exp(KAPPA0 * ds)

    def model(p):
        kap, k, A, B = p[:4]
        f = np.exp(-kap * ds) * (A * np.cos(k * ds) + B * np.sin(k * ds))
        if with_real:
            f = f + p[4] * np.exp(-KAPPA_REAL * ds)
        return f

    A0 = np.max(np.abs(v * w))
    p0 = [KAPPA0, K0, A0, 0.0] + ([0.0] if with_real else [])
    best = None
    for k_init in (K0 * 0.95, K0, K0 * 1.05):
        for ph in range(4):
            p0[1] = k_init
            p0[2], p0[3] = A0 * np.cos(ph * np.pi / 2), A0 * np.sin(ph * np.pi / 2)
            r = least_squares(lambda p: (model(p) - v) * w, p0, x_scale="jac")
            if best is None or r.cost < best.cost:
                best = r
    p = best.x
    out = dict(kappa=float(p[0]), k=float(p[1]), A=float(p[2]), B=float(p[3]), cost=float(best.cost),
               rel_rms=float(np.sqrt(np.mean(((model(p) - v) * w) ** 2)) / np.sqrt(np.mean((v * w) ** 2))))
    if with_real:
        out["C"] = float(p[4])
        out["C_over_amp"] = float(p[4] / np.hypot(p[2], p[3]))
    return out


# ---------------- main ----------------

def main():
    os.makedirs(OUT, exist_ok=True)
    t_start = time.time()
    logf = open(os.path.join(OUT, "run_log.txt"), "w")

    def log(s):
        print(s, flush=True)
        logf.write(s + "\n")
        logf.flush()

    rec = dict(spec="spec/MG194_two_particle_gravity.md", b=BSTAR, a=14.0,
               platform=dict(platform=platform.platform(), python=sys.version.split()[0], numpy=np.__version__,
                             scipy=scipy.__version__), workers=WORKERS)

    def save():
        json.dump(rec, open(os.path.join(OUT, "results.json"), "w"), indent=1, default=float)

    log(f"MG194 start {time.strftime('%Y-%m-%d %H:%M:%S UTC', time.gmtime())}")
    g = m.Grid(100.0, 1024)
    U, res = soliton(g)
    np.savez_compressed(os.path.join(OUT, "soliton_bstar.npz"), U=U)
    lam0, others, psi, Up, J = m.spectrum_and_adjoint(g, U)
    mob = m.mobilities(g, U, psi)
    np.savez_compressed(os.path.join(OUT, "adjoint_bstar.npz"), psi=psi)
    rec["soliton"] = dict(residual=res, depth=float(U[1].min() - m.hom()[1]), translation_eig=float(np.real(lam0)),
                          least_stable_other=float(np.max(np.real(others))), mobilities=mob)
    log(f"soliton residual {res:.2e}; mobilities {mob}")
    g2 = m.Grid(200.0, 2048)
    U2, res2 = embed(g, U, g2)
    np.savez_compressed(os.path.join(OUT, "soliton_bstar_L200.npz"), U=U2)
    rec["soliton_L200_residual"] = res2
    save()

    # ---- Part A continuation recomputation (separate method, src/mg192_check.py logic) ----
    import mg192_check as chk
    chk.B = BSTAR
    zU = np.concatenate([U[s] for s in range(3)])
    chk.Zref = zU.copy()
    U0p = np.concatenate([chk.d1(zU[s * 1024:(s + 1) * 1024]) for s in range(3)])
    z0, r0, _ = chk.newton(np.concatenate([zU, [0.0]]), 0.0, "G", U0p)
    chk.Zref = z0[:3 * 1024].copy()
    U0p = np.concatenate([chk.d1(chk.Zref[s * 1024:(s + 1) * 1024]) for s in range(3)])
    cont = {}
    for kap in (1e-4, 2e-4):
        zp, rp, _ = chk.newton(z0.copy(), kap, "G", U0p)
        zm, rm, _ = chk.newton(z0.copy(), -kap, "G", U0p)
        cont[str(kap)] = dict(mu_a=float((zp[-1] - zm[-1]) / (2 * kap)), residuals=[rp, rm])
    rec["partA_continuation"] = cont
    log(f"continuation mu_a {cont}")
    save()

    # ---- Part A dynamics and Part B dynamics, in a pool ----
    L = 100.0
    a_jobs = [(0.002, L / 16), (0.002, -L / 16), (0.002, 3 * L / 16), (0.002, -3 * L / 16), (0.001, 3 * L / 16),
              (0.0, L / 16), (0.0, 3 * L / 16)]
    ds = [round(8.0 + 0.25 * i, 2) for i in range(89)]
    b_jobs = [(d, 100.0, 1024) for d in ds] + [(0.0, 100.0, 1024), (30.0, 200.0, 2048)]
    with Pool(WORKERS) as pool:
        ra = pool.map_async(part_a_run, a_jobs)
        rb = pool.map_async(pair_run, b_jobs)
        a_res, b_res = ra.get(), rb.get()
    series = {}
    for r in a_res:
        series[f"A_eps{r['eps']}_x0{r['x0']}"] = np.array(r.pop("series"))
        if r["eps"] > 0:
            r["predicted_v"] = mob["mu_a"] * r["dadx"]
            r["rel_err"] = (r["v"] - r["predicted_v"]) / r["predicted_v"]
        log(f"A eps={r['eps']} x0={r['x0']:+.3f}: v={r['v']:.5e} pred={r.get('predicted_v', 0):.5e} "
            f"rel={r.get('rel_err')} depth dev {r['depth_maxdev']:.4f}")
    for r in b_res:
        series[f"B_d{r['d']}_L{int(r['L'])}"] = np.array(r.pop("series"))
    rec["partA_runs"] = a_res
    rec["partB_runs"] = b_res
    np.savez_compressed(os.path.join(OUT, "trajectories.npz"), **series)
    save()

    # ---- gates A ----
    act = [r for r in a_res if r["eps"] > 0]
    nul = [r for r in a_res if r["eps"] == 0]
    A0 = (res <= 1e-10 and all(abs(c["mu_a"] - mob["mu_a"]) <= 1e-3 * abs(mob["mu_a"]) for c in cont.values())
          and all(abs(r["v"]) < 1e-12 for r in nul))
    A1 = all(abs(r["rel_err"]) <= 0.10 and r["depth_maxdev"] <= 0.05 for r in act) if A0 else None
    rec["gates"] = dict(A0=A0, A1=A1, partA="CONFIRMED" if (A0 and A1) else ("NOT CONFIRMED" if A0 else "A0 FAILED"))
    log(f"Gates A: {rec['gates']}")

    # ---- Part B floor, classification, fit ----
    main_runs = sorted([r for r in b_res if r["L"] == 100.0 and r["d"] > 0], key=lambda r: r["d"])
    null = next(r for r in b_res if r["d"] == 0.0)
    r200 = next(r for r in b_res if r["L"] == 200.0)
    r100_30 = next(r for r in main_runs if r["d"] == 30.0)
    floor = max(10 * abs(null["v_rel"]), 2 * abs(r200["v_rel"] - r100_30["v_rel"]))
    dd = np.array([r["d"] for r in main_runs])
    vv = np.array([r["v_rel"] for r in main_runs])
    resolved = np.array([(not r["merged"]) and abs(r["v_rel"]) > floor for r in main_runs])
    B0_pre = abs(null["v_rel"]) < 1e-12
    label, z, gaps = classify(dd, vv, resolved)
    fit = fit_decay(dd[resolved], vv[resolved]) if resolved.sum() >= 4 else None
    fit5 = fit_decay(dd[resolved], vv[resolved], with_real=True) if resolved.sum() >= 5 else None
    B2 = None if fit is None else bool(abs(fit["kappa"] - KAPPA0) <= 0.10 * KAPPA0 and abs(fit["k"] - K0) <= 0.05 * K0)
    rec["partB"] = dict(floor=floor, null_v=null["v_rel"], L200_vs_L100_at30=[r200["v_rel"], r100_30["v_rel"]],
                        n_resolved=int(resolved.sum()), n_merged=int(sum(r["merged"] for r in main_runs)),
                        n_moving=int(sum(r["moving"] for r in main_runs)), B1=label, zeros=z, zero_gaps=gaps,
                        B2_fit=fit, B2=("SUPPORTED" if B2 else "NOT SUPPORTED") if B2 is not None else "NOT SCORED",
                        B5_fit=fit5)
    log(f"Part B: floor {floor:.3e}, resolved {resolved.sum()}/89, B1 {label}, zeros {np.round(z, 3).tolist()}, "
        f"B2 {rec['partB']['B2']} fit {fit}")
    save()

    # ---- B-lin ----
    dl = np.round(np.arange(6.0, 40.0001, 0.05), 3)
    vlin = lin_curve(g, U, psi, dl)
    zl = [zz for zz in zeros_of(dl, vlin) if 8.0 <= zz <= 30.0]
    vlin_at = lin_curve(g, U, psi, dd)
    rel = np.abs(vlin_at[resolved] - vv[resolved]) / np.abs(vv[resolved])
    rec["partB_lin"] = dict(d=dl.tolist(), v_rel=vlin.tolist(), zeros=zl,
                            median_rel_err_resolved=float(np.median(rel)) if len(rel) else None)
    log(f"B-lin zeros {np.round(zl, 3).tolist()}; median rel err vs B-dyn {rec['partB_lin']['median_rel_err_resolved']}")
    save()

    # ---- B-eq ----
    with Pool(WORKERS) as pool:
        eq = pool.map(pair_equilibrium, zl)
    rec["partB_eq"] = eq
    for e in eq:
        log(f"B-eq from d0={e['d0']:.3f}: {e}")
    rec["partB"]["B0"] = bool(B0_pre and all(e.get("residual", 0) <= 1e-10 for e in eq if e["converged"]))
    rec["wall_seconds"] = time.time() - t_start
    save()
    log(f"done, wall {rec['wall_seconds']:.0f} s")


if __name__ == "__main__":
    main()

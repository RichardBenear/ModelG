"""MG8 Part 0 (spec/MG8_proposal.md draft 4 §5; registered at 3ab4ec8). Unregistered streams seed_for(8, 9, ...) only.
Writes results/MG8/part0/decisions.json (proceed flag) and part0_report.json. Usage: python src/part0_mg8.py"""
import os, sys, json, time, itertools, platform
import numpy as np
import numba as nb
from scipy import stats
sys.path.insert(0, os.path.dirname(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "tests"))
import mg8
from mg8 import J, R, L, P, TAU0, H
from runner import seed_for
from lattice import tau_leap_euler_1d, euler_mean_map

ROOT = os.path.join(os.path.dirname(__file__), "..")
OUT = os.path.join(ROOT, "results", "MG8", "part0")
CUTOFF = 20.0
NREP, N0 = 2000, 10 ** 6
RS = (0.01, 0.05, 0.2)


def orders(sched, N):
    if sched == R:
        return [[(s + k) % N for k in range(N)] for s in range(N)]
    if sched == P:
        return [list(o) for o in itertools.permutations(range(N))]
    return [list(range(N))]


# ------------------------------------------------------------------ 1(a) matrices
def bench_matrices():
    res = {}
    for r in RS:
        dR = max(np.abs(mg8.map_matrix(np.array(o), r) - mg8.chain_matrix(o, r)).max() for o in orders(R, 12))
        dP = max(np.abs(mg8.map_matrix(np.array(o), r) - mg8.chain_matrix(o, r)).max() for o in orders(P, 6))
        res[str(r)] = dict(R12=float(dR), P6=float(dP))
    ok = all(v["R12"] <= 1e-12 and v["P6"] <= 1e-12 for v in res.values())
    return dict(maxdiff=res, pass_=ok)


# ------------------------------------------------------------------ 1(b) stochastic benchmark (conditional construction)
def raw_moments(A, d, kmax=6):
    return [float((d ** k * A).sum()) for k in range(kmax + 1)]


def bench_stochastic():
    out = {}
    ok = True
    for sched, N in ((R, 12), (P, 6), (J, 12)):
        dmax = N
        d = np.arange(-dmax, dmax + 1)
        cond_cache = {}
        for r in RS:
            Abar = np.mean([mg8.displacement_dist(o, r, 0, dmax) if sched != J else jdist(r, N, dmax)
                            for o in orders(sched, N)], axis=0)
            ps, z = [], {1: [], 2: [], 3: []}
            summ = {1: [], 2: [], 3: []}
            cons = True
            for rep in range(NREP):
                cnt, order = mg8.bench_one_step(N0, N, r, sched, seed_for(8, 9, 2, 10 * sched + RS.index(r), rep), dmax)
                cons &= int(cnt.sum()) == N0
                ok_key = (r, tuple(order))
                if ok_key not in cond_cache:
                    cond_cache[ok_key] = jdist(r, N, dmax) if sched == J else mg8.displacement_dist(list(order), r, 0, dmax)
                A = cond_cache[ok_key]
                ps.append(pooled_chisq(cnt, N0 * A))
                m = raw_moments(A, d)
                for k in (1, 2, 3):
                    sk = float((d ** k * cnt).sum() / N0)
                    summ[k].append(sk)
                    se = np.sqrt(max(m[2 * k] - m[k] ** 2, 1e-300) / N0)
                    z[k].append((sk - m[k]) / se)
            mb = raw_moments(Abar, d)
            ks = float(stats.kstest(ps, "uniform").pvalue)
            zc = {k: dict(mean=float(np.mean(z[k])), sd=float(np.std(z[k], ddof=1))) for k in z}
            across = {}
            for k in (1, 2, 3):
                s = np.array(summ[k]); se = s.std(ddof=1) / np.sqrt(NREP)
                across[k] = dict(avg=float(s.mean()), target=mb[k], se=float(se),
                                 ok=bool(abs(s.mean() - mb[k]) <= 4 * se + 1e-15))
            mean0 = abs(np.mean(summ[1])) <= 4 * np.std(summ[1], ddof=1) / np.sqrt(NREP) + 1e-15
            c_ok = (ks >= 1e-3 and min(ps) >= 1e-6 and
                    all(abs(zc[k]["mean"]) <= 4 / np.sqrt(NREP) and 0.9 <= zc[k]["sd"] <= 1.1 for k in zc))
            a_ok = all(v["ok"] for v in across.values()) and bool(mean0)
            this = dict(ks_p=ks, min_p=float(min(ps)), z=zc, across=across, mean_zero=bool(mean0), conservation=bool(cons),
                        pass_=bool(c_ok and a_ok and cons))
            ok &= this["pass_"]
            out[f"{mg8.NAMES[sched]}_N{N}_r{r}"] = this
    return dict(cases=out, pass_=bool(ok))


def pooled_chisq(obs, exp, emin=5.0):
    """Chi-square p-value with tail pooling: cells are merged from each end inward until every pooled cell has
    expected count >= emin (attempt 2; attempt 1 pooled all low cells into one cell that could itself be < 5)."""
    obs = np.asarray(obs, float); exp = np.asarray(exp, float)
    nz = np.where(exp > 0)[0]
    lo, hi = nz[0], nz[-1]
    groups = []
    i = lo
    acc_o = acc_e = 0.0
    while i <= hi:                                  # left tail
        acc_o += obs[i]; acc_e += exp[i]; i += 1
        if acc_e >= emin:
            break
    groups.append((acc_o, acc_e))
    j = hi
    racc_o = racc_e = 0.0
    while j >= i:                                   # right tail
        racc_o += obs[j]; racc_e += exp[j]; j -= 1
        if racc_e >= emin:
            break
    for k in range(i, j + 1):
        groups.append((obs[k], exp[k]))
    if racc_e > 0:
        if racc_e >= emin:
            groups.append((racc_o, racc_e))
        else:
            o_, e_ = groups[-1]; groups[-1] = (o_ + racc_o, e_ + racc_e)
    outside = obs.sum() - sum(g[0] for g in groups)  # tokens in zero-probability cells
    if outside > 0:
        return 0.0
    o = np.array([g[0] for g in groups]); e = np.array([g[1] for g in groups])
    return float(stats.chisquare(o, e * o.sum() / e.sum()).pvalue)


def jdist(r, N, dmax):
    A = np.zeros(2 * dmax + 1)
    A[dmax - 1] = r; A[dmax] = 1 - 2 * r; A[dmax + 1] = r
    return A


# ------------------------------------------------------------------ 1(c) refinement, 1(d) formulas
def bench_refinement():
    res = {}
    ok = True
    for sched, N in ((R, 12), (P, 6)):
        D, Drand = [], []
        for tau in (0.2, 0.1, 0.05, 0.025, 0.0125):
            n = int(round(1 / tau))
            Tbar = np.mean([mg8.map_matrix(np.array(o), tau) for o in orders(sched, N)], axis=0)
            TJ = mg8.map_matrix(np.arange(N), tau, sync=True)
            TJn = np.linalg.matrix_power(TJ, n)
            D.append(float(np.linalg.norm(np.linalg.matrix_power(Tbar, n) - TJn)))
            if sched == R:
                rng = np.random.default_rng(seed_for(8, 9, 5, 0, 0))
                M = np.eye(N)
                for _ in range(n):
                    s = int(rng.integers(N))
                    M = mg8.map_matrix(np.array([(s + k) % N for k in range(N)]), tau) @ M
                Drand.append(float(np.linalg.norm(M - TJn)))
        rat = [D[i + 1] / D[i] for i in range(4)]
        p = all(x < 1 for x in rat) and all(0.4 <= x <= 0.6 for x in rat[-2:])
        ok &= p
        res[mg8.NAMES[sched]] = dict(delta=D, ratios=rat, random_seam_delta=Drand or None, pass_=bool(p))
    return dict(cases=res, pass_=bool(ok))


def bench_formulas():
    out = {}
    for r in RS:
        N = 12
        d = np.arange(-N, N + 1)
        Abar = np.mean([mg8.displacement_dist(o, r, 0, N) for o in orders(R, N)], axis=0)
        out[str(r)] = dict(var_avg=float((d ** 2 * Abar).sum()), var_formula=2 * r + 2 * r * r / (1 - r),
                           m3_avg=float((d ** 3 * Abar).sum()), m3_formula=6 * r * r / (1 - r) ** 2)
    return out


# ------------------------------------------------------------------ 2, 3, 5: regression, mirror, clipping, timing
def regression_mirror():
    import test_mg8 as T
    T.test_J_regression_stochastic(); T.test_det_J_equals_euler_mean_map(); T.test_mirror_stochastic_and_deterministic()
    return dict(pass_=True)


def clipping_and_timing():
    Pm = mg8.params()
    ref, c = mg8.reference_start()
    args = (Pm["a"], Pm["b"], Pm["g"], Pm["p"], Pm["q"], Pm["dx"], Pm["dy"])
    spu = mg8.steps_per_unit(TAU0)
    out = {}
    for sch in (J, R, P):
        n = [np.rint(1000.0 * c[i]).astype(np.int64) for i in range(3)]
        t0 = time.time()
        *_, clips = mg8.run_stoch(n[0], n[1], n[2], sch, 1000.0, *args, H, TAU0, 50 * spu, spu, seed_for(8, 9, 3, sch, 0))
        w = time.time() - t0
        out[f"stoch_{mg8.NAMES[sch]}"] = dict(clips=clips.tolist(), clip_frac=(clips / (500 * 50 * spu)).tolist(),
                                              core_s_per_unit=w / 50)
    for sch in (J, R):
        G, X, Y = (c[i].copy() for i in range(3))
        t0 = time.time()
        clips = mg8.run_det(G, X, Y, sch, *args, H, TAU0, 50 * spu, seed_for(8, 9, 3, 10 + sch, 0))
        out[f"det_{mg8.NAMES[sch]}_k0"] = dict(clips=clips.tolist(), core_s_per_unit=(time.time() - t0) / 50)
    tau6 = TAU0 / 64
    G, X, Y = (c[i].copy() for i in range(3))
    t0 = time.time()
    mg8.run_det(G, X, Y, R, *args, H, tau6, 2 * mg8.steps_per_unit(tau6), seed_for(8, 9, 3, 20, 0))
    det_tau0_eq = (time.time() - t0) / (2 * 64)               # core-s per tau0-equivalent time unit
    sto = max(out[f"stoch_{n}"]["core_s_per_unit"] for n in "JRP")
    det_units = 4 * 127 * 300 + 2 * 64 * 300 + 300 + 64 * 50
    proj = (25 * 500 * sto + det_units * det_tau0_eq) / 3600
    out.update(det_core_s_per_tau0_unit=det_tau0_eq, projected_core_h=proj)
    return out


def synthetic():
    import test_mg8 as T
    for n in ("test_shape_error_shift_and_deformation", "test_centre_velocity_with_wrap", "test_tracking_and_fate_synthetic",
              "test_labels_production", "test_labels_refinement", "test_b2_statistics_synthetic"):
        getattr(T, n)()
    return dict(pass_=True)


def brute_alignment():
    """§5.4(a) with the registered brute-force grid (step 1e-5 h over +-h around the true shift)."""
    ref, c = mg8.reference_start()
    phi, x = ref["phi"], ref["x"]
    from measure import fourier_shift
    out = []
    ok = True
    for amp in (1e-3, 1e-2):
        odd = x * np.exp(-x ** 2 / 4)
        A = fourier_shift(phi + amp * np.outer(np.abs(phi).max(axis=1), odd / np.abs(odd).max()), H, -0.05)
        E, sg, _ = mg8.shape_error(A, phi, H)
        An, Bn = mg8.nyquist_free(A), mg8.nyquist_free(phi)
        nB = np.sqrt((Bn ** 2).sum(axis=1))
        best = (np.inf, 0.0)
        for g in 0.05 + np.arange(-H, H, 1e-5 * H):
            e = float(np.max(np.sqrt(((fourier_shift(An, H, g) - Bn) ** 2).sum(axis=1)) / nB))
            if e < best[0]:
                best = (e, g)
        p = abs(E - best[0]) <= 0.01 * best[0] and abs(sg - best[1]) <= 1e-4 * H
        ok &= p
        out.append(dict(amp=amp, E=E, sigma=sg, E_bf=best[0], sigma_bf=best[1], pass_=bool(p)))
    return dict(cases=out, pass_=bool(ok))


# ------------------------------------------------------------------ supplement (review MG037 item 6)
def _first_order(sched, N, seed):
    """The order a production kernel draws first after seeding (the order draw comes first in every step)."""
    o = np.zeros(N, dtype=np.int64)
    _first_order_jit(sched, N, seed, o)
    return o


@nb.njit(cache=True)
def _first_order_jit(sched, N, seed, o):
    np.random.seed(seed)
    mg8._order(sched, N, o)


def production_kernel_crosscheck(nrep=500, n0=10 ** 4):
    """Diffusion-only cross-check of the PRODUCTION kernels against the conditional maps (G only, all reactions off,
    h = 1 so the hop probability is r = tau). (a) run_det over 20 steps vs the product of map_matrix for the realised
    orders (seams regenerated with mg8.seam_chunk): 1e-12. (b) run_stoch, one step from n0 tokens at node 0: given the
    order, the counts are Multinomial(n0, T_o e_0) exactly (tokens move independently); per replicate a pooled
    chi-square p-value (cells merged to expected >= 5); pass if the KS test of the p-values against U(0,1) has
    p > 1e-3, min p > 1e-6, and every node's aggregate z of (observed - expected) is within 4."""
    Z = (0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0)                 # a, b, g, p, q, dx, dy
    det, ok = [], True
    for sched in (J, R, L):
        for r in (0.05, 0.2):
            N, nst, sd = 12, 20, seed_for(8, 9, 8, sched, int(r * 100))
            rng = np.random.default_rng(sd)
            G0 = rng.random(N) + 0.5
            G, X, Y = G0.copy(), np.zeros(N), np.zeros(N)
            mg8.run_det(G, X, Y, sched, *Z, 1.0, r, nst, sd)
            seams = mg8.seam_chunk(sched, N, nst, sd)
            v = G0.copy()
            for s in seams:
                if sched == J:
                    v = mg8.map_matrix(np.arange(N), r, sync=True) @ v
                elif sched == R:
                    v = mg8.map_matrix(np.array([(s + k) % N for k in range(N)]), r) @ v
                else:
                    sp_ = N - 1 - s
                    v = mg8.map_matrix(np.array([(sp_ - k) % N for k in range(N)]), r) @ v
            d = float(np.max(np.abs(G - v)) / np.max(np.abs(v)))
            p_ = d <= 1e-12 and not X.any() and not Y.any()
            ok &= p_
            det.append(dict(sched=mg8.NAMES[sched], r=r, max_rel=d, pass_=bool(p_)))
    sto = []
    for sched, N in ((J, 12), (R, 12), (P, 6)):
        for r in (0.05, 0.2):
            ps, dev, var = [], np.zeros(N), np.zeros(N)
            for rep in range(nrep):
                sd = seed_for(8, 9, 9, 10 * sched + int(r * 100), rep)
                o = _first_order(sched, N, sd)
                pr = mg8.map_matrix(o, r, sync=(sched == J))[:, 0]
                nG = np.zeros(N, dtype=np.int64); nG[0] = n0
                nX = np.zeros(N, dtype=np.int64); nY = np.zeros(N, dtype=np.int64)
                rG, *_ = mg8.run_stoch(nG, nX, nY, sched, 1.0, *Z, 1.0, r, 1, 1, sd)
                obs = rG[0].astype(float)
                ex = n0 * pr
                cells = np.array([k % N for k in range(-(N // 2), N - N // 2)])   # ring order centred on node 0
                oc, ec = obs[cells], ex[cells]
                if np.any(oc[ec == 0] > 0):                  # tokens on a node the order cannot reach
                    ps.append(0.0)
                else:                                        # unreachable nodes (exp = obs = 0) are dropped
                    ps.append(pooled_chisq(oc[ec > 0], ec[ec > 0]))
                dev += obs - ex; var += n0 * pr * (1 - pr)
            z = dev / np.sqrt(np.where(var > 0, var, 1.0))
            ks = float(stats.kstest(ps, "uniform").pvalue)
            p_ = ks > 1e-3 and min(ps) > 1e-6 and float(np.max(np.abs(z))) <= 4.0
            ok &= p_
            sto.append(dict(sched=mg8.NAMES[sched], N=N, r=r, ks_p=ks, min_p=float(min(ps)),
                            max_abs_z=float(np.max(np.abs(z))), pass_=bool(p_)))
    return dict(deterministic=det, stochastic=sto, pass_=bool(ok))


def main_supplement():
    """Supplementary Part 0 check added after review MG037 (not a rerun of Part 0; decisions.json is unchanged)."""
    t0 = time.time()
    rep = dict(platform=platform.platform(), python=platform.python_version(), numpy=np.__version__,
               production_kernel_crosscheck=production_kernel_crosscheck())
    rep["wall_s"] = time.time() - t0
    json.dump(rep, open(os.path.join(OUT, "supplement_MG037.json"), "w"), indent=1, default=float)
    print(json.dumps(rep, indent=1, default=float))


def main():
    os.makedirs(OUT, exist_ok=True)
    t0 = time.time()
    rep = dict(platform=platform.platform(), python=platform.python_version(), numpy=np.__version__)
    rep["1a_matrices"] = bench_matrices(); print("1a", rep["1a_matrices"]["pass_"], flush=True)
    rep["1b_stochastic"] = bench_stochastic(); print("1b", rep["1b_stochastic"]["pass_"], flush=True)
    rep["1c_refinement"] = bench_refinement(); print("1c", rep["1c_refinement"]["pass_"], flush=True)
    rep["1d_formulas"] = bench_formulas()
    rep["2_regression_mirror"] = regression_mirror(); print("2 ok", flush=True)
    rep["4_synthetic"] = synthetic(); rep["4a_alignment_bruteforce"] = brute_alignment()
    print("4", rep["4a_alignment_bruteforce"]["pass_"], flush=True)
    rep["3_5_clipping_timing"] = clipping_and_timing()
    rep["wall_s"] = time.time() - t0
    passes = [rep[k]["pass_"] for k in ("1a_matrices", "1b_stochastic", "1c_refinement", "2_regression_mirror",
                                         "4_synthetic", "4a_alignment_bruteforce")]
    proj = rep["3_5_clipping_timing"]["projected_core_h"]
    dec = dict(all_checks_pass=all(passes), projected_core_h=proj, cutoff_core_h=CUTOFF,
               proceed=bool(all(passes) and proj <= CUTOFF), platform=rep["platform"])
    json.dump(rep, open(os.path.join(OUT, "part0_report.json"), "w"), indent=1, default=float)
    json.dump(dec, open(os.path.join(OUT, "decisions.json"), "w"), indent=1)
    print(json.dumps(dec, indent=1))


if __name__ == "__main__":
    main_supplement() if "--supplement" in sys.argv else main()

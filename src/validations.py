"""MG004: pre-gate validations of protocol v1.3 (frozen). Usage: python src/validations.py {v3,v1,v2,v4,v6}
Results in results/MG004/<v>.json and <v>.txt."""
import json, os, sys, time
import numpy as np
from multiprocessing import Pool
sys.path.insert(0, os.path.dirname(__file__))
from runner import one_run, run_cell, seed_for, reference
from lattice import tau_leap_euler_1d, steady_state, P17
from stats_mg2 import welch_diff_bounds, bias_bounds, decide, a1_diagnostics, _studentized_loo
from scipy import stats

OUT = "results/MG004"
os.makedirs(OUT, exist_ok=True)
H, TAU0 = 0.2, 0.4 * 0.2 * 0.2 / 12 / 2
G_V = 0                     # gate id for validations (protocol §7)
LOG = []


def log(s):
    print(s, flush=True); LOG.append(s)


def save(name, obj):
    json.dump(obj, open(f"{OUT}/{name}.json", "w"), indent=1, default=float)
    open(f"{OUT}/{name}.txt", "w").write("\n".join(LOG) + "\n")


def tau_table():
    p = f"{OUT}/v3.json"
    return {float(k): v for k, v in json.load(open(p))["tau_by_omega"].items()} if os.path.exists(p) else {}


# ------------------------------------------------------------------------------------------------ V3
def v3(pool):
    oms = sorted(set([10, 30, 100, 300, 1000, 3000, 10000] + [300, 100, 50, 30, 20, 15, 10, 7, 5, 3, 2, 1] + [3, 30] + [300, 3000]))
    log(f"V3 clip pilot: b=29, homogeneous start, T=20, one run per Omega; threshold clip fraction <= 1e-4 per node-step")
    table, rows = {}, []
    for k, om in enumerate(oms):
        tau = TAU0
        for halv in range(3):
            r = one_run(dict(b=29.0, Om=float(om), h=H, tau=tau, T=20.0, start="homog", kind="none",
                             rng=seed_for(G_V, 0, k, 3, 0, halv)))
            cf = r["clips"] / (500 * round(20.0 / tau))
            rows.append(dict(Omega=om, tau=tau, clip_frac_node_step=cf))
            log(f"  Omega={om:6g} tau={tau:.3e}: clips={r['clips']} fraction/node-step={cf:.2e}")
            if cf <= 1e-4:
                table[om] = tau
                break
            tau /= 2
        else:
            table[om] = None
            log(f"  Omega={om}: EXCLUDED (clip fraction > 1e-4 at tau/4)")
    save("v3", dict(rows=rows, tau_by_omega=table))


# ------------------------------------------------------------------------------------------------ V1
def _v1_job(args):
    from ssa import ssa_ring
    kind, Om, tau, seed = args
    P = {**P17, "b": 29.0}
    S0 = steady_state(P)
    N = 30
    rs = np.random.default_rng(seed ^ 0x1234)
    nG, nX, nY = (rs.poisson(Om * S0[k], N).astype(np.int64) for k in range(3))
    T, dtr = 500.0, 0.5
    t0 = time.time()
    if kind == "ssa":
        rec, nev = ssa_ring(nG, nX, nY, float(Om), P["a"], P["b"], P["g"], P["p"], P["q"], P["dx"], P["dy"], H, T, dtr, seed)
        rec = rec.astype(float)
    else:
        xs = np.zeros(N)
        rG, rX, rY, rt, clips = tau_leap_euler_1d(nG, nX, nY, float(Om), P["a"], P["b"], P["g"], P["p"], P["q"],
                                                  P["dx"], P["dy"], H, tau, int(round(T / tau)), 0.0, 0.0, xs,
                                                  int(round(dtr / tau)), seed)
        rec = np.stack([rG, rX, rY], 1).astype(float); nev = clips
    return kind, Om, tau, rec / Om, time.time() - t0, nev


def _v1_stats(c):
    """c: (T_snap, 3, N) concentrations after burn-in -> per-batch statistics (20 batches)."""
    out = []
    for blk in np.array_split(c, 20, axis=0):
        m = blk.mean((0, 2)); v = blk.var((0, 2))
        X = blk[:, 1, :]
        cov = np.mean((X - X.mean()) * (np.roll(X, 1, axis=1) - X.mean()))
        out.append(np.concatenate([m, v, [cov]]))
    out = np.array(out)
    return out.mean(0), out.std(0, ddof=1) / np.sqrt(out.shape[0])


def v1(pool):
    tt = tau_table()
    log("V1: exact SSA vs tau_leap_euler_1d; ring N=30, h=0.2, b=29, T=500, discard t<50, 20 batch means")
    names = ["mean G", "mean X", "mean Y", "var G", "var X", "var Y", "nn-cov X"]
    jobs = []
    for k, Om in enumerate((3.0, 30.0)):
        tau = tt.get(Om, TAU0) or TAU0
        jobs += [("ssa", Om, tau, seed_for(G_V, 0, k, 1, 0)), ("tau", Om, tau, seed_for(G_V, 0, k, 1, 1))]
    res = pool.map(_v1_job, jobs)
    allpass, rows = True, []
    for Om in (3.0, 30.0):
        r = {kind: (rec, wall, nev, tau) for kind, O, tau, rec, wall, nev in res if O == Om}
        s_m, s_se = _v1_stats(r["ssa"][0][100:]); t_m, t_se = _v1_stats(r["tau"][0][100:])
        log(f" Omega={Om:g} (tau={r['tau'][3]:.3e}; SSA events {r['ssa'][2]:.3g}, {r['ssa'][1]:.0f}s; tau clips {r['tau'][2]})")
        for i, nm in enumerate(names):
            tol = max(0.03 * abs(s_m[i]), 3 * np.hypot(s_se[i], t_se[i]))
            ok = abs(t_m[i] - s_m[i]) <= tol
            allpass &= ok
            rows.append(dict(Omega=Om, stat=nm, ssa=s_m[i], ssa_se=s_se[i], tau=t_m[i], tau_se=t_se[i], tol=tol, ok=bool(ok)))
            log(f"   {nm:9s}: SSA {s_m[i]:.5f}±{s_se[i]:.5f}  tau {t_m[i]:.5f}±{t_se[i]:.5f}  |diff| {abs(t_m[i]-s_m[i]):.5f} tol {tol:.5f} {'ok' if ok else 'FAIL'}")
    log(f"V1 {'PASS' if allpass else 'FAIL'}")
    save("v1", dict(rows=rows, verdict="PASS" if allpass else "FAIL"))


# ------------------------------------------------------------------------------------------------ helpers
def run_records(runs):
    """Audit record of every attempted run (no profiles)."""
    return [dict(replicate=r.get("replicate"), rng=r["rng"], eligible=bool(r.get("eligible")), events=r.get("events"),
                 n_per_block=r.get("n_per_block"), farvar=r.get("farvar"), clip_frac=r.get("clip_frac")) for r in runs]


def run_profiles(elig, ref):
    """eligible runs -> run-level profiles (M,3,n) and block profiles (M,K,3,n)."""
    B = np.array([r["blocks"] for r in elig])          # (M, K, 3, n)
    return B.mean(1), B


def window_ref(ref):
    win = np.abs(ref["x"]) < 10.0
    return ref["phi"][:, win]


# ------------------------------------------------------------------------------------------------ V2
def v2(pool):
    log("V2 (protocol v1.3.1, Amendment 1): tau-halving equivalence, b=29.08, Omega=300, 12 eligible runs each at tau and tau/2 (cap 24)")
    ref = reference(29.08, H); phis = window_ref(ref)
    tt = tau_table(); tau = tt.get(300.0, TAU0) or TAU0
    res = {}
    for j, (lab, tv) in enumerate((("tau", tau), ("tau/2", tau / 2))):
        base = dict(b=29.08, Om=300.0, h=H, tau=tv, T=100.0, start="soliton", kind="mg2")
        el, al = run_cell(pool, base, (G_V, 1, 0, 30 + j), 12, 24, log=log)
        res[lab] = (el, al)
        log(f"  {lab}: {len(el)} eligible of {len(al)} attempted")
    out = evaluate_v2(res, phis, ref)
    if out["verdict"] == "INCONCLUSIVE":
        log("  V2 INCONCLUSIVE -> one extension of 12 + 12 eligible runs")
        for j, (lab, tv) in enumerate((("tau", tau), ("tau/2", tau / 2))):
            base = dict(b=29.08, Om=300.0, h=H, tau=tv, T=100.0, start="soliton", kind="mg2")
            el, al = run_cell(pool, base, (G_V, 1, 0, 30 + j), 12, 24, start_rep=len(res[lab][1]), log=log)
            res[lab] = (res[lab][0] + el, res[lab][1] + al)
        out = evaluate_v2(res, phis, ref)
        if out["verdict"] == "INCONCLUSIVE":
            out["verdict"] = "FAIL (INCONCLUSIVE after extension counts as FAIL)"
    log(f"V2 {out['verdict']}")
    out["runs"] = {lab: run_records(al) for lab, (el, al) in res.items()}
    save("v2", out)


def evaluate_v2(res, phis, ref):
    (eA, _), (eB, _) = res["tau"], res["tau/2"]
    if len(eA) < 12 or len(eB) < 12:
        return dict(verdict="INCONCLUSIVE", reason="insufficient eligible runs")
    RA, _ = run_profiles(eA, ref); RB, _ = run_profiles(eB, ref)
    verdict, rows = "PASS", []
    for s, nm in enumerate("GXY"):
        nrm = float(np.linalg.norm(phis[s]))
        b = welch_diff_bounds(RA[:, s], RB[:, s], nrm)
        prof_ok = b["upper"] < 0.0025
        prof_fail = b["diff"] > 0.0025
        vA = np.array([r["farvar"][s] for r in eA]); vB = np.array([r["farvar"][s] for r in eB])
        d = vB.mean() - vA.mean()
        se = np.sqrt(vA.var(ddof=1) / len(vA) + vB.var(ddof=1) / len(vB))
        dfw = se**4 / ((vA.var(ddof=1) / len(vA)) ** 2 / (len(vA) - 1) + (vB.var(ddof=1) / len(vB)) ** 2 / (len(vB) - 1))
        tq = stats.t.ppf(0.95, dfw)
        lo, hi = d - tq * se, d + tq * se
        m = 0.05 * vA.mean()
        var_ok = (lo > -m) and (hi < m)
        var_fail = abs(d) > m
        rows.append(dict(field=nm, prof_diff=b["diff"], prof_upper=b["upper"], var_tau=vA.mean(), var_tau2=vB.mean(),
                         var_diff=d, var_ci90=[lo, hi], var_margin=m))
        log(f"   {nm}: profile |mu_tau - mu_tau/2|/|phi*| = {100*b['diff']:.3f}% (upper {100*b['upper']:.3f}%, margin 0.25%); "
            f"far-field Omega*Var {vA.mean():.3f} vs {vB.mean():.3f}, diff 90% CI [{lo:+.3f}, {hi:+.3f}] (margin ±{m:.3f})")
        if prof_fail or var_fail:
            verdict = "FAIL"
        elif not (prof_ok and var_ok) and verdict != "FAIL":
            verdict = "INCONCLUSIVE"
    return dict(verdict=verdict, rows=rows, n=[len(eA), len(eB)])


# ------------------------------------------------------------------------------------------------ V4
def v4(pool):
    from measure import analyze_run, fourier_shift
    log("V4 (protocol v1.3.1): alignment validation, b=29; synthetic phi* + homogeneous-run noise, random circular shifts; 12 runs per Omega")
    ref = reference(29.0, H); phis = window_ref(ref); S0 = ref["S0"]
    tt = tau_table()
    out = {}
    for k, Om in enumerate((100.0, 1000.0, 10000.0)):
        tau = tt.get(Om, TAU0) or TAU0
        jobs = [dict(b=29.0, Om=Om, h=H, tau=tau, T=100.0, start="homog", kind="homog", keep_fields=True,
                     rng=seed_for(G_V, 0, k, 4, rep)) for rep in range(12)]
        runs = pool.map(one_run, jobs)
        rs = np.random.default_rng(seed_for(G_V, 0, k, 40, 0))
        L = 100.0
        elig = []
        for r in runs:
            s = rs.uniform(-L / 2, L / 2)
            sol = fourier_shift(ref["phi"], H, -s)            # soliton centred at x = s
            phi = r["phi"].astype(float) + sol[None, :, :]
            n = Om * (S0[None, :, None] + phi)
            a = analyze_run(n[:, 0], n[:, 1], n[:, 2], r["times"], Om, H, S0, ref["A_s"], ref["lam"], x0=s)
            if a["eligible"]:
                elig.append(a)
        R, Bk = run_profiles(elig, ref) if elig else (None, None)
        cell = dict(n_eligible=len(elig))
        for s_, nm in enumerate("GXY"):
            if len(elig) >= 2:
                bb = bias_bounds(R[:, s_], phis[s_])
                diag = a1_diagnostics(Bk[:, :, s_], phis[s_])
                cell[nm] = dict(Eb_point=bb["Eb_point"], E_total=bb["E_total"], lower=bb["lower"], upper=bb["upper"],
                                decision=decide(bb, 0.003), a1_ok=diag["ok"])
                log(f"  Omega={Om:g} {nm}: alignment-induced bias point {100*bb['Eb_point']:.4f}%  bounds [{100*bb['lower']:.4f}, {100*bb['upper']:.4f}]%  "
                    f"(eps 0.3%: {decide(bb, 0.003)}; A1 diag {'ok' if diag['ok'] else 'flag'}); eligible {len(elig)}/12")
        out[str(Om)] = cell
    c = out["10000.0"]
    if c["n_eligible"] < 12:
        v = "INCONCLUSIVE (blocking): insufficient eligible synthetic runs at Omega=1e4"
    else:
        ds = [c[nm]["decision"] for nm in "GXY"]
        v = "PASS" if all(d == "PASS" for d in ds) else ("FAIL" if any(d == "FAIL" for d in ds) else "INCONCLUSIVE (blocking)")
    out["verdict"] = v
    log(f"V4 {v}")
    save("v4", out)


# ------------------------------------------------------------------------------------------------ V6
def v6(pool):
    log("V6 (protocol v1.3.1): A1 pilot on real survival-conditioned runs, b=29.08, Omega in {300, 3000}, 12 eligible each (cap 24)")
    ref = reference(29.08, H); phis = window_ref(ref)
    tt = tau_table()
    cells, out = {}, {}
    for k, Om in enumerate((300.0, 3000.0)):
        tau = tt.get(Om, TAU0) or TAU0
        base = dict(b=29.08, Om=Om, h=H, tau=tau, T=100.0, start="soliton", kind="mg2")
        el, al = run_cell(pool, base, (G_V, 1, k, 6), 12, 24, log=log)
        log(f"  Omega={Om:g}: {len(el)} eligible of {len(al)} attempted; events in ineligible runs: "
            f"{[r['events'] for r in al if not r.get('eligible')]}")
        cells[Om] = el
        out.setdefault("runs", {})[str(Om)] = run_records(al)
    ok_all = True
    pooled = {nm: [] for nm in "GXY"}
    for Om, el in cells.items():
        if len(el) < 12:
            ok_all = False
            log(f"  Omega={Om:g}: insufficient eligible runs -> V6 FAIL")
            continue
        R, Bk = run_profiles(el, ref)
        out[str(Om)] = {}
        for s, nm in enumerate("GXY"):
            blocks = Bk[:, :, s]
            dg = a1_diagnostics(blocks, phis[s])
            ok_all &= dg["ok"]
            M, K, d = blocks.shape
            runs = blocks.mean(1); mbar = runs.mean(0)
            W = (blocks - runs[:, None, :]).reshape(M * K, d)
            _, _, Vt = np.linalg.svd(W - W.mean(0), full_matrices=False)
            bdir = (mbar - phis[s]) / np.linalg.norm(mbar - phis[s])
            ratios = []
            for v in (bdir, Vt[0], Vt[1]):
                pr = (runs - mbar) @ v
                pb = W @ v
                ratios.append(float(pr.var(ddof=1) / (pb.var(ddof=1) * (K / (K - 1)) / K)))
                z = (pr - pr.mean()) / pr.std(ddof=1)
                pooled[nm].append(z)
            out[str(Om)][nm] = dict(a1_ok=dg["ok"], details=dg["details"], var_ratio_run_vs_block=ratios,
                                    block_dev_quantiles=np.quantile(W @ Vt[0], [0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99]).tolist())
            log(f"  Omega={Om:g} {nm}: A1 diagnostics {'pass' if dg['ok'] else 'FLAG'} (Shapiro p min {min(dg['details']['shapiro_p']):.3g}, "
                f"max block |z| {max(dg['details']['max_block_z']):.2f}, max run t {max(dg['details']['max_run_t']):.2f} vs crit {dg['details']['run_t_crit']:.2f}); "
                f"run/block variance ratios {np.round(ratios, 2).tolist()}")
    # pooled audit: per field, per direction index, 24 standardized values; Bonferroni over all runs and directions
    audit_ok = True
    for nm in "GXY":
        Z = pooled[nm]
        if not Z:
            continue
        ndir = 3
        groups = [np.concatenate([Z[c * ndir + j] for c in range(len(Z) // ndir)]) for j in range(ndir)]
        for j, zz in enumerate(groups):
            tcrit = stats.t.ppf(1 - 0.01 / (2 * zz.size * ndir * 3), zz.size - 2)
            mx = float(np.max(np.abs(_studentized_loo(zz))))
            audit_ok &= mx <= tcrit
            log(f"  pooled audit {nm} dir{j}: max |t| {mx:.2f} vs crit {tcrit:.2f} ({'ok' if mx <= tcrit else 'FLAG'})")
    v = "PASS" if (ok_all and audit_ok) else "FAIL"
    out["verdict"] = v
    log(f"V6 {v}")
    save("v6", out)


# ================================================================================ protocol v1.3.2 (Amendment 2)
def evaluate_v2b(res, phis):
    """A2 rules. Profiles: PASS upper<0.25%, FAIL lower>0.25%, else INCONCLUSIVE.
    Variances (signed, 90% CI, margin m=5% of tau value): PASS -m<lo & hi<m; FAIL lo>m or hi<-m; else INC."""
    (eA, _), (eB, _) = res["tau"], res["tau/2"]
    if len(eA) < 12 or len(eB) < 12:
        return dict(verdict="INCONCLUSIVE", reason="insufficient eligible runs", n=[len(eA), len(eB)])
    RA = np.array([r["blocks"] for r in eA]).mean(1); RB = np.array([r["blocks"] for r in eB]).mean(1)
    dec, rows = [], []
    for s, nm in enumerate("GXY"):
        b = welch_diff_bounds(RA[:, s], RB[:, s], float(np.linalg.norm(phis[s])))
        pd = "PASS" if b["upper"] < 0.0025 else ("FAIL" if b["lower"] > 0.0025 else "INCONCLUSIVE")
        vA = np.array([r["farvar"][s] for r in eA]); vB = np.array([r["farvar"][s] for r in eB])
        d = vB.mean() - vA.mean()
        sa, sb = vA.var(ddof=1) / len(vA), vB.var(ddof=1) / len(vB)
        se = np.sqrt(sa + sb); dfw = (sa + sb) ** 2 / (sa**2 / (len(vA) - 1) + sb**2 / (len(vB) - 1))
        tq = stats.t.ppf(0.95, dfw); lo, hi = d - tq * se, d + tq * se; m = 0.05 * vA.mean()
        vd = "PASS" if (lo > -m and hi < m) else ("FAIL" if (lo > m or hi < -m) else "INCONCLUSIVE")
        dec += [pd, vd]
        rows.append(dict(field=nm, prof_diff=b["diff"], prof_lower=b["lower"], prof_upper=b["upper"], prof_decision=pd,
                         var_tau=vA.mean(), var_tau2=vB.mean(), var_diff=d, var_ci90=[lo, hi], var_margin=m, var_decision=vd))
        log(f"   {nm}: profile diff {100*b['diff']:.4f}% bounds [{100*b['lower']:.4f}, {100*b['upper']:.4f}]% (margin 0.25%) -> {pd}; "
            f"far-field Omega*Var {vA.mean():.3f} vs {vB.mean():.3f}, diff 90% CI [{lo:+.3f}, {hi:+.3f}] (margin ±{m:.3f}) -> {vd}")
    v = "FAIL" if "FAIL" in dec else ("PASS" if all(x == "PASS" for x in dec) else "INCONCLUSIVE")
    return dict(verdict=v, rows=rows, n=[len(eA), len(eB)])


def v2b(pool):
    LOG.clear()
    log("V2 (protocol v1.3.2, Amendment 2): tau-halving equivalence, b=29.08, Omega=1e8, 12 eligible each at tau and tau/2 (cap 24, one extension)")
    ref = reference(29.08, H); phis = window_ref(ref)
    Om = 1e8; tau = TAU0
    res = {}
    for j, (lab, tv) in enumerate((("tau", tau), ("tau/2", tau / 2))):
        base = dict(b=29.08, Om=Om, h=H, tau=tv, T=100.0, start="soliton", kind="mg2")
        el, al = run_cell(pool, base, (G_V, 1, 8, 50 + j), 12, 24, log=log)
        res[lab] = (el, al); log(f"  {lab}: {len(el)} eligible of {len(al)} attempted")
    out = evaluate_v2b(res, phis)
    if out["verdict"] == "INCONCLUSIVE":
        log("  V2 INCONCLUSIVE -> one extension of 12 + 12 eligible runs")
        for j, (lab, tv) in enumerate((("tau", tau), ("tau/2", tau / 2))):
            base = dict(b=29.08, Om=Om, h=H, tau=tv, T=100.0, start="soliton", kind="mg2")
            el, al = run_cell(pool, base, (G_V, 1, 8, 50 + j), 12, 24, start_rep=len(res[lab][1]), log=log)
            res[lab] = (res[lab][0] + el, res[lab][1] + al)
        out = evaluate_v2b(res, phis)
        if out["verdict"] == "INCONCLUSIVE":
            out["verdict"] = "FAIL (INCONCLUSIVE after extension counts as FAIL)"
    log(f"V2 {out['verdict']}")
    out["runs"] = {lab: run_records(al) for lab, (el, al) in res.items()}
    out["protocol"] = "v1.3.2"
    save("v2b", out)


def v4b(pool):
    from measure import analyze_run, fourier_shift
    LOG.clear()
    log("V4 (protocol v1.3.2, Amendment 2): alignment validation, b=29; synthetic phi* + homogeneous-run noise; 12 runs per Omega in {1e4, 1e6, 1e8}; decided at 1e8")
    ref = reference(29.0, H); phis = window_ref(ref); S0 = ref["S0"]
    out = {}
    for k, Om in enumerate((1e4, 1e6, 1e8)):
        jobs = [dict(b=29.0, Om=Om, h=H, tau=TAU0, T=100.0, start="homog", kind="homog", keep_fields=True,
                     rng=seed_for(G_V, 0, 10 + k, 42, rep)) for rep in range(12)]
        runs = pool.map(one_run, jobs)
        rs = np.random.default_rng(seed_for(G_V, 0, 10 + k, 43, 0))
        elig = []
        for r in runs:
            s = rs.uniform(-50.0, 50.0)
            phi = r["phi"].astype(float) + fourier_shift(ref["phi"], H, -s)[None, :, :]
            n = Om * (S0[None, :, None] + phi)
            a = analyze_run(n[:, 0], n[:, 1], n[:, 2], r["times"], Om, H, S0, ref["A_s"], ref["lam"], x0=s)
            if a["eligible"]:
                elig.append(a)
        cell = dict(n_eligible=len(elig))
        if len(elig) >= 2:
            Bk = np.array([a["blocks"] for a in elig]); R = Bk.mean(1)
            for s_, nm in enumerate("GXY"):
                bb = bias_bounds(R[:, s_], phis[s_]); diag = a1_diagnostics(Bk[:, :, s_], phis[s_])
                cell[nm] = dict(Eb_point=bb["Eb_point"], E_total=bb["E_total"], lower=bb["lower"], upper=bb["upper"],
                                decision=decide(bb, 0.003), a1_ok=diag["ok"])
                log(f"  Omega={Om:g} {nm}: alignment-induced bias point {100*bb['Eb_point']:.5f}% bounds [{100*bb['lower']:.5f}, {100*bb['upper']:.5f}]% "
                    f"(eps 0.3%: {decide(bb, 0.003)}; A1 {'ok' if diag['ok'] else 'flag'}); eligible {len(elig)}/12")
        else:
            log(f"  Omega={Om:g}: eligible {len(elig)}/12 -- insufficient")
        out[f"{Om:g}"] = cell
    c = out["1e+08"]
    if c["n_eligible"] < 12:
        v = "INCONCLUSIVE (blocking): insufficient eligible synthetic runs at Omega=1e8"
    else:
        ds = [c[nm]["decision"] for nm in "GXY"]; a1 = all(c[nm]["a1_ok"] for nm in "GXY")
        v = ("INCONCLUSIVE (blocking): A1 not supported" if not a1 else
             "PASS" if all(d == "PASS" for d in ds) else ("FAIL" if any(d == "FAIL" for d in ds) else "INCONCLUSIVE (blocking)"))
    out["verdict"] = v; out["protocol"] = "v1.3.2"
    log(f"V4 {v}")
    save("v4b", out)


def _v7_job(args):
    """One homogeneous run at b=29; detections, tracks, and deterministic confirmation of persistent ones."""
    from measure import smooth_periodic, detect, track, circ
    from lattice import lattice_ode_1d
    r = one_run(dict(args, kind="homog", keep_fields=True))
    ref = reference(29.0, H)
    phi, times = r["phi"].astype(float), r["times"]
    x = ref["x"]; L = 100.0
    ps = smooth_periodic(phi[:, 1, :], H)
    dets = [detect(ps[s], x, H, ref["thr"], ref["lam"]) for s in range(len(times))]
    tr = track(dets, times, L, x0=None)
    durs = [float(q["t_last"] - q["t_start"]) for q in tr["tracks"]]
    pers = [q for q in tr["tracks"] if q["t_last"] - q["t_start"] >= 20.0]
    cls = []
    for q in pers[:10]:
        tc = q["t_start"] + 20.0
        s = int(np.argmin(np.abs(times - tc)))
        y0 = np.concatenate([phi[s, 0], phi[s, 1], phi[s, 2]])
        xx, pG, pX, pY, mdt = lattice_ode_1d({**P17, "b": 29.0}, 500, H, 300.0, seed_scale=0.0, y0=y0)
        if pX.min() < -0.5 * abs(ref["A"]) and mdt < 1e-6:
            c = "confirmed"
        elif np.max(np.abs(np.concatenate([pG, pX, pY]))) < 1e-3:
            c = "excursion"
        else:
            c = "unresolved"
        cls.append(dict(t_start=float(q["t_start"]), duration=float(q["t_last"] - q["t_start"]), cls=c, ode_min_phiX=float(pX.min()), ode_maxdt=float(mdt)))
    return dict(Om=args["Om"], rng=args["rng"], n_dets=[len(d) for d in dets], track_durations=durs,
                n_persistent=len(pers), classified=cls, clip_frac=r["clip_frac"])


def v7(pool):
    LOG.clear()
    log("V7 (protocol v1.3.2, descriptive): detector characterization, b=29, homogeneous start, MG3 Omega grid, 4 runs, T=100; "
        "persistent detections classified by deterministic confirmation (max 10 per run). PI checkpoint before MG3.")
    tt = tau_table()
    grid = [1, 2, 3, 5, 7, 10, 15, 20, 30, 50, 100, 300]
    jobs = [dict(b=29.0, Om=float(om), h=H, tau=(tt.get(float(om)) or TAU0), T=100.0, start="homog",
                 rng=seed_for(G_V, 0, k, 7, rep)) for k, om in enumerate(grid) for rep in range(4)]
    res = pool.map(_v7_job, jobs)
    out = {"protocol": "v1.3.2", "cells": {}}
    for om in grid:
        rr = [r for r in res if r["Om"] == float(om)]
        nd = np.concatenate([r["n_dets"] for r in rr])
        durs = np.concatenate([r["track_durations"] for r in rr]) if any(len(r["track_durations"]) for r in rr) else np.zeros(0)
        cls = [c for r in rr for c in r["classified"]]
        cnt = {k: sum(c["cls"] == k for c in cls) for k in ("confirmed", "excursion", "unresolved")}
        npers = sum(r["n_persistent"] for r in rr)
        out["cells"][str(om)] = dict(mean_dets_per_snapshot=float(nd.mean()), frac_snapshots_with_dets=float((nd > 0).mean()),
                                     n_tracks=int(durs.size), track_duration_quantiles=(np.quantile(durs, [0.5, 0.9, 0.99]).tolist() if durs.size else []),
                                     n_persistent=int(npers), classified=cnt, n_unclassified=int(npers - len(cls)), runs=rr)
        log(f"  Omega={om:4d}: dets/snapshot {nd.mean():.3f} (snapshots with any {100*(nd>0).mean():.1f}%); tracks {durs.size}, "
            f"duration median/90%/99% {np.round(np.quantile(durs,[.5,.9,.99]),1).tolist() if durs.size else '-'}; persistent>=20: {npers} "
            f"-> confirmed {cnt['confirmed']}, excursion {cnt['excursion']}, unresolved {cnt['unresolved']}, unclassified {npers-len(cls)}")
    log("V7 complete (descriptive). Awaiting PI checkpoint before MG3.")
    save("v7", out)


if __name__ == "__main__":
    which = sys.argv[1]
    with Pool(int(os.environ.get("MODELG_WORKERS", "2"))) as pool:
        dict(v3=v3, v1=v1, v2=v2, v4=v4, v6=v6, v2b=v2b, v4b=v4b, v7=v7)[which](pool)

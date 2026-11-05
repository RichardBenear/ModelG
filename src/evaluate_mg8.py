"""MG8 evaluator (spec/MG8_proposal.md draft 4, §3.3 and §4). Reads results/MG8/runs/*.json and results/MG8/fields/*.npz,
writes results/MG8/evaluation.json and results/MG8/MG8_evaluation.md. Usage: python src/evaluate_mg8.py [outdir]"""
import os, sys, json
import numpy as np
from scipy import stats
sys.path.insert(0, os.path.dirname(__file__))
import mg8
from mg8 import TAU0, H

ROOT = os.path.join(os.path.dirname(__file__), "..")
BOOT_N = 2000


def load(out):
    runs = {}
    rd = os.path.join(out, "runs")
    if not os.path.isdir(rd):
        return runs
    for f in sorted(os.listdir(rd)):
        if f.endswith(".json"):
            runs[f[:-5]] = json.load(open(os.path.join(rd, f)))
    return runs


def fields(out, key):
    z = np.load(os.path.join(out, "fields", key + ".npz"))
    return {int(k[1:]): z[k] for k in z.files}


def units(rec):
    """Nominal unit index of every record entry. Recorded physical times are u * 1.00005 (spec/MG8_prerun.md §1), so
    windows and checkpoints are selected by the unit index, never by a physical-time bound (review MG037 item 1)."""
    if "unit" in rec:
        return np.asarray(rec["unit"], int)
    return np.rint(np.asarray(rec["t"], float)).astype(int)


def finite_units(rec):
    return np.asarray(rec.get("finite_units", [bool(rec.get("finite", True))] * len(rec["t"])), bool)


def window_v(rec, u0, u1):
    """Least-squares centre velocity over the units u0..u1 inclusive (slope against physical time)."""
    t = np.array(rec["t"], float); c = mg8.unwrap(rec["centre"], 100.0); u = units(rec)
    w = (u >= u0) & (u <= u1)
    if w.sum() < u1 - u0 + 1 or not np.all(finite_units(rec)[u <= u1]):
        return None
    return mg8.slope(t[w], c[w])


def single_through(rec, U):
    """SINGLE (§3.3): exactly one detection at every unit 50..U inclusive, every unit present, and the fields finite at
    every unit up to U. A failure after U does not affect the assessment at U."""
    u, n, fin = units(rec), np.asarray(rec["n_dets"]), finite_units(rec)
    have = set(u[u <= U].tolist())
    if any(k not in have for k in range(1, U + 1)):
        return False
    w = (u >= 50) & (u <= U)
    return bool(np.all(n[w] == 1) and np.all(fin[u <= U]))


def det_stats(out, runs, k, rep, T, Jf, Efloor):
    key = f"det_R_k{k}_r{rep}"
    rec = runs[key]
    Rf = fields(out, key)
    single = single_through(rec, T)
    E = E1 = v = v1 = None
    rep_ = {}
    if all(u in Rf and np.all(np.isfinite(Rf[u])) for u in (T, T - 50)):
        rep_ = mg8.shape_report(Rf[T], Jf[T], H)
        E = rep_["E"]
        E1 = mg8.shape_error(Rf[T - 50], Jf[T - 50], H)[0]
        v = window_v(rec, T - 150, T); v1 = window_v(rec, T - 200, T - 50)
    st = E is not None and v is not None and v1 is not None and mg8.settled(E, E1, v, v1)
    return dict(single=single, E=E, E1=E1, v=v, v1=v1, settled=bool(st), sigma=rep_.get("sigma"),
                E_per_species=rep_.get("per"), E_raw=rep_.get("E_raw"), E_nyq=rep_.get("E_nyq"))


# ------------------------------------------------------------------ implementation checks (§4.1; review MG037 item 2)
def mirror_det(out, runs, k):
    kL, kR = f"det_L_k{k}_r0", f"det_R_k{k}_r0"
    if kL not in runs or kR not in runs:
        return dict(status="INCOMPLETE")
    Lf, Rf = fields(out, kL), fields(out, kR)
    TL = runs[kL]["job"]["T"]
    want = list(range(50, TL + 1, 50))
    if any(u not in Lf or u not in Rf for u in want):
        return dict(status="FAIL", reason="missing checkpoint", units=want)
    d = max(float(np.max(np.abs(Lf[u][:, ::-1] - Rf[u])) / max(1e-300, float(np.max(np.abs(Rf[u]))))) for u in want)
    return dict(status="PASS" if d <= 1e-12 else "FAIL", max_rel=d, units=want)


def mirror_stoch(runs):
    kL, kR = "sto_L_k0_r0", "sto_R_k0_r0"
    if kL not in runs or kR not in runs:
        return dict(status="INCOMPLETE")
    ok = runs[kL]["counts_sha256"] == runs[kR]["counts_sha256"]
    return dict(status="PASS" if ok else "FAIL", hash_L=runs[kL]["counts_sha256"], hash_R=runs[kR]["counts_sha256"])


def impl_checks(out, runs):
    """Mandatory implementation checks: PASS, FAIL, or INCOMPLETE (records not there yet)."""
    ch = {}
    dec = os.path.join(out, "part0", "decisions.json")
    if os.path.exists(dec):
        d = json.load(open(dec))
        ch["part0"] = dict(status="PASS" if d.get("proceed") and d.get("all_checks_pass") else "FAIL")
    else:
        ch["part0"] = dict(status="INCOMPLETE")
    ch["mirror_det_k0"] = mirror_det(out, runs, 0)
    ch["mirror_det_k6"] = mirror_det(out, runs, 6)
    ch["mirror_stoch"] = mirror_stoch(runs)
    return ch


def run_status(out):
    p = os.path.join(out, "status.json")
    return json.load(open(p)) if os.path.exists(p) else dict(state="none")


WITHHELD = "WITHHELD: implementation check failed (STOP, PI decides)"


def b1_keys(k=None):
    import run_mg8
    return [run_mg8.key(j) for j in run_mg8.jobs() if j["kind"] == "det" and (k is None or j["k"] == k)]


def evaluate_b1(out, runs, checks=None, state=None):
    checks = checks if checks is not None else impl_checks(out, runs)
    state = state if state is not None else run_status(out)["state"]
    res = dict(levels={}, checks=checks, state=state)
    for k in range(7):
        Jrec = runs.get(f"det_J_k{k}_r0")
        if Jrec is None or any(f"det_R_k{k}_r{j}" not in runs for j in range(3)):
            res["levels"][k] = None
            continue
        Jf = fields(out, f"det_J_k{k}_r0")
        T = 300
        Ef = mg8.shape_error(Jf[T - 50], Jf[T], H)[0] + 1e-9
        vJ = window_v(Jrec, T - 150, T)
        vf = abs(vJ if vJ is not None else float("nan")) + 1e-8
        reps = [det_stats(out, runs, k, j, T, Jf, Ef) for j in range(3)]
        Es = [r["E"] for r in reps]
        lev = dict(Efloor=Ef, vJ=vJ, vfloor=vf, reps=reps,
                   all_single=all(r["single"] for r in reps), all_settled=all(r["settled"] for r in reps))
        if all(e is not None for e in Es):
            Ek = float(np.median(Es))
            lev.update(Ek=Ek, S=(max(Es) - min(Es)) / max(Ek, 10 * Ef))
        else:
            lev.update(Ek=float("nan"), S=float("inf"))
        lev["labels"] = [mg8.label_production(r["single"], r["E"] if r["E"] is not None else float("nan"), Ef,
                                              r["v"] if r["v"] is not None else float("nan"), vf, r["settled"])
                         for r in reps]
        lev["summary"] = mg8.combine_replicates(lev["labels"])
        res["levels"][k] = lev
    # horizon (k = 6, J and R rep 0 to t = 600)
    try:
        Jf = fields(out, "det_J_k6_r0"); Rf = fields(out, "det_R_k6_r0"); rec = runs["det_R_k6_r0"]
        hz = mg8.shape_report(Rf[600], Jf[600], H)
        E550 = mg8.shape_error(Rf[550], Jf[550], H)[0]
        v, v1 = window_v(rec, 450, 600), window_v(rec, 400, 550)
        res["horizon"] = dict(single=single_through(rec, 600),
                              settled=bool(v is not None and v1 is not None and mg8.settled(hz["E"], E550, v, v1)),
                              E300=res["levels"][6]["reps"][0]["E"], E600=hz["E"], E600_raw=hz["E_raw"],
                              E600_nyq=hz["E_nyq"], Efloor=res["levels"][6]["Efloor"])
    except (KeyError, TypeError):
        res["horizon"] = None
    res["mirror"] = {k: checks[f"mirror_det_k{k}"] for k in (0, 6)}
    # final labels only when every applicable implementation check passed (§4.1) and the run state allows it (§6)
    lv = res["levels"]
    failed = [n for n, c in checks.items() if c["status"] == "FAIL"]
    k0_done = all(kk in runs for kk in b1_keys(0))
    b1_done = all(kk in runs for kk in b1_keys())
    final = state in ("complete", "BUDGET-SHORT")
    if failed:
        res["production"] = res["refinement"] = WITHHELD
        res["failed_checks"] = failed
        return res
    if not final:
        res["production"] = res["refinement"] = f"NOT FINAL (production state: {state})"
    else:
        if not k0_done:
            res["production"] = "NOT SCORED (BUDGET-SHORT)" if state == "BUDGET-SHORT" else "NOT SCORED (k = 0 incomplete)"
        elif checks["part0"]["status"] != "PASS" or checks["mirror_det_k0"]["status"] != "PASS":
            res["production"] = "NOT FINAL (implementation check incomplete)"
        else:
            res["production"] = lv[0]["summary"]
        if not b1_done:
            res["refinement"] = ("NOT RESOLVED BY THIS LADDER (BUDGET-SHORT)" if state == "BUDGET-SHORT"
                                 else "NOT SCORED (B1 incomplete)")
        elif any(checks[n]["status"] != "PASS" for n in ("part0", "mirror_det_k0", "mirror_det_k6")):
            res["refinement"] = "NOT FINAL (implementation check incomplete)"
        elif all(lv.get(k) for k in range(7)) and res["horizon"]:
            res["refinement"], res["refinement_info"] = mg8.label_refinement(lv, res["horizon"])
        if state == "BUDGET-SHORT":
            res["budget_short"] = True
    return res


# ------------------------------------------------------------------ B2 (descriptive, §4.4)
def welch(a, b):
    a, b = np.asarray(a, float), np.asarray(b, float)
    va, vb = a.var(ddof=1) / a.size, b.var(ddof=1) / b.size
    d = a.mean() - b.mean()
    df = (va + vb) ** 2 / (va ** 2 / (a.size - 1) + vb ** 2 / (b.size - 1))
    hw = stats.t.ppf(0.975, df) * np.sqrt(va + vb)
    return float(d), float(d - hw), float(d + hw)


def stoch_E(A, B):
    return mg8.shape_error(np.asarray(A), np.asarray(B), H)[0]


def boot_E(profA, profJ, seed):
    rng = np.random.default_rng(seed)
    vals = []
    for _ in range(BOOT_N):
        a = [profA[i] for i in rng.integers(0, len(profA), len(profA))]
        j = [profJ[i] for i in rng.integers(0, len(profJ), len(profJ))]
        vals.append(stoch_E(np.mean(a, axis=0), np.mean(j, axis=0)))
    return float(np.percentile(vals, 2.5)), float(np.percentile(vals, 97.5))


def dp(rec):
    p = np.array(rec["positions"], float)
    p = np.unwrap(p * 2 * np.pi / 100.0) * 100.0 / (2 * np.pi)
    dd = p[10:] - p[:-10]
    dd = dd[np.isfinite(dd)]
    return float(np.mean(dd ** 2) / (2 * 10.0)) if dd.size else None


def cp(k, n):
    """Exact (Clopper-Pearson) 95% interval for k of n."""
    lo = 0.0 if k == 0 else float(stats.beta.ppf(0.025, k, n - k + 1))
    hi = 1.0 if k == n else float(stats.beta.ppf(0.975, k + 1, n - k))
    return [lo, hi]


def evaluate_b2(runs):
    from runner import seed_for
    arms = {}
    for a in ("J", "R", "P"):
        rs = [runs[k] for k in sorted(runs) if k.startswith(f"sto_{a}_")]
        if not rs:
            continue
        ev = {t: sum(1 for r in rs if any(e[0] == t for e in r["events"])) for t in ("loss", "replication", "merge")}
        el = [r for r in rs if r["eligible"]]
        arms[a] = dict(n=len(rs), runs_with_event=ev, track_fail=sum(r["track_fail"] for r in rs),
                       n_eligible=len(el), v=[r["v"] for r in el], profiles=[r["profile"] for r in el],
                       Dp=[dp(r) for r in el], farvarX=[r["farvar"][1] for r in rs],
                       clip_frac_max=max(max(r["clip_frac"]) for r in rs))
        arms[a]["event_ci"] = {t: cp(ev[t], len(rs)) for t in ev}
    res = dict(arms={a: {k: v for k, v in d.items() if k != "profiles"} for a, d in arms.items()})
    if "J" in arms:
        J = arms["J"]
        pj = J["profiles"]
        if len(pj) >= 2:
            h1 = (len(pj) + 1) // 2
            res["J_split_half"] = stoch_E(np.mean(pj[:h1], axis=0), np.mean(pj[h1:], axis=0))
        else:
            res["J_split_half"] = None
        for a in ("R", "P"):
            if a not in arms:
                continue
            A = arms[a]
            comp = {}
            if A["n_eligible"] >= 5 and J["n_eligible"] >= 5:
                d, lo, hi = welch(A["v"], J["v"])
                comp["drift_diff"] = dict(d=d, lo=lo, hi=hi,
                                          verdict="equivalent within margin" if (lo > -1e-3 and hi < 1e-3) else "imprecise or different")
                rp = mg8.shape_report(np.mean(A["profiles"], axis=0), np.mean(pj, axis=0), H)
                E = rp["E"]
                blo, bhi = boot_E(A["profiles"], pj, seed_for(8, 2, 0, 0, 0))
                comp["E_stoch"] = dict(E=E, lo=blo, hi=bhi, E_raw=rp["E_raw"], E_nyq=rp["E_nyq"],
                                       verdict="equivalent within margin" if bhi < 1e-2 else "imprecise or different")
            else:
                comp["note"] = "NOT ESTIMATED (fewer than 5 eligible runs)"
            res[f"{a}_vs_J"] = comp
    res["mirror_stoch"] = mirror_stoch(runs)
    return res


def main():
    out = sys.argv[1] if len(sys.argv) > 1 else os.path.join(ROOT, "results", "MG8")
    runs = load(out)
    checks, st = impl_checks(out, runs), run_status(out)
    ev = dict(state=st, checks=checks, B1=evaluate_b1(out, runs, checks, st["state"]))
    if any(c["status"] == "FAIL" for c in checks.values()):
        ev["B2"] = WITHHELD                                   # §4.1: STOP, nothing is reported as a result
    else:
        ev["B2"] = evaluate_b2(runs)
        if checks["mirror_stoch"]["status"] != "PASS" or st["state"] not in ("complete", "BUDGET-SHORT"):
            ev["B2"]["note"] = "NOT FINAL (stochastic mirror check or production incomplete)"
        elif st["state"] == "BUDGET-SHORT":
            ev["B2"]["note"] = "BUDGET-SHORT: descriptive on the completed runs"
    json.dump(ev, open(os.path.join(out, "evaluation.json"), "w"), indent=1, default=float)
    print(json.dumps(dict(state=st["state"], checks={n: c["status"] for n, c in checks.items()},
                          production=ev["B1"].get("production"), refinement=ev["B1"].get("refinement")), indent=1))


if __name__ == "__main__":
    main()

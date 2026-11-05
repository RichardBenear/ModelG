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
    for f in sorted(os.listdir(rd)):
        if f.endswith(".json"):
            runs[f[:-5]] = json.load(open(os.path.join(rd, f)))
    return runs


def fields(out, key):
    z = np.load(os.path.join(out, "fields", key + ".npz"))
    return {int(k[1:]): z[k] for k in z.files}


def window_v(rec, t0, t1):
    t = np.array(rec["t"]); c = mg8.unwrap(rec["centre"], 100.0)
    w = (t >= t0 - 1e-6) & (t <= t1 + 1e-6)
    return mg8.slope(t[w], c[w])


def single_through(rec, t1):
    t = np.array(rec["t"]); n = np.array(rec["n_dets"])
    w = (t >= 50 - 1e-6) & (t <= t1 + 1e-6)
    return bool(len(t) and t[-1] >= t1 - 1e-3 and np.all(n[w] == 1) and rec["finite"])


def det_stats(out, runs, k, rep, T, Jf, Efloor):
    key = f"det_R_k{k}_r{rep}"
    rec = runs[key]
    Rf = fields(out, key)
    single = single_through(rec, T)
    E = E1 = v = v1 = None
    if T in Rf and (T - 50) in Rf and np.all(np.isfinite(Rf[T])):
        E = mg8.shape_error(Rf[T], Jf[T], H)[0]
        E1 = mg8.shape_error(Rf[T - 50], Jf[T - 50], H)[0]
        v = window_v(rec, T - 150, T); v1 = window_v(rec, T - 200, T - 50)
    st = E is not None and mg8.settled(E, E1, v, v1)
    return dict(single=single, E=E, E1=E1, v=v, v1=v1, settled=bool(st))


def evaluate_b1(out, runs):
    res = dict(levels={}, mirror={})
    for k in range(7):
        Jrec = runs.get(f"det_J_k{k}_r0")
        if Jrec is None or any(f"det_R_k{k}_r{j}" not in runs for j in range(3)):
            res["levels"][k] = None
            continue
        Jf = fields(out, f"det_J_k{k}_r0")
        T = 300
        Ef = mg8.shape_error(Jf[T - 50], Jf[T], H)[0] + 1e-9
        vJ = window_v(Jrec, T - 150, T)
        vf = abs(vJ) + 1e-8
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
        E600 = mg8.shape_error(Rf[600], Jf[600], H)[0]; E550 = mg8.shape_error(Rf[550], Jf[550], H)[0]
        v, v1 = window_v(rec, 450, 600), window_v(rec, 400, 550)
        res["horizon"] = dict(single=single_through(rec, 600), settled=bool(mg8.settled(E600, E550, v, v1)),
                              E300=res["levels"][6]["reps"][0]["E"], E600=E600, Efloor=res["levels"][6]["Efloor"])
    except (KeyError, TypeError):
        res["horizon"] = None
    # mirror checks (deterministic)
    for k in (0, 6):
        kL, kR = f"det_L_k{k}_r0", f"det_R_k{k}_r0"
        if kL in runs and kR in runs:
            Lf, Rf = fields(out, kL), fields(out, kR)
            d = max(float(np.max(np.abs(Lf[u][:, ::-1] - Rf[u])) / max(1e-300, float(np.max(np.abs(Rf[u])))))
                    for u in Lf)
            res["mirror"][k] = dict(max_rel=d, pass_=d <= 1e-12)
    lv = res["levels"]
    if lv.get(0):
        res["production"] = lv[0]["summary"]
    if all(lv.get(k) for k in range(7)) and res["horizon"]:
        res["refinement"], res["refinement_info"] = mg8.label_refinement(lv, res["horizon"])
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
                E = stoch_E(np.mean(A["profiles"], axis=0), np.mean(pj, axis=0))
                blo, bhi = boot_E(A["profiles"], pj, seed_for(8, 2, 0, 0, 0))
                comp["E_stoch"] = dict(E=E, lo=blo, hi=bhi,
                                       verdict="equivalent within margin" if bhi < 1e-2 else "imprecise or different")
            else:
                comp["note"] = "NOT ESTIMATED (fewer than 5 eligible runs)"
            res[f"{a}_vs_J"] = comp
    if "sto_L_k0_r0" in runs and "sto_R_k0_r0" in runs:
        res["mirror_stoch"] = runs["sto_L_k0_r0"]["counts_sha256"] == runs["sto_R_k0_r0"]["counts_sha256"]
    return res


def main():
    out = sys.argv[1] if len(sys.argv) > 1 else os.path.join(ROOT, "results", "MG8")
    runs = load(out)
    ev = dict(B1=evaluate_b1(out, runs), B2=evaluate_b2(runs))
    json.dump(ev, open(os.path.join(out, "evaluation.json"), "w"), indent=1, default=float)
    print(json.dumps(dict(production=ev["B1"].get("production"), refinement=ev["B1"].get("refinement")), indent=1))


if __name__ == "__main__":
    main()

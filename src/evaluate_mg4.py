"""MG4 evaluator (proposal §5.3 and §7; GATES.md MG4). Reads Part 0 decisions, Part A run records and saved fields;
no integration. Writes results/MG4/evaluation.json and prints the verdict with every criterion's inputs.
    python src/evaluate_mg4.py
Scope as re-stated by the PI on 2026-10-31: the soliton runs on T, D and R decide C1 to C3. The h = 0.3 check is reported
only. N3 enters only through C3(a).
"""
import os, sys, json, glob
import numpy as np
sys.path.insert(0, os.path.dirname(__file__))
import graphs as gr
import graph_ode as go
import measure_graph as mg
import run_mg4_partA as pa

ROOT = os.environ.get('MG4_EVAL_OUT') or pa.ROOT
TOL = dict(a=0.10, b=0.10, c=0.05, frac=0.80, c3a=0.05, c3b_pass=0.10, c3b_fail=0.20)


def load_runs():
    return {os.path.basename(p)[:-5]: json.load(open(p)) for p in sorted(glob.glob(os.path.join(pa.OUT, "runs", "*.json")))}


def c1(cls_by_key, fam):
    per = {s: sum(cls_by_key[f"{fam}{s}_p{i}_sol"] == "single" for i in range(8)) for s in range(4)}
    tot = sum(per.values())
    if tot <= 20 or min(per.values()) <= 3:
        v = "FAIL"
    elif tot >= 28 and min(per.values()) >= 6:
        v = "PASS"
    else:
        v = "INCONCLUSIVE"
    return dict(verdict=v, total=tot, per_graph=per)


def c2(an, fam, ref):
    rows = [a for k, a in an.items() if k.startswith(fam) and not k.startswith("h03") and a["cls"] == "single"]
    if len(rows) < 8:
        return dict(verdict="INCONCLUSIVE", reason="fewer than 8 single-particle runs", n=len(rows))
    rc = np.asarray(ref["rc"])
    ea = np.array([mg.profile_l2(rc, np.asarray(a["prof"]["X"], float), np.asarray(ref["X"])) for a in rows])
    ebs = np.array([abs(a["depth_s"] / ref["depth_s"] - 1) for a in rows])
    ebu = np.array([abs(a["depth_u"] / ref["depth_u"] - 1) for a in rows])
    ec = np.array([abs(a["lam"] / ref["lam"] - 1) if a["lam"] is not None else np.inf for a in rows])
    m = dict(a=float(np.median(ea)), b_s=float(np.median(ebs)), b_u=float(np.median(ebu)), c=float(np.median(ec)))
    frac = float(np.mean(ea < TOL["a"]))
    meds_ok = m["a"] < TOL["a"] and m["b_s"] <= TOL["b"] and m["b_u"] <= TOL["b"] and m["c"] <= TOL["c"]
    v = ("PASS" if frac >= TOL["frac"] else "INCONCLUSIVE") if meds_ok else "FAIL"
    return dict(verdict=v, n=len(rows), medians=m, frac_a=frac,
                per_run=dict(a=ea.tolist(), b_s=ebs.tolist(), b_u=ebu.tolist(), c=[float(x) for x in ec]))


def c3a(dec, fam):
    sp = {k: v for k, v in dec["spectra"].items() if v["fam"] == fam}
    means = [v["mean_ratio"] for v in sp.values()]
    if fam == "N3" or fam == "T":
        ok = all(abs(x - 1) <= TOL["c3a"] for x in means)
        return dict(passes=ok, means=means)
    if all(abs(x - 1) <= TOL["c3a"] for x in means):
        v = "PASS"
    elif abs(np.mean(means) - 1) > TOL["c3a"]:
        v = "FAIL"
    else:
        v = "INCONCLUSIVE"
    return dict(verdict=v, means=means, panel_mean=float(np.mean(means)))


def c3b_family(runs, dec, fam, cont):
    keys = [k for k in runs if runs[k]["fam"] == fam and runs[k]["torus"] == "main" and runs[k]["ic"] == "soliton"]
    vals = {k: pa.c3b(fam, runs[k]["slot"], np.asarray(runs[k]["centre"]), dec, cont) for k in sorted(keys)}
    med = {S: float(np.median([v[S] for v in vals.values()])) for S in "GXY"}
    if all(x <= TOL["c3b_pass"] for x in med.values()):
        v = "PASS"
    elif any(x > TOL["c3b_fail"] for x in med.values()):
        v = "FAIL"
    else:
        v = "INCONCLUSIVE"
    return dict(verdict=v, medians=med, per_run=vals)


def combine_c3(a, b, n3_passes):
    if n3_passes:
        return "INCONCLUSIVE"
    if "FAIL" in (a, b):
        return "FAIL"
    return "PASS" if a == b == "PASS" else "INCONCLUSIVE"


def main():
    dec = json.load(open(pa.DEC))
    theta = json.load(open(os.path.join(pa.OUT, "threshold.json")))["theta"]
    runs = load_runs()
    cont = go.continuum2d()
    out = dict(theta=theta, n_runs=len(runs))
    # 1. Part 0 stop
    if dec.get("stop_for_PI"):
        out["verdict"] = "NO VERDICT (Part 0 stop for PI)"
        return finish(out)
    # 2. Budget stop: every registered gate run present?
    need = [j["key"] for j in pa.jobs_for(dec)["T"] + pa.jobs_for(dec)["main"]]
    missing = [k for k in need if k not in runs]
    out["missing_runs"] = missing
    # 3. Implementation checks
    bad = [k for k, r in runs.items() if not r["balances"]["ok"]]
    out["balance_failures"] = bad
    an = {k: pa.analyse(r, theta, dec) for k, r in runs.items() if r["ic"] == "soliton"}
    out["classes"] = {k: a["cls"] for k, a in an.items()}
    out["stationary"] = {k: runs[k]["stationary"] for k in an}
    out["drift"] = {k: a.get("drift") for k, a in an.items()}
    # T reference particle: mean over its 4 placements at the same bins
    Tk = [f"T_p{i}_sol" for i in range(4)]
    Tsingle = all(an[k]["cls"] == "single" for k in Tk if k in an) and all(k in an for k in Tk)
    if Tsingle:
        rc = np.asarray(an[Tk[0]]["rc"])
        X = np.nanmean([an[k]["prof"]["X"] for k in Tk], axis=0)
        ref = dict(rc=rc.tolist(), X=X.tolist(), depth_u=float(X[0]),
                   depth_s=float(np.mean([an[k]["depth_s"] for k in Tk])), lam=mg.wavelength(rc, X))
        out["T_reference"] = ref
    fams = {}
    n3 = c3a(dec, "N3")
    t3a = c3a(dec, "T")
    t3b = c3b_family(runs, dec, "T", cont)
    for fam in ("D", "R"):
        f = dict(C1=c1(out["classes"], fam) if not missing else None,
                 C2=c2(an, fam, ref) if Tsingle else dict(verdict="INCONCLUSIVE", reason="no T reference"),
                 C3a=c3a(dec, fam), C3b=c3b_family(runs, dec, fam, cont))
        f["C3"] = combine_c3(f["C3a"]["verdict"], f["C3b"]["verdict"], n3["passes"])
        fams[fam] = f
    out.update(families=fams, N3_C3a=n3, T_C3a=t3a, T_C3b=t3b, T_single_all=Tsingle)
    # h = 0.3 check (reported only)
    out["h03"] = {k: dict(cls=a["cls"], drift=a.get("drift"), depth_u=a.get("depth_u"), lam=a.get("lam"),
                          h=runs[k]["h"]) for k, a in an.items() if k.startswith("h03")}
    # Precedence (§7)
    if missing:
        v = "INCONCLUSIVE (budget stop)"
    elif bad:
        v = "NOT SCORED (implementation failure: investigate balance residuals)"
    else:
        pos_fail = (not Tsingle) or (not t3a["passes"]) or t3b["verdict"] != "PASS"
        crit = [fams[f][c] if c == "C3" else fams[f][c]["verdict"] for f in ("D", "R") for c in ("C1", "C2", "C3")]
        if pos_fail:
            v = "INCONCLUSIVE (positive-control failure)" + (
                "; FAIL against a failed comparator among D/R" if "FAIL" in crit else "")
        elif "FAIL" in crit:
            v = "FAIL"
        elif "INCONCLUSIVE" in crit:
            v = "INCONCLUSIVE"
        else:
            v = "PASS"
    out["verdict"] = v
    return finish(out)


def finish(out):
    path = os.path.join(ROOT, "evaluation.json")
    json.dump(out, open(path, "w"), indent=1, default=float)
    print("MG4 verdict:", out["verdict"])
    for fam, f in out.get("families", {}).items():
        print(fam, {c: (f[c] if isinstance(f[c], str) else (f[c] or {}).get("verdict")) for c in ("C1", "C2", "C3a", "C3b", "C3")})
    print("written", path)
    return out


if __name__ == "__main__":
    main()

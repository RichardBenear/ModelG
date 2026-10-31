"""MG4c evaluator (spec/MG4c_proposal.md draft 3 §5; GATES.md MG4c and its pre-run implementation specification).
Reads results/MG4c/ (driver_log.json, runs/, fields/, resolution.json, part3_checks.json); no integration. Writes
results/MG4c/evaluation.json and prints the window status (W1/W2/W3/INCOMPLETE/NOT SCORED) and the operator status
(O1/O2/O3/O-invalid) with their inputs.
    python src/evaluate_mg4c.py
Classes are recomputed here from the saved fields with the registered thresholds (not copied from the driver log), and
compared with the driver's; a disagreement is reported and stops the evaluation.
"""
import os, sys, json
import numpy as np
sys.path.insert(0, os.path.dirname(__file__))
import graphs as gr
import measure_graph as mg
import run_mg4c as rm

OUT = rm.OUT
CTRL_TOL = 1e-6


def records():
    return rm.all_records()


def single_shape(rec, theta):
    """Depth (smoothed profile at the particle centre) and lambda of a single particle, as MG4's analyse()."""
    dec = json.load(open(rm.DEC))
    g = rm.graph_for(rec, dec)[0]
    f = np.load(os.path.join(OUT, "fields", rec["key"] + ".npz"))
    checks = [mg.particles(g, f["checks_sX"][k].astype(float), theta) for k in range(len(rm.CHECKS))]
    cls, track = mg.classify(g, checks)
    if cls != "single":
        return None
    cen = track[-1]
    fin = f["final"].astype(float)
    rc, pX, _ = mg.radial_profile(g, fin[1], cen)
    rcs, ps, _ = mg.radial_profile(g, mg.smoother(g) @ fin[1], cen)
    return dict(depth_s=float(ps[0]), depth_u=float(pX[0]), lam=mg.wavelength(rc, pX))


def window_status(tables, part1_done, unresolved):
    """Proposal §5 window statuses from the class tables {b: table} of every tested b (b = 29 included)."""
    if unresolved:
        return "NOT SCORED"
    if not part1_done:
        return "INCOMPLETE"
    if any(rm.all_single(t) for t in tables.values()):
        return "W1"
    if not any(c == "single" for t in tables.values() for c in t["D"]):
        return "W2"
    return "W3"


def operator_status(checks_ok, part3):
    if not checks_ok or part3 is None:
        return "O-invalid" if not checks_ok else "not run"
    cT, cD = part3.get("w_b29.000000_T_p0"), part3.get("w_b29.000000_D0_p0")
    if cT is None or cD is None:
        return "not run"
    if cT != "single" or not (part3.get("T_control_maxdiff_X", np.inf) <= CTRL_TOL):
        return "O-invalid"
    return {"single": "O1", "multiplied": "O2", "decayed": "O2", "unresolved": "O3"}[cD]


def main():
    log = json.load(open(os.path.join(OUT, "driver_log.json")))
    recs = records()
    out = dict(driver_status=log.get("status"), core_h=log.get("core_h"), runs=len(recs), platform=log.get("platform"))
    # balance flags: original flags kept; resolution file read before any class
    res_path = os.path.join(OUT, "resolution.json")
    res = json.load(open(res_path)) if os.path.exists(res_path) else {}
    flagged = [k for k, r in recs.items() if not r["balances"]["ok"]]
    unresolved = [k for k in flagged if not res.get(k, {}).get("resolved", False)]
    out["balance"] = dict(flagged=flagged, resolved=[k for k in flagged if k not in unresolved], unresolved=unresolved)
    # tested b values: grid points whose five runs exist, plus refinement points
    tested = sorted({r["b"] for r in recs.values() if r["op"] is None}, reverse=True)
    complete = [b for b in tested if all(j["key"] in recs for j in rm.jobs_at(b))]
    part1_done = all(b in complete for b in rm.GRID)
    tables = {29.0: rm.table_at(29.0, recs)}
    for b in complete:
        tables[b] = rm.table_at(b, recs)
    # agreement with the driver's tables
    drv = {**log.get("grid", {}), **log.get("refined", {})}
    mism = [rm.bkey(b) for b, t in tables.items() if rm.bkey(b) in drv and drv[rm.bkey(b)]["D"] != t["D"]]
    if mism:
        out["error"] = f"classes differ from the driver log at {mism}"
        return finish(out)
    rows = {}
    for b, t in sorted(tables.items(), reverse=True):
        row = dict(theta=t["theta"], T_eligible=t["eligible"], D=dict(zip(["D0_p0", "D0_p1", "D1_p0", "D1_p1"], t["D"])),
                   all_single=rm.all_single(t), any_multiplied=rm.any_mult(t), n_single=sum(c == "single" for c in t["D"]),
                   source="grid" if b in [29.0] + rm.GRID else "refinement")
        for k in ("reason", "T_min", "T_class_theta29", "D_theta29", "D_min", "source"):
            if k in t:
                row[k if k != "source" else "records"] = t[k]
        if row["all_single"] and b != 29.0:                       # depth and lambda against T only where T is eligible
            if t["eligible"]:
                Tsh = single_shape(recs[f"{rm.bkey(b)}_T_p0"], t["theta"])
                Dsh = [single_shape(recs[f"{rm.bkey(b)}_D{s}_p{p}"], t["theta"]) for s, p in rm.DRUNS]
                row["shape_vs_T"] = dict(T=Tsh, D=Dsh,
                                         depth_s_ratio=[d["depth_s"] / Tsh["depth_s"] for d in Dsh],
                                         lam_ratio=[(d["lam"] / Tsh["lam"]) if d["lam"] and Tsh["lam"] else None
                                                    for d in Dsh])
            else:
                row["shape_vs_T"] = "unavailable (T not an eligible reference at this b)"
        rows[rm.bkey(b)] = row
    out["tables"] = rows
    out["case"] = log.get("case")
    out["edges"] = log.get("edges")
    out["case_B"] = log.get("case_B")
    out["all_single_b"] = [b for b, t in sorted(tables.items(), reverse=True) if rm.all_single(t)]
    out["window_status"] = window_status(tables, part1_done, unresolved)
    if not part1_done:
        out["missing_part1"] = [b for b in rm.GRID if b not in complete]
    # Part 3
    chk_path = os.path.join(OUT, "part3_checks.json")
    checks_ok = json.load(open(chk_path))["ok"] if os.path.exists(chk_path) else False
    p3 = {k: rm.classify_run(recs[k], rm.THETA_29) for k in ("w_b29.000000_D0_p0", "w_b29.000000_T_p0") if k in recs}
    for k in ("T_control_maxdiff_X", "T_control_maxdiff_all"):
        if k in log.get("part3", {}):
            p3[k] = log["part3"][k]
    out["part3"] = dict(checks_ok=checks_ok, **p3)
    out["operator_status"] = operator_status(checks_ok, p3 if p3 else None)
    return finish(out)


def finish(out):
    path = os.path.join(OUT, "evaluation.json")
    json.dump(out, open(path, "w"), indent=1, default=float)
    if "error" in out:
        print("MG4c evaluation stopped:", out["error"])
        return out
    print("MG4c window status:", out["window_status"], "| operator status:", out["operator_status"])
    for b, r in out["tables"].items():
        print(f"  {b}: theta {r['theta']:.3f} (T eligible {r['T_eligible']}) D {list(r['D'].values())}")
    print("written", path)
    return out


if __name__ == "__main__":
    main()

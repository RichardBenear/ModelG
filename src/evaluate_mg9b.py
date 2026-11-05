"""MG9b evaluation from the saved records (spec/MG9b_proposal.md draft 3 §4-§5; registered, frozen at 185aceb).

Reads results/MG9b/{partA_*, partB_*, vb1_*, vb2_*, vb3_*, ledger}.json, runs/*/*.json, the MG3b b = 29.00 cells
(results/MG3b/evaluation.json, used as recorded) and the calibration labels (results/MG9b/calibration/
calibration.json). Writes results/MG9b/evaluation.json. Usage: python src/evaluate_mg9b.py
"""
import os, sys, json, glob
import numpy as np
sys.path.insert(0, os.path.dirname(__file__))
import mg9b
from mg3 import garwood

ROOT = os.path.join(os.path.dirname(__file__), "..")
OUT = os.environ.get("MG9B_OUT", os.path.join(ROOT, "results", "MG9b"))
P = lambda *a: os.path.join(OUT, *a)


def _load(fn):
    return json.load(open(fn)) if os.path.exists(fn) else None


def runs_of(cell):
    return sorted((json.load(open(f)) for f in glob.glob(P("runs", cell, "*.json"))), key=lambda r: r["rep"])


def cell_stats(recs, n_used=None):
    """Per-cell rate with risk sets and landmark sensitivities (§4.1). Only the first n_used reps (the batching rule)
    and complete runs count."""
    rr = [r for r in recs[: n_used or len(recs)] if r["complete"]]
    n = sum(r["event"] for r in rr)
    E = sum(r["exposure"] for r in rr)
    lo, hi = garwood(n)
    out = dict(runs=len(rr), events=int(n), exposure=float(E), startup=int(sum(r["startup"] for r in rr)),
               k=n / (mg9b.L * E) if E > 0 else None, k_ci=[lo / (mg9b.L * E), hi / (mg9b.L * E)] if E > 0 else None)
    out["risk_set"] = {str(t): sum(1 for r in rr if not r["startup"] and (r["t_first"] is None or r["t_first"] >= t))
                       for t in (50, 100, 150)}
    lm = {}
    for t in (100.0, 150.0):
        at = [r for r in rr if not r["startup"] and (r["t_first"] is None or r["t_first"] >= t)]
        ev = [r for r in at if r["event"] and r["t_first"] is not None and r["t_first"] <= mg9b.TC]
        ex = sum(min(r["t_first"] if (r["event"] and r["t_first"] is not None) else mg9b.TC, mg9b.TC) - t for r in at)
        l2, h2 = garwood(len(ev))
        lm[str(int(t))] = dict(at_risk=len(at), events=len(ev), exposure=float(ex),
                              k=len(ev) / (mg9b.L * ex) if ex > 0 else None,
                              k_ci=[l2 / (mg9b.L * ex), h2 / (mg9b.L * ex)] if ex > 0 else None)
    out["landmarks"] = lm
    return out


def historical_cells():
    ev = _load(os.path.join(ROOT, "results", "MG3b", "evaluation.json"))
    return [(float(Om), int(c["n"]), float(c["exposure"])) for Om, c in ev["part2"]["cells"].items()]


def event_shapes(recs):
    fe = [r for r in recs if r["event"] and r["first"] is not None]
    conf = {(r["rep"], c["sc"]): c for r in recs for c in r["confirmations"]}
    shapes = [conf.get((r["rep"], r["first"]["s_conf"])) for r in fe]
    shapes = [s for s in shapes if s is not None]
    return dict(events=len(shapes), single_minimum=sum(1 for s in shapes if s["n_deep_end"] == 1),
                n_persistent_hist={str(k): sum(1 for s in shapes if s["n_persistent"] == k)
                                   for k in sorted({s["n_persistent"] for s in shapes})})


def evaluate():
    cal = _load(P("calibration", "calibration.json"))
    labels = (cal or {}).get("labels") or {}
    by_b, per_b = {29.0: historical_cells()}, {}
    per_b["29.00"] = dict(source="MG3b (used as recorded)", cells=by_b[29.0], fit=mg9b.fit_b(by_b[29.0]))
    lg = _load(P("ledger.json")) or {}
    bstop = lg.get("budget_stop")
    # amendment 1: study-level flag (any budget stop that truncated rate-data collection), not only short cells
    budget_truncated = bool(bstop) and not bstop.get("rate_complete", False)
    for b in mg9b.B_NEW:
        key = f"{b:.2f}"
        v1, pa, pb = _load(P(f"vb1_{key}.json")), _load(P(f"partA_{key}.json")), _load(P(f"partB_{key}.json"))
        v2, v3 = _load(P(f"vb2_{key}.json")), _load(P(f"vb3_{key}.json"))
        d = dict(vb1=v1 and v1["result"], vb2=v2 and v2["result"], vb3=v3 and v3["result"])
        if v1 is not None and not v1["pass_"]:
            d["status"] = "REFERENCE-INVALID"
        elif pa is None:
            d["status"] = "NOT RUN"
        elif pa["Om_T"] is None:
            d["status"] = "NOT PLACED"
        else:
            d["Om_T"] = pa["Om_T"]
            cells, stats, recs_all = [], {}, []
            for c in (pb or {}).get("cells", []):
                recs = runs_of(f"b{b:.2f}_B_{c['Om']:g}")
                s = cell_stats(recs, None if c.get("budget_short") else c["runs"])
                s.update(capped=c["capped"], budget_short=c.get("budget_short", False))
                stats[f"{c['Om']:g}"] = s
                cells.append((c["Om"], s["events"], s["exposure"]))
                recs_all += recs
            d.update(cells=stats, identifiable=mg9b.identifiable(cells), fit=mg9b.fit_b(cells),
                     event_shape=event_shapes(recs_all),
                     capped=any(s["capped"] for s in stats.values()),
                     budget_short=any(s["budget_short"] for s in stats.values()))
            if d["identifiable"]:
                by_b[b] = cells
            d["status"] = "SLOPE-ELIGIBLE" if d["identifiable"] else "NOT SLOPE-ELIGIBLE"
        flags = [x for x in (d["vb2"], d["vb3"]) if x not in (None, "PASS")]
        if d.get("status") == "SLOPE-ELIGIBLE" and (d["vb2"] is None or d["vb3"] is None):
            flags.append("CONTROL-INSUFFICIENT (control not run)")
        d["lattice_flags"] = flags
        per_b[key] = d
    tr = mg9b.fit_trend(by_b)
    br = mg9b.branch_of(tr.get("b_used", []))
    cal_ok = labels.get(br, {}).get("label") == "calibrated"
    prim = mg9b.primary(tr, cal_ok, budget_truncated)
    q = []
    if tr.get("lr") and tr["lr"]["p"] < 0.05:
        q.append("LINEAR MODEL MISFIT")
    nm = mg9b.nonmonotone({b: per_b[f"{b:.2f}"]["fit"] for b in by_b})
    if nm:
        q.append("NON-MONOTONE")
    for b in mg9b.B_NEW:
        d = per_b[f"{b:.2f}"]
        if d["lattice_flags"]:
            q.append(f"LATTICE-FLAGGED ({b:.2f})")
        if d.get("capped"):
            q.append(f"CAPPED ({b:.2f})")
        if d.get("budget_short"):
            q.append(f"BUDGET-SHORT ({b:.2f})")
    if budget_truncated:
        q.append("BUDGET-TRUNCATED (study)")
    elif bstop:
        q.append("BUDGET STOP AFTER RATE COLLECTION (controls only)")
    if br in ("B2-20", "B2-35"):
        q.append(f"TWO-B (missing {'29.35' if br == 'B2-20' else '29.20'})")
    out = dict(spec="spec/MG9b_proposal.md draft 3 (frozen 185aceb) with amendment 1", per_b=per_b, trend=tr, branch=br,
               budget_stop=bstop, budget_truncated=budget_truncated,
               branch_label=labels.get(br), primary=prim, qualifiers=q,
               fitted_range=[min(tr["b_used"]), max(tr["b_used"])] if tr.get("b_used") else None,
               ledger_core_h=sum(e["wall_s"] for e in (_load(P("ledger.json")) or {"entries": []})["entries"]) / 3600)
    json.dump(out, open(P("evaluation.json"), "w"), indent=1, default=lambda o: o.tolist() if hasattr(o, "tolist") else float(o))
    return out


if __name__ == "__main__":
    o = evaluate()
    print(o["primary"], o["qualifiers"], o["branch"])

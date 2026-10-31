"""MG5 evaluator (proposal §7-§8; GATES.md MG5). Recomputes every class from the committed sparse check records
(results/MG5/checks/*.npz) with the registered thresholds, stops if any differs from the driver's record, then applies
run statuses, validity rules and the registered outcomes. Writes results/MG5/evaluation.json.
    python src/evaluate_mg5.py
"""
import os, sys, json
import numpy as np
sys.path.insert(0, os.path.dirname(__file__))
import evolve as ev
import run_mg5 as R

PINNED_CONTROL_D = 0.1


def classes_from_checks(path, theta):
    recs, _ = R.load_checks(path)
    parts = [ev.detect_sparse(r, theta) for r in recs]
    Lx, Ly = recs[0]["Lx"], recs[0]["Ly"]
    out = dict(n0=len(parts[0]), counts=[len(p) for p in parts[1:]],
               depths=[[q["depth"] for q in p] for p in parts[1:]])
    if len(parts[0]) != 1:
        out["cls"] = dict(cls="unresolved", sub="start not single")
        return out
    r = ev.classify_run(parts[1:], parts[0][0]["centre"], Lx, Ly)
    if r.get("cls") == "single":
        r["D_p"] = ev.msd_dp(r["track"])
    out["cls"] = r
    out["c0"] = np.asarray(parts[0][0]["centre"]).tolist()
    return out


def majority(classes, want, k=3):
    return sum(c == want for c in classes) >= k


def score(table):
    """table: {key: dict(part, nu, status, cls (dict or None))}. Pure function of the per-run results (proposal §8)."""
    out = dict(partA={}, partB={})
    grp = lambda part, nu: [v for v in table.values() if v["part"] == part and v["nu"] == nu]
    # validity: Part A frozen controls
    ctrl = grp("A", 0.0)
    ctrl_ok = (len(ctrl) == 4 and all(v["status"] == R.COMPLETE for v in ctrl)
               and all(v["cls"]["cls"] == "single" and v["cls"].get("max_d", np.inf) < PINNED_CONTROL_D for v in ctrl))
    out["partA"]["controls"] = [dict(status=v["status"], cls=(v["cls"] or {}).get("cls"),
                                     max_d=(v["cls"] or {}).get("max_d")) for v in ctrl]
    if len(ctrl) < 4 or any(v["status"] != R.COMPLETE for v in ctrl):
        out["partA"]["validity"] = "NOT SCORED (controls incomplete)"
    elif not ctrl_ok:
        out["partA"]["validity"] = "INVALID-CONTROL"
    else:
        out["partA"]["validity"] = "valid"
    rates = {}
    for nu in (0.002, 0.02, 0.2):
        g = grp("A", nu)
        if len(g) < 4 or any(v["status"] != R.COMPLETE for v in g):
            rates[nu] = dict(scored=False, statuses=[v["status"] for v in g])
            continue
        cls = [v["cls"]["cls"] for v in g]
        r = dict(scored=True, classes=cls, subs=[v["cls"].get("sub") for v in g], survives=majority(cls, "single"))
        if r["survives"]:
            mob = [v["cls"].get("mobility") for v in g if v["cls"]["cls"] == "single"]
            r["mobility"] = "PINNED" if majority(mob, "pinned") else "MOBILE" if majority(mob, "mobile") else "MIXED"
            r["max_d"] = [v["cls"].get("max_d") for v in g]
            r["D_p"] = [v["cls"].get("D_p") for v in g]
        rates[nu] = r
    out["partA"]["rates"] = {f"{k:g}": v for k, v in rates.items()}
    if out["partA"]["validity"] != "valid":
        out["gate"] = out["partA"]["validity"]
    elif not all(r["scored"] for r in rates.values()):
        out["gate"] = "INCOMPLETE"
    else:
        surv = [nu for nu, r in rates.items() if r["survives"]]
        out["gate"] = "PASS" if len(surv) == 3 else "FAIL" if not surv else "PARTIAL"
        if out["gate"] == "PARTIAL":
            out["largest_surviving_rate"] = max(surv)
    # Part B
    fz, fa = grp("B", 0.0), grp("B", 0.2)
    if len(fz) < 4 or len(fa) < 4 or any(v["status"] != R.COMPLETE for v in fz + fa):
        out["partB"]["outcome"] = "NOT SCORED"
    else:
        cz, ca = [v["cls"]["cls"] for v in fz], [v["cls"]["cls"] for v in fa]
        out["partB"].update(frozen=cz, fast=ca)
        if not majority(cz, "multiplied"):
            out["partB"]["outcome"] = "A-void"
        elif majority(ca, "single"):
            out["partB"]["outcome"] = "A1"
        elif majority(ca, "multiplied"):
            out["partB"]["outcome"] = "A2"
        else:
            out["partB"]["outcome"] = "A3"
    return out


def main():
    rd = os.path.join(R.OUT, "runs")
    table, mismatches = {}, []
    for f in sorted(os.listdir(rd)):
        if not f.endswith(".json"):
            continue
        rec = json.load(open(os.path.join(rd, f)))
        spec = rec["spec"]
        cls = None
        if rec["status"] == R.COMPLETE:
            mine = classes_from_checks(os.path.join(R.OUT, "checks", rec["key"] + ".npz"), rec["theta"])
            theirs = rec["summary"]
            if mine["counts"] != theirs["counts"] or mine["cls"].get("cls") != theirs["cls"].get("cls"):
                mismatches.append(rec["key"])
            cls = mine["cls"]
        table[rec["key"]] = dict(part=spec["part"], nu=spec["nu"], status=rec["status"], cls=cls,
                                 flags=len(rec.get("flags", [])))
    if mismatches:
        raise SystemExit(f"class recomputed from sparse records differs from the driver for: {mismatches}")
    expected = {s["key"] for s in R.run_specs()}
    for k in expected - set(table):
        s = next(x for x in R.run_specs() if x["key"] == k)
        table[k] = dict(part=s["part"], nu=s["nu"], status=R.INCOMPLETE, cls=None, flags=0)
    out = score(table)
    out["runs"] = table
    ch, done = R.spent()
    out["core_h"] = ch
    json.dump(out, open(os.path.join(R.OUT, "evaluation.json"), "w"), indent=1, default=float)
    print("Gate MG5:", out["gate"])
    for nu, r in out["partA"]["rates"].items():
        print(f"  nu = {nu}: {r}")
    print("Part B:", out["partB"].get("outcome"), out["partB"])
    print(f"core-h {ch:.2f}")


if __name__ == "__main__":
    main()

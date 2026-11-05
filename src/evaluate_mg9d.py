"""MG9d evaluator (spec/MG9d_proposal.md draft 3 §4). Re-reads results/MG9d/runs/*.json and floor.json and
rewrites results/MG9d/evaluation.json; prints the reading and the component labels.
    python src/evaluate_mg9d.py
"""
import os, sys, json
sys.path.insert(0, os.path.dirname(__file__))
import run_mg9d as D


def main():
    with open(os.path.join(D.OUT, "floor.json")) as f:
        vf = json.load(f)["v_f"]
    recs = {k: D.load_rec(k) for k in D.ORDER}
    ev = D.evaluate(recs, vf)
    lg = D.ledger()
    ev["integrations"] = len(lg["integrations"])
    ev["core_h"] = sum(r["wall_s"] or 0.0 for r in lg["integrations"]) / 3600.0
    ev["platform"] = D.platform()
    with open(os.path.join(D.OUT, "evaluation.json"), "w") as f:
        json.dump(ev, f, indent=1, default=float)
    print(json.dumps(dict(reading=ev["reading"], components=ev["components"],
                          valid={k: v["valid"] for k, v in ev["validity"].items()}), indent=1))
    return ev


if __name__ == "__main__":
    main()

"""MG9c evaluator (spec/MG9c_proposal.md draft 2, §§4-7). Reads results/MG9c/ and writes results/MG9c/evaluation.json.
    python src/evaluate_mg9c.py
"""
import os, sys, json
import numpy as np
sys.path.insert(0, os.path.dirname(__file__))
import run_mg9c as C


def main():
    pts, desc = {}, {}
    for b in C.Q1_B + (C.DESC_B,):
        rec = C.load_rec(C.main_key(b))
        if rec is None:
            entry = dict(label="UNRESOLVED-NOT-RUN", replay_ok=None, status=None)
        else:
            rb = C.load_rec(f"robust_{b:.5f}")
            entry = dict(label=rec["label"], replay_ok=rec.get("replay", {}).get("ok"),
                         replay=rec.get("replay"), status=rb["status"] if rb else None,
                         trials=rb["trials"] if rb else None, T=rec["T"], extended=rec["extended"],
                         final_edges=[rec["track"][-1]["xL"], rec["track"][-1]["xR"]],
                         residual=rec["residual"], info={k: v for k, v in rec["info"].items()})
            if rec["label"] == "SETTLED" and rb is None:
                entry["status"] = "ROBUSTNESS-UNRESOLVED (not run)"
        (desc if b == C.DESC_B else pts)[f"{b:.2f}"] = entry
    q1 = C.q1_outcome({k: dict(label=v["label"], replay_ok=v["replay_ok"],
                               status=("ROBUSTNESS-UNRESOLVED" if (v["status"] or "").startswith("ROBUSTNESS-UNRESOLVED")
                                       else v["status"]))
                       for k, v in pts.items()})
    q2 = C.load_rec("q2")
    width = {}
    for b in C.Q1_B + (C.DESC_B,):
        for k in (f"front_D02L50_{b:.5f}", f"front_D02w_{b:.5f}"):
            r = C.load_rec(k)
            if r:
                width[k] = dict(label=r["label"], final_edges=[r["track"][-1]["xL"], r["track"][-1]["xR"]])
    lg = C.ledger()
    out = dict(Q1=q1, points=pts, descriptive_29_22=desc, Q2=q2, width_trials=width,
               integrations=len(lg["integrations"]),
               core_h=sum(r["wall_s"] for r in lg["integrations"]) / 3600.0, platform=C.platform())
    with open(os.path.join(C.OUT, "evaluation.json"), "w") as f:
        json.dump(out, f, indent=1, default=float)
    print(json.dumps(dict(Q1=q1, Q2=(q2 or {}).get("comparison"), integrations=out["integrations"],
                          core_h=round(out["core_h"], 3)), indent=1))
    return out


if __name__ == "__main__":
    main()

"""MG8 post-hoc diagnostic (unregistered; review MG038 item 3). Homogeneous start at J's steady state S0 (no particle),
deterministic J and R, k = 0 and k = 3, 60 time units, seam seed 123 (not a registered stream). Reports the offset of the
mean field from S0 and the residual spatial range. Writes results/MG8/diag_homogeneous.json.
Usage: python src/diag_mg8_homogeneous.py"""
import os, sys, json, platform
import numpy as np
sys.path.insert(0, os.path.dirname(__file__))
import mg8
from mg8 import J, R, TAU0, H

ROOT = os.path.join(os.path.dirname(__file__), "..")


def main():
    Pm = mg8.params()
    ref, _ = mg8.reference_start()
    S0 = ref["S0"]
    out = dict(note="post-hoc, unregistered; homogeneous start at S0, 60 units, seam seed 123", S0=list(map(float, S0)),
               platform=platform.platform(), numpy=np.__version__, cases=[])
    for k in (0, 3):
        tau = TAU0 / 2 ** k
        spu = mg8.steps_per_unit(tau)
        for sch in (J, R):
            G, X, Y = (np.full(500, S0[i]) for i in range(3))
            for u in range(60):
                mg8.run_det(G, X, Y, sch, Pm["a"], Pm["b"], Pm["g"], Pm["p"], Pm["q"], Pm["dx"], Pm["dy"], H, tau, spu,
                            123 if u == 0 else -1)
            out["cases"].append(dict(k=k, sched=mg8.NAMES[sch],
                                     mean_offset=[float(G.mean() - S0[0]), float(X.mean() - S0[1]), float(Y.mean() - S0[2])],
                                     range=[float(np.ptp(G)), float(np.ptp(X)), float(np.ptp(Y))]))
    json.dump(out, open(os.path.join(ROOT, "results", "MG8", "diag_homogeneous.json"), "w"), indent=1)
    for c in out["cases"]:
        print(c)


if __name__ == "__main__":
    main()

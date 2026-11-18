"""EX009 part B: EX008 test model at small consumption rate k (G feeds X), radial 3D, R = 60."""
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ex008_longrange as e  # noqa: E402

OUT = os.path.join(e.ROOT, "results", "EX009_followups")
out = []
prev = None
for k in (0.05, 0.02, 0.01, 0.005, 0.002, 0.001):
    M = e.Model(60.0, 0.05, k, 1, 0)
    if prev is None:
        d = np.load(os.path.join(e.OUT, "particle_GX_k0.05_R60.npz"))
        u = d["u"]
    else:
        u = prev
    u, res, ok = M.newton(u)
    run = dict(k=k, ok=bool(ok), residual=res)
    if ok and u[2, 0] > 1:
        a = e.analyse(M, u)
        dr = M.drift(u)
        run.update(M=a["consumption_M"], core_G=a["core_G"], far_dev=a["far"]["rel_dev_from_point_sink"],
                   eig_l0=e.leading_eigs(M, u, 0)[0], eig_l1=e.leading_eigs(M, u, 1)[:2], mu=dr["mu"],
                   pair_coeff=2 * dr["mu"] * a["consumption_M"] / (4 * np.pi))
        prev = u
    print(run, flush=True)
    out.append(run)
ks = np.array([r["k"] for r in out if "mu" in r])
for key in ("M", "mu", "pair_coeff"):
    v = np.array([abs(r[key]) for r in out if "mu" in r])
    m = ks <= 0.01
    print(key, "log-log slope (k <= 0.01):", np.polyfit(np.log(ks[m]), np.log(v[m]), 1)[0])
json.dump(out, open(os.path.join(OUT, "partB_smallk.json"), "w"), indent=1, default=float)

"""EX004 post-hoc: budget ODE check, G inside/ahead of the front by band, front speed, comoving decay-rate prediction."""
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ex003_screening as x3  # noqa: E402

OUT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "results", "EX004_growth")
rec = json.load(open(os.path.join(OUT, "results.json")))
post = []
for g, tag in zip(rec["armG"], ["", "_fine", "", ""]):
    d = np.load(os.path.join(OUT, f"growth_eps{g['eps']:.3g}{tag}.npz"))
    r, A = d["r"], d["series"]
    q = g["q"]
    dr = r[1] - r[0]
    V = 4 * np.pi * ((r + dr / 2) ** 3 - (r - dr / 2) ** 3) / 3
    t, IG, IX = A[:, 0], A[:, 2], A[:, 3]
    # exact budget: dIG/dt = -q IG + 0.1 IX (central differences, interior points)
    lhs = np.gradient(IG, t)[5:-5]
    rhs = (-q * IG + 0.1 * IX)[5:-5]
    ode_err = float(np.max(np.abs(lhs - rhs)) / np.max(np.abs(rhs)))
    m = t >= 100
    v = float(np.polyfit(t[m], A[m, 1], 1)[0])
    ks = g["kappa_slow_linear"] if g["eps"] < 1 else x3.rates(q, g["b"])[1]["decay"]   # monotone G mode
    kap_comoving = 0.5 * (v + np.sqrt(v * v + 4 * ks * ks))
    bands = {}
    for name in [k for k in d.files if k.startswith("t")] + ["final"]:
        phi = d[name]
        tt = float(name[1:]) if name != "final" else g["T_end"]
        i = int(np.argmin(np.abs(t - tt)))
        fr = A[i, 1]
        if fr < 20:
            continue
        b = {}
        for lo, hi in ((0, fr - 6), (fr - 6, fr + 2), (fr + 2, fr + 5), (fr + 5, fr + 10), (fr + 10, fr + 20), (fr + 20, r[-1] + 1)):
            mm = (r >= lo) & (r < hi)
            b[f"{lo - fr:+.0f}..{hi - fr:+.0f}" if lo > 0 else "interior"] = float(np.sum(V[mm] * phi[0][mm]))
        b["interior_mean_phiG"] = float(np.sum(V * phi[0] * (r < fr - 6)) / np.sum(V * (r < fr - 6)))
        bands[name] = dict(t=tt, front=float(fr), **b)
    # decay rate ahead of the final front, from log|phiG| on [front+3, front+12]
    fr = A[-1, 1]
    mm = (r > fr + 3) & (r < fr + 12)
    y = d["final"][0][mm] * r[mm]
    kap_fit = float(-np.polyfit(r[mm], np.log(np.abs(y)), 1)[0]) if np.all(y > 0) or np.all(y < 0) else None
    post.append(dict(eps=g["eps"], tag=tag, ode_rel_err=ode_err, front_speed=v, kappa_static=ks,
                     kappa_comoving_pred=float(kap_comoving), kappa_ahead_fit=kap_fit,
                     sign_ahead=float(np.sign(np.mean(y))), bands=bands))
    print(g["eps"], tag, {k: v_ for k, v_ in post[-1].items() if k != "bands"})
    for k, b in bands.items():
        print("   ", k, {kk: round(vv, 3) for kk, vv in b.items()})
json.dump(post, open(os.path.join(OUT, "posthoc.json"), "w"), indent=1, default=float)

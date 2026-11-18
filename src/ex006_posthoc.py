"""EX006 post-hoc: far-field decay of v_rel(d), and a separate prediction of the long-range part from the single
particle's G tail and the adjoint mobility: v_rel ~ 2 mu_G * d(phi_G tail)/dx at the partner's position."""
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ex006_g2_1d as w  # noqa: E402

rec = json.load(open(os.path.join(w.OUT, "results.json")))
out = []
for a in rec["partA"]:
    s = a["sigma"]
    km = a["rates"][0]["decay"]
    B = sorted([b for b in rec["partB"] if b["sigma"] == s and b["dt"] == w.DT and b["L"] == w.CASES[s]["L"]
                and b["d"] > 0], key=lambda b: b["d"])
    d = np.array([b["d"] for b in B])
    v = np.array([b["v_rel"] for b in B])
    far = d >= 45
    c = np.polyfit(d[far], np.log(v[far]), 1)
    # single-particle G tail fit: phi_G = A e^{-km |x|} for |x| in [30, L/2 - 30]
    p = np.load(os.path.join(w.OUT, f"particle_s{s:g}.npz"))
    x, phiG = p["x"], p["U"][0] - p["h0"][0]
    msk = (x > 30) & (x < w.CASES[s]["L"] / 2 - 30)
    Lh = w.CASES[s]["L"]
    # periodic: phi = A (e^{-km x} + e^{-km (L - x)})
    basis = np.exp(-km * x[msk]) + np.exp(-km * (Lh - x[msk]))
    A = float(np.sum(basis * phiG[msk]) / np.sum(basis * basis))
    pred = 2 * a["mu_G"] * A * (-km) * np.exp(-km * d[far])     # d/dx of A e^{-km x} is -km A e^{-km x}
    first_pos = float(d[np.argmax(np.cumsum((v < 0)[::-1])[::-1] == 0)]) if np.any(v < 0) else float(d[0])
    o = dict(sigma=s, kappa_linear=km, kappa_from_vrel=float(-c[0]), amp_at_d0=float(np.exp(c[1])),
             G_tail_A=A, mu_G=a["mu_G"], pred_over_measured_median=float(np.median(pred / v[far])),
             last_attractive_d=float(d[v < 0].max()), v_at_40=float(v[d == 40][0]))
    print(o)
    out.append(o)
json.dump(out, open(os.path.join(w.OUT, "posthoc.json"), "w"), indent=1)

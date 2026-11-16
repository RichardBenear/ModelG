"""MG193 figure: existence window of the 1D continuum soliton in (a, b). Reads window.json and posthoc.json."""
import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
R = os.path.join(ROOT, "results", "MG193_existence_window")
w = json.load(open(os.path.join(R, "window.json")))
p = json.load(open(os.path.join(R, "posthoc.json")))
fig, ax = plt.subplots(figsize=(6.4, 4.4), dpi=150)
pts = [(14.0, w["a14"]["b_low"], w["a14"]["b_high"], "registered")]
for k, v in p["strip"].items():
    if v.get("reached"):
        pts.append((float(k), v["b_low"], v["b_high"], "post hoc"))
for a, lo, hi, lab in pts:
    ax.plot([a, a], [lo, hi], lw=4, color="#1f6fb4" if lab == "registered" else "#9ab7d6", solid_capstyle="butt")
for key, v in w["a_lines"].items():
    ax.plot([v["a_low"], v["a_high"]], [v["b"], v["b"]], lw=4, color="#c2571a", solid_capstyle="butt")
bc = sorted([(14.0, w["a14"]["b_c"])] + [(float(k), v["b_c"]) for k, v in p["strip"].items() if v.get("reached")])
ax.plot(*zip(*bc), "x--", color="0.45", lw=1)
ax.plot(14.0, w["a14"]["b_star"], "o", color="k", ms=4, label="recommended point (14, 29.205)")
ax.plot([], [], lw=4, color="#1f6fb4", label="b-window at a = 14 (registered)")
ax.plot([], [], lw=4, color="#9ab7d6", label="b-window at a = 13, 15, 16 (post hoc)")
ax.plot([], [], lw=4, color="#c2571a", label="a-windows at b = 29, 29.205 (registered)")
ax.plot([], [], "x--", color="0.45", lw=1, label="Turing threshold b_c(a)")
ax.set_xlabel("a")
ax.set_ylabel("b")
ax.set_title("MG193: 1D continuum soliton existence window (all edges are folds)")
ax.legend(fontsize=7, loc="upper left")
fig.tight_layout()
fig.savefig(os.path.join(ROOT, "paper", "figures", "MG193_existence_window.png"))
print("ok")

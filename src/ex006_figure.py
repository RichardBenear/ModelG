"""EX006 figure: G-2 two-particle relative velocity vs separation (1D)."""
import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
R = os.path.join(ROOT, "results", "EX006_g2_drift_pair")
rec = json.load(open(os.path.join(R, "results.json")))
post = {p["sigma"]: p for p in json.load(open(os.path.join(R, "posthoc.json")))}
cols = {0.0003: "#3b6fb6", 0.001: "#c0504d", 0.003: "#2e9e6b"}
Ls = {0.0003: 300.0, 0.001: 200.0, 0.003: 200.0}
fig, ax = plt.subplots(figsize=(8, 4.8))
for s, c in cols.items():
    B = sorted([b for b in rec["partB"] if b["sigma"] == s and b["dt"] == 0.005 and b["L"] == Ls[s] and b["d"] > 0],
               key=lambda b: b["d"])
    d = np.array([b["d"] for b in B])
    v = np.array([b["v_rel"] for b in B])
    ax.semilogy(d, np.where(v > 0, v, np.nan), "o-", ms=3, color=c, label=f"sigma = {s}: repulsive (v_rel > 0)")
    ax.semilogy(d, np.where(v < 0, -v, np.nan), "x", ms=4, color=c, alpha=0.6)
    p = post[s]
    dd = np.linspace(30, 70, 50)
    ax.semilogy(dd, p["amp_at_d0"] * np.exp(-p["kappa_linear"] * dd), ":", color=c, lw=1)
ax.set(xlabel="separation d", ylabel="|v_rel|  (o: apart, x: together)", ylim=(1e-12, 1e-1),
       title="EX006 (EXPLORATORY): Model G-2 pair force in 1D\n"
             "short range: oscillating; long range: repulsive, decay = monotone G rate (dotted)")
ax.legend(fontsize=8, frameon=False)
fig.tight_layout()
fig.savefig(os.path.join(ROOT, "paper", "figures", "EX006_g2_pair.png"), dpi=130)

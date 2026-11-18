"""EX007 figure: drift mobility and long-range pair velocity vs f (share of X supplied by G)."""
import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
r = json.load(open(os.path.join(ROOT, "results", "EX007_g2_gx", "results.json")))
f = np.array([s["f"] for s in r["scan"]])
mu = np.array([s["mu_G"] for s in r["scan"]])
A = np.array([s["tail_A"] for s in r["scan"]])
fig, ax = plt.subplots(1, 2, figsize=(12, 4.2))
ax[0].plot(f, mu * 1e6, "o-", color="#3b6fb6", label="drift mobility mu_G (x 1e-6): > 0 climbs G, < 0 descends")
ax[0].plot(f, A, "s--", color="#c0504d", label="far G-tail amplitude A (< 0 = well)")
ax[0].axhline(0, color="k", lw=0.5)
ax[0].set(xlabel="f = share of X supply coming from G (f = 0 is pure G-2)", title="Drift flips toward lower G near f = 0.28; well stays")
ax[0].legend(fontsize=8, frameon=False)
cols = {40.0: "#3b6fb6", 50.0: "#2e9e6b", 60.0: "#8e5ea2"}
for d, c in cols.items():
    v = np.array([[p["f"], p["v_rel"]] for p in r["pairs"] if p["d"] == d])
    v = v[np.argsort(v[:, 0])]
    pred = np.array([s["pred_far_force_amp"] * np.exp(-s["rates"][0]["decay"] * d) for s in r["scan"]])
    ax[1].plot(v[:, 0], v[:, 1] / np.abs(v[:, 1]).max(), "o-", color=c, label=f"d = {d:g}: measured")
    ax[1].plot(f, pred / np.abs(v[:, 1]).max(), "x:", color=c, label=f"d = {d:g}: 2 mu_G A kappa e^(-kappa d)")
ax[1].axhline(0, color="k", lw=0.5)
ax[1].set(xlabel="f", ylabel="v_rel (scaled per d)", title="Long-range pair force: repulsive -> ATTRACTIVE for f >= 0.4")
ax[1].legend(fontsize=7, frameon=False, ncol=2)
fig.suptitle("EX007 (EXPLORATORY): hybrid Model G-2 + G -> X, sigma = 0.001, 1D", fontsize=10)
fig.tight_layout()
fig.savefig(os.path.join(ROOT, "paper", "figures", "EX007_g2_gx.png"), dpi=130)

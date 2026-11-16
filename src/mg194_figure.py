"""MG194 figure: two-particle relative speed vs separation (scaled by exp(0.4608 d)), linear response, equilibria."""
import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
R = os.path.join(ROOT, "results", "MG194_two_particle")
r = json.load(open(os.path.join(R, "results.json")))
K = 0.4608152900342353
runs = sorted([x for x in r["partB_runs"] if x["L"] == 100 and x["d"] > 0], key=lambda x: x["d"])
d = np.array([x["d"] for x in runs])
v = np.array([x["v_rel"] for x in runs])
dl, vl = np.array(r["partB_lin"]["d"]), np.array(r["partB_lin"]["v_rel"])
f = r["partB"]["B2_fit"]
dd = np.linspace(8, 30, 800)
fit = np.exp(-f["kappa"] * dd) * (f["A"] * np.cos(f["k"] * dd) + f["B"] * np.sin(f["k"] * dd))
fig, axs = plt.subplots(2, 1, figsize=(7, 6.2), dpi=150, sharex=True)
ax = axs[0]
ax.plot(d, v, "o", ms=3, color="#1f6fb4", label="time integration (B-dyn)")
ax.axhline(0, color="0.6", lw=0.8)
ax.set_yscale("symlog", linthresh=1e-9)
ax.set_ylabel("v_rel = dd/dt  (< 0: attraction)")
ax.legend(fontsize=8)
ax.set_title("MG194: two-particle interaction, 1D continuum, a = 14, b = 29.205")
ax = axs[1]
ax.plot(d, v * np.exp(K * d), "o", ms=3, color="#1f6fb4", label="time integration")
ax.plot(dd, fit * np.exp(K * dd), "-", color="#1f6fb4", lw=1, label=f"fit: kappa {f['kappa']:.4f}, k {f['k']:.4f}")
m = (dl >= 8) & (dl <= 30)
ax.plot(dl[m], vl[m] * np.exp(K * dl[m]), "--", color="#c2571a", lw=1, label="linear response (tail overlap)")
for e in r["partB_eq"]:
    if e["converged"]:
        ax.plot(e["d_star"], 0, "o" if e["stable"] else "o", mfc="k" if e["stable"] else "w", mec="k", ms=5)
ax.plot([], [], "o", mfc="k", mec="k", label="stable pair equilibrium")
ax.plot([], [], "o", mfc="w", mec="k", label="unstable pair equilibrium")
ax.axhline(0, color="0.6", lw=0.8)
ax.set_xlabel("separation d")
ax.set_ylabel("v_rel × exp(0.4608 d)")
ax.legend(fontsize=7, ncol=2)
fig.tight_layout()
fig.savefig(os.path.join(ROOT, "paper", "figures", "MG194_two_particle.png"))
print("ok")

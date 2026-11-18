"""EX008 figure: 1/r G field of the test-model particle and its drift mobility."""
import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
R = os.path.join(ROOT, "results", "EX008_longrange")
rec = json.load(open(os.path.join(R, "results.json")))
fig, ax = plt.subplots(1, 3, figsize=(15, 4.2))
for Rr, c in ((60, "#3b6fb6"), (120, "#c0504d")):
    d = np.load(os.path.join(R, f"particle_GX_k0.5_R{Rr}.npz"))
    r, u = d["r"], d["u"]
    run = [x for x in rec["runs"] if x["variant"] == "G feeds X" and x["k"] == 0.5 and x["R"] == Rr][0]
    Mc = run["consumption_M"]
    ax[0].plot(r, r * (1 - u[3]), color=c, label=f"R = {Rr}: r (G_inf - G)")
    ax[0].plot(r, Mc / (4 * np.pi) * (1 - r / Rr), ":", color="k", lw=1)
ax[0].set(xlabel="r", ylabel="r x G deficit", title="G deficit = M/(4 pi r) (dotted: point sink, Dirichlet at R)")
ax[0].legend(fontsize=8, frameon=False)
d = np.load(os.path.join(R, "particle_GX_k0.5_R60.npz"))
r, u = d["r"], d["u"]
ax[1].plot(r, u[0] - 140 / 9, label="phi_X")
ax[1].plot(r, u[1] - 29 / (140 / 9), label="phi_Y")
ax[1].plot(r, u[2], label="W (core only)")
ax[1].plot(r, 10 * (u[3] - 1), label="10 x (G - G_inf)")
ax[1].set(xlim=(0, 10), xlabel="r", title="Particle (G feeds X, k = 0.5)")
ax[1].legend(fontsize=8, frameon=False)
for name, c in (("G feeds X", "#3b6fb6"), ("G feeds Y", "#2e9e6b")):
    xs = [x for x in rec["runs"] if x["variant"] == name and x["R"] == 60 and "drift" in x]
    ax[2].plot([x["k"] for x in xs], [x["drift"]["mu"] for x in xs], "o-", color=c, label=f"{name}: mobility mu")
    ax[2].plot([x["k"] for x in xs], [2 * x["pair_force_coeff"] for x in xs], "s--", color=c,
               label=f"{name}: pair v_rel x d^2 = 2 mu M / 4 pi")
ax[2].axhline(0, color="k", lw=0.5)
ax[2].set(xlabel="k (consumption rate)", title="mu < 0: moves toward lower G -> 1/d^2 attraction")
ax[2].legend(fontsize=7, frameon=False)
fig.suptitle("EX008 (EXPLORATORY): designed long-range test model, 3D", fontsize=10)
fig.tight_layout()
fig.savefig(os.path.join(ROOT, "paper", "figures", "EX008_longrange.png"), dpi=130)

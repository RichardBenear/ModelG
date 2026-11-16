"""EX003 figure: G field of the 3D particle as q -> g p (eps -> 0)."""
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
R = os.path.join(ROOT, "results", "EX003_screening")
fig, axs = plt.subplots(1, 2, figsize=(9, 3.6), dpi=150)
cols = ["#1f6fb4", "#2a9d8f", "#e9a23b", "#c2571a", "#7b3294"]
for eps, c in zip(("1", "0.3", "0.1", "0.03", "0.01"), cols):
    d = np.load(os.path.join(R, f"profile_eps{eps}.npz"))
    r, phiG = d["r"], d["u"][0] - 140 / 9
    axs[0].plot(r, phiG, color=c, lw=1, label=f"eps = {eps}")
    axs[1].semilogy(r, np.abs(r * phiG) + 1e-30, color=c, lw=0.8, label=f"eps = {eps}")
axs[0].set_xlim(0, 8)
axs[0].axhline(0, color="0.6", lw=0.6)
axs[0].set_xlabel("r")
axs[0].set_ylabel("phi_G")
axs[0].legend(fontsize=7)
axs[1].set_xlim(0, 120)
axs[1].set_ylim(1e-10, 1)
axs[1].set_xlabel("r")
axs[1].set_ylabel("|r * phi_G|  (long-range halo is G EXCESS)")
axs[1].legend(fontsize=7)
fig.suptitle("EX003: 3D particle G field as q -> g p  (q = 0.1 + 0.9 eps, a = 14 eps)", fontsize=10)
fig.tight_layout()
fig.savefig(os.path.join(ROOT, "paper", "figures", "EX003_screening.png"))
print("ok")

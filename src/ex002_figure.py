"""EX002 figure: radial G perturbation of the 3D particle, standard and broadened-core variant."""
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
R = os.path.join(ROOT, "results", "EX002_far_field_3d")
G0 = 14.0 / 0.9
fig, axs = plt.subplots(1, 2, figsize=(9, 3.6), dpi=150)
for name, lab, col in (("std", "standard Model G", "#1f6fb4"), ("delta", "broadened core (P&L eq. 24-25)", "#c2571a")):
    d = np.load(os.path.join(R, f"profile_R60_neumann_{name}.npz"))
    r, phiG = d["r"], d["u"][0] - G0
    axs[0].plot(r, phiG, color=col, lw=1, label=lab)
    axs[1].semilogy(r, np.abs(r * phiG) + 1e-30, color=col, lw=0.8, label=lab)
rr = np.linspace(2, 40, 100)
axs[1].semilogy(rr, 0.5 * np.exp(-0.5141 * (rr - 2)), "k--", lw=0.8, label="exp(-0.514 r) (linear theory)")
axs[1].semilogy(rr, np.full_like(rr, 0.3), "k:", lw=0.8, label="1/r field would be flat here")
axs[0].set_xlim(0, 12)
axs[0].axhline(0, color="0.6", lw=0.6)
axs[0].set_xlabel("r")
axs[0].set_ylabel("phi_G = G - G0")
axs[0].legend(fontsize=7)
axs[1].set_xlim(0, 40)
axs[1].set_ylim(1e-12, 10)
axs[1].set_xlabel("r")
axs[1].set_ylabel("|r * phi_G|")
axs[1].legend(fontsize=6.5, loc="lower left")
fig.suptitle("EX002: 3D Model G particle (a = 14, b = 29), G field", fontsize=10)
fig.tight_layout()
fig.savefig(os.path.join(ROOT, "paper", "figures", "EX002_far_field_3d.png"))
print("ok")

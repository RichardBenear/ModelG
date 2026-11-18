"""EX005 figure: Model G-2 3D particle (kappa = 0): G hill in the core, G deficit in the far tail."""
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
R = os.path.join(ROOT, "results", "EX005_model_g2")
cols = {0.001: "#3b6fb6", 0.005: "#2e9e6b", 0.01: "#c0504d", 0.05: "#8e5ea2"}
fig, ax = plt.subplots(1, 3, figsize=(14, 4))
for s, c in cols.items():
    d = np.load(os.path.join(R, f"profile_k0_s{s:g}.npz"))
    r, phi = d["r"], d["u"] - d["h0"][:, None]
    ax[0].plot(r, phi[0], color=c, label=f"sigma = {s}")
    ax[1].plot(r, phi[2] / np.abs(phi[2]).max(), color=c, lw=1, ls="--")
    ax[1].plot(r, phi[0] / np.abs(phi[0]).max(), color=c, lw=1.5)
    y = r * phi[0]
    ax[2].semilogy(r, np.where(y < 0, -y, np.nan), color=c, label=f"sigma = {s}: r·phi_G < 0 (deficit)")
    ax[2].semilogy(r, np.where(y > 0, y, np.nan), color=c, lw=0.6, alpha=0.5)
ax[0].set(xlim=(0, 10), xlabel="r", ylabel="phi_G", title="Core: a G hill (net G > 0)")
ax[1].set(xlim=(0, 10), xlabel="r", ylabel="normalized", title="phi_G (solid) tracks phi_Y (dashed)")
ax[2].set(xlim=(0, 85), ylim=(1e-6, 1e3), xlabel="r", ylabel="|r phi_G|",
          title="Far tail: deficit, decay rate ~ sqrt(sigma)\n(thin lines: positive parts)")
ax[0].legend(fontsize=8, frameon=False)
ax[2].legend(fontsize=7, frameon=False)
for a in ax[:2]:
    a.axhline(0, color="k", lw=0.5)
fig.suptitle("EX005 (EXPLORATORY): Model G-2 particle in 3D, B -> G irreversible (kappa = 0)", fontsize=10)
fig.tight_layout()
fig.savefig(os.path.join(ROOT, "paper", "figures", "EX005_model_g2.png"), dpi=130)

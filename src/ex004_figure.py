"""EX004 figure: G field of a growing 3D structure (eps = 0.03, 0.01) and its net G over time."""
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
R = os.path.join(ROOT, "results", "EX004_growth")
fig, ax = plt.subplots(1, 2, figsize=(11, 4))
for eps, c in ((0.03, "#3b6fb6"), (0.01, "#c0504d")):
    d = np.load(os.path.join(R, f"growth_eps{eps:.3g}.npz"))
    r, phi, A = d["r"], d["final"], d["series"]
    w = int(round(3.066 / (r[1] - r[0])))
    sm = np.convolve(phi[0], np.ones(w) / w, mode="same")
    ax[0].plot(r, sm, color=c, label=f"eps = {eps}: phi_G averaged over one wavelength")
    ax[0].axvline(A[-1, 1], color=c, ls=":", lw=1)
    ax[1].plot(A[:, 0], A[:, 2], color=c, label=f"eps = {eps}: net G (whole domain)")
    ax[1].plot(A[:, 0], A[:, 5], color=c, ls="--", label=f"eps = {eps}: G beyond front + 2")
ax[0].axhline(0, color="k", lw=0.5)
ax[0].set(xlim=(0, 130), xlabel="r", ylabel="phi_G (smoothed)", title="t = 400: deficit inside, excess just ahead of front (dotted)")
ax[1].axhline(0, color="k", lw=0.5)
ax[1].set(xlabel="t", ylabel="integral of phi_G", title="Net G deficit grows; surroundings carry an excess")
for a in ax:
    a.legend(fontsize=7, frameon=False)
fig.suptitle("EX004 (EXPLORATORY): growing Model G structure in 3D, near-marginal ether", fontsize=10)
fig.tight_layout()
fig.savefig(os.path.join(ROOT, "paper", "figures", "EX004_growth.png"), dpi=130)

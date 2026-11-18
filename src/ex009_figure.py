"""EX009 figure: direct two-particle approach speed (axisymmetric) vs the 1/d^2 prediction; small-k scaling."""
import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
R = os.path.join(ROOT, "results", "EX009_followups")
a = json.load(open(os.path.join(R, "partA_pair2d.json")))
b = json.load(open(os.path.join(R, "partB_smallk.json")))
fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
d = np.array([p["d"] for p in a["pairs"]])
v = np.array([p["v"] for p in a["pairs"]])
dd = np.linspace(7, 32, 100)
ax[0].loglog(d[v < 0], -v[v < 0], "o", color="#3b6fb6", label="direct (r, z) solve: approach speed")
ax[0].loglog(d[v > 0], v[v > 0], "x", color="#c0504d", ms=8, label="d = 8: repelled (short-range Turing force)")
ax[0].loglog(dd, -a["mu_radial"] * a["M_radial"] / (4 * np.pi * dd ** 2), "--", color="k",
             label="prediction |mu| M / (4 pi d^2) from EX008")
ax[0].set(xlabel="separation d", ylabel="speed of each particle", title="Two particles: 1/d^2 attraction confirmed (d >= 15)")
ax[0].legend(fontsize=8, frameon=False)
k = np.array([x["k"] for x in b])
ax[1].loglog(k, [x["M"] for x in b], "o-", label="M (consumption) ~ k")
ax[1].loglog(k, [-x["mu"] for x in b], "s-", label="|mu| (drift) ~ k")
ax[1].loglog(k, [-x["pair_coeff"] for x in b], "^-", label="pair force coeff ~ k^2")
ax[1].set(xlabel="consumption rate k", title="Smaller k: weaker force, range stays infinite (1/r)")
ax[1].legend(fontsize=8, frameon=False)
fig.suptitle("EX009 (EXPLORATORY): EX008 test model follow-ups", fontsize=10)
fig.tight_layout()
fig.savefig(os.path.join(ROOT, "paper", "figures", "EX009_followups.png"), dpi=130)

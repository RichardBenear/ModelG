"""MG195 figure: |v| vs d for P1-P3, C1, C2 with predictions; ratio panel."""
import json, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
o = json.load(open(os.path.join(ROOT, "results", "MG195_inverse_square", "results.json")))
P = {(r["point"], r["d"]): r for r in o["pairs"] if r["h"] == 0.1 and r["far"] == 150.0}
DS = np.array([18, 22, 27, 33, 40.])
dd = np.linspace(16, 42, 100)
fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
cols = dict(P1="#1f77b4", P2="#2ca02c", P3="#9467bd", C2="#d62728", C1="#7f7f7f")
lab = dict(P1="P1 G feeds X, k=0.2", P2="P2 G feeds Y, k=0.2", P3="P3 G feeds X, k=0.05",
           C2="C2 screened, lambda=0.01", C1="C1 null k=0 (|v|)")
for n in ("P1", "P2", "P3", "C2"):
    rad = o["radial"][n]
    kap = np.sqrt(0.01) if n == "C2" else 0.0
    ax[0].plot(dd, abs(rad["mu"] * rad["M"] / (4 * np.pi)) * np.exp(-kap * dd) * (1 / dd ** 2 + kap / dd), color=cols[n], lw=1)
    ax[0].plot(DS, [abs(P[(n, d)]["v"]) for d in DS], "o", color=cols[n], label=lab[n])
    ax[1].plot(DS, [P[(n, d)]["v"] / P[(n, d)]["prediction"] for d in DS], "o-", color=cols[n], label=n)
ax[0].plot(DS, [abs(P[("C1", d)]["v"]) for d in DS], "s", color=cols["C1"], label=lab["C1"])
ax[0].axhline(o["verdicts"]["floor"], color="k", ls=":", lw=1, label="null floor 10 max|v_C1|")
ax[0].set(xscale="log", yscale="log", xlabel="separation d", ylabel="|v| (approach speed)", title="Pair speed vs prediction (lines)")
ax[0].legend(fontsize=7)
ax[1].axhspan(0.95, 1.05, color="0.9")
ax[1].axhline(1, color="k", lw=0.8)
ax[1].set(xlabel="separation d", ylabel="v / prediction", title="Ratio (band = registered 5%, d >= 22)")
ax[1].legend(fontsize=8)
fig.tight_layout()
fig.savefig(os.path.join(ROOT, "paper", "figures", "MG195_inverse_square.png"), dpi=130)

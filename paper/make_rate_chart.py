"""Chart for the plain-language summary: spontaneous creation rate vs Omega (MG3 and MG3b). Data from the records:
MG3 (GATES.md, MG3 result) and MG3b (results/MG3b/MG3b_evaluation.md). Run: python paper/make_rate_chart.py"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import chi2
from scipy.optimize import minimize

L = 100.0
mg3 = dict(om=[150, 120, 100], k=[4.05e-6, 4.63e-5, 1.58e-4], n=[5, 10, 11])
mg3b = dict(om=[130, 120, 110, 100], n=[25, 44, 28, 21], E=[15461, 12892, 4464, 1866])
mg3b["k"] = [n / (L * e) for n, e in zip(mg3b["n"], mg3b["E"])]

def garwood(n, k):
    lo = chi2.ppf(0.025, 2 * n) / 2 / n * k
    hi = chi2.ppf(0.975, 2 * n + 2) / 2 / n * k
    return lo, hi

# Poisson GLM log k = a - S (Omega - mean), as the MG3b slope fit
om = np.array(mg3b["om"], float); n = np.array(mg3b["n"]); E = np.array(mg3b["E"]) * L; c = om.mean()
nll = lambda p: -(n * (p[0] - p[1] * (om - c)) - E * np.exp(p[0] - p[1] * (om - c))).sum()
a, S = minimize(nll, [np.log(n.sum() / E.sum()), 0.07]).x
print("S =", S)

BLUE, ORANGE, INK, MUTED, SURF = "#2a78d6", "#eb6834", "#0b0b0b", "#52514e", "#fcfcfb"
fig, ax = plt.subplots(figsize=(7.2, 4.4), dpi=150)
fig.patch.set_facecolor(SURF); ax.set_facecolor(SURF)
x = np.linspace(95, 155, 100)
ax.plot(x, np.exp(a - S * (x - c)), color=MUTED, lw=1.5, ls="--", zorder=1)
ax.text(153, 1.3e-4, f"fit: rate halves every {np.log(2)/S:.0f} units of Ω\n(S = {S:.3f})",
        color=MUTED, fontsize=8.5, ha="right")
for d, col, mk, lab, dx in ((mg3, BLUE, "o", "MG3 (5–11 events per point)", -1.2),
                            (mg3b, ORANGE, "s", "MG3b (21–44 events per point)", 1.2)):
    lo, hi = zip(*[garwood(nn, kk) for nn, kk in zip(d["n"], d["k"])])
    xs = np.array(d["om"]) + dx
    ax.errorbar(xs, d["k"], yerr=[np.array(d["k"]) - lo, np.array(hi) - d["k"]], fmt=mk, color=col, ms=7,
                mec=SURF, mew=1.5, elinewidth=2, capsize=0, label=lab, zorder=3)
ax.set_yscale("log")
ax.set_xlabel("Ω  (tokens per unit concentration: larger Ω = weaker noise)", color=INK, fontsize=9.5)
ax.set_ylabel("creation rate\n(per unit length per unit time)", color=INK, fontsize=9.5)
ax.set_title("Particles appearing on their own: the rate falls steeply as noise weakens", color=INK, fontsize=11,
             loc="left")
ax.grid(True, which="major", color="#e4e3df", lw=0.8); ax.set_axisbelow(True)
for s in ("top", "right"): ax.spines[s].set_visible(False)
for s in ("left", "bottom"): ax.spines[s].set_color("#b8b7b0")
ax.tick_params(colors=MUTED, labelsize=8.5)
ax.legend(frameon=False, fontsize=8.5, loc="lower left", labelcolor=INK)
ax.text(0.0, -0.2, "b = 29, 1D ring, h = 0.2. Bars: 95% Poisson intervals (Garwood). Fit: MG3b's four points. Points nudged ±1 in Ω so they do not overlap.",
        transform=ax.transAxes, fontsize=7.5, color=MUTED)
fig.tight_layout()
fig.savefig("paper/creation_rate_vs_omega.png", facecolor=SURF)

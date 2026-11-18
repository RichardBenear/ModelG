"""Illustration: Model G reaction scheme and the 3D particle it forms (our EX002 simulation, a = 14, b = 29).
Panels: (a) reaction network, (b) radial concentration profiles, (c) cross-sections through the particle.
Drawn in the spirit of LaViolette's particle illustrations; all data are from our own solve, not from the book."""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
from matplotlib.colors import TwoSlopeNorm

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
d = np.load(os.path.join(ROOT, "results", "EX002_far_field_3d", "profile_R60_neumann_std.npz"))
r, (G, X, Y) = d["r"], d["u"]
X0 = 14.0 / 0.9
h0 = dict(G=X0, X=X0, Y=29.0 / X0)
phi = dict(G=G - h0["G"], X=X - h0["X"], Y=Y - h0["Y"])
col = dict(G="#2a9d5c", X="#c0392b", Y="#2c6fbb")

fig = plt.figure(figsize=(13, 8.6), dpi=150)
gs = fig.add_gridspec(2, 3, height_ratios=[1, 1.05], hspace=0.32, wspace=0.28)

# (a) reaction network
ax = fig.add_subplot(gs[0, 0:2]); ax.set_xlim(0, 10); ax.set_ylim(0, 4.6); ax.axis("off")
def node(x, y, s, c, big=True):
    ax.add_patch(FancyBboxPatch((x - .38, y - .3), .76, .6, boxstyle="round,pad=0.05", fc=c if big else "white",
                                ec=c, lw=1.6, alpha=0.9 if big else 1))
    ax.text(x, y, s, ha="center", va="center", fontsize=15, weight="bold", color="white" if big else c)
def arrow(p, q, c="0.25", rad=0.0, lab=None, lp=None, two=False):
    ax.add_patch(FancyArrowPatch(p, q, arrowstyle="<|-|>" if two else "-|>", mutation_scale=14, lw=1.5, color=c,
                                 connectionstyle=f"arc3,rad={rad}", shrinkA=20, shrinkB=20))
    if lab: ax.text(*lp, lab, fontsize=9, color=c, ha="center", va="center", style="italic")
node(1.0, 3.6, "A", "0.45", False); node(3.2, 3.6, "G", col["G"]); node(5.6, 3.6, "X", col["X"])
node(8.6, 3.6, "Ω", "0.45", False); node(5.6, 1.2, "Y", col["Y"]); node(3.2, 1.2, "B", "0.45", False)
node(8.6, 1.2, "Z", "0.45", False)
arrow((1.0, 3.6), (3.2, 3.6), two=True, lab="k1", lp=(2.1, 3.95))
arrow((3.2, 3.6), (5.6, 3.6), two=True, lab="k2", lp=(4.4, 3.95))
arrow((5.6, 3.6), (8.6, 3.6), two=True, lab="k5", lp=(7.1, 3.95))
arrow((5.6, 3.6), (5.6, 1.2), rad=0.35, c=col["Y"], lab="B + X → Y + Z  (k3)", lp=(4.15, 2.4))
arrow((5.6, 1.2), (5.6, 3.6), rad=0.35, c=col["X"], lab="2X + Y → 3X  (k4)\nautocatalysis", lp=(7.6, 2.75))
arrow((3.2, 1.2), (5.15, 2.3), c="0.45", rad=-0.15); arrow((5.25, 2.0), (8.6, 1.2), c="0.45", rad=0.25)
ax.text(0.05, 4.45, "(a) Model G reaction network (LaViolette, Subquantum Kinetics)", fontsize=11, weight="bold")
ax.text(0.05, 0.25, "A, B, Ω, Z: reservoirs held fixed.   G, X, Y: the three reacting ether variables, which diffuse "
        "(D_G = D_X = 1, D_Y = 12).\nX feeds Y, Y feeds back to X: an activator-inhibitor loop. G is the upstream "
        "supply of X.", fontsize=8.5, color="0.25", va="bottom")

# (b) radial profiles
axb = fig.add_subplot(gs[0, 2])
for s in ("X", "Y", "G"):
    axb.plot(r, phi[s], color=col[s], lw=1.6, label=f"φ_{s} = {s} − {s}₀")
axb.axhline(0, color="0.6", lw=0.6); axb.set_xlim(0, 14)
axb.set_xlabel("distance from particle centre r"); axb.set_ylabel("concentration departure from ether")
axb.set_title("(b) radial profiles of the particle", fontsize=11, weight="bold", loc="left")
axb.legend(fontsize=8.5, frameon=False)
axb.annotate("core: X well,\nY peak, G well", xy=(0.3, phi["X"][6]), xytext=(3.0, -9.0), fontsize=8.5,
             arrowprops=dict(arrowstyle="->", color="0.4"))
axb.annotate("X-rich shell, then\ndamped Turing shells", xy=(2.2, phi["X"].max()), xytext=(5.5, -4.5), fontsize=8.5,
             arrowprops=dict(arrowstyle="->", color="0.4"))

# (c) cross-sections: slice through centre, z = 0 plane
L, n = 14.0, 401
xx = np.linspace(-L, L, n); XX, ZZ = np.meshgrid(xx, xx); RR = np.hypot(XX, ZZ)
for k, s in enumerate(("X", "Y", "G")):
    axc = fig.add_subplot(gs[1, k])
    f = np.interp(RR, r, phi[s])
    f = np.sign(f) * np.abs(f) ** 0.3          # compressed scale so the faint outer shells show
    m = np.abs(f).max()
    im = axc.imshow(f, extent=(-L, L, -L, L), origin="lower", cmap="RdBu_r", norm=TwoSlopeNorm(0, -m, m))
    axc.set_title(f"(c{k+1}) slice through particle: φ_{s}  (scale |φ|^0.3)", fontsize=10.5, weight="bold", loc="left")
    axc.set_xlabel("x"); axc.set_ylabel("z") if k == 0 else axc.set_yticklabels([])
    cb = fig.colorbar(im, ax=axc, fraction=0.046, pad=0.03); cb.ax.tick_params(labelsize=7)
    cb.set_ticks([]); cb.set_label("below ether  ←  →  above ether", fontsize=7)
fig.text(0.5, 0.015, "Data: our steady 3D radially symmetric Model G particle (EX002; a = 14, b = 29, g = 0.1, "
         "p = q = 1, D_Y = 12). Red = above the ether value, blue = below. Net G of the particle is exactly zero "
         "(G-deficit core balanced by a G-excess shell).\nColour in (c) uses sign(φ)·|φ|^0.3 so the weak outer shells are visible.", ha="center", fontsize=8.5, color="0.3")
fig.savefig(os.path.join(ROOT, "paper", "figures", "ModelG_particle_dynamics.png"), bbox_inches="tight")
print("ok", {s: (float(phi[s].min()), float(phi[s].max())) for s in phi})

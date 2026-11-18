"""Side-by-side: original Model G particle (EX002) and the conserved-G 'Model Gc' particle (EX008, G feeds X, k = 0.2).
All data from our own solves."""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
col = dict(G="#2a9d5c", X="#c0392b", Y="#2c6fbb", W="#8e44ad")

d = np.load(os.path.join(ROOT, "results", "EX002_far_field_3d", "profile_R60_neumann_std.npz"))
r1 = d["r"]; G, X, Y = d["u"]; X0 = 14.0 / 0.9
m1 = dict(G=G - X0, X=X - X0, Y=Y - 29.0 / X0)
d = np.load(os.path.join(ROOT, "results", "EX008_longrange", "particle_GX_k0.2_R120.npz"))
r2 = d["r"]; X, Y, W, G = d["u"]; A, B = 140 / 9, 29.0
m2 = dict(X=X - A, Y=Y - B / A, W=W, G=G - 1.0)

fig = plt.figure(figsize=(13, 15), dpi=140)
gs = fig.add_gridspec(4, 4, height_ratios=[0.55, 1, 1, 1], hspace=0.42, wspace=0.35)
titles = ("Original Model G (LaViolette)", "Model Gc (our conserved-G test model, EX008)")
rx = ("A ⇌ G ⇌ X ⇌ Ω\nB + X → Y + Z\n2X + Y → 3X\n\nG is made and lost everywhere in the ether\n"
      "(−qG + gX), so the ether rebalances G\naround the particle: net G = 0, screened field.",
      "A → X,  X → Ω\nB + X → Y + Z\n2X + Y → 3X\nW: core-only, grows on Y (Y W² gain;\n   W, W³ losses), W = 0 in the ether\nG + W → W + X   (rate k G W)\n\n"
      "G has no source or loss in the ether, only\nsupply from far away: consumption in the core\n"
      "leaves an unscreened 1/r G well.")
for c in range(2):
    ax = fig.add_subplot(gs[0, 2 * c:2 * c + 2]); ax.axis("off")
    ax.text(0.0, 1.0, titles[c], fontsize=13, weight="bold", va="top")
    ax.text(0.0, 0.72, rx[c], fontsize=10, va="top", family="DejaVu Sans")

# row 1: radial profiles
for c, (r, m, sp) in enumerate(((r1, m1, ("X", "Y", "G")), (r2, m2, ("X", "Y", "W", "G")))):
    ax = fig.add_subplot(gs[1, 2 * c:2 * c + 2])
    for s in sp:
        lab = "W (core species)" if s == "W" else f"φ_{s} = {s} − {s}₀" + (" (×10)" if s == "G" else "")
        ax.plot(r, m[s] * (10 if s == "G" else 1), color=col[s], lw=1.6, label=lab)
    ax.axhline(0, color="0.6", lw=0.6); ax.set_xlim(0, 14)
    ax.set_xlabel("distance from particle centre r"); ax.set_ylabel("departure from ether value")
    ax.set_title("radial profiles", fontsize=10.5, loc="left", weight="bold")
    ax.legend(fontsize=8, frameon=False, loc="lower right")

# row 2: slices of X and G (compressed colour scale)
L, n = 14.0, 301
xx = np.linspace(-L, L, n); XX, ZZ = np.meshgrid(xx, xx); RR = np.hypot(XX, ZZ)
for c, (r, m) in enumerate(((r1, m1), (r2, m2))):
    for k, s in enumerate(("X", "G")):
        ax = fig.add_subplot(gs[2, 2 * c + k])
        f = np.interp(RR, r, m[s]); f = np.sign(f) * np.abs(f) ** 0.3; mm = np.abs(f).max()
        ax.imshow(f, extent=(-L, L, -L, L), origin="lower", cmap="RdBu_r", norm=TwoSlopeNorm(0, -mm, mm))
        ax.set_title(f"slice: φ_{s}", fontsize=10.5, loc="left", weight="bold")
        ax.set_xticks([-10, 0, 10]); ax.set_yticks([-10, 0, 10])

# row 3: far G field, r * |G - G_inf|
for c, (r, m, R) in enumerate(((r1, m1, 60), (r2, m2, 120))):
    ax = fig.add_subplot(gs[3, 2 * c:2 * c + 2])
    y = np.abs(r * m["G"]); ax.semilogy(r, y + 1e-30, color=col["G"], lw=1.5, label="r · |G − G₀|" if c == 0 else "r · |G − G∞|")
    if c == 0:
        rr = np.linspace(2, 40, 50); ax.semilogy(rr, 0.5 * np.exp(-0.514 * (rr - 2)), "k--", lw=0.8,
                                                 label="exp(−0.514 r): screened")
        ax.set_xlim(0, 45); ax.set_ylim(1e-12, 10)
        ax.set_title("far G field: decays exponentially (short range)", fontsize=10.5, loc="left", weight="bold")
    else:
        M = 5.3178
        rr = np.linspace(5, 100, 50); ax.semilogy(rr, M / (4 * np.pi) * (1 - rr / R), "k--", lw=0.8,
                                                  label="M/4π · (1 − r/R): pure 1/r well (R = 120 boundary)")
        ax.set_xlim(0, 110); ax.set_ylim(1e-3, 3)
        ax.set_title("far G field: follows 1/r out to the boundary (long range)", fontsize=10.5, loc="left", weight="bold")
    ax.set_xlabel("r"); ax.legend(fontsize=8, frameon=False, loc="lower left")

fig.text(0.5, 0.04, "Data: our steady 3D radial solves. Left: Model G, a = 14, b = 29 (EX002). Right: Model Gc, Brusselator X/Y "
         "(A = 140/9, B = 29), core species W, conserved G with G feeding X at k = 0.2 (EX008/MG195 point P1; M = 5.32).\n"
         "Slices use sign(φ)·|φ|^0.3 so faint shells show; red above ether, blue below. Model Gc is our construction, "
         "not LaViolette's.", ha="center", fontsize=8.5, color="0.3")
fig.savefig(os.path.join(ROOT, "paper", "figures", "ModelG_vs_Gc_particles.png"), bbox_inches="tight")
print("ok", {s: (float(m2[s].min()), float(m2[s].max())) for s in m2})

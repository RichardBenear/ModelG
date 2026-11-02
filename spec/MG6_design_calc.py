"""MG6 design calculation (descriptive; no dynamics, no MG6 run). Uses MG4c's stored b = 28 start states.
Expected node-count excess within radius R of the particle centre if the local node density were rho * exp(kappa * phi_G),
against the Poisson counting sd of a disc. Run from the repo root: python spec/MG6_design_calc.py
"""
import json, os, sys
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
import run_mg4c as r4

LX, LY = 19.92, 20.08
dec = json.load(open(r4.DEC))
for slot, pl in [(0, 0), (0, 1), (1, 0), (1, 1)]:
    g, h, c = r4.graph_for(dict(kind="D", slot=slot, pl=pl), dec)
    G, X, Y = np.load(f"results/MG4c/fields/b28.000000_D{slot}_p{pl}.npz")["final"]
    pos = np.asarray(g["pos"]); cen = pos[np.argmin(X)]
    d = pos - cen
    d[:, 0] -= LX * np.round(d[:, 0] / LX); d[:, 1] -= LY * np.round(d[:, 1] / LY)
    rr = np.hypot(*d.T); rho = g["N"] / (LX * LY); dA = LX * LY / g["N"]
    print(f"D{slot}_p{pl}: rho {rho:.2f}, phi_G min {G.min():.3f} max {G.max():.3f}")
    for kap in (1, 2, 4):
        row = []
        for R in (1.5, 3, 5, 8):
            m = rr < R
            row.append(f"R={R}: {rho*dA*np.sum(np.exp(kap*G[m])-1):+.1f} (sd {np.sqrt(rho*np.pi*R*R):.1f})")
        print(f"  kappa {kap}: " + "; ".join(row))

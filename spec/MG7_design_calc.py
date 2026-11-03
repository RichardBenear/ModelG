"""MG7 design calculation (descriptive; no dynamics, no MG7 run). Static prediction of the coupled density response from
MG4c's stored settled states at several b on the four D start states: P0(1.5) at kappa = 2 (nodes), the G-well integral
M- = integral of max(-phi_G, 0), and the first zero crossing of the azimuthal mean phi_G. Node-sum approximation (each node
carries area A/N), centred on the minimum of phi_X. Run from the repo root: python spec/MG7_design_calc.py
"""
import json, os, sys
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
import run_mg4c as r4

dec = json.load(open(r4.DEC))
BS = ("27.500000", "27.750000", "28.000000", "28.250000", "28.500000")
print("b       start   P0(1.5)k=2  P0(3)k=2   M-      r0(phi_G=0)  stationary")
tab = {}
for b in BS:
    for slot, pl in [(0, 0), (0, 1), (1, 0), (1, 1)]:
        g, h, c = r4.graph_for(dict(kind="D", slot=slot, pl=pl), dec)
        G, X, Y = np.load(f"results/MG4c/fields/b{b}_D{slot}_p{pl}.npz")["final"]
        LX, LY = g["Lx"], g["Ly"]
        pos = np.asarray(g["pos"]); cen = pos[np.argmin(X)]
        d = pos - cen
        d[:, 0] -= LX * np.round(d[:, 0] / LX); d[:, 1] -= LY * np.round(d[:, 1] / LY)
        rr = np.hypot(*d.T); dA = LX * LY / g["N"]; rho = g["N"] / (LX * LY)
        P = {R: rho * dA * np.sum(np.expm1(2 * G[rr < R])) for R in (1.5, 3.0)}
        Mm = dA * np.sum(np.maximum(-G, 0))
        edges = np.arange(0, 4.01, 0.25); mid = 0.5 * (edges[1:] + edges[:-1])
        prof = np.array([G[(rr >= a) & (rr < e)].mean() for a, e in zip(edges[:-1], edges[1:])])
        i = np.flatnonzero((prof[:-1] < 0) & (prof[1:] >= 0))
        r0 = float(mid[i[0]] + (mid[i[0] + 1] - mid[i[0]]) * (-prof[i[0]]) / (prof[i[0] + 1] - prof[i[0]])) if len(i) else None
        st = json.load(open(f"results/MG4c/runs/b{b}_D{slot}_p{pl}.json")).get("stationary")
        tab.setdefault(b, []).append((P[1.5], Mm))
        print(f"{float(b):.3f}  D{slot}_p{pl}   {P[1.5]:+8.1f}   {P[3.0]:+8.1f}   {Mm:6.3f}   {r0 if r0 is None else round(r0, 2)}   {st}")
for b, v in tab.items():
    v = np.array(v)
    print(f"b = {float(b):.3f}: mean P0(1.5) {v[:, 0].mean():+.1f} (sd {v[:, 0].std(ddof=1):.1f}), mean M- {v[:, 1].mean():.3f}")

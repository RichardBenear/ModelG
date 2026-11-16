"""MG192 post-hoc diagnostic (not a registered test): core depth of the soliton in the dynamic runs.
Re-runs the four eps>0 dynamic cases to T=60 and records the X-core depth (min X minus local background X)."""
import json, os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import mg192_mobility as m
OUT = m.OUT
d = np.load(os.path.join(OUT, "soliton_adjoint_L100_N1024.npz")); U = d["U"]
g = m.Grid(100.0, 1024); h0 = m.hom(); res = {}
for eps in (0.001, 0.002):
    for x0 in (0.0, 12.5):
        bg, a, _ = m.modulated_background(g, eps)
        u0 = bg + np.array([g.shift(U[s], x0) for s in range(3)]) - h0[:, None]
        st = m.ETDRK4(g, 0.002); vh = np.fft.fft(u0, axis=1); rows = []
        for n in range(0, 30001):
            if n % 1250 == 0:
                u = np.real(np.fft.ifft(vh, axis=1)); j = int(np.argmin(u[1]))
                rows.append([n * 0.002, float(g.x[j]), float(u[1, j] - bg[1, j]), float(a[j])])
            if n < 30000: vh = st.step(vh, a)
        res[f"eps{eps}_x0{x0}"] = rows
        print(eps, x0, rows[0], rows[-1], flush=True)
json.dump(res, open(os.path.join(OUT, "posthoc_core_depth.json"), "w"), indent=1)

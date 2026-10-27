"""MG002 driver: MG1 null controls N1/N2 and integral-balance diagnostic. Writes results/MG002/controls.json."""
import json, sys
import numpy as np
sys.path.insert(0, "src")
from mg_continuum import run, analyze

H = 0.025
out = {"N1": [], "N2": [], "N2_threshold": {}, "diag": []}
for n in (1, 2, 3):
    g, P, s, rhs = run(n, h=H, seed_scale=0.0)
    m = float(np.abs(s.y[:, -1]).max())
    out["N1"].append(dict(n=n, max_abs_phi=m, first_step_retry=s.first_step_used))
    print(f"N1 n={n}: max|phi|={m:.2e} (first_step={s.first_step_used})", flush=True)
    g, P, s, rhs = run(n, h=H, seed_scale=1e-3)
    m = float(np.abs(s.y[:, -1]).max())
    out["N1"].append(dict(n=n, seed_scale=1e-3, max_abs_phi=m))
    print(f"N1b n={n} seed 1e-3: max|phi|(t=100)={m:.2e}", flush=True)
    def nucleates(sc):
        g, P, s, rhs = run(n, h=H, seed_scale=sc, t_end=100)
        return float(s.y[g.N, -1])  # phi_X(0)
    for sc in (0.25, 0.5, 0.75):
        x0 = nucleates(sc)
        out["N2"].append(dict(n=n, scale=sc, phiX0=x0, nucleated=bool(x0 < -1.0)))
        print(f"N2 n={n} scale={sc}: phiX(0,t=100)={x0:.4f}", flush=True)
    lo, hi = 0.0, 1.0
    for _ in range(12):
        mid = 0.5 * (lo + hi)
        if nucleates(mid) < -1.0: hi = mid
        else: lo = mid
    out["N2_threshold"][n] = [lo, hi]; print(f"N2 n={n}: threshold scale in [{lo:.4f}, {hi:.4f}] (judged at t=100)", flush=True)
    # near-threshold check: does a sub-threshold seed nucleate later than t=100?
    g, P, s, rhs = run(n, h=H, seed_scale=lo * 0.99, t_end=400)
    x0 = float(s.y[g.N, -1]); out["N2"].append(dict(n=n, scale=lo * 0.99, t=400, phiX0=x0))
    print(f"N2 n={n} scale={lo*0.99:.4f} at t=400: phiX(0)={x0:.4f}", flush=True)

for n in (1, 2, 3):
    for R in (50.0, 100.0, 200.0):
        for t_end in (100.0, 300.0):
            g, P, s, rhs = run(n, R=R, h=H, t_end=t_end)
            r, _ = analyze(g, P, s, rhs)
            it = r["integrals"]
            out["diag"].append(dict(n=n, R=R, t=t_end, **it))
            print(f"diag n={n} R={R:.0f} t={t_end:.0f}: I_G={it['I']['G']:+.3e} I_X={it['I']['X']:+.3e} "
                  f"I_Y={it['I']['Y']:+.4f} F_G={it['flux']['G']:+.2e} F_X={it['flux']['X']:+.2e} "
                  f"dIG={it['dIdt']['G']:+.2e} dIsum={sum(it['dIdt'].values()):+.2e} "
                  f"resG={it['balance_G_residual']:+.1e} resS={it['balance_sum_residual']:+.1e}", flush=True)
json.dump(out, open("results/MG002/controls.json", "w"), indent=1)

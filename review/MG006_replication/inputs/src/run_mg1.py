"""MG002 driver: MG1 main runs, grid ladder, tolerance check. Writes results/MG002/main.json."""
import json, sys, time
sys.path.insert(0, "src")
from mg_continuum import run, analyze

cfg = json.load(open("configs/MG1_main.json"))
out = []
for n in cfg["dims"]:
    for h in cfg["h_ladder"]:
        for t_end in (cfg["t_end"], cfg["t_supplementary"]):
            t0 = time.time()
            g, P, s, rhs = run(n, R=cfg["R"], h=h, t_end=t_end)
            r, _ = analyze(g, P, s, rhs)
            r["wall_s"] = time.time() - t0; r["rtol"] = 1e-8
            out.append(r)
            f = r["fields"]
            print(f"n={n} h={h} t={t_end}: amp G {f['G']['amp']:.4f} X {f['X']['amp']:.4f} Y {f['Y']['amp']:.4f} "
                  f"r0 {f['G']['r0']:.4f}/{f['X']['r0']:.4f}/{f['Y']['r0']:.4f} lam {r['lambda']:.4f} "
                  f"maxdt {r['max_abs_dphidt']:.1e} ({r['wall_s']:.0f}s)", flush=True)
    g, P, s, rhs = run(n, R=cfg["R"], h=cfg["h_ladder"][-1], t_end=cfg["t_end"], rtol=1e-10, atol=1e-12)
    r, _ = analyze(g, P, s, rhs); r["rtol"] = 1e-10; out.append(r)
    f = r["fields"]
    print(f"n={n} TOLCHECK: amp X {f['X']['amp']:.6f} lam {r['lambda']:.5f}", flush=True)
json.dump(out, open("results/MG002/main.json", "w"), indent=1)

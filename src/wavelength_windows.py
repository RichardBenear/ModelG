"""MG002a: sensitivity of the pooled zero-spacing wavelength to the radial window.
Writes results/MG002/wavelength_windows.txt."""
import sys
sys.path.insert(0, "src")
from mg_continuum import run, wavelength

W = [None, (5, 15), (10, 25), (15, 30), (20, 30)]
lines = ["h=0.0125, t=100 (1D also t=300). 'registered' = MG002 estimator (r<0.6R, core spacing dropped).",
         "dim  t    " + "  ".join(f"{('registered' if w is None else f'{w[0]}<r<{w[1]}'):>11}" for w in W)]
for n, t in ((1, 100), (1, 300), (2, 100), (3, 100)):
    g, P, s, rhs = run(n, h=0.0125, t_end=t)
    y = s.y[:, -1]; N = g.N
    F = dict(G=y[:N], X=y[N:2 * N], Y=y[2 * N:])
    vals = [wavelength(g, F, window=w)[0] for w in W]
    lines.append(f"{n}D  {t:<4d} " + "  ".join(f"{v:11.5f}" for v in vals))
    lines[-1] += f"   spread {100*(max(vals)-min(vals))/3.08:.2f}% of 3.08"
open("results/MG002/wavelength_windows.txt", "w").write("\n".join(lines) + "\n")
print("\n".join(lines))

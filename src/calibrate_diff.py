"""MG003b Part C: calibration of the two-sample difference bound (stats_mg2.diff_bounds) used by V2
and MG3-C2b, at M = 12 and M = 8 per side, with the Part A noise families. Reports coverage of
[lower, upper] for |muA - muB|/|phi*| and the false-equivalence rate (upper < margin while the true
difference >= margin). Writes results/MG003b/calibration_diff.txt."""
import sys
import numpy as np
sys.path.insert(0, "src")
from stats_mg2 import diff_bounds

rng = np.random.default_rng(20261028)
d = 100
u = rng.standard_normal(d); u /= np.linalg.norm(u)
Q, _ = np.linalg.qr(np.column_stack([u, rng.standard_normal((d, 2))]))

def noise(kind, M, sigma):
    if kind == "iso": return rng.standard_normal((M, d)) * sigma / np.sqrt(d)
    if kind == "rank1b": return np.outer(rng.standard_normal(M) * sigma, u)
    if kind == "rank3": return (rng.standard_normal((M, 3)) * sigma / np.sqrt(3)) @ Q.T
    if kind in ("twoP+", "twoP-"):
        z = np.where(rng.random(M) < 0.75, 1 / np.sqrt(3), -np.sqrt(3))
        return np.outer((z if kind == "twoP+" else -z) * sigma, u)
    if kind == "t3":
        g = rng.standard_normal((M, 3)); w = rng.chisquare(3, M) / 3
        return (g / np.sqrt(w)[:, None] / np.sqrt(3.0) * sigma / np.sqrt(3)) @ Q.T

R = 3000; lines = []
for M, margin, deltas, sigmas in ((12, 0.0025, (0.0, 0.001, 0.0025, 0.005), (0.0003, 0.001, 0.003)),
                                  (8, 0.10, (0.0, 0.05, 0.10, 0.20), (0.01, 0.03, 0.1))):
    worst_cov, worst_fe = 1.0, 0.0
    for kind in ("iso", "rank1b", "rank3", "twoP+", "twoP-", "t3"):
        for delta in deltas:
            for s in sigmas:
                cov = fe = 0
                for _ in range(R):
                    A = delta * u + noise(kind, M, s); B = noise(kind, M, s)
                    b = diff_bounds(A, B, 1.0)
                    cov += b["lower"] - 1e-15 <= delta <= b["upper"] + 1e-15
                    fe += (b["upper"] < margin) and (delta >= margin)
                worst_cov = min(worst_cov, cov / R); worst_fe = max(worst_fe, fe / R)
    lines.append(f"M={M} per side, margin {100*margin:.2f}%: worst coverage {worst_cov:.4f}; worst false-equivalence rate {worst_fe:.4f}")
    print(lines[-1], flush=True)
open("results/MG003b/calibration_diff.txt", "w").write("\n".join(lines) + "\n")

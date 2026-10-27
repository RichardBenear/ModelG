"""MG003b: calibration of the MG2 inference layer (src/stats_mg2.py) on synthetic run ensembles with
known bias. Profile dimension d = 100 (|x|<10 at h=0.2). Reference phi* has unit norm, so all errors are
relative. Run profiles m_i = phi* + beta*u + noise_i, noise types:
  iso    isotropic Gaussian, total RMS norm sigma
  rank1b Gaussian along the bias direction u (worst case for norm bounds)
  rank3  Gaussian in a 3-dim subspace containing u
  twoP+  two-point along u: z = +1/sqrt3 (p=3/4), -sqrt3 (p=1/4)  (reviewer's counterexample family)
  twoP-  mirrored two-point
  t3     multivariate Student-t (3 dof) in the rank-3 subspace, RMS sigma
Part A: coverage of [lower, upper] for |bias|, and false-PASS / false-FAIL rates at eps = 1%.
Part B: coverage of the run-bootstrap alpha CI (true alpha = 1); calibrated widening factor.
Writes results/MG003b/calibration.txt, calibration.json."""
import json, sys
import numpy as np
sys.path.insert(0, "src")
from stats_mg2 import bias_bounds, decide, alpha_ci

rng = np.random.default_rng(20261027)
d = 100
phi = rng.standard_normal(d); phi /= np.linalg.norm(phi)
u = rng.standard_normal(d); u -= (u @ phi) * phi; u /= np.linalg.norm(u)
Q, _ = np.linalg.qr(np.column_stack([u, rng.standard_normal((d, 2))]))
KINDS = ["iso", "rank1b", "rank3", "twoP+", "twoP-", "t3"]

def noise(kind, M, sigma):
    if kind == "iso":
        return rng.standard_normal((M, d)) * sigma / np.sqrt(d)
    if kind == "rank1b":
        return np.outer(rng.standard_normal(M) * sigma, u)
    if kind == "rank3":
        return (rng.standard_normal((M, 3)) * sigma / np.sqrt(3)) @ Q.T
    if kind in ("twoP+", "twoP-"):
        z = np.where(rng.random(M) < 0.75, 1 / np.sqrt(3), -np.sqrt(3))
        return np.outer((z if kind == "twoP+" else -z) * sigma, u)
    if kind == "t3":
        g = rng.standard_normal((M, 3)); w = rng.chisquare(3, M) / 3
        return (g / np.sqrt(w)[:, None] / np.sqrt(3.0) * sigma / np.sqrt(3)) @ Q.T
    raise ValueError(kind)

lines, out = [], {"A": [], "B": []}
EPS = 0.01
# ---------------- Part A
R = 3000
worst = dict(cov=1.0, fp=0.0, ff=0.0)
for M in (4, 6, 12):
    for kind in KINDS:
        for beta in (0.0, 0.005, 0.01, 0.02):
            for sigma in (0.001, 0.003, 0.01, 0.03):
                cov = fp = ff = inc = 0
                for _ in range(R):
                    b = bias_bounds(phi + beta * u + noise(kind, M, sigma), phi)
                    cov += b["lower"] - 1e-15 <= beta <= b["upper"] + 1e-15
                    dcs = decide(b, EPS)
                    fp += (dcs == "PASS") and beta >= EPS
                    ff += (dcs == "FAIL") and beta <= EPS
                    inc += dcs == "INCONCLUSIVE"
                rec = dict(M=M, kind=kind, beta=beta, sigma=sigma, coverage=cov / R,
                           false_pass=fp / R, false_fail=ff / R, inconclusive=inc / R)
                out["A"].append(rec)
                worst["cov"] = min(worst["cov"], rec["coverage"])
                worst["fp"] = max(worst["fp"], rec["false_pass"]); worst["ff"] = max(worst["ff"], rec["false_fail"])
lines.append(f"Part A ({R} reps/cell, M in 4/6/12, 6 noise types, bias 0-2%, sigma 0.1-3%):")
lines.append(f"  worst coverage of [lower,upper] = {worst['cov']:.4f}; worst false-PASS rate (true bias >= 1%) = {worst['fp']:.4f}; "
             f"worst false-FAIL rate (true bias <= 1%) = {worst['ff']:.4f}")
for M in (4, 6, 12):
    sub = [r for r in out["A"] if r["M"] == M]
    lines.append(f"  M={M}: min coverage {min(r['coverage'] for r in sub):.4f}; max false-PASS {max(r['false_pass'] for r in sub):.4f}; "
                 f"max false-FAIL {max(r['false_fail'] for r in sub):.4f}")
    # power: P(PASS) when true bias 0 or 0.5% and noise small
    for beta in (0.0, 0.005):
        for sigma in (0.001, 0.003):
            pw = [r["inconclusive"] for r in sub if r["beta"] == beta and r["sigma"] == sigma]
            lines.append(f"    bias {100*beta:.1f}% sigma {100*sigma:.1f}%: P(INCONCLUSIVE) worst over noise types {max(pw):.3f}")
print("\n".join(lines), flush=True)

# ---------------- Part B: alpha CI calibration
OM = np.array([100, 300, 1000, 3000, 10000.0])
RB, B = 300, 400
lines.append(f"Part B (alpha CI; true alpha = 1; Omega = {OM.astype(int).tolist()}; {RB} reps x {B} bootstrap):")
best = None
for widen in (1.0, 1.25, 1.5, 2.0):
    worst_cov_b, rows = 1.0, []
    for M in (6, 12):
        for kind in ("iso", "rank3", "twoP+", "t3"):
            for ratio in (0.1, 0.3, 1.0):          # run-level noise / bias at Omega = 100
                beta100 = 0.05
                cov = n_def = 0
                for _ in range(RB):
                    runs = [phi + beta100 * (100 / o) * u + noise(kind, M, ratio * beta100 * np.sqrt(100 / o)) for o in OM]
                    a, (lo, hi), fv = alpha_ci({"ref": phi, "runs": runs}, OM, B=B, widen=widen, rng=rng)
                    if np.isfinite(lo):
                        n_def += 1; cov += lo <= 1.0 <= hi
                c = cov / max(n_def, 1)
                rows.append(dict(widen=widen, M=M, kind=kind, ratio=ratio, coverage=c, defined=n_def / RB))
                worst_cov_b = min(worst_cov_b, c)
    out["B"].extend(rows)
    mc = 2 * np.sqrt(0.95 * 0.05 / RB)
    lines.append(f"  widen={widen}: worst coverage {worst_cov_b:.3f} (target >= 0.95 - 2 MC SE = {0.95 - mc:.3f})")
    print(lines[-1], flush=True)
    if worst_cov_b >= 0.95 - mc:
        best = widen
        break
lines.append(f"  calibrated widening factor: {best}")
for r in [r for r in out["B"] if r["widen"] == best] if best else []:
    lines.append(f"    M={r['M']} {r['kind']:6s} noise/bias@100={r['ratio']}: coverage {r['coverage']:.3f}, CI defined {r['defined']:.2f}")
out["widen"] = best
open("results/MG003b/calibration.txt", "w").write("\n".join(lines) + "\n")
json.dump(out, open("results/MG003b/calibration.json", "w"), indent=1)
print("\n".join(lines[-30:]))

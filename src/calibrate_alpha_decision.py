"""MG003b Part B2: unconditional calibration of the MG2-C2 DECISION rule (protocol §5.3) at widening
w = 2.0, M = 12. True alpha in {1.0 (correct O(1/Omega) bias), 0.5 (noise-like scaling)}; six run
families; run-noise/bias ratio at Omega=100 in {0.1, 0.3, 1}. Reports P(PASS), P(FAIL), P(INCONCLUSIVE).
Acceptance: under alpha=1, P(FAIL) <= 0.05 in every scenario; under alpha=0.5, P(PASS) <= 0.05.
Writes results/MG003b/calibration_alpha_decision.txt."""
import sys
import numpy as np
sys.path.insert(0, "src")
from stats_mg2 import alpha_ci

rng = np.random.default_rng(20261029)
d = 100
phi = rng.standard_normal(d); phi /= np.linalg.norm(phi)
u = rng.standard_normal(d); u -= (u @ phi) * phi; u /= np.linalg.norm(u)
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

def decide_c2(a, lo, hi):
    if not np.isfinite(lo): return "INCONCLUSIVE"
    if 0.7 <= a <= 1.3 and not (lo <= 0.5 <= hi): return "PASS"
    if hi < 0.7 or lo > 1.3: return "FAIL"
    return "INCONCLUSIVE"

OM = np.array([100, 300, 1000, 3000, 10000.0]); W = 2.0; RB, B, M = 300, 400, 12
lines = [f"Part B2: C2 decision calibration, w={W}, M={M}, {RB} reps x {B} bootstrap"]
worst_ff, worst_fp = 0.0, 0.0
for alpha in (1.0, 0.5):
    for kind in ("iso", "rank1b", "rank3", "twoP+", "twoP-", "t3"):
        for ratio in (0.1, 0.3, 1.0):
            beta100 = 0.05; cnt = {"PASS": 0, "FAIL": 0, "INCONCLUSIVE": 0}
            for _ in range(RB):
                runs = [phi + beta100 * (100 / o) ** alpha * u + noise(kind, M, ratio * beta100 * np.sqrt(100 / o)) for o in OM]
                a, (lo, hi), fv, keep = alpha_ci({"ref": phi, "runs": runs}, OM, B=B, widen=W, rng=rng)
                cnt[decide_c2(a, lo, hi)] += 1
            p = {k: v / RB for k, v in cnt.items()}
            if alpha == 1.0: worst_ff = max(worst_ff, p["FAIL"])
            else: worst_fp = max(worst_fp, p["PASS"])
            lines.append(f"  true alpha {alpha} {kind:6s} noise/bias@100={ratio}: P(PASS) {p['PASS']:.3f} P(FAIL) {p['FAIL']:.3f} P(INC) {p['INCONCLUSIVE']:.3f}")
            print(lines[-1], flush=True)
lines.append(f"worst false-FAIL (true alpha=1): {worst_ff:.3f}; worst false-PASS (true alpha=0.5): {worst_fp:.3f}; "
             f"acceptance (<= 0.05 each): {'MET' if worst_ff <= 0.05 and worst_fp <= 0.05 else 'NOT MET'}")
print(lines[-1])
open("results/MG003b/calibration_alpha_decision.txt", "w").write("\n".join(lines) + "\n")

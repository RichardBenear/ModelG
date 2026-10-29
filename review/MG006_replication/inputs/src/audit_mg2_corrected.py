"""MG005: corrected gate-level A1 audit applied to the MG2 data (descriptive; NOT a retroactive registered
replacement). Same projections as the registered pooled audit in evaluate_mg2.py (run means projected on the
bias direction and the first two within-run principal components, per field), but each run is studentised
against the other runs of its own cell, Bonferroni over all run x field x direction tests at 1%.
Usage: python src/audit_mg2_corrected.py > results/MG005/audit_corrected.txt"""
import os, sys, pickle
import numpy as np
sys.path.insert(0, os.path.dirname(__file__))
from stats_mg2 import audit_within_cell
from runner import reference

CELLS = "results/MG2/cells"
projs, labels = [], []
for fn in sorted(os.listdir(CELLS)):
    c = pickle.load(open(f"{CELLS}/{fn}", "rb"))
    if "n_eligible" not in c or c["n_eligible"] < 12:
        continue
    ref = reference(c["b"], c["h"])
    win = np.abs(ref["x"]) < 10.0
    el = [r for r in c["runs"] if r.get("eligible")]
    Bk = np.array([r["blocks"] for r in el])
    for s, nm in enumerate("GXY"):
        blocks = Bk[:, :, s]; M, K, d = blocks.shape
        runs = blocks.mean(1); mbar = runs.mean(0)
        W = (blocks - runs[:, None, :]).reshape(M * K, d)
        _, _, Vt = np.linalg.svd(W - W.mean(0), full_matrices=False)
        bd = mbar - ref["phi"][s][win]
        bd = bd / np.linalg.norm(bd) if np.linalg.norm(bd) > 0 else Vt[2]
        for j, v in enumerate((bd, Vt[0], Vt[1])):
            projs.append((runs - mbar) @ v)
            labels.append((c["name"], nm, j, [r["replicate"] for r in el]))
res = audit_within_cell(projs)
cells = sorted({l[0] for l in labels})
print(f"Corrected within-cell audit (descriptive), {len(cells)} decision-eligible cells: {', '.join(cells)}")
print(f"tests {res['n_tests']}; max studentised |t| {res['max_stat']:.4f}; Bonferroni critical value {res['crit']:.4f} (1%)")
if res["flagged"]:
    for i, r, t in res["flagged"]:
        nm, f, j, reps = labels[i]
        print(f"FLAG {nm} field {f} dir {j} replicate {reps[r]}: |t| = {t:.3f}")
else:
    print("No run flagged.")

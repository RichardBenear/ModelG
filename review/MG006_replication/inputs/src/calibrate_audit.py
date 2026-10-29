"""MG005: calibration of the gate-level audits (null false-alarm rate and power against one outlying run).
Design matches MG2: 13 cells x 12 runs x 9 (field, direction) projections. Gaussian null; alternative: one run in
one cell shifted by delta (in units of the run-level SD). Compares the registered pooled audit (cell-standardised
before pooling) with the corrected within-cell audit (stats_mg2.audit_within_cell).
Usage: python src/calibrate_audit.py > results/MG005/audit_calibration.txt"""
import os, sys
import numpy as np
from scipy import stats
sys.path.insert(0, os.path.dirname(__file__))
from stats_mg2 import _studentized_loo, audit_within_cell

C, M, D, NSIM = 13, 12, 9, 2000
rng = np.random.default_rng(20261028)


def registered(projs):
    flag = False
    for d in range(D):
        z = np.concatenate([(p - p.mean()) / p.std(ddof=1) for p in projs[d::D]])
        tcrit = stats.t.ppf(1 - 0.01 / (2 * z.size * 9), z.size - 2)
        flag |= np.max(np.abs(_studentized_loo(z))) > tcrit
    return flag


def sim(delta):
    fr = fc = 0
    for _ in range(NSIM):
        projs = [rng.standard_normal(M) for _ in range(C * D)]
        if delta:
            projs[0][0] += delta
        fr += registered(projs)
        fc += bool(audit_within_cell(projs)["flagged"])
    return fr / NSIM, fc / NSIM


print(f"Audit calibration: {C} cells x {M} runs x {D} projections, {NSIM} simulations per row (seed 20261028).")
print("delta (run-level SDs) | registered pooled audit flag rate | corrected within-cell audit flag rate")
for delta in (0, 4, 6, 8, 10, 15, 20):
    a, b = sim(delta)
    print(f"{delta:>5} | {a:.4f} | {b:.4f}", flush=True)

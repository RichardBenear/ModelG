"""MG9c Draft 1 checks: existing MG9 tracks and static decision logic only.

Usage: python design_checks.py CHECKOUT
No Model G integration or production-module imports.
"""
import json
import sys
from pathlib import Path
import numpy as np
from scipy.optimize import curve_fit

root = Path(sys.argv[1]) / 'results/MG9/runs'
lam = 100 / 34
out = {'reviewed_revision': '4466be3471056b6e2a9507fe6d98498b17de61ba', 'fits': {}}
for b in [29.25, 29.30, 29.35]:
    r = json.loads((root / f'front_grid_{b:.4f}.json').read_text())
    tr = r['track']
    t = np.array([v['t'] for v in tr])
    mask = t >= 600
    sides = {}
    for side, sign in [('xL', -1), ('xR', 1)]:
        y = sign * np.array([v[side] for v in tr])
        def model(t, limit, amplitude, tau):
            return limit - amplitude * np.exp(-(t - 600) / tau)
        p, _ = curve_fit(model, t[mask], y[mask], p0=[y[-1] + .22, 1.2, 720])
        exp_rms = np.sqrt(np.mean((y[mask] - model(t[mask], *p))**2))
        straight = np.polyval(np.polyfit(t[mask], y[mask], 1), t[mask])
        line_rms = np.sqrt(np.mean((y[mask] - straight)**2))
        sides[side] = dict(tau=float(p[2]), remaining_at_1800=float(p[0]-model(1800,*p)),
                          exponential_rms=float(exp_rms), linear_rms=float(line_rms),
                          predicted_D3=float(model(12000,*p)-model(6000,*p)))
    changes = [dict(t=v['t'], n=v['n']) for i,v in enumerate(tr)
               if i == 0 or v['n'] != tr[i-1]['n']]
    out['fits'][str(b)] = dict(sides=sides, cell_counts=sorted(set(v['n'] for v in tr)),
                             cell_count_changes=changes)

# A midpoint sequence allowed by the proposed three-step bisection.
lo, hi = 29.35, 29.40
steps = []
for _ in range(3):
    mid = (lo + hi) / 2
    steps.append(mid)
    lo = mid  # SETTLED at each midpoint, not a prediction.
br = [29.395313, 29.396094]
out['q2_counterexample'] = dict(midpoints_settled=steps, bracket=[lo,hi],
    draft_label='AT b_r' if lo >= br[0]-.01 else 'other',
    hypothetical_edge_above_br=29.398,
    edge_inside_bracket=lo < 29.398 < hi,
    edge_above_entire_br_benchmark=29.398 > br[1])
assert lo < 29.398 < hi and 29.398 > br[1]
out['thresholds'] = dict(wavelength=lam, settled_displacement=lam/40,
    settled_speed=lam/40000, initial_window_displacement_speed=lam/40/6000,
    early_stop_displacement=lam/400)
# Both large, opposite outward displacements miss all three clauses of rule 3.
dL, dR = .3*lam, -.3*lam
out['opposite_motion'] = dict(DL=dL, DR=dR,
    rule3_advances=dL>=lam/4 and dR>=lam/4,
    rule3_recedes=dL<=-lam/4 and dR<=-lam/4,
    rule3_asymmetric=(abs(dL)>=lam/4) != (abs(dR)>=lam/4))
dest = Path(__file__).with_name('design_checks.json')
dest.write_text(json.dumps(out, indent=2)+'\n')
print(json.dumps(out, indent=2))

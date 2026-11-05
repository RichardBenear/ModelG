"""MG9b Draft 3 budget-stop rule examples. No simulations or production imports."""
import json
from pathlib import Path

cases=[
    dict(name='budget_before_29_20_starts',fitted_b=[29.00,29.35],short_cell_in_fit=False,
         budget_truncated_rate_design=True),
    dict(name='budget_after_two_completed_cells_before_third',fitted_b=[29.00,29.20,29.35],
         short_cell_in_fit=False,budget_truncated_rate_design=True),
    dict(name='budget_inside_production_cell',fitted_b=[29.00,29.20,29.35],
         short_cell_in_fit=True,budget_truncated_rate_design=True),
    dict(name='ordinary_completed_design',fitted_b=[29.00,29.20,29.35],
         short_cell_in_fit=False,budget_truncated_rate_design=False)]
for c in cases:
    # Assume a directional interval and a calibrated nominal available-b branch.
    c['draft3_budget_fallback']=c['short_cell_in_fit']
    c['recommended_budget_fallback']=c['budget_truncated_rate_design']
assert [c['name'] for c in cases if c['draft3_budget_fallback']!=c['recommended_budget_fallback']]==[
    'budget_before_29_20_starts','budget_after_two_completed_cells_before_third']
out=dict(revision='111a07a43ef9d5e27d0543a55132f94adcaedfd0',cases=cases,
         scope='Logical examples only; assumes a directional interval and a calibrated nominal b branch.')
Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))

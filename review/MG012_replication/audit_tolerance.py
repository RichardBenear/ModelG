"""Measure the registered-vs-implemented T operator tolerance; no dynamics."""
from pathlib import Path
import sys,json
import numpy as np
ROOT=Path(__file__).resolve().parent
IN=ROOT/'inputs'
sys.path.insert(0,str(IN/'src'))
import graphs,run_mg4c
dec=json.loads((IN/'results/MG4/part0/decisions.json').read_text())
g,h,_=run_mg4c.graph_for(run_mg4c.job(29,'T',op='cotan'),dec)
a=float(np.linalg.norm(g['dx'][0]))
reference=graphs.laplacian(g)/(1.5*a*a)
weighted=graphs.cotan(g)['L']
error=float(abs(weighted-reference).max())
absolute_limit=1e-12
scaled_limit=1e-12*float(abs(reference).max())
y=np.load(IN/'results/MG4c/fields/w_b29.000000_T_p0.npz')['final']
old=np.load(IN/'results/MG4/part0/steps/Tstate_122x142.npy').reshape(y.shape)
out=dict(reviewed_commit=json.loads((ROOT/'download_manifest.json').read_text())['commit'],
    max_entry_difference=error,registered_entrywise_absolute_limit=absolute_limit,
    implemented_scaled_limit=scaled_limit,passes_literal_absolute_check=error<absolute_limit,
    passes_implemented_check=error<scaled_limit,
    max_entry_difference_using_MG4_calibrated_h=float(abs(weighted-graphs.laplacian(g)/h**2).max()),
    T_final_field_max_difference=float(abs(y-old).max()),T_field_limit=1e-6,
    recommendation='Record and explicitly adjudicate roundoff-level tolerance deviation; no new trajectories indicated')
assert not out['passes_literal_absolute_check'] and out['passes_implemented_check']
(ROOT/'tolerance_audit.json').write_text(json.dumps(out,indent=2))
print(json.dumps(out,indent=2))

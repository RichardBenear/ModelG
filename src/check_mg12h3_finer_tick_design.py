"""MG145 independent saved-layout and precision arithmetic checks."""
from pathlib import Path
import numpy as np,json,hashlib
from run_mg12h3_sparse_combined import load,ROOT
BASE=ROOT/'results/MG12H3/finer_tick_design'
def run():
 q=json.loads((BASE/'planning.json').read_text())
 for rel,sha in q['source_hashes'].items():assert hashlib.sha256((ROOT/rel).read_bytes()).hexdigest()==sha
 for row in q['precision']:
  assert row['nominal_halfwidth']>row['transport_baseline_standard_error']>=0 and 0<=row['approximate_zero_difference_equivalence_power']<=1
 directory=ROOT/'results/MG12H3/multitick_pilot/production/o0_a0_Omega10_d1';initial,model=load(ROOT/'results/MG12H3/ring_feasibility/o0_a0_Omega10.npz');pop=initial.copy()
 with np.load(BASE/'saved_pilot_compact_layout.npz') as z:hist=z['channel_histograms'];final=z['final'];available=z['available']
 for step,path in enumerate(sorted(directory.glob('tick_*.npz'))):
  with np.load(path) as z:unique,n=np.unique(z['channel_ids'],return_counts=True);expected=np.zeros(9000,np.uint32);expected[unique]=n;np.testing.assert_array_equal(expected,hist[step]);np.testing.assert_array_equal(pop,z['initial']);np.testing.assert_array_equal(final[step],z['final']);np.testing.assert_array_equal(available[step],z['available'])
  reserve=np.zeros(pop.size,np.int64);change=reserve.copy()
  for c,count in zip(model.channels,hist[step]):
   if count:
    for i,v in c.reserve:reserve[i]+=int(count)*v
    for i,v in c.change:change[i]+=int(count)*v
  np.testing.assert_array_equal(pop.ravel()-reserve,available[step].ravel());np.testing.assert_array_equal(pop.ravel()+change,final[step].ravel());assert available[step].min()>=0;pop=final[step]
 summary=[]
 for r in (32,48,64):
  rows=[x for x in q['precision'] if x['replicates']==r];assert len(rows)==24;summary.append(dict(replicates=r,transport_baseline_worst_contrast_halfwidth=max(x['nominal_halfwidth'] for x in rows if x['observable']=='contrast'),transport_baseline_worst_power_halfwidth=max(x['nominal_halfwidth'] for x in rows if x['observable']=='high_power'),approximate_zero_effect_family_union_bound=max(0,1-sum(1-x['approximate_zero_difference_equivalence_power'] for x in rows))))
 (BASE/'validation.json').write_text(json.dumps(dict(passed=True,saved_ticks_reconstructed=len(hist),all_old_files_unchanged=True,precision_summary=summary,limitations='Saved-only layout validation, not production writer qualification or exact coupled-model power.'),indent=2)+'\n');print(json.dumps(summary))
if __name__=='__main__':run()

"""MG157 separate saved-data analysis recomputation; no new random draws."""
import os
for k in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):os.environ[k]='1'
import numpy as np,time,json,math
from mg12h3_coupling_long import BASE,read,hashes
from mg12h3_finer_precision import setup
from mg10e2c1_pool import atomic_json
start=time.perf_counter();plan=read(BASE/'plan.json');source=read(BASE/'analysis_results.json');initial,contrast=setup();rows=[]
with np.load(BASE/'bootstrap_indices.npz') as z:indices=z['indices'];assert int(z['seed'])==plan['seeds']['bootstrap']
for c in range(2):
 groups={};costs={}
 for kind in ('paired','independent'):
  delta=[];cost=[]
  for r in range(8):
   path=BASE/'production'/f'c{c}_{kind}_r{r:02d}'
   with np.load(path/'observations.npz') as z:
    left=contrast.weights@z['counts'][-1,0].ravel();right=contrast.weights@z['counts'][-1,1].ravel();delta.append(left-right)
   cost.append(read(path.with_suffix('.json'))['elapsed_seconds'])
  groups[kind]=np.asarray(delta);costs[kind]=np.asarray(cost)
 def variance(a,axis=0):return np.sum((a-np.mean(a,axis=axis,keepdims=True))**2,axis=axis)/7
 raw=variance(groups['paired'])/variance(groups['independent']);cost=float(sum(costs['paired'])/sum(costs['independent']));ip,ii=indices[c];vp=variance(groups['paired'][ip,2],axis=1);vi=variance(groups['independent'][ii,2],axis=1);boot=vp/vi*costs['paired'][ip].mean(axis=1)/costs['independent'][ii].mean(axis=1);finite=bool(np.isfinite(boot).all());upper=float(np.percentile(boot,97.5)) if finite else None;expected=source['comparisons'][c];np.testing.assert_allclose(raw,expected['variance_ratio'][-1],rtol=1e-10,atol=1e-12);assert (upper is None)==(expected['primary_Y_bootstrap_upper_97_5'] is None)
 if upper is not None:np.testing.assert_allclose(upper,expected['primary_Y_bootstrap_upper_97_5'],rtol=1e-10,atol=1e-12)
 supported=finite and upper<1 and float(raw[2])<1;assert supported==expected['pilot_supported'];rows.append(dict(comparison=c,raw_ratio=raw.tolist(),cost_ratio=cost,Y_cost_adjusted_ratio=float(raw[2])*cost,Y_upper=upper,pilot_supported=supported))
assert all(r['pilot_supported'] for r in rows)==source['pilot_benefit_supported'];atomic_json(BASE/'separate_analysis_recomputation.json',dict(passed=True,comparisons=rows,elapsed_seconds=time.perf_counter()-start,source_hashes=hashes(),method='Recomputed projections from raw endpoint counts, centered sums of squares and stored bootstrap indices; same author/shared weights and seeds, no external replication.'));print(json.dumps(rows))

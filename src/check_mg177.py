import os
for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS'):os.environ[k]='1'
import json,gzip,math,time
import numpy as np
from scipy.special import logsumexp
from run_mg176 import ROOT,RAW,OUT as PRIOR,models,read,atomic_json
from audit_mg177 import OUT
from mg174_modes import ModeTracker

def check():
 start=time.perf_counter();assert not(OUT/'check_attempt.json').exists();atomic_json(OUT/'check_attempt.json',dict(status='RUNNING'))
 try:
  a=read(OUT/'audit.json');model,L=next((m,L) for ri,m,L in models() if ri==2);base=ModeTracker(model,L,0,100);real=np.fft.irfft(base.ft,n=500,axis=2);templates=np.array([np.roll(t,k,axis=1) for t in real for k in range(500)]).reshape(6500,1000);distance=np.minimum(np.arange(500),500-np.arange(500))*.2;eligible=np.broadcast_to(distance[None,:]>6.4,(13,500));N=eligible.sum();frames=0;maxdiff=0.;filters=0
  for path in a['conditional_paths']:
   tag=f"r2_rho0_scale{path['scale']:g}_{path['case']}";d=json.loads(gzip.decompress((PRIOR/(tag+'.json.gz')).read_bytes()));rep=path['replicate'];rec=d['paths'][rep]
   with np.load(RAW/(tag+'.npz')) as z:fields=z['fields'][rep]
   P=base.P/path['estimated_variance'];filtered_templates=np.fft.irfft(np.einsum('fab,kbf->kaf',P,base.ft),n=500,axis=2);norms=np.sum(real*filtered_templates,axis=(1,2));rs=rm=ps=.5;pg=np.zeros((13,500));pg[eligible]=.5/N;previous={k:.2 for k in rec['records'][0]['class_weights']}
   for field,truth,row,original in zip(fields,rec['truth'],path['rows'],rec['records']):
    k,ss=truth;yf=np.fft.rfft(field-field.mean(axis=1)[:,None],axis=1);phase=np.exp(-2j*np.pi*np.arange(251)*k/500);e=yf-base.ft[4]*phase;filtered=np.roll(np.fft.irfft(np.einsum('fab,bf->af',P,e),n=500,axis=1),-k,axis=1);grid=(templates@filtered.ravel()).reshape(13,500)-.5*norms[:,None];grid=np.where(eligible,grid,-np.inf);bf=float(logsumexp(grid)-math.log(N));np.testing.assert_allclose(bf,row['oracle_primary_uniform_logBF'],rtol=1e-7,atol=4e-6);np.testing.assert_allclose(grid[4,250],row['oracle_pair_logBF'],rtol=1e-7,atol=4e-6);maxdiff=max(maxdiff,abs(bf-row['oracle_primary_uniform_logBF']));rS=.98*rs+.005*rm;rM=(.98*rm+.005*rs)*math.exp(bf);rs,rm=rS/(rS+rM),rM/(rS+rM)
    pS=.98*ps+.005*pg.sum();pG=.98*pg+.005*ps/N;mx=max(0.,float(grid[eligible].max()));pS*=math.exp(-mx);pG=np.where(eligible,pG*np.exp(grid-mx),0);z=pS+pG.sum();ps=pS/z;pg=pG/z;np.testing.assert_allclose([rm,pg.sum()],[row['oracle_primary_redrawn_multiple_weight'],row['oracle_primary_persistent_multiple_weight']],rtol=1e-5,atol=2e-7);filters+=2
    prior=math.log((.98*previous['multiple']+.005*(1-previous['multiple']))/(.98*previous['single']+.005*(1-previous['single'])));post=math.log(original['class_weights']['multiple']/original['class_weights']['single']);np.testing.assert_allclose([prior,post-prior],[row['prior_log_multiple_single'],row['integrated_log_evidence_ratio']],rtol=1e-10,atol=1e-10);previous=original['class_weights'];frames+=1
  assert frames==288;atomic_json(OUT/'separate_check.json',dict(passed=True,full_secondary_grids_checked=frames,eligible_grid_points_per_frame=int(N),temporal_filter_steps_checked=filters,prior_evidence_decompositions=frames,max_logBF_difference=maxdiff,elapsed_seconds=time.perf_counter()-start,limits='Real-space dense template products and probability-domain filters; same templates/precision/author, not external replication.'));atomic_json(OUT/'check_attempt.json',dict(status='COMPLETE',elapsed_seconds=time.perf_counter()-start))
 except BaseException:
  import traceback
  atomic_json(OUT/'check_attempt.json',dict(status='FAILED',error=traceback.format_exc(),elapsed_seconds=time.perf_counter()-start));raise
if __name__=='__main__':check()

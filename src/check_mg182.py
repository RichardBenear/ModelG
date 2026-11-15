import os
for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS'):os.environ[k]='1'
import time,json,hashlib
import numpy as np
from scipy.integrate import solve_ivp
from scipy.sparse import coo_matrix
from mg182_gscale import setup,ROOT,PATH
from prepare_mg182 import OUT,RAW,read,atomic_json
from mg11n1 import continuum_rhs

def check():
 start=time.perf_counter();assert not(OUT/'check_attempt_v2.json').exists();atomic_json(OUT/'check_attempt_v2.json',dict(status='RUNNING'));traces=events=drifts=0
 try:
  results=read(OUT/'scaling_checks.json');benchmark=read(OUT/'benchmark.json')['rows']
  for row in results['rows']:
   c=setup(row['omega']);z=c['pop']/c['omega'];linksA=[(i,j,w) for kind,i,j,w in c['links'] if kind=='ALink'];linksB=[(i,j,w) for kind,i,j,w in c['links'] if kind=='BLink'];assert linksA==linksB;g=dict(n=500,src=np.array([x[0] for x in linksA]),dst=np.array([x[1] for x in linksA]),weight=np.array([x[2] for x in linksA]));case=dict(omega=c['omega'],a0=c['a']*c['omega'],b0=c['b']*c['omega'],diff=c['diff']);expected=continuum_rhs(0,z.ravel(),case,g,np.array(c['coupling'])).reshape(7,500);expected[1]-=z[1]*z[2]/c['omega'];expected[2]+=z[1]*z[2]/c['omega'];rows=[];cols=[];values=[]
   for i,ch in enumerate(c['model'].channels):
    for pos,delta in ch.change:rows.append(pos);cols.append(i);values.append(delta)
    Gchanges=[(pos,delta) for pos,delta in ch.change if pos<500]
    if Gchanges:assert ch.rule in ['R1','R2','R2r','H_G']
   S=coo_matrix((values,(rows,cols)),shape=(3500,len(c['model'].channels))).tocsr();rates=np.array([ch.rate(c['pop'].ravel()) for ch in c['model'].channels]);actual=(S@rates/c['omega']).reshape(7,500);np.testing.assert_allclose(actual,expected,rtol=1e-8,atol=4e-8);assert np.all(np.asarray(S.sum(axis=0))==0);drifts+=1
   for b in benchmark:
    if b['omega']!=row['omega'] or b['paired']:continue
    with np.load(RAW/f'benchmark_{b["index"]}.npz') as z:
     av=z['initial'].ravel().copy();hist=np.zeros(len(c['model'].channels),np.int64);last=0.
     for index,mark in zip(z['ids'],z['marks']):
      ch=c['model'].channels[index];assert last<mark<b['duration'];last=mark;assert ch.rate(av)>0
      for pos,n in ch.reserve:assert av[pos]>=n;av[pos]-=n
      hist[index]+=1
     np.testing.assert_array_equal(hist,z['histogram']);np.testing.assert_array_equal(av,z['available'].ravel());np.testing.assert_array_equal(z['initial'].ravel()+S@hist,z['final'].ravel());traces+=1;events+=len(z['ids'])
   if row['omega']==160:
    with np.load(RAW/'full_tick_probe.npz') as z:
     for side in (0,1):np.testing.assert_array_equal(z['initial'].ravel()+S@(z['histogram'][0]+z['histogram'][side+1]),z['final'][side].ravel())
  c=setup(10);times=np.array([0,.005,.01,.1]);g=dict(n=500,src=np.array([i for kind,i,j,w in c['links'] if kind=='ALink']),dst=np.array([j for kind,i,j,w in c['links'] if kind=='ALink']),weight=np.array([w for kind,i,j,w in c['links'] if kind=='ALink']));case=dict(omega=10,a0=c['a']*10,b0=c['b']*10,diff=c['diff']);solution=solve_ivp(lambda t,y:continuum_rhs(t,y,case,g,np.array(c['coupling'])),(0,.1),c['reference'].ravel(),method='RK45',t_eval=times,max_step=.0005,rtol=1e-10,atol=1e-12);assert solution.success
  with np.load(RAW/'continuum_reference.npz') as z:err=float(np.max(abs(solution.y.T.reshape(-1,7,500)-z['fields'])));gerr=float(np.max(abs(solution.y.T.reshape(-1,7,500)[:,0]-z['fields'][:,0])));np.testing.assert_allclose(solution.y.T.reshape(-1,7,500),z['fields'],rtol=2e-8,atol=2e-8);assert gerr<1e-8
  import importlib
  from mg12h3_contrast_estimator import Contrast
  from mg12h3_pattern_persistence import homogenize
  for omega in(10,40,160):
   c=setup(omega);u=c['reference'][0]-c['reference'][0].mean();u-=u.mean();u/=np.linalg.norm(u);u-=u.mean();weights=np.zeros((1,3500));weights[0,:500]=u/omega;contrast=Contrast(c['model'],weights)
   for pop in(c['pop'],homogenize(c['pop'])):np.testing.assert_allclose(contrast.direct(pop)[0],((pop[0]/omega)-(pop[0]/omega).mean())@u,rtol=1e-10,atol=1e-12)
  atomic_json(OUT/'separate_check.json',dict(passed=True,independent_mean_field_drift_cases=drifts,chronological_single_traces=traces,chronological_single_events=events,large_pair_endpoint_reconstructions=2,independent_RK45_vs_BDF_max_error=err,G_max_error=gerr,proposed_projection_identity_checks=6,production_started=False,elapsed_seconds=time.perf_counter()-start,limits='Independent older continuum RHS/sparse stoichiometry/RK45 and trace arithmetic; same author and underlying Model G equations, not external replication.'));atomic_json(OUT/'check_attempt_v2.json',dict(status='COMPLETE',elapsed_seconds=time.perf_counter()-start));atomic_json(OUT/'local_raw_manifest.json',dict(directory=str(RAW),files={p.name:dict(bytes=p.stat().st_size,sha256=hashlib.sha256(p.read_bytes()).hexdigest()) for p in RAW.iterdir()},full_raw_published=False))
 except BaseException:
  import traceback
  atomic_json(OUT/'check_attempt_v2.json',dict(status='FAILED',elapsed_seconds=time.perf_counter()-start,error=traceback.format_exc()));raise
if __name__=='__main__':check()

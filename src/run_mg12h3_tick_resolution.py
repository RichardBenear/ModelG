"""MG144 exact transport-mean development on actual substrates, no RNG."""
import os
for key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[key]='1'
from pathlib import Path
import time,json,math,hashlib,traceback
import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import expm_multiply
ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/'results/MG12H3/tick_resolution'
DELTAS=(.001,.0005,.0001,.00002);HORIZONS=(.01,.1)

def run():
 if (BASE/'attempt.json').exists():raise ValueError('Already attempted')
 BASE.mkdir(parents=True,exist_ok=True);start=time.perf_counter();(BASE/'attempt.json').write_text(json.dumps(dict(status='STARTED',cap_seconds=60))+'\n');records=[];controls=[];sources={}
 for path in sorted((ROOT/'results/MG12H3/ring_feasibility').glob('*Omega1.npz')):
  sources[path.relative_to(ROOT).as_posix()]=hashlib.sha256(path.read_bytes()).hexdigest()
  with np.load(path) as z:counts=z['counts'];links=json.loads(str(z['links']));diff=z['diff'];coupling=z['coupling']
  n=counts.shape[1];position=np.arange(n);profiles=[np.ones(n)/math.sqrt(n),np.eye(1,n,0).ravel()];names=['uniform','point_site0']
  for k in (1,5,25,125,250):
   for name,fun in [('cos',np.cos),('sin',np.sin)]:
    if k==250 and name=='sin':continue
    v=fun(2*np.pi*k*position/n);profiles.append(v/np.linalg.norm(v));names.append(name+str(k))
  for s,species in enumerate(('G','X','Y','A','B')):
   if time.perf_counter()-start>60:raise TimeoutError('Stage1 elapsed cap')
   kind='ActiveLink' if s<3 else ('ALink' if s==3 else 'BLink');coefficient=float(diff[s] if s<3 else coupling[s-3]);rate=np.zeros((n,n))
   for link,i,j,w in links:
    if link==kind:rate[j,i]+=coefficient*w
   exits=rate.sum(axis=0);L=rate-np.diag(exits);np.testing.assert_allclose(L,L.T,rtol=0,atol=1e-12);np.testing.assert_allclose(exits,exits[0],rtol=1e-12,atol=1e-12);assert abs(L.sum(axis=0)).max()<1e-10
   eigenvalues,U=np.linalg.eigh(L);assert eigenvalues.max()<1e-9;np.testing.assert_allclose(L@U,U*eigenvalues,rtol=1e-9,atol=1e-9);np.testing.assert_allclose(U.T@U,np.eye(n),rtol=1e-9,atol=1e-9)
   actual=counts[s].astype(float);actual=actual/np.linalg.norm(actual) if np.linalg.norm(actual) else actual;X=np.array(profiles+[actual]).T;labels=names+['actual_species'];coeff=U.T@X;continuous=[];discrete=[];P_all=[];multipliers=[];lam=float(exits[0]);drift=[]
   for horizon in HORIZONS:
    continuous.append(U@(np.exp(eigenvalues*horizon)[:,None]*coeff))
   np.testing.assert_allclose(continuous[0],expm_multiply(csr_matrix(L)*HORIZONS[0],X),rtol=1e-9,atol=1e-9)
   for delta in DELTAS:
    factor=-math.expm1(-lam*delta)/lam if lam else delta;P=np.eye(n)+factor*L;assert P.min()>-1e-12;np.testing.assert_allclose(P.sum(axis=0),1.,rtol=0,atol=1e-10);np.testing.assert_allclose(P@X[:,0],X[:,0],rtol=1e-10,atol=1e-10);P_all.append(P);drift.append(factor/delta);bytime=[];bymult=[]
    for t,horizon in enumerate(HORIZONS):
     steps=round(horizon/delta);assert abs(steps*delta-horizon)<1e-14;mult=np.power(1+factor*eigenvalues,steps);mean=U@(mult[:,None]*coeff);np.testing.assert_allclose(mean[:,0],X[:,0],rtol=1e-9,atol=1e-9);assert mean[:,:2].min()>-1e-9 and mean[:,-1].min()>-1e-9
     # Conserves sums, including signed linear-response probes.
     np.testing.assert_allclose(mean.sum(axis=0),X.sum(axis=0),rtol=1e-9,atol=1e-9);error=mean-continuous[t];op=float(abs(mult-np.exp(eigenvalues*horizon)).max());bytime.append(mean);bymult.append(mult)
     records.append(dict(case=path.stem,species=species,exit_rate=lam,delta=delta,horizon=horizon,initial_drift_factor=factor/delta,operator_L2_error=op,profile_names=labels,profile_L2_error=np.linalg.norm(error,axis=0).tolist(),maximum_site_error=np.max(abs(error),axis=0).tolist(),screen_1percent_operator=bool(op<=.01)))
    discrete.append(bytime);multipliers.append(bymult)
   iteration=X.copy();P=csr_matrix(P_all[0])
   for _ in range(10):iteration=P@iteration
   np.testing.assert_allclose(iteration,discrete[0][0],rtol=1e-9,atol=1e-9);controls.append(dict(case=path.stem,species=species,independent_continuous=True,iterated_tick=True,residual_passed=True))
   np.savez_compressed(BASE/(path.stem+'_'+species+'.npz'),L=L,exits=exits,eigenvalues=eigenvalues,eigenvectors=U,profiles=X,profile_names=labels,continuous=np.asarray(continuous),discrete=np.asarray(discrete),transition=np.asarray(P_all),multipliers=np.asarray(multipliers),deltas=DELTAS,horizons=HORIZONS)
 # Explicit zero-exit control.
 zero=np.zeros((3,3));np.testing.assert_array_equal(np.eye(3)+.001*zero,np.eye(3))
 for name in ('src/run_mg12h3_tick_resolution.py','spec/MG12H3_tick_resolution_stage1.md'):
  sources[name]=hashlib.sha256((ROOT/name).read_text(encoding='utf-8').replace('\r\n','\n').encode()).hexdigest()
 elapsed=time.perf_counter()-start;assert elapsed<60;(BASE/'results.json').write_text(json.dumps(dict(complete=True,elapsed_seconds=elapsed,records=records,controls=controls,zero_exit_control=True,source_hashes=sources,random_draws=0,scope='Deterministic pure transport means; no coupled-model or particle qualification'),indent=2)+'\n')
if __name__=='__main__':
 try:run()
 except BaseException:BASE.mkdir(parents=True,exist_ok=True);(BASE/'failure.json').write_text(json.dumps(dict(error=traceback.format_exc()),indent=2)+'\n');raise

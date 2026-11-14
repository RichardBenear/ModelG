import os
for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS'):os.environ[k]='1'
import json,gzip,time,traceback,math
import numpy as np
from scipy.optimize import minimize_scalar
from scipy.stats import chi2,binomtest
from scipy.linalg import solve_triangular
from run_mg179 import OUT,PRIOR,INPUT,SEED,read,atomic_json,models
from mg174_modes import ModeTracker

def independent_fit(stats):
 qs,A,B,C,first=stats;D=998;M=15968;n=63872
 def score(r,joint):
  qt=first+(A-2*r*B+r*r*C)/(1-r*r);df=n+M if joint else n;scale=(qt+(qs if joint else 0))/df;v=-.5*(df*(math.log(scale)+1)+62874*math.log(1-r*r))
  if not joint:v+=-.5*M*(math.log(qs/M)+1)
  return v
 maxima=[]
 for joint in (False,True):
  grid=np.linspace(-.99,.99,401);values=np.array([score(r,joint) for r in grid]);i=int(values.argmax());lo=grid[max(0,i-1)];hi=grid[min(400,i+1)];opt=minimize_scalar(lambda r:-score(r,joint),bounds=(lo,hi),method='bounded',options={'xatol':1e-12});maxima.append(score(opt.x,joint))
 return max(0.,2*(maxima[0]-maxima[1]))
def check():
 start=time.perf_counter();assert not(OUT/'check_attempt.json').exists();atomic_json(OUT/'check_attempt.json',dict(status='RUNNING'))
 try:
  a=read(OUT/'audit.json');stats=json.loads(gzip.decompress((OUT/'calibration_statistics.json.gz').read_bytes()));profiles=0;maximum=0.;draws=0
  selected=[r for i,r in enumerate(stats) if i%17==0 or r['phase']==1]
  for r in a['reference_profiles']+selected:
   stored=r.get('profile',r.get('fit'));lr=independent_fit(r['stats']);np.testing.assert_allclose(lr,stored['LR'],rtol=1e-5,atol=3e-6);maximum=max(maximum,abs(lr-stored['LR']));profiles+=1
  for r in stats:
   fit=r['fit'];np.testing.assert_allclose(chi2.sf(fit['LR'],1),fit['p'],rtol=1e-10,atol=1e-12)
  for r in stats:
   if r['replicate']>=3:continue
   rng=np.random.default_rng(np.random.SeedSequence([SEED,r['phase'],r['reg'],r['replicate']]));qs=float(rng.chisquare(15968));y=rng.normal(size=(64,998));amp=1.25 if r['phase'] else 1.;rho=r['rho'];y[0]*=amp
   for t in range(1,64):y[t]=rho*y[t-1]+np.sqrt(1-rho*rho)*amp*y[t]
   np.testing.assert_allclose(qs,r['stats'][0]);rr=r['fit']['shared']['rho'];K=rr**abs(np.arange(64)[:,None]-np.arange(64)[None,:]);chol=np.linalg.cholesky(K);z=solve_triangular(chol,y,lower=True);qt=float(np.sum(z*z));_,A,B,C,first=r['stats'];np.testing.assert_allclose(qt,first+(A-2*rr*B+rr*rr*C)/(1-rr*rr),rtol=1e-10,atol=1e-6);draws+=1
  for row in a['null_tails']:
   sample=[r for r in stats if r['phase']==0 and r['rho']==row['rho']];count=sum(r['fit']['p']<row['tail'] for r in sample);pv=float(binomtest(count,len(sample),row['tail']).pvalue);assert count==row['count'] and (pv>=.01/9)==row['screen_passed'];np.testing.assert_allclose(pv,row['binomial_p'])
  null=[r for r in stats if r['phase']==0];defects=[r for r in stats if r['phase']==1];assert sum(r['fit']['p']<.001/12 for r in null)==a['rare_tail_count'];assert sum(not r['fit']['valid'] or r['fit']['p']<.001/12 for r in defects)==a['defects_detected'];cache={ri:(m,L) for ri,m,L in models()};grids=steps=labels=0
  for path in a['component_paths']:
   ri=0 if path['regime']==.1 else 2;model,L=cache[ri];base=ModeTracker(model,L,0,100);tag=f"r{ri}_rho0_scale{path['noise']:g}_{path['case']}";d=json.loads(gzip.decompress((PRIOR/(tag+'.json.gz')).read_bytes()));old=d['paths'][path['replicate']]
   with np.load(INPUT/(tag+'.npz')) as z:seq=z['fields'][path['replicate']]
   P=base.P/old['calibration']['scale'];rs=rm=ps=.5;pg=np.full(48,.5/48);shifts=list(range(-24,0))+list(range(1,25))
   for field,truth,row,r in zip(seq,old['truth'],path['rows'],old['records']):
    lab=r['winner'] if r['class_weights'][r['winner']]>=.95 and r['innovation_q']<=r['qcut'] else 'uncertain'
    if lab=='single' and (r['single_center_spread']>4 or r['single_width_spread']>.25):lab='uncertain'
    assert lab==r['label']==row['old_label'];labels+=1;k,width=truth;template=model.signal(width,k);yf=np.fft.rfft(field-field.mean(axis=1)[:,None],axis=1);mu=np.fft.rfft(template,axis=1);e0=yf-mu
    def q(e):return float(np.sum(np.fft.irfft(e,n=500,axis=1)*np.fft.irfft(np.einsum('fab,bf->af',P,e),n=500,axis=1)))
    q0=q(e0);values=[]
    for shift in shifts:
     shifted=template.copy();shifted[1]=np.roll(shifted[1],shift);e=yf-np.fft.rfft(shifted,axis=1);values.append(-.5*(q(e)-q0))
    np.testing.assert_allclose(values,row['shift_logBF'],rtol=1e-7,atol=5e-6);grids+=1;v=np.array(values);bf=float(np.log(np.exp(v-v.max()).sum())+v.max()-math.log(48));rS=.98*rs+.005*rm;rM=(.98*rm+.005*rs)*math.exp(bf);rs,rm=rS/(rS+rM),rM/(rS+rM);pS=.98*ps+.005*pg.sum();pG=(.98*pg+.005*ps/48)*np.exp(v);z=pS+pG.sum();ps,pg=pS/z,pG/z;np.testing.assert_allclose([rm,pg.sum()],[row['redrawn_displaced_weight'],row['persistent_displaced_weight']],rtol=1e-5,atol=3e-7);steps+=2
   assert time.perf_counter()-start<120
  assert labels==576;atomic_json(OUT/'separate_check.json',dict(passed=True,independent_profile_fits=profiles,max_LR_difference=maximum,fit_pvalues=len(stats),selected_generator_covariance_checks=draws,full_shift_grids=grids,temporal_filter_steps=steps,unchanged_labels=labels,elapsed_seconds=time.perf_counter()-start,limits='Independent optimizer bracketing, dense AR covariance Cholesky, real-space shifted templates; shared input/template/noise model and same author, not external replication.'));atomic_json(OUT/'check_attempt.json',dict(status='COMPLETE',elapsed_seconds=time.perf_counter()-start))
 except BaseException:atomic_json(OUT/'check_attempt.json',dict(status='FAILED',error=traceback.format_exc(),elapsed_seconds=time.perf_counter()-start));raise
if __name__=='__main__':check()

import os
for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS'):os.environ[k]='1'
import time,json,gzip,hashlib,traceback,math
import numpy as np
from scipy.stats import binomtest,chi2
from scipy.special import logsumexp
from run_mg178 import ROOT,OUT as PRIOR,RAW as INPUT,read,atomic_json,models
from mg174_modes import ModeTracker
from check_mg174 import inner,mu
from mg179_profile import fit_reference,D,N,M
OUT=ROOT/'results/MG12H3/profile_component_audit';SEED=17900
FILES=['src/mg179_profile.py','src/run_mg179.py','spec/MG179_profile_component_audit.md','src/run_mg178.py','src/mg178_background.py','src/mg178_persistent.py','src/run_mg176.py','src/analyze_mg170.py','src/mg167_xy.py','src/mg161_local_shape.py','src/check_mg174.py']
def register():
 OUT.mkdir(exist_ok=True);assert not(OUT/'plan.json').exists();paths=sorted(PRIOR.glob('r*.json.gz'))+sorted(INPUT.glob('r*.npz'))+[PRIOR/'plan.json',PRIOR/'analysis.json'];atomic_json(OUT/'plan.json',dict(files={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths+[ROOT/x for x in FILES]},seed=SEED,cap_seconds=600,workers=1,null_replicates=3000,rhos=[0.,.8,.95],defect_replicates=256,no_ModelG_draws=True))
def realinner(a,b,P):return np.sum(np.fft.irfft(a,n=500,axis=-1)*np.fft.irfft(np.einsum('fab,...bf->...af',P,b),n=500,axis=-1),axis=(-1,-2))
def generated(phase,reg,rep,rho,amplitude):
 rng=np.random.default_rng(np.random.SeedSequence([SEED,phase,reg,rep]));qspace=float(rng.chisquare(M));y=rng.normal(size=(N,D));y[0]*=amplitude
 for t in range(1,N):y[t]=rho*y[t-1]+np.sqrt(1-rho*rho)*amplitude*y[t]
 A=float(np.sum(y[1:]**2));B=float(np.sum(y[1:]*y[:-1]));C=float(np.sum(y[:-1]**2));qfirst=float(np.sum(y[0]**2));return [qspace,A,B,C,qfirst]
def run():
 start=time.perf_counter();assert not(OUT/'attempt.json').exists();atomic_json(OUT/'attempt.json',dict(status='RUNNING'));references=[];paths=[];grids=0;labels=0;maximum=0.
 try:
  plan=read(OUT/'plan.json');assert all(hashlib.sha256((ROOT/k).read_bytes()).hexdigest()==v for k,v in plan['files'].items());cache={ri:(m,L) for ri,m,L in models()}
  for file in sorted(PRIOR.glob('r*.json.gz')):
   d= json.loads(gzip.decompress(file.read_bytes()));ri=0 if d['regime']==.1 else 2;model,L=cache[ri];base=ModeTracker(model,L,0,100);P=base.P;ft=base.ft;w=base.freqweights
   with np.load(INPUT/(file.name[:-8]+'.npz')) as z:spatial=z['spatial'];temporal=z['temporal'];fields=z['fields']
   for old,sp,tm,seq in zip(d['paths'],spatial,temporal,fields):
    ys=np.fft.rfft(sp-sp.mean(axis=2)[:,:,None],axis=2);yt=np.fft.rfft(tm-tm.mean(axis=2)[:,:,None],axis=2);stats=[float(realinner(ys,ys,P).sum()),float(realinner(yt[1:],yt[1:],P).sum()),float(realinner(yt[1:],yt[:-1],P).sum()),float(realinner(yt[:-1],yt[:-1],P).sum()),float(realinner(yt[0],yt[0],P))];cal=old['calibration'];np.testing.assert_allclose(stats[0],sum(cal['q']),rtol=1e-8,atol=3e-6);fit=fit_reference(*stats);references.append(dict(regime=d['regime'],rho=d['rho'],noise=d['noise_amplitude'],case=d['case'],primary=d['primary'],replicate=old['replicate'],stats=stats,old_adequate=cal['adequate'],old_temporal_p=cal['checks']['temporal_variance_p'],profile=fit,profile_temporal_tail=not fit['valid'] or fit['p']<.001/12))
    if d['rho']!=0 or d['case'] not in ('displacedY','moving'):continue
    Pe=P/cal['scale'];shifts=np.r_[np.arange(-24,0),np.arange(1,25)];ls=lm=-math.log(2);ps=-math.log(2);pg=np.full(48,-math.log(2)-math.log(48));rows=[]
    for t,(field,truth,r) in enumerate(zip(seq,old['truth'],old['records'])):
     lab=r['winner'] if max(r['class_weights'].values())>=.95 and r['innovation_q']<=r['qcut'] else 'uncertain'
     if lab=='single' and(r['single_center_spread']>4 or r['single_width_spread']>.25):lab='uncertain'
     assert lab==r['label'] and r['emitted']==(lab if cal['adequate'] else 'noise_unqualified');labels+=1;k,s=truth;yf=np.fft.rfft(field-field.mean(axis=1)[:,None],axis=1);shared=mu('single',k,round((s-.5)/.125),None,ft);err=yf-shared;q0=inner(err,err,Pe,w);g=[]
     for shift in shifts:
      shifted=shared.copy();shifted[1]*=np.exp(-2j*np.pi*np.arange(251)*shift/500);e=yf-shifted;value=-.5*(inner(e,e,Pe,w)-q0);direct=-.5*(float(realinner(e,e,Pe))-float(realinner(err,err,Pe)));np.testing.assert_allclose(value,direct,rtol=1e-8,atol=4e-6);maximum=max(maximum,abs(value-direct));grids+=1;g.append(value)
     g=np.array(g);bf=float(logsumexp(g)-math.log(48));known=float(g[np.where(shifts==16)[0][0]]);difference=shared.copy();difference[0]=0;difference[1]*=np.exp(-2j*np.pi*np.arange(251)*16/500)-1;separation=float(inner(difference,difference,Pe,w));rS=float(np.logaddexp(ls+math.log(.98),lm+math.log(.005)));rM=float(np.logaddexp(lm+math.log(.98),ls+math.log(.005)))+bf;z=float(np.logaddexp(rS,rM));ls,lm=rS-z,rM-z;pS=float(np.logaddexp(ps+math.log(.98),logsumexp(pg)+math.log(.005)));pG=np.logaddexp(pg+math.log(.98),ps+math.log(.005)-math.log(48))+g;z=float(np.logaddexp(pS,logsumexp(pG)));ps,pg=pS-z,pG-z;rows.append(dict(frame=t,old_label=r['label'],old_emitted=r['emitted'],class_weights=r['class_weights'],selected_center=r['center'],selected_width=r['width'],truth_center=k,rho_estimate=cal['rho'],innovation_q=r['innovation_q'],qcut=r['qcut'],center_spread=r['single_center_spread'],width_spread=r['single_width_spread'],oracle_shift16_logBF=known,uniform_shift_logBF=bf,oracle_squared_separation=separation,redrawn_displaced_weight=float(np.exp(lm)),persistent_displaced_weight=float(np.exp(logsumexp(pg))),shift_logBF=g.tolist()))
    paths.append(dict(regime=d['regime'],noise=d['noise_amplitude'],case=d['case'],replicate=old['replicate'],rows=rows))
   if time.perf_counter()-start>600:raise TimeoutError('MG179cap')
  results=[];tails=[]
  for reg,rho in enumerate(plan['rhos']):
   values=[]
   for rep in range(plan['null_replicates']):
    stats=generated(0,reg,rep,rho,1.);fit=fit_reference(*stats);values.append(fit['p']);results.append(dict(phase=0,reg=reg,rho=rho,replicate=rep,stats=stats,fit=fit))
    if time.perf_counter()-start>600:raise TimeoutError('MG179cap')
   for tail in(.05,.01,.001):
    count=sum(p<tail for p in values);pv=float(binomtest(count,len(values),tail).pvalue);tails.append(dict(rho=rho,tail=tail,count=count,n=len(values),binomial_p=pv,screen_passed=pv>=.01/9))
  defects=[]
  for rep in range(plan['defect_replicates']):
   stats=generated(1,2,rep,.95,1.25);fit=fit_reference(*stats);defects.append(dict(phase=1,reg=2,rho=.95,replicate=rep,stats=stats,fit=fit))
  nullvalid=all(r['fit']['valid'] for r in results);detected=sum(not r['fit']['valid'] or r['fit']['p']<.001/12 for r in defects);atomic_json(OUT/'audit.json',dict(reference_profiles=references,component_paths=paths,unchanged_labels=labels,realspace_shift_checks=grids,max_shift_difference=maximum,null_tails=tails,null_screen_passed=nullvalid and all(r['screen_passed'] for r in tails),rare_tail_count=sum(r['fit']['p']<.001/12 for r in results),null_n=len(results),defects_detected=detected,defect_n=len(defects),defect_screen_passed=detected>=.9*len(defects),elapsed_seconds=time.perf_counter()-start,no_new_tracker_verdict=True));(OUT/'calibration_statistics.json.gz').write_bytes(gzip.compress(json.dumps(results+defects,separators=(',',':')).encode(),mtime=0));atomic_json(OUT/'attempt.json',dict(status='COMPLETE',elapsed_seconds=time.perf_counter()-start))
 except BaseException:atomic_json(OUT/'attempt.json',dict(status='FAILED',error=traceback.format_exc(),elapsed_seconds=time.perf_counter()-start));raise
if __name__=='__main__':register() if __import__('sys').argv[1]=='register' else run()

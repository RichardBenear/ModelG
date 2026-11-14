import os
for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS'):os.environ[k]='1'
import time,hashlib,gzip,json,traceback
import numpy as np
from run_mg174 import ROOT,OLD,BENCH,read,atomic_json,geometry,POSITIVES,NEGATIVES,FILES as PRIORFILES,verify
from mg167_xy import XYCore
from mg176_acquisition import estimate,tracker,outcome
OUT=ROOT/'results/MG12H3/acquisition_noise';RAW=ROOT/'local_raw/MG176_acquisition'
FILES=PRIORFILES+['src/mg176_acquisition.py','src/run_mg176.py','spec/MG176_acquisition_noise.md']
def binding():return {p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in FILES}
def register():
 OUT.mkdir(exist_ok=True);assert not(OUT/'plan.json').exists();inputs=[OLD/'prepared_inputs.npz',BENCH]+[OLD/'production'/f'homogeneous_r{r}'/'observations.npz' for r in (0,1)];atomic_json(OUT/'plan.json',dict(source=binding(),inputs={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs},seed=17600,cap_seconds=1200,workers=1,threads=1,replicates=2,frames=24,calibration_fields=16,primary_cases=POSITIVES+NEGATIVES,sensitivities=['weak_static','weak_moving','wrong_rho','shape_X','contaminated_reference'],new_ModelG_draws=False))
def models():
 with np.load(OLD/'prepared_inputs.npz') as z:reference=z['reference'];initial=z['initial'];omega=float(z['omega'])
 inputs=[]
 for r in (0,1):
  with np.load(OLD/'production'/f'homogeneous_r{r}'/'observations.npz') as z:inputs.append(z['counts'][:,0,:3]/omega)
 cov=read(BENCH)['regime_inputs']
 for ri in (0,2):
  when=[.1,.5,1.][ri];model=XYCore(reference,initial,omega,np.array(cov[ri]['covariance'])[1:,1:]);samples=np.array([x[round(when*1000)] for x in inputs]);samples-=samples.mean(axis=2)[:,:,None];ft=np.fft.rfft(samples,axis=2);S=np.einsum('rsf,rtf->fst',ft,np.conj(ft))/(2*500);S=np.stack([S[max(0,i-2):min(251,i+3)].mean(axis=0) for i in range(251)]);diag=np.zeros_like(S)
  for sp in range(3):diag[:,sp,sp]=S[:,sp,sp].real
  S=.9*S+.1*diag;ev,V=np.linalg.eigh(S);ev=np.maximum(ev,1e-8*ev.max(axis=1)[:,None]);L=V*np.sqrt(ev)[:,None,:];L[0]=0;yield ri,model,L

def run():
 p=read(OUT/'plan.json');assert p['source']==binding() and not(OUT/'attempt.json').exists();assert all(hashlib.sha256((ROOT/k).read_bytes()).hexdigest()==h for k,h in p['inputs'].items());RAW.mkdir(exist_ok=False);start=time.perf_counter();atomic_json(OUT/'attempt.json',dict(status='RUNNING'));rng=np.random.default_rng(p['seed']);rows=[];checks=0
 try:
  for ri,model,L in models():
   np.savez_compressed(RAW/f'model_r{ri}.npz',L=L,FT=np.array([np.fft.rfft(model.signal(s,0),axis=1) for s in np.arange(.5,2.0001,.125)]))
   def noise():return np.fft.irfft(np.einsum('fij,jf->if',L,np.fft.rfft(rng.normal(size=(3,500)),axis=1))*model.x.smoother,n=500,axis=1)[1:]
   settings=[(rho,scale,case,True) for rho in (0.,.95) for scale in (.8,1.,1.25) for case in p['primary_cases']]
   if ri==2:settings +=[(.95,1.,case,False) for case in p['sensitivities']]
   for rho,scale,case,primary in settings:
    paths=[];allfields=[];allcal=[]
    for rep in range(2):
     seed=(100+71*rep)%500;calfields=np.array([noise()*scale for _ in range(16)])
     if case=='shape_X':calfields[:,0]*=1.5
     if case=='contaminated_reference':calfields+=model.signal(1,(seed+250)%500)
     cal=estimate(model,L,calfields,seed);obj=tracker(model,L,rho,seed,cal);records=[];truths=[];fields=[];n=noise();actualrho=.8 if case=='wrong_rho' else rho
     for t in range(24):
      if t:n=actualrho*n+np.sqrt(1-actualrho**2)*noise()
      mapped=case if primary or case.startswith('weak_') else 'static';k,s=geometry(mapped,t,seed);signal=model.signal(s,k)*(.8 if case.startswith('weak_') else 1.)
      if case=='absent' or(case=='vanish' and t>=10):signal[:]=0
      if case=='missingX' and t>=10:signal[0]=0
      if case=='missingY' and t>=10:signal[1]=0
      if case=='displacedY' and t>=10:signal[1]=np.roll(signal[1],16)
      if case=='split' and t>=10:signal+=model.signal(s,(k+250)%500)
      nn=n*scale
      if case=='shape_X':nn=nn.copy();nn[0]*=1.5
      field=signal+nn;r=obj.step(field);r['emitted']=r['label'] if cal['adequate'] else 'noise_unqualified';records.append(r);truths.append([k,s]);fields.append(field)
      if rep==0 and t==10:checks+=verify(obj,False)
      if time.perf_counter()-start>1200:raise TimeoutError('MG176 cap')
     positive=case in POSITIVES or not primary;result=outcome(records,truths,positive,case=='absent');paths.append(dict(replicate=rep,seed=seed,calibration=cal,truth=truths,records=records,outcome=result));allfields.append(fields);allcal.append(calfields)
    tag=f'r{ri}_rho{rho:g}_scale{scale:g}_{case}';np.savez_compressed(RAW/(tag+'.npz'),fields=allfields,calibration=allcal);data=dict(regime=[.1,.5,1.][ri],rho=rho,noise_amplitude=scale,case=case,primary=primary,paths=paths);payload=json.dumps(data,separators=(',',':')).encode();(OUT/(tag+'.json.gz')).write_bytes(gzip.compress(payload,mtime=0));rows.append(dict(regime=data['regime'],rho=rho,noise_amplitude=scale,case=case,primary=primary,success=sum(x['outcome']['success'] for x in paths),reference_adequate=sum(x['calibration']['adequate'] for x in paths),acquired=sum(x['outcome']['acquisition'] is not None for x in paths)))
    if sum(x.stat().st_size for x in RAW.iterdir())>250000000:raise RuntimeError('raw250MBcap')
    atomic_json(OUT/'progress.json',dict(completed_cases=len(rows),elapsed_seconds=time.perf_counter()-start))
  atomic_json(OUT/'analysis.json',dict(rows=rows,primary_screen_passed=all(r['success']==2 for r in rows if r['primary']),endpoint_checks=checks,ready_for_real_persistence=False,G_status='UNRESOLVED',elapsed_seconds=time.perf_counter()-start));atomic_json(OUT/'attempt.json',dict(status='COMPLETE',elapsed_seconds=time.perf_counter()-start))
 except BaseException:atomic_json(OUT/'attempt.json',dict(status='FAILED',error=traceback.format_exc(),elapsed_seconds=time.perf_counter()-start));raise
if __name__=='__main__':register() if __import__('sys').argv[1]=='register' else run()

import os
for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS'):os.environ[k]='1'
import time,hashlib,gzip,json,traceback,math
import numpy as np
from scipy.special import logsumexp
from run_mg176 import ROOT,OLD,BENCH,models,read,atomic_json
from run_mg174 import geometry,POSITIVES,NEGATIVES,FILES as OLD_FILES
from mg178_background import calibrate,make_tracker
from mg176_acquisition import outcome
OUT=ROOT/'results/MG12H3/persistent_secondary';RAW=ROOT/'local_raw/MG178_persistent'
EXTRA=['secondary_static','secondary_counter','secondary_initial','secondary_vanish'];SENS=['weak_static','weak_moving','actual_rho','shape_X','contaminated_reference','temporal_shift']
FILES=OLD_FILES+['src/run_mg176.py','src/mg176_acquisition.py','src/mg178_persistent.py','src/mg178_background.py','src/run_mg178.py','spec/MG178_persistent_secondary.md']
def binding():return {p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in FILES}
def register():
 OUT.mkdir(exist_ok=True);assert not(OUT/'plan.json').exists();inputs=[OLD/'prepared_inputs.npz',BENCH]+[OLD/'production'/f'homogeneous_r{r}'/'observations.npz' for r in(0,1)];atomic_json(OUT/'plan.json',dict(source=binding(),inputs={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs},seed=17800,cap_seconds=2000,workers=1,threads=1,primary=POSITIVES+NEGATIVES,extra=EXTRA,sensitivities=SENS,new_ModelG_draws=False))
def result(records,truth,case,primary):
 if case not in EXTRA:return outcome(records,truth,case in POSITIVES or not primary,case=='absent')
 target='multiple' if case=='secondary_initial' else 'single';times=[t for t in range(1,6) if records[t-1]['emitted']==records[t]['emitted']==target];acq=min(times) if times else None;checks={12:'multiple',20:'single'} if case=='secondary_vanish' else {15:'multiple'};okay=acq is not None and all(records[t]['emitted']==v for t,v in checks.items());return dict(acquisition=acq,acquisition_target=target,deadline_labels={str(t):records[t]['emitted'] for t in checks},success=bool(okay))
def directchecks(obj):
 s=obj.audit_snapshot;checks=0
 for motion in ('S','M'):
  ends=s['ends'][('multiple',motion)];(k,si),grid=next(iter(ends.items()));flat=int(np.argmax(grid));ss,kk=np.unravel_index(flat,grid.shape);vals=[];mu=obj.mean('multiple',k,si,(kk,ss));rho=0 if obj.frame_index==0 else obj.rho;distance=np.minimum((np.arange(500)-k)%500,(k-np.arange(500))%500)*.2;eligible=np.broadcast_to(distance[None,:]>6.4,(13,500))
  for p in s['old']:
   allowed=k==p['center'] and si==p['si'] if motion=='S' else min((k-p['center'])%500,(p['center']-k)%500)<=2 and abs(si-p['si'])<=1
   if not allowed:continue
   lp=-math.log(5655)
   if p['kind']=='multiple' and p['second'] is not None:
    oldk,oldsi=p['second'];options=[(sj,kj%500) for sj in range(max(0,oldsi-1),min(13,oldsi+2)) for kj in range(oldk-2,oldk+3) if eligible[sj,kj%500]]
    if options:
     z={q:.5/len(options) for q in options}
     if eligible[oldsi,oldk]:z[(oldsi,oldk)]=z.get((oldsi,oldk),0)+.5
     if (ss,kk) not in z:continue
     lp=math.log(z[(ss,kk)]/sum(z.values()))
   e=s['innovation']-(mu-rho*p['mean']);q=float(obj.inner(e,e)/s['variance']);nloc=1 if motion=='S' else 5*(min(13,p['si']+2)-max(0,p['si']-1));vals.append(p['log']+math.log(.98 if p['kind']=='multiple' else .005)+math.log(.98 if p['motion']==motion else .02)-math.log(nloc)+lp-.5*q)
  np.testing.assert_allclose(logsumexp(vals),grid[ss,kk],rtol=1e-8,atol=4e-6);checks+=1
 return checks

def run():
 plan=read(OUT/'plan.json');assert plan['source']==binding() and not(OUT/'attempt.json').exists();assert all(hashlib.sha256((ROOT/k).read_bytes()).hexdigest()==v for k,v in plan['inputs'].items());RAW.mkdir(exist_ok=False);start=time.perf_counter();atomic_json(OUT/'attempt.json',dict(status='RUNNING'));rng=np.random.default_rng(plan['seed']);rows=[];checks=0
 try:
  for ri,model,L in models():
   def noise():return np.fft.irfft(np.einsum('fij,jf->if',L,np.fft.rfft(rng.normal(size=(3,500)),axis=1))*model.x.smoother,n=500,axis=1)[1:]
   settings=[(rho,scale,case,True) for rho in(0.,.95) for scale in(.8,1.,1.25) for case in plan['primary']]
   if ri==2:settings += [(rho,1.25,case,True) for rho in(0.,.95) for case in EXTRA]+[(.95,1.,case,False) for case in SENS]
   for rho,scale,case,primary in settings:
    paths=[];allfields=[];allspatial=[];alltemporal=[]
    for rep in range(2):
     seed=(100+71*rep)%500;spatial=np.array([noise()*scale for _ in range(16)]);ref_rho=.8 if case=='actual_rho' else rho;n=noise();temporal=[]
     for t in range(64):
      if t:n=ref_rho*n+np.sqrt(1-ref_rho**2)*noise()
      temporal.append(n*scale)
     temporal=np.array(temporal)
     if case=='shape_X':spatial[:,0]*=1.5;temporal[:,0]*=1.5
     if case=='contaminated_reference':spatial+=model.signal(1,(seed+250)%500);temporal+=model.signal(1,(seed+250)%500)
     cal=calibrate(model,L,spatial,temporal,seed);obj=make_tracker(model,L,seed,cal);records=[];fields=[];truths=[];n=noise();actualrho=.8 if case in ('actual_rho','temporal_shift') else rho
     for t in range(24):
      if t:n=actualrho*n+np.sqrt(1-actualrho**2)*noise()
      mapped=case if case in plan['primary'] or case.startswith('weak_') else 'moving' if case in EXTRA else 'static';k,s=geometry(mapped,t,seed);signal=model.signal(s,k)*(.8 if case.startswith('weak_') else 1.)
      if case=='absent' or(case=='vanish' and t>=10):signal[:]=0
      if case=='missingX' and t>=10:signal[0]=0
      if case=='missingY' and t>=10:signal[1]=0
      if case=='displacedY' and t>=10:signal[1]=np.roll(signal[1],16)
      present=(case in ['split','secondary_static','secondary_counter'] and t>=10) or case=='secondary_initial' or(case=='secondary_vanish' and 6<=t<14)
      if present:
       second=(seed+250)%500 if case=='secondary_static' else (seed+250-t)%500 if case=='secondary_counter' else (k+250)%500;signal+=model.signal(s,second)
      nn=n*scale
      if case=='shape_X':nn=nn.copy();nn[0]*=1.5
      field=signal+nn;obj.frame_index=t;r=obj.step(field);r['emitted']=r['label'] if cal['adequate'] else 'noise_unqualified';records.append(r);fields.append(field);truths.append([k,s])
      if rep==0 and t in(0,10):checks+=directchecks(obj)
      if time.perf_counter()-start>2000:raise TimeoutError('2000s cap')
     paths.append(dict(replicate=rep,seed=seed,calibration=cal,truth=truths,records=records,outcome=result(records,truths,case,primary)));allfields.append(fields);allspatial.append(spatial);alltemporal.append(temporal)
    tag=f'r{ri}_rho{rho:g}_scale{scale:g}_{case}';np.savez_compressed(RAW/(tag+'.npz'),fields=allfields,spatial=allspatial,temporal=alltemporal);d=dict(regime=[.1,.5,1.][ri],rho=rho,noise_amplitude=scale,case=case,primary=primary,paths=paths);(OUT/(tag+'.json.gz')).write_bytes(gzip.compress(json.dumps(d,separators=(',',':')).encode(),mtime=0));rows.append(dict(regime=d['regime'],rho=rho,noise_amplitude=scale,case=case,primary=primary,success=sum(p['outcome']['success'] for p in paths),reference_adequate=sum(p['calibration']['adequate'] for p in paths),acquired=sum(p['outcome']['acquisition'] is not None for p in paths)))
    if sum(p.stat().st_size for p in RAW.iterdir())>400000000:raise RuntimeError('400MBcap')
    atomic_json(OUT/'progress.json',dict(completed_cases=len(rows),elapsed_seconds=time.perf_counter()-start))
  atomic_json(OUT/'analysis.json',dict(rows=rows,primary_screen_passed=all(r['success']==2 for r in rows if r['primary']),endpoint_checks=checks,ready_for_real_persistence=False,G_status='UNRESOLVED',elapsed_seconds=time.perf_counter()-start));atomic_json(OUT/'attempt.json',dict(status='COMPLETE',elapsed_seconds=time.perf_counter()-start))
 except BaseException:atomic_json(OUT/'attempt.json',dict(status='FAILED',error=traceback.format_exc(),elapsed_seconds=time.perf_counter()-start));raise
if __name__=='__main__':register() if __import__('sys').argv[1]=='register' else run()

import os
for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS'):os.environ[k]='1'
import time,hashlib,gzip,json,traceback
import numpy as np
from run_mg180 import ROOT,OLD,BENCH,models,read,atomic_json,directchecks,FILES as OLD_FILES
from run_mg174 import geometry
from mg180_background import calibrate,make_tracker
from mg181_acquisition import metrics,episodes
OUT=ROOT/'results/MG12H3/acquisition_feasibility';RAW=ROOT/'local_raw/MG181_acquisition';PRIOR=ROOT/'results/MG12H3/component_alignment'
FILES=OLD_FILES+['src/mg181_acquisition.py','src/run_mg181.py','spec/MG181_acquisition_feasibility.md']
def binding():return {p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in FILES}
def register():
 OUT.mkdir(exist_ok=True);assert not(OUT/'plan.json').exists();inputs=[OLD/'prepared_inputs.npz',BENCH]+[OLD/'production'/f'homogeneous_r{r}'/'observations.npz' for r in(0,1)]+list(PRIOR.glob('r*.json.gz'));atomic_json(OUT/'plan.json',dict(source=binding(),inputs={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs},seed=18100,cap_seconds=2000,workers=1,threads=1,new_ModelG_draws=False,frames=40,paths=200))
def settings():
 out=[]
 for rho in(0.,.8,.95):
  for scale in(.8,1.25):
   for case in ['static','moving','stop_start','birth_moving','weak_static','absent','Xonly','Yonly','displaced_initial','split_initial']:out.append((rho,scale,case,0,4 if case in ['static','moving','stop_start','birth_moving'] else 2))
 for rho in(0.,.95):
  for onset in(3,8):
   for case in ['damage_X','damage_Y','damage_displaced','damage_split']:out.append((rho,1.25,case,onset,2))
 return out

def run():
 plan=read(OUT/'plan.json');assert plan['source']==binding() and not(OUT/'attempt.json').exists();assert all(hashlib.sha256((ROOT/k).read_bytes()).hexdigest()==v for k,v in plan['inputs'].items());RAW.mkdir(exist_ok=False);start=time.perf_counter();atomic_json(OUT/'attempt.json',dict(status='RUNNING'));rng=np.random.default_rng(plan['seed']);saved=[];rows=[];checks=0
 try:
  for file in sorted(PRIOR.glob('r*.json.gz')):
   d=json.loads(gzip.decompress(file.read_bytes()))
   if not d['primary']:continue
   for p in d['paths']:
    hits=episodes(p['records'],2);saved.append(dict(regime=d['regime'],rho=d['rho'],noise=d['noise_amplitude'],case=d['case'],replicate=p['replicate'],old_success=p['outcome']['success'],old_acquisition=p['outcome']['acquisition'],first_anytime_hit=hits[0] if hits else None,episode_times=hits,not_a_rescore=True))
  model,L=next((m,L) for ri,m,L in models() if ri==2)
  def noise():return np.fft.irfft(np.einsum('fij,jf->if',L,np.fft.rfft(rng.normal(size=(3,500)),axis=1))*model.x.smoother,n=500,axis=1)[1:]
  for rho,scale,case,onset,reps in settings():
   paths=[];allfields=[];allspatial=[];alltemporal=[]
   for rep in range(reps):
    seed=(100+71*rep)%500;spatial=np.array([noise()*scale for _ in range(16)]);n=noise();temporal=[]
    for t in range(64):
     if t:n=rho*n+np.sqrt(1-rho*rho)*noise()
     temporal.append(n*scale)
    temporal=np.array(temporal);cal=calibrate(model,L,spatial,temporal,seed);obj=make_tracker(model,L,seed,cal);records=[];truths=[];fields=[];n=noise()
    for t in range(40):
     if t:n=rho*n+np.sqrt(1-rho*rho)*noise()
     if case=='birth_moving':k=(seed+max(0,t-8))%500;s=1.
     else:k,s=geometry('static' if case=='weak_static' else case if case in ['static','moving','stop_start'] else 'moving',t,seed)
     signal=model.signal(s,k)*(.8 if case=='weak_static' else 1.)
     if case=='absent' or(case=='birth_moving' and t<8):signal[:]=0
     if case=='Xonly' or(case=='damage_Y' and t>=onset):signal[1]=0
     if case=='Yonly' or(case=='damage_X' and t>=onset):signal[0]=0
     if case=='displaced_initial' or(case=='damage_displaced' and t>=onset):signal[1]=np.roll(signal[1],16)
     if case=='split_initial' or(case=='damage_split' and t>=onset):signal+=model.signal(s,(k+250)%500)
     field=signal+n*scale;obj.frame_index=t;r=obj.step(field);r['emitted']=r['label'] if cal['adequate'] else 'noise_unqualified';records.append(r);truths.append([k,s]);fields.append(field)
     if rep==0 and t in(0,10):checks+=directchecks(obj)
     if time.perf_counter()-start>2000:raise TimeoutError('MG181cap')
    score={str(h):metrics(records,truths,case,onset,h) for h in(2,3)};paths.append(dict(replicate=rep,seed=seed,calibration=cal,truth=truths,records=records,metrics=score));allfields.append(fields);allspatial.append(spatial);alltemporal.append(temporal)
   tag=f'rho{rho:g}_scale{scale:g}_{case}_onset{onset}';np.savez_compressed(RAW/(tag+'.npz'),fields=allfields,spatial=allspatial,temporal=alltemporal);d=dict(rho=rho,noise=scale,case=case,onset=onset,paths=paths);(RAW/(tag+'_records.json.gz')).write_bytes(gzip.compress(json.dumps(d,separators=(',',':')).encode(),mtime=0));public=json.loads(json.dumps(d));[r.pop('beam') for p in public['paths'] for r in p['records']];(OUT/(tag+'.json.gz')).write_bytes(gzip.compress(json.dumps(public,separators=(',',':')).encode(),mtime=0));rows.append(dict(rho=rho,noise=scale,case=case,onset=onset,n=reps,reference_adequate=sum(p['calibration']['adequate'] for p in paths),paths=[dict(replicate=p['replicate'],metrics=p['metrics']) for p in paths]));atomic_json(OUT/'progress.json',dict(completed_cases=len(rows),elapsed_seconds=time.perf_counter()-start))
   if sum(p.stat().st_size for p in RAW.iterdir())>400000000:raise RuntimeError('400MBrawcap')
  atomic_json(OUT/'analysis.json',dict(saved_development=saved,rows=rows,endpoint_checks=checks,no_qualification_verdict=True,old_MG180_unchanged=True,ready_for_real_persistence=False,elapsed_seconds=time.perf_counter()-start));atomic_json(OUT/'attempt.json',dict(status='COMPLETE',elapsed_seconds=time.perf_counter()-start))
 except BaseException:atomic_json(OUT/'attempt.json',dict(status='FAILED',error=traceback.format_exc(),elapsed_seconds=time.perf_counter()-start));raise
if __name__=='__main__':register() if __import__('sys').argv[1]=='register' else run()

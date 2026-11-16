"""Approved MG191 deterministic qualification then conditional gradient production."""
import os
for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[key]='1'
import sys,json,time,hashlib,traceback,shutil,threading,_thread
from pathlib import Path
import numpy as np
from mg191_gravity import EPS,FINE,phase,model_for,reference_data,initialize,evolve_chunk,stability,uniform_stability,measure_series,velocity_blocks,wrap
from mg10e2c1_pool import atomic_json,run as pool_run
from run_mg176 import ROOT
PUBLIC=ROOT/'results/MG191_gravity_gradient';RAW=ROOT/'local_raw/MG191_gravity_gradient'
FILES=('src/run_mg191.py','src/mg191_gravity.py','src/gravity_mg190.py','src/lattice.py','src/runner.py','src/measure.py','src/mg10e2c1_pool.py','spec/MG191_gravity_gradient_draft.md','spec/MG191_implementation_contract.md','results/MG003/ref_profiles.npz')
def read(path):return json.loads(Path(path).read_text(encoding='utf-8'))
def binding():return {p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in FILES}
def bgid(job):return f"{job['arm']}_h{job['h']}_e{job.get('eps',0.):g}"
def storage_bytes():
 total=0
 for directory,_,files in os.walk(RAW):
  for name in files:
   if name.endswith('.tmp'):continue
   try:total+=os.stat(os.path.join(directory,name)).st_size
   except FileNotFoundError:pass
 return total

def register():
 assert read(PUBLIC/'approval.json')['approved_worker_seconds']==7200 and read(PUBLIC/'implementation_checks.json')['passed'];assert not(PUBLIC/'plan.json').exists() and not RAW.exists();RAW.mkdir();(RAW/'inputs').mkdir();jobs=[]
 for h in (.2,.1):
  for j in range(8):jobs.append(dict(id=f'N_h{h}_p{j}',arm='N',h=h,eps=0.,position=phase(j),phase=j,kind='production',duration=2000))
 for arm in ('A','R','C'):
  for eps in EPS:
   for j in range(8):jobs.append(dict(id=f'{arm}_e{eps}_p{j}',arm=arm,h=.2,eps=-eps if arm=='R' else eps,position=-phase(j) if arm=='R' else phase(j),phase=j,kind='production',duration=2000))
 for sign in (-1,1):
  for j in range(8):jobs.append(dict(id=f'L_s{sign}_p{j}',arm='L',h=.2,eps=sign*.001,position=phase(j),phase=j,kind='production',duration=2000))
 for eps in FINE:
  for j in range(8):jobs.append(dict(id=f'H_e{eps}_p{j}',arm='H',h=.1,eps=eps,position=phase(j),phase=j,kind='production',duration=2000))
 assert len(jobs)==120;models={bgid(j):{k:j[k] for k in ('arm','h','eps')} for j in jobs};qualification=[]
 for arm in ('qa','qb'):
  for h in ((.2,.1) if arm=='qa' else (.2,)):
   for eps in (EPS if h==.2 else FINE):
    for sign in (-1,1):
     job=dict(arm=arm,h=h,eps=sign*eps,kind='viability',duration=400,position=0.,id=f'{arm}_h{h}_e{sign*eps:g}');qualification.append(job);models[bgid(job)]={k:job[k] for k in ('arm','h','eps')}
 inputs={}
 for key,config in models.items():
  model=model_for(config);background,residual,iterations=model.steady();gradient=np.fft.ifft(1j*2*np.pi*np.fft.fftfreq(model.n,d=model.h)*np.fft.fft(background.reshape(3,-1)[0])).real;path=RAW/'inputs'/(key+'.npz');np.savez_compressed(path,background=background.reshape(3,-1),a=model.a,b=model.b,gradient=gradient,x=model.x);inputs[key]=dict(config=config,sha256=hashlib.sha256(path.read_bytes()).hexdigest(),residual_inf=residual,iterations=iterations)
  if config['arm'] not in ('qa','qb'):qualification.append(dict(id='bg_'+key,kind='background',duration=0,**config))
 assert len(qualification)==31
 for job in jobs+qualification:job['background_key']=bgid(job)
 atomic_json(PUBLIC/'plan.json',dict(source=binding(),approval=read(PUBLIC/'approval.json'),jobs=jobs,qualification_jobs=qualification,inputs=inputs,qualification_cap=600,production_cap=6000,pre_run_charge=60,analysis_reserve=540,wall_cap=3600,job_cap=600,workers=8,raw_cap=10000000000,forecast_coarse=25,forecast_fine=45,no_stochastic_seeds=True))

def input_for(job,plan):
 key=job['background_key'];path=RAW/'inputs'/(key+'.npz');assert hashlib.sha256(path.read_bytes()).hexdigest()==plan['inputs'][key]['sha256']
 with np.load(path) as f:return {k:f[k].copy() for k in f.files}

def worker(job,directory):
 plan=read(PUBLIC/'plan.json');assert plan['source']==binding() and job in plan['jobs']+plan['qualification_jobs'];path=Path(directory)/job['id'];path.mkdir(exist_ok=False);start=time.perf_counter();clips=0
 try:
  inp=input_for(job,plan);model=model_for(job);_,ref=reference_data(job['h']);background=inp['background'];atomic_json(path/'metadata.json',dict(job=job,reference=ref,source=plan['source'],input_sha256=plan['inputs'][job['background_key']]['sha256']))
  if job['kind']=='background':
   result=stability(model,background.ravel());atomic_json(path/'stability.json',result);atomic_json(Path(directory)/(job['id']+'.json'),dict(status='COMPLETE',job=job,passed=result['passed'],stability=result,seconds=time.perf_counter()-start));return
  uniform=uniform_stability(model) if job['kind']=='viability' else None;current=initialize(job,background);np.savez_compressed(path/'initial.npz',fields=current);previous=job.get('position',0.);lost=False;observations=[];valid=[];times=[]
  for chunk in range(job['duration']//50):
   if time.perf_counter()-start>plan['job_cap']-5:raise TimeoutError('MG191 per-job cap')
   fields=evolve_chunk(model,current,50.);initial=current.copy();current=fields[-1].copy();tt=chunk*50+np.arange(101)*.5;np.savez_compressed(path/f'chunk_{chunk:03d}.npz',fields=fields,times=tt,initial=initial,final=current)
   use=fields if chunk==0 else fields[1:];obs,ok,previous,lost=measure_series(use,background,job['h'],ref,previous,lost);observations.extend(obs);valid.extend(ok.tolist());times.extend((tt if chunk==0 else tt[1:]).tolist());atomic_json(path/'progress.json',dict(model_time=(chunk+1)*50,target=job['duration'],lost_original=bool(lost)))
   if storage_bytes()>plan['raw_cap']:raise RuntimeError('MG191 10GB raw cap')
  observations=np.array(observations);valid=np.array(valid);times=np.array(times);assert len(times)==job['duration']*2+1;np.savez_compressed(path/'observations.npz',times=times,core=observations,original_valid=valid)
  if job['kind']=='viability':
   keep=times>=100;shape=valid&(observations[:,1]>=.6*ref['depth'])&(observations[:,2]>=.75*ref['width'])&(observations[:,2]<=1.25*ref['width']);passed=bool(shape[keep].all() and uniform['passed']);result=dict(status='COMPLETE',job=job,passed=passed,uniform_stability=uniform,valid_shape_frames=int(shape[keep].sum()),required_shape_frames=int(keep.sum()),minimum_depth_ratio=float(np.nanmin(observations[keep,1])/ref['depth']),width_ratio_range=[float(np.nanmin(observations[keep,2])/ref['width']),float(np.nanmax(observations[keep,2])/ref['width'])],seconds=time.perf_counter()-start)
  else:
   result=velocity_blocks(observations,times,valid,ref);np.savez_compressed(path/'track.npz',center_unwrapped=result.pop('center_unwrapped'),shape_valid=result.pop('shape_valid'));gradient=np.interp(observations[:,0],inp['x'],inp['gradient'],period=100);minimum_x=float(inp['x'][np.argmin(background[0])]);distance=abs(wrap(observations[:,0]-minimum_x));near=(distance<=ref['wavelength'])&valid;first=np.flatnonzero(near)
   for block in result['blocks']:
    i=block['left']*2;j=block['right']*2;block['minus_G_gradient']=float(-np.mean(gradient[i:j+1])) if block['valid'] else None
   good=[b for b in result['blocks'] if b['valid']];force=np.array([b['minus_G_gradient'] for b in good]);velocity=np.array([b['velocity'] for b in good]);den=float(force@force);result.update(status='COMPLETE',job=job,mobility=None if den<1e-24 else float(force@velocity/den),median_abs_late_speed=None if not good else float(np.median(abs(velocity))),G_minimum_x=minimum_x,first_near_minimum=None if not len(first) else float(times[first[0]]),near_minimum_censored=not bool(len(first)),final_distance_to_minimum=None if not valid[-1] else float(distance[-1]),seconds=time.perf_counter()-start)
  atomic_json(path/'summary.json',result);atomic_json(Path(directory)/(job['id']+'.json'),result)
 except BaseException:atomic_json(Path(directory)/(job['id']+'.json'),dict(status='FAILED',job=job,seconds=time.perf_counter()-start,error=traceback.format_exc()));raise

def admission(plan,ledger):
 records={r['id']:r for r in ledger['records']};assert len(records)==31 and all(r['status']=='COMPLETE' for r in records.values());admitted=[];excluded=[]
 def ok(arm,h,eps):return records[f'{arm}_h{h}_e{eps:g}']['passed']
 coarse={eps:all(ok(arm,.2,sign*eps) for arm in ('qa','qb') for sign in (-1,1)) for eps in EPS}
 for job in plan['jobs']:
  reasons=[];eps=abs(job['eps'])
  if not records['bg_'+job['background_key']]['passed']:reasons.append('background stability')
  if job['arm'] in ('A','R','C','H') and not coarse[eps]:reasons.append('uniform coarse viability')
  if job['arm']=='H' and not all(ok('qa',.1,s*eps) for s in (-1,1)):reasons.append('uniform fine viability')
  (excluded if reasons else admitted).append(dict(job=job,reasons=reasons) if reasons else job)
 available=[eps for eps in EPS if all(any(j['arm']==arm and abs(j['eps'])==eps and j['phase']==p for j in admitted) for arm in ('A','R','C') for p in range(8))]
 return dict(admitted_jobs=admitted,excluded=excluded,coarse_eligible_eps=available,proceed=len(available)>=2)

def run():
 plan=read(PUBLIC/'plan.json');assert plan['source']==binding() and shutil.disk_usage(ROOT).free>20000000000 and not(PUBLIC/'attempt.json').exists();start=time.perf_counter();atomic_json(PUBLIC/'attempt.json',dict(status='RUNNING',phase='qualification'));timer=threading.Timer(plan['wall_cap'],_thread.interrupt_main);timer.daemon=True;timer.start()
 try:
  q=pool_run(plan['qualification_jobs'],RAW/'qualification',Path(__file__).resolve(),ROOT,plan['qualification_cap'],workers=8,timeout=plan['job_cap'],forecast={j['id']:20 if j['h']==.1 else 5 for j in plan['qualification_jobs']});assert q['charged_worker_seconds']<=600
  if len(q['records'])!=31 or any(r['status']!='COMPLETE' for r in q['records']):atomic_json(PUBLIC/'attempt.json',dict(status='STOPPED_QUALIFICATION',seconds=time.perf_counter()-start,charged_worker_seconds=q['charged_worker_seconds']+plan['pre_run_charge'],stop_reason=q['stop_reason']));return
  allowed=admission(plan,q);atomic_json(PUBLIC/'qualification.json',allowed)
  if not allowed['proceed']:atomic_json(PUBLIC/'attempt.json',dict(status='INSUFFICIENT_ELIGIBILITY',seconds=time.perf_counter()-start,charged_worker_seconds=q['charged_worker_seconds']+plan['pre_run_charge']));return
  atomic_json(PUBLIC/'attempt.json',dict(status='RUNNING',phase='production',qualified_eps=allowed['coarse_eligible_eps']));jobs=allowed['admitted_jobs'];p=pool_run(jobs,RAW/'production',Path(__file__).resolve(),ROOT,plan['production_cap'],workers=8,timeout=plan['job_cap'],forecast={j['id']:plan['forecast_coarse'] if j['h']==.2 else plan['forecast_fine'] for j in jobs});complete=len(p['records'])==len(jobs) and all(r['status']=='COMPLETE' for r in p['records']);atomic_json(PUBLIC/'attempt.json',dict(status='COMPLETE' if complete else 'STOPPED_PRODUCTION',seconds=time.perf_counter()-start,qualified_eps=allowed['coarse_eligible_eps'],planned_production=len(jobs),production_records=len(p['records']),charged_worker_seconds=q['charged_worker_seconds']+p['charged_worker_seconds']+plan['pre_run_charge'],analysis_reserve=plan['analysis_reserve'],stop_reason=p['stop_reason']))
 except BaseException:atomic_json(PUBLIC/'attempt.json',dict(status='FAILED',seconds=time.perf_counter()-start,error=traceback.format_exc()));raise
 finally:timer.cancel()
if __name__=='__main__':
 if sys.argv[1]=='register':register()
 elif sys.argv[1]=='run':run()
 elif sys.argv[1]=='--worker':worker(json.loads(sys.argv[2]),sys.argv[3])
 else:raise ValueError('register/run/--worker')

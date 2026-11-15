"""MG185 proposed Phase 2 study; approval-gated, immutable legacy kernel."""
import os
for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[k]='1'
import sys,time,json,hashlib,traceback,shutil,threading,_thread
from pathlib import Path
import numpy as np
from scipy.stats import t
from mg184_phase2 import ref,initial,advance,summarize,TAU,OMEGAS
from runner import seed_for
from mg10e2c1_pool import atomic_json,run as pool_run
from run_mg176 import ROOT
PUBLIC=ROOT/'results/MG185_Phase2_G';RAW=ROOT/'local_raw/MG185_Phase2_G';DEV=ROOT/'results/MG12H3/phase2_G_development'
FILES=['src/run_mg185_phase2_G.py','src/score_mg185_phase2_G.py','src/mg184_phase2.py','src/lattice.py','src/measure.py','src/runner.py','src/mg10e2c1_pool.py','spec/MG185_phase2_G_draft.md','results/MG003/ref_profiles.npz']
def read(p):return json.loads(Path(p).read_text(encoding='utf-8'))
def binding():return {p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in FILES}
def register():
 if not(PUBLIC/'approval.json').exists():raise RuntimeError('PI approval of the concrete MG185 criteria/budget is required before scientific registration.')
 a=read(PUBLIC/'approval.json');assert a['prepared_per_omega']==30 and a['controls_per_omega']==6 and a['duration']==1000 and a['total_worker_cap']==360000;assert read(DEV/'separate_check.json')['passed'] and not(PUBLIC/'plan.json').exists();jobs=[];used=set()
 def seeds(value):
  if type(value)is int:used.add(value)
  elif isinstance(value,(dict,list)):
   for x in (value.values() if isinstance(value,dict) else value):seeds(x)
 def visit(value):
  if isinstance(value,dict):
   for k,x in value.items():seeds(x) if 'seed' in k.lower() else visit(x)
  elif isinstance(value,list):
   for x in value:visit(x)
 for p in (ROOT/'results').rglob('*.json'):visit(read(p))
 for phase in('qualification','main'):
  for omega in ((3000,) if phase=='qualification' else OMEGAS):
   for di,tau in enumerate((TAU,TAU/2) if phase=='qualification' else (TAU,)):
    for ci,condition in enumerate(('prepared','homogeneous')):
     for r in range(12 if phase=='qualification' else 30 if ci==0 else 6):
      label=f'{phase}_Omega{omega}_d{di}_{condition}_r{r:02d}';T=100 if phase=='qualification' else 1000;seeds_out=[]
      for segment in range(T//10+1):
       nonce=0;s=seed_for(185,omega,ci,di+(0 if phase=='qualification' else 10),r,segment)
       while s in used:nonce+=1;s=seed_for(185,omega,ci,100+di+nonce*1000,r,segment)
       used.add(s);seeds_out.append(s)
      jobs.append(dict(id=label,phase=phase,omega=omega,tau=tau,condition=condition,replicate=r,duration=T,seed=seeds_out[0],segment_seeds=seeds_out[1:]))
 assert len(jobs)==228;atomic_json(PUBLIC/'plan.json',dict(jobs=jobs,source=binding(),approval=a,workers=8,qualification_cap=28800,main_cap=331200,qualification_wall_cap=7200,main_wall_cap=50400,job_cap=7200,storage_cap=10000000000,forecast=read(DEV/'proposal.json')['job_forecasts'],no_bridge_dependency=True))
def jsonsafe(value):return json.loads(json.dumps(value,default=lambda x:x.tolist() if isinstance(x,np.ndarray) else x.item()))
def worker(job,directory):
 plan=read(PUBLIC/'plan.json');assert plan['source']==binding() and job in plan['jobs'];path=Path(directory)/job['id'];path.mkdir(exist_ok=False);start=time.perf_counter();r=ref();pop=initial(job['omega'],job['condition'],job['seed'],r);frames=[pop[None,:,:].copy()];alltimes=[np.array([0.])];clips=0
 try:
  np.savez_compressed(path/'initial.npz',counts=pop,reference_phi=r['phi'],reference_S0=r['S0'],omega=job['omega']);atomic_json(path/'metadata.json',dict(job=job,source=binding(),kernel='unchanged tau_leap_euler_1d',stream_layout='independent registered seeds per10modeltimeunits'))
  for chunk,seed in enumerate(job['segment_seeds']):
   if time.perf_counter()-start>plan['job_cap']-30:raise TimeoutError('Per-job cap')
   previous=pop.copy();out=advance(pop,job['omega'],10.,chunk*10.,seed,tau=job['tau']);logical=chunk*10+np.arange(1,21)*.5;np.testing.assert_allclose(out['times'],logical,rtol=0,atol=1e-7);assert out['counts'].shape==(20,3,500) and out['counts'].min()>=0 and np.array_equal(out['final'],out['counts'][-1]);pop=out['final'];clips+=out['clips'];np.savez_compressed(path/f'chunk_{chunk:03d}.npz',initial=previous,counts=out['counts'],times=logical,kernel_times=out['times'],final=pop,seed=seed,clips=out['clips']);frames.append(out['counts']);alltimes.append(logical);atomic_json(path/'progress.json',dict(model_time=(chunk+1)*10,target=job['duration'],chunks=chunk+1,clips=clips))
   if sum(p.stat().st_size for p in RAW.rglob('*') if p.is_file())>plan['storage_cap']:raise RuntimeError('Local10GBcap')
   if clips:
    atomic_json(Path(directory)/(job['id']+'.json'),dict(status='NUMERICAL_UNQUALIFIED',job=job,seconds=time.perf_counter()-start,clips=clips,completed_model_time=(chunk+1)*10,physical_failure_claim=False));return
  counts=np.concatenate(frames);times=np.concatenate(alltimes);summary,arrays=summarize(counts,times,job['omega'],job['condition'],job['duration'],r);np.savez_compressed(path/'GXY_measurements.npz',**arrays);atomic_json(path/'summary.json',jsonsafe(summary));numeric=float(arrays['fixed_projection'][(times>=50)&(times<=100),0].mean()) if job['phase']=='qualification' else None;atomic_json(Path(directory)/(job['id']+'.json'),dict(status='COMPLETE',job=job,seconds=time.perf_counter()-start,clips=clips,historical_X_persistence=summary['historical_X_persistence'],strict_sampled_unique_X_track=summary['strict_sampled_unique_X_track'],end_ambiguous=summary['end_unmatched_or_ambiguous'],numeric_G_block_mean=numeric))
 except BaseException:atomic_json(Path(directory)/(job['id']+'.json'),dict(status='FAILED',job=job,seconds=time.perf_counter()-start,clips=clips,error=traceback.format_exc()));raise

def qualification_gate(directory):
 records=read(Path(directory)/'ledger.json')['records'];assert len(records)==48 and all(r['status']=='COMPLETE' and r['clips']==0 for r in records);margin=.05*np.linalg.norm(ref()['phi'][0]-ref()['phi'][0].mean());out=[]
 for condition in('prepared','homogeneous'):
  a=np.array([r['numeric_G_block_mean'] for r in records if r['job']['condition']==condition and r['job']['tau']==TAU]);b=np.array([r['numeric_G_block_mean'] for r in records if r['job']['condition']==condition and r['job']['tau']==TAU/2]);assert len(a)==len(b)==12;va=a.var(ddof=1)/12;vb=b.var(ddof=1)/12;df=(va+vb)**2/(va*va/11+vb*vb/11) if va+vb else 22;mean=float(b.mean()-a.mean());half=float(t.ppf(.975,df)*np.sqrt(va+vb));ci=[mean-half,mean+half];out.append(dict(condition=condition,mean=mean,CI95=ci,margin=float(margin),classification='EQUIVALENT_IN_SCOPE' if ci[0]>=-margin and ci[1]<=margin else 'INCONCLUSIVE_OR_DIFFERENT'))
 return dict(passed=all(r['classification']=='EQUIVALENT_IN_SCOPE' for r in out),checks=out,does_not_validate_T1000_all_observables=True)
def run_phase(phase):
 plan=read(PUBLIC/'plan.json');assert plan['source']==binding();flag=PUBLIC/(phase+'_attempt.json');assert not flag.exists()
 if phase=='main':assert read(PUBLIC/'qualification.json')['passed']
 directory=RAW/phase;directory.mkdir(parents=True,exist_ok=True);cap=plan['qualification_cap' if phase=='qualification' else 'main_cap'];wall=plan['qualification_wall_cap' if phase=='qualification' else 'main_wall_cap'];jobs=[j for j in plan['jobs'] if j['phase']==phase];forecasts={j['id']:plan['forecast'][phase][str(j['omega'])]*j['duration']*(TAU/j['tau']) for j in jobs};start=time.perf_counter();atomic_json(flag,dict(status='RUNNING'));stopped=threading.Event()
 def interrupt():stopped.set();_thread.interrupt_main()
 timer=threading.Timer(wall,interrupt);timer.daemon=True;timer.start()
 try:
  ledger=pool_run(jobs,directory,Path(__file__).resolve(),ROOT,cap,workers=8,timeout=plan['job_cap'],forecast=forecasts);done=len(ledger['records'])==len(jobs) and all(r['status']=='COMPLETE' for r in ledger['records']);atomic_json(flag,dict(status='COMPLETE' if done else 'STOPPED_OR_UNQUALIFIED',wall_seconds=time.perf_counter()-start,charged_worker_seconds=ledger['charged_worker_seconds'],stop_reason=ledger['stop_reason']))
  if phase=='qualification' and done:atomic_json(PUBLIC/'qualification.json',qualification_gate(directory))
 except BaseException:atomic_json(flag,dict(status='STOPPED' if stopped.is_set() else 'FAILED',wall_seconds=time.perf_counter()-start,error=traceback.format_exc()));raise
 finally:timer.cancel()
def run():
 assert not RAW.exists() and __import__('shutil').disk_usage(ROOT).free>20000000000;run_phase('qualification')
 if (PUBLIC/'qualification.json').exists() and read(PUBLIC/'qualification.json')['passed']:run_phase('main')
if __name__=='__main__':
 if sys.argv[1]=='register':register()
 elif sys.argv[1]=='run':run()
 elif sys.argv[1]=='--worker':worker(json.loads(sys.argv[2]),sys.argv[3])
 else:raise ValueError('register/run/--worker')

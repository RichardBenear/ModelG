"""MG185 authorized infrastructure continuation; original artifacts stay immutable."""
import os
for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[k]='1'
import sys,json,time,hashlib,traceback,threading,_thread
from pathlib import Path
import numpy as np
from run_mg185_phase2_G import ROOT,PUBLIC,RAW,read,binding,jsonsafe
from mg184_phase2 import initial,ref,advance,summarize
from mg10e2c1_pool import atomic_json,run as pool_run
PLAN=PUBLIC/'continuation_plan.json'

def storage_bytes(root):
 total=0
 for directory,_,files in os.walk(root):
  for name in files:
   if name.endswith('.tmp'):continue
   try:total+=os.stat(os.path.join(directory,name)).st_size
   except FileNotFoundError:pass # Atomic replacement can remove an enumerated name.
 return total

def prefix(job):
 directory=RAW/'main'/job['id'];pop=initial(job['omega'],job['condition'],job['seed']);chunks=[]
 paths=sorted(directory.glob('chunk_*.npz'))
 for i,path in enumerate(paths):
  assert path.name==f'chunk_{i:03d}.npz'
  with np.load(path) as a:
   assert np.array_equal(pop,a['initial']) and int(a['seed'])==job['segment_seeds'][i] and int(a['clips'])==0
   assert a['counts'].shape==(20,3,500) and a['counts'].min()>=0 and np.array_equal(a['final'],a['counts'][-1])
   np.testing.assert_allclose(a['times'],i*10+np.arange(1,21)*.5,atol=1e-7,rtol=0)
   np.testing.assert_allclose(a['kernel_times'],a['times'],atol=1e-7,rtol=0)
   pop=a['final'].copy()
  chunks.append(dict(path=str(path.relative_to(ROOT)),sha256=hashlib.sha256(path.read_bytes()).hexdigest()))
 return chunks

def register():
 assert not PLAN.exists();original=read(PUBLIC/'plan.json');assert binding()==original['source'];ledger=read(RAW/'main/ledger.json');complete={r['id'] for r in ledger['records'] if r['status']=='COMPLETE'}
 jobs=[j for j in original['jobs'] if j['phase']=='main' and j['id'] not in complete];prefixes={j['id']:prefix(j) for j in jobs}
 assert len(jobs)==144 and len(complete)==36
 remaining=original['main_cap']-ledger['charged_worker_seconds'];wall=original['main_wall_cap']-read(PUBLIC/'main_attempt.json')['wall_seconds']
 atomic_json(PLAN,dict(approval='PI: Do the fix and restart',original_source=original['source'],resume_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),jobs=jobs,prefixes=prefixes,prior_charge=ledger['charged_worker_seconds'],remaining_worker_seconds=remaining,remaining_wall_seconds=wall,completed_ids=sorted(complete),original_ledger_sha256=hashlib.sha256((RAW/'main/ledger.json').read_bytes()).hexdigest(),criteria_unchanged=True))

def verify(plan):
 assert plan['original_source']==binding() and plan['resume_sha256']==hashlib.sha256(Path(__file__).read_bytes()).hexdigest()

def worker(job,directory):
 plan=read(PLAN);verify(plan);assert job in plan['jobs'];start=time.perf_counter();path=Path(directory)/job['id'];path.mkdir(exist_ok=False);r=ref();pop=initial(job['omega'],job['condition'],job['seed'],r);frames=[pop[None].copy()];times=[np.array([0.])];clips=0
 try:
  for entry in plan['prefixes'][job['id']]:
   source=ROOT/entry['path'];assert hashlib.sha256(source.read_bytes()).hexdigest()==entry['sha256']
   with np.load(source) as a:frames.append(a['counts'].copy());times.append(a['times'].copy());pop=a['final'].copy()
  first=len(plan['prefixes'][job['id']]);atomic_json(path/'progress.json',dict(model_time=first*10,target=1000,resumed_chunks=first,clips=0))
  prior=next((v['charged_worker_seconds'] for v in read(RAW/'main/ledger.json')['records'] if v['id']==job['id']),0)
  for chunk in range(first,100):
   if time.perf_counter()-start+prior>7170:raise TimeoutError('Cumulative per-job cap')
   seed=job['segment_seeds'][chunk];previous=pop.copy();a=advance(pop,job['omega'],10.,chunk*10.,seed,tau=job['tau']);logical=chunk*10+np.arange(1,21)*.5
   np.testing.assert_allclose(a['times'],logical,atol=1e-7,rtol=0);assert a['counts'].shape==(20,3,500) and a['counts'].min()>=0 and np.array_equal(a['final'],a['counts'][-1]);pop=a['final'];clips+=a['clips']
   np.savez_compressed(path/f'chunk_{chunk:03d}.npz',initial=previous,counts=a['counts'],times=logical,kernel_times=a['times'],final=pop,seed=seed,clips=a['clips']);frames.append(a['counts']);times.append(logical);atomic_json(path/'progress.json',dict(model_time=(chunk+1)*10,target=1000,resumed_chunks=first,clips=clips))
   if storage_bytes(RAW)>10000000000:raise RuntimeError('Local10GBcap')
   if clips:
    atomic_json(Path(directory)/(job['id']+'.json'),dict(status='NUMERICAL_UNQUALIFIED',job=job,clips=clips,physical_failure_claim=False));return
  summary,arrays=summarize(np.concatenate(frames),np.concatenate(times),job['omega'],job['condition'],1000,r);np.savez_compressed(path/'GXY_measurements.npz',**arrays);atomic_json(path/'summary.json',jsonsafe(summary));atomic_json(Path(directory)/(job['id']+'.json'),dict(status='COMPLETE',job=job,clips=clips,seconds=time.perf_counter()-start,historical_X_persistence=summary['historical_X_persistence'],strict_sampled_unique_X_track=summary['strict_sampled_unique_X_track'],end_ambiguous=summary['end_unmatched_or_ambiguous'],summary_path=str((path/'summary.json').relative_to(ROOT))))
 except BaseException:atomic_json(Path(directory)/(job['id']+'.json'),dict(status='FAILED',job=job,clips=clips,error=traceback.format_exc()));raise

def run():
 plan=read(PLAN);verify(plan);directory=RAW/'continuation';assert not directory.exists();assert read(PUBLIC/'qualification.json')['passed'];original=read(PUBLIC/'plan.json');forecasts={j['id']:original['forecast']['main'][str(j['omega'])]*(1000-10*len(plan['prefixes'][j['id']])) for j in plan['jobs']};flag=PUBLIC/'continuation_attempt.json';assert not flag.exists();start=time.perf_counter();atomic_json(flag,dict(status='RUNNING'));timer=threading.Timer(plan['remaining_wall_seconds'],_thread.interrupt_main);timer.daemon=True;timer.start()
 try:
  ledger=pool_run(plan['jobs'],directory,Path(__file__).resolve(),ROOT,original['main_cap'],workers=8,prior=plan['prior_charge'],timeout=7200,forecast=forecasts)
  atomic_json(flag,dict(status='COMPLETE' if len(ledger['records'])==144 and all(r['status']=='COMPLETE' for r in ledger['records']) else 'STOPPED_OR_UNQUALIFIED',charged_worker_seconds_including_original=ledger['charged_worker_seconds'],wall_seconds=time.perf_counter()-start,stop_reason=ledger['stop_reason']))
 finally:timer.cancel()
if __name__=='__main__':
 if sys.argv[1]=='register':register()
 elif sys.argv[1]=='run':run()
 elif sys.argv[1]=='--worker':worker(json.loads(sys.argv[2]),sys.argv[3])
 else:raise ValueError('register/run/--worker')

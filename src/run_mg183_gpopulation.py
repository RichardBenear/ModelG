"""MG183 proposed G population pilot; approval-gated, no automatic launch."""
import os
for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[k]='1'
import sys,time,json,hashlib,traceback,shutil
from pathlib import Path
import numpy as np
from mg182_gscale import setup,projection,ROOT
from mg12h3_pattern_persistence import homogenize
from mg12h3_contrast_estimator import Contrast,tick,aggregate_audit
from mg12h3_tick_coupling import pair,replay
from mg12h3_finer_precision import ordered_audit
from mg10e2c1_pool import atomic_json,run as pool_run
from runner import seed_for
from prepare_mg182 import OUT as DEVELOPMENT
PUBLIC=ROOT/'results/MG12H3/G_population_pilot';RAW=ROOT/'local_raw/MG183_G_population'
FILES=['src/run_mg183_gpopulation.py','src/score_mg183_gpopulation.py','spec/MG183_G_population_pilot_draft.md','src/mg182_gscale.py','src/mg12h3_sparse.py','src/mg12h3_contrast_estimator.py','src/mg12h3_tick_coupling.py','src/mg12h3_pattern_persistence.py','src/mg12h3_finer_precision.py','src/mg10e2c1_pool.py']
def binding():return {p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in FILES}
def read(p):return json.loads(Path(p).read_text(encoding='utf-8'))
def register():
 approval=PUBLIC/'approval.json'
 if not approval.exists():raise RuntimeError('Scientific production requires PI approval of the concrete readiness criteria and budget; preparation does not authorize this run.')
 p=read(approval);assert p['duration']==.01 and p['worker_cap']==86400 and p['workers']==8 and p['replicates_per_condition']==4;assert read(DEVELOPMENT/'separate_check.json')['passed'];assert not(PUBLIC/'plan.json').exists();jobs=[]
 for omega in(160,40,10):
  for ci,condition in enumerate(('prepared','homogeneous')):
   for r in range(4):jobs.append(dict(id=f'Omega{omega}_{condition}_r{r}',omega=omega,condition=condition,replicate=r,seed=seed_for(12,omega,ci,18300,r),paired=omega==160,sentinel=r==0))
 assert len({j['seed'] for j in jobs})==24;atomic_json(PUBLIC/'plan.json',dict(jobs=jobs,source=binding(),approval=p,duration=.01,delta=.0001,fine_delta=.00002,observations=.001,worker_cap=86400,job_cap=10800,workers=8,storage_cap=2000000000,forecast=read(DEVELOPMENT/'proposal.json')['job_forecast_seconds']))
def worker(job,directory):
 directory=Path(directory);plan=read(PUBLIC/'plan.json');assert plan['source']==binding() and job in plan['jobs'];path=directory/job['id'];path.mkdir(exist_ok=False);start=time.perf_counter();deadline=start+plan['job_cap']-10
 try:
  c=setup(job['omega']);initial=c['pop'] if job['condition']=='prepared' else homogenize(c['pop']);model=c['model'];u=c['reference'][0]-c['reference'][0].mean();u-=u.mean();u/=np.linalg.norm(u);u-=u.mean();weights=np.zeros((1,3500));weights[0,:500]=u/job['omega'];contrast=Contrast(model,weights);rng=np.random.default_rng(job['seed']);pop=np.stack([initial,initial]) if job['paired'] else initial[None,:,:];frames=[pop.copy()];events=0
  np.savez_compressed(path/'initial.npz',counts=initial,reference=c['reference'],omega=job['omega']);atomic_json(path/'metadata.json',dict(job=job,source=binding(),condition_mass_per_species=initial.sum(axis=1).tolist()))
  for step in range(100):
   previous=pop.copy();before=json.dumps(rng.bit_generator.state);sentinel=job['sentinel'] and step==0
   if job['paired']:
    out=pair(model,*pop,[.0001,.00002],[0,.0001],job['seed'],deadline=deadline,rng=rng,sentinel=sentinel,max_events=2000000);pop=out['counts'][:,-1];hist=out['histogram']
    for side in(0,1):
     flat=previous[side].ravel().copy()
     for i,ch in enumerate(model.channels):
      for coordinate,v in ch.change:flat[coordinate]+=int(hist[0,i]+hist[side+1,i])*v
     np.testing.assert_array_equal(flat,pop[side].ravel())
    if sentinel:
     replay(model,*previous,[.0001,.00002],[0,.0001],out);import gzip;(path/'sentinel.json.gz').write_bytes(gzip.compress(json.dumps(out,default=lambda x:x.tolist() if isinstance(x,np.ndarray) else x.item()).encode(),mtime=0))
    available=out['final_available_before_reset']
   else:
    out=tick(contrast,pop[0],.0001,rng,deadline=deadline,max_events=2000000,sentinel=sentinel);aggregate_audit(contrast,pop[0],out)
    if sentinel:
     ordered_audit(contrast,pop[0],out,.0001);np.savez_compressed(path/'sentinel.npz',ids=[e['channel'] for e in out['trace']],marks=[e['time'] for e in out['trace']])
    pop=out['counts'][None,:,:];hist=out['histogram'];available=out['available'][None,:,:]
   assert all(int(x.sum())==int(initial.sum()) for x in pop) and pop.min()>=0;events+=out['events'];np.savez_compressed(path/f'checkpoint_{step:03d}.npz',initial=previous,final=pop,available=available,histogram=hist,rng_before=before,rng_after=json.dumps(rng.bit_generator.state),events=out['events'])
   if (step+1)%10==0:frames.append(pop.copy())
   if sum(p.stat().st_size for p in RAW.rglob('*') if p.is_file())>plan['storage_cap']:raise RuntimeError('2GBlocalrawcap')
   atomic_json(path/'progress.json',dict(ticks=step+1,required=100,events=events))
  values=np.array(frames);np.savez_compressed(path/'observations.npz',counts=values,times=np.arange(11)*.001,omega=job['omega']);atomic_json(directory/(job['id']+'.json'),dict(status='COMPLETE',job=job,seconds=time.perf_counter()-start,events=events))
 except BaseException:atomic_json(directory/(job['id']+'.json'),dict(status='FAILED',job=job,seconds=time.perf_counter()-start,error=traceback.format_exc()));raise

def run():
 p=read(PUBLIC/'plan.json');assert p['source']==binding() and not(PUBLIC/'attempt.json').exists() and not RAW.exists();assert shutil.disk_usage(ROOT).free>8000000000;RAW.mkdir();atomic_json(PUBLIC/'attempt.json',dict(status='RUNNING'));start=time.perf_counter()
 import threading,_thread
 wall_stop=threading.Event()
 def interrupt_wall():wall_stop.set();_thread.interrupt_main()
 timer=threading.Timer(10800,interrupt_wall);timer.daemon=True;timer.start()
 try:ledger=pool_run(p['jobs'],RAW/'production',Path(__file__).resolve(),ROOT,p['worker_cap'],workers=8,timeout=p['job_cap'],forecast={j['id']:p['forecast'][str(j['omega'])] for j in p['jobs']})
 except KeyboardInterrupt:
  ledger=read(RAW/'production'/'ledger.json');atomic_json(PUBLIC/'attempt.json',dict(status='STOPPED',elapsed_seconds=time.perf_counter()-start,charged_worker_seconds=ledger['charged_worker_seconds'],stop_reason='WALL CAP' if wall_stop.is_set() else 'INTERRUPTED'));raise
 finally:timer.cancel()
 done=len(ledger['records'])==24 and all(x['status']=='COMPLETE' and x['returncode']==0 for x in ledger['records']);atomic_json(PUBLIC/'attempt.json',dict(status='COMPLETE' if done else 'STOPPED',elapsed_seconds=time.perf_counter()-start,charged_worker_seconds=ledger['charged_worker_seconds'],stop_reason=ledger['stop_reason']))
if __name__=='__main__':
 if sys.argv[1]=='--worker':worker(json.loads(sys.argv[2]),sys.argv[3])
 elif sys.argv[1]=='register':register()
 elif sys.argv[1]=='run':run()
 else:raise ValueError('Choose register/run')

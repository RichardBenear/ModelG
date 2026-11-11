"""MG140 auditable multi-tick wrapper of frozen sparse Workspace."""
import os
for key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[key]='1'
from pathlib import Path
import numpy as np,json,time,hashlib,math,sys,traceback
from datetime import datetime,timezone
from run_mg12h3_sparse_combined import ROOT,load
from mg10e2c1_pool import atomic_json
from runner import seed_for
BASE=ROOT/'results/MG12H3/multitick_pilot'
FILES=['src/mg12h3_multitick.py','src/prepare_mg12h3_multitick.py','src/run_mg12h3_multitick.py','spec/MG12H3_multitick_pilot.md','src/mg12h3_sparse.py','src/run_mg12h3_sparse_combined.py','src/mg10e2c1_pool.py']
def hashes():
 h={p:hashlib.sha256((ROOT/p).read_text(encoding='utf-8').replace('\r\n','\n').encode()).hexdigest() for p in FILES}
 for p in (ROOT/'results/MG12H3/ring_feasibility').glob('*.npz'):h[p.relative_to(ROOT).as_posix()]=hashlib.sha256(p.read_bytes()).hexdigest()
 return h

def jobs():
 cases=sorted(p.stem for p in (ROOT/'results/MG12H3/ring_feasibility').glob('*.npz'))
 return [dict(id=case+'_d'+str(level),case=case,level=level,delta=.001/(level+1),duration=.01) for case in cases for level in range(2)]

def audit_tick(model,initial,final,available,ids,marks,delta):
 assert len(ids)==len(marks) and np.all(np.diff(marks)>0) and np.all(marks>0) and np.all(marks<delta)
 counts=np.bincount(ids,minlength=len(model.channels));reserve=np.zeros(initial.size,np.int64);change=reserve.copy()
 for c,n in zip(model.channels,counts):
  if n:
   for i,v in c.reserve:reserve[i]+=int(n)*v
   for i,v in c.change:change[i]+=int(n)*v
 np.testing.assert_array_equal(initial.ravel()-reserve,available.ravel());np.testing.assert_array_equal(initial.ravel()+change,final.ravel());assert available.min()>=0 and final.min()>=0 and final.sum()==initial.sum()

def run_job(job,seed,directory,deadline):
 directory=Path(directory);directory.mkdir(parents=True,exist_ok=True);path=directory/job['id'];path.mkdir(exist_ok=False);clock=time.perf_counter();started=datetime.now(timezone.utc).isoformat();pop,model=load(ROOT/f'results/MG12H3/ring_feasibility/{job["case"]}.npz');nsteps=round(job['duration']/job['delta']);assert abs(nsteps*job['delta']-job['duration'])<1e-14;rng=np.random.default_rng(seed);summaries=[];frames=[pop.copy()];times=[0.];events=0
 for step in range(nsteps):
  tickstart=time.perf_counter();initial=pop.copy();work=model.workspace(initial);ids=[];marks=[];checks=[];local=0.
  def check_rates():
   fresh=np.array([c.rate(work.available) for c in model.channels]);total=math.fsum(fresh);np.testing.assert_allclose(work.tree.values,fresh,rtol=1e-10,atol=1e-6);np.testing.assert_allclose(work.tree.total(),total,rtol=1e-10,atol=1e-6)
   if work.positive:
    prefix=np.cumsum(fresh)
    for u in (1e-9,.5,1-1e-9):assert work.tree.select(u)==int(np.searchsorted(prefix,u*total,side='right'))
   checks.append(dict(events=len(ids),total=total,tree_total=work.tree.total(),max_leaf_error=float(abs(fresh-work.tree.values).max())))
  try:
   while work.positive:
    if time.perf_counter()>deadline:raise TimeoutError('Multi-tick job cap')
    total=work.tree.total();assert total>0 and math.isfinite(total);wait=float(rng.exponential(1/total));assert wait>0 and math.isfinite(wait);local+=wait
    if local>=job['delta']:break
    if len(ids)>=1000000:raise ValueError('Tick event cap')
    index=work.tree.select(float(rng.random()));work.accept(index);ids.append(index);marks.append(local)
    if len(ids)%10000==0:check_rates()
   check_rates()
  finally:
   final=work.output.reshape(initial.shape);available=work.available.reshape(initial.shape);np.savez_compressed(path/f'tick_{step:04d}.npz',initial=initial,final=final,available=available,channel_ids=np.asarray(ids,np.int32),marks=np.asarray(marks),step=step,delta=job['delta'],seed=seed,job=json.dumps(job));atomic_json(path/f'tick_{step:04d}_checks.json',dict(records=checks))
  audit_tick(model,initial,final,available,np.asarray(ids,np.int32),np.asarray(marks),job['delta']);pop=final;events+=len(ids);summaries.append(dict(step=step,events=len(ids),wall_seconds=time.perf_counter()-tickstart));atomic_json(path/'progress.json',dict(ticks_complete=len(summaries),ticks_required=nsteps,events=events))
  if (step+1)%(job['level']+1)==0:frames.append(pop.copy());times.append((step+1)*job['delta'])
 np.savez_compressed(path/'frames.npz',counts=np.asarray(frames),times=times,job=json.dumps(job),seed=seed);meta=dict(status='COMPLETE',job=job,seed=seed,ticks=nsteps,events=events,elapsed_seconds=time.perf_counter()-clock,started=started,finished=datetime.now(timezone.utc).isoformat(),tick_summaries=summaries);atomic_json(directory/(job['id']+'.json'),meta);return meta

def audit_job(job,seed,directory):
 directory=Path(directory);path=directory/job['id'];meta=json.loads((directory/(job['id']+'.json')).read_text());assert meta['status']=='COMPLETE' and meta['job']==job and meta['seed']==seed;pop,model=load(ROOT/f'results/MG12H3/ring_feasibility/{job["case"]}.npz');frames=[pop.copy()];events=0;nsteps=round(job['duration']/job['delta'])
 for step in range(nsteps):
  with np.load(path/f'tick_{step:04d}.npz') as z:
   assert int(z['seed'])==seed and int(z['step'])==step and float(z['delta'])==job['delta'] and json.loads(str(z['job']))==job;np.testing.assert_array_equal(pop,z['initial']);audit_tick(model,pop,z['final'],z['available'],z['channel_ids'],z['marks'],job['delta']);pop=z['final'].copy();events+=len(z['channel_ids'])
  records=json.loads((path/f'tick_{step:04d}_checks.json').read_text())['records'];assert records and records[-1]['events']==meta['tick_summaries'][step]['events'];assert all(abs(x['tree_total']-x['total'])<=1e-6+1e-10*abs(x['total']) for x in records)
  if (step+1)%(job['level']+1)==0:frames.append(pop.copy())
 with np.load(path/'frames.npz') as z:np.testing.assert_array_equal(z['counts'],frames);np.testing.assert_allclose(z['times'],np.arange(len(frames))*.001,rtol=0,atol=1e-14);assert int(z['seed'])==seed
 assert events==meta['events'];return dict(passed=True,ticks=nsteps,events=events)

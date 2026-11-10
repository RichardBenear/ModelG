"""MG136 combined bounded development; no formal qualification."""
import os
for key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[key]='1'
from pathlib import Path
import time,json,math,hashlib,traceback,sys
import numpy as np
from mg12h3_sparse import Model
from prepare_mg12h3_bridge import fixture,config,ROOT,seed_for
BASE=ROOT/'results/MG12H3/sparse_combined'
FILES=['src/run_mg12h3_sparse_combined.py','src/mg12h3_sparse.py','src/prepare_mg12h3_bridge.py','spec/MG12H3_sparse_combined_development.md']
def hashes():
 h={p:hashlib.sha256((ROOT/p).read_text(encoding='utf-8').replace('\r\n','\n').encode()).hexdigest() for p in FILES}
 for p in (ROOT/'results/MG12H3/ring_feasibility').glob('*.npz'):h[p.relative_to(ROOT).as_posix()]=hashlib.sha256(p.read_bytes()).hexdigest()
 for c in range(4):
  p=ROOT/f'results/MG12H3/mixed_qualification/references/c{c}_reference.npz';h[p.relative_to(ROOT).as_posix()]=hashlib.sha256(p.read_bytes()).hexdigest()
 return h

def load(p):
 with np.load(p) as z:pop=z['counts'];cfg=(float(z['omega']),list(zip(z['a0'],z['b0'])),json.loads(str(z['links'])),tuple(z['diff']),tuple(z['coupling']),(True,)*6)
 return pop,Model(pop.shape[1],cfg)

def save(name,value):
 BASE.mkdir(parents=True,exist_ok=True);(BASE/name).write_text(json.dumps(value,indent=2)+'\n')
def plan():
 assert not (BASE/'plan.json').exists();cases=sorted(p.stem for p in (ROOT/'results/MG12H3/ring_feasibility').glob('*.npz'));tiny=[[seed_for(12,c,0,13600,r) for r in range(256)] for c in range(4)];actual=[seed_for(12,c,0,13601,0) for c in range(8)];values=[v for a in tiny for v in a]+actual;assert len(values)==len(set(values))
 prior=set()
 def collect(x):
  if isinstance(x,dict):
   for v in x.values():collect(v)
  elif isinstance(x,list):
   for v in x:collect(v)
  elif type(x) is int and 0<=x<2**32:prior.add(x)
 for p in (ROOT/'results/MG12H3').rglob('registration.json'):
  q=json.loads(p.read_text());collect(q.get('seeds',{}))
 assert not set(values)&prior;save('plan.json',dict(source_hashes=hashes(),tiny_seeds=tiny,actual_seeds=actual,cases=cases,cap=180,actual_cap=20,prior_seed_count=len(prior),authorization='PI go ahead and combine batches and test; bounded exploratory development, not formal qualification'))

def execute_tick(model,pop,delta,rng,deadline):
 work=model.workspace(pop);clock=0.;ids=[];times=[]
 while work.positive:
  if time.perf_counter()>deadline:raise TimeoutError('Tick deadline')
  total=work.tree.total();assert total>0 and math.isfinite(total);wait=float(rng.exponential(1/total));assert wait>0 and math.isfinite(wait);clock+=wait
  if clock>=delta:break
  if len(ids)>=1000000:raise ValueError('Event cap')
  index=work.tree.select(float(rng.random()));work.accept(index);ids.append(index);times.append(clock)
 return work,np.asarray(ids,np.int32),np.asarray(times,float)

def audit(model,pop,ids,times,final,available,delta,deadline):
 a=pop.copy().ravel();out=pop.copy().ravel();assert len(ids)==len(times) and np.all(np.diff(times)>0) and np.all(times<delta) and np.all(times>0)
 for step,index in enumerate(ids):
  if step%1024==0 and time.perf_counter()>deadline:raise TimeoutError('Replay deadline')
  c=model.channels[int(index)];assert c.rate(a)>0
  for i,n in c.reserve:assert a[i]>=n;a[i]-=n
  for i,v in c.change:out[i]+=v
 np.testing.assert_array_equal(a,available.ravel());np.testing.assert_array_equal(out,final.ravel());assert out.min()>=0 and out.sum()==pop.sum()

def run():
 p=json.loads((BASE/'plan.json').read_text());assert p['source_hashes']==hashes() and not (BASE/'attempt.json').exists();save('attempt.json',dict(status='STARTED'));start=time.perf_counter();deadline=start+180;stress=[];tiny=[];actual=[]
 for case in p['cases']:
  pop,model=load(ROOT/f'results/MG12H3/ring_feasibility/{case}.npz');work=model.workspace(pop)
  for step in range(10000):
   if time.perf_counter()>deadline:raise TimeoutError('Combined cap')
   index=work.tree.select(((step+1)*.6180339887498949)%1);work.accept(index)
   if step%1000==999:
    fresh=np.array([c.rate(work.available) for c in model.channels]);np.testing.assert_allclose(work.tree.values,fresh,rtol=1e-10,atol=1e-6);np.testing.assert_allclose(work.tree.total(),math.fsum(fresh),rtol=1e-10,atol=1e-6)
    for u in (1e-9,.5,1-1e-9):assert work.tree.select(u)==int(np.searchsorted(np.cumsum(fresh),u*work.tree.total(),side='right'))
  stress.append(dict(case=case,updates=10000,passed=True))
 save('stress.json',dict(passed=True,records=stress))
 for c in range(4):
  state=fixture(c);model=Model(len(state.sites),config(state));arrays=[];events=0
  with (BASE/f'tiny_c{c}.jsonl').open('x') as f:
   for seed in p['tiny_seeds'][c]:
    rng=np.random.default_rng(seed);pop=state.counts();saved=[pop.copy()];traces=[]
    for step in range(100):
     work,ids,times=execute_tick(model,pop,.001,rng,deadline);out=work.output.reshape(pop.shape);available=work.available.reshape(pop.shape);audit(model,pop,ids,times,out,available,.001,deadline);events+=len(ids);traces.extend([dict(step=step,channel=int(i),time=float(t)) for i,t in zip(ids,times)]);pop=out
     if step in (49,99):saved.append(pop.copy())
    arrays.append(saved);f.write(json.dumps(dict(seed=seed,events=traces))+'\n')
  a=np.asarray(arrays);np.savez_compressed(BASE/f'tiny_c{c}.npz',counts=a,seeds=p['tiny_seeds'][c]);dist=[]
  with np.load(ROOT/f'results/MG12H3/mixed_qualification/references/c{c}_reference.npz') as z:support=z['states'];prob=z['probabilities']
  for t in range(2):
   index={tuple(v):i for i,v in enumerate(support)};freq=np.bincount([index[tuple(v)] for v in a[:,t+1].reshape(256,-1)],minlength=len(support))/256;dist.append(float(abs(freq-prob[t]).sum()/2))
  tiny.append(dict(fixture=c,paths=256,accepted_events=events,observed_TV=dist,formal_verdict=None))
 save('tiny_results.json',dict(complete=True,records=tiny,scope='Exploratory endpoint samples, not powered formal qualification'))
 for case,seed in zip(p['cases'],p['actual_seeds']):
  clock=time.perf_counter();end=min(deadline,clock+20);pop,model=load(ROOT/f'results/MG12H3/ring_feasibility/{case}.npz');work,ids,times=execute_tick(model,pop,.001,np.random.default_rng(seed),end);sampling=time.perf_counter()-clock;out=work.output.reshape(pop.shape);available=work.available.reshape(pop.shape);np.savez_compressed(BASE/f'actual_{case}.npz',initial=pop,final=out,available=available,channel_ids=ids,times=times,seed=seed);auditstart=time.perf_counter();audit(model,pop,ids,times,out,available,.001,end);actual.append(dict(case=case,seed=seed,events=len(ids),sampling_seconds=sampling,replay_seconds=time.perf_counter()-auditstart,total_seconds=time.perf_counter()-clock));save('actual_progress.json',dict(records=actual))
 assert hashes()==p['source_hashes'];save('results.json',dict(complete=True,development_checks_passed=True,elapsed_seconds=time.perf_counter()-start,stress=stress,tiny=tiny,actual=actual,formal_qualification=False,source_hashes=hashes()))
if __name__=='__main__':
 try:plan() if '--plan' in sys.argv else run()
 except BaseException:save('failure.json',dict(error=traceback.format_exc()));raise

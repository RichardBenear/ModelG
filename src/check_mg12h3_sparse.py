"""MG135 deterministic equivalence and forced-work benchmarks; no RNG."""
import os
for key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[key]='1'
from pathlib import Path
import json,time,math,numpy as np
from collections import defaultdict
from mg12h3_sparse import Model,Tree,tick
from mg12h3_exchangeable import channels
from prepare_mg12h3_bridge import fixture,config,ROOT
BASE=ROOT/'results/MG12H3/sparse_development'
def signature(rule,reserve,change):return rule,tuple(sorted(reserve)),tuple(sorted(change))
def groups_dense(a,cfg):
 result=defaultdict(float)
 for rule,change,reserve,rate in channels(a,cfg):result[signature(rule,[(int(i),int(v)) for i,v in enumerate(reserve.ravel()) if v],[(int(i),int(v)) for i,v in enumerate(change.ravel()) if v])]+=rate
 return result
def groups_sparse(model,work):
 result=defaultdict(float)
 for c,rate in zip(model.channels,work.tree.values):
  if rate:result[signature(c.rule,c.reserve,c.change)]+=rate
 return result
def compare(model,work,cfg):
 a=groups_dense(work.available.reshape(7,model.n),cfg);b=groups_sparse(model,work);assert a.keys()==b.keys()
 for key in a:np.testing.assert_allclose(a[key],b[key],rtol=1e-11,atol=1e-11)
 np.testing.assert_allclose(work.tree.total(),math.fsum(work.tree.values),rtol=1e-11,atol=1e-11)
 if work.positive:
  for u in (.1234567,.5432109,.9234567):
   idx=work.tree.select(u);cum=np.cumsum(work.tree.values);expected=int(np.searchsorted(cum,u*work.tree.total(),side='right'));assert idx==expected

def check():
 start=time.perf_counter();states=0;updates=0
 for c in range(4):
  state=fixture(c);cfg=config(state);model=Model(len(state.sites),cfg)
  with np.load(ROOT/f'results/MG12H3/bridge/c{c}_reference.npz') as z:support=z['states']
  for flat in support:
   if time.perf_counter()-start>60:raise TimeoutError('Algebra check cap')
   work=model.workspace(flat.reshape(7,model.n));compare(model,work,cfg)
   for step in range(3):
    if not work.positive:break
    before=work.output.copy();idx=work.tree.select(.4123456);channel=work.accept(idx);expected=before.copy()
    for i,value in channel.change:expected[i]+=value
    np.testing.assert_array_equal(work.output,expected);assert work.available.min()>=0 and work.output.min()>=0 and sum(map(int,work.output))==sum(map(int,flat));compare(model,work,cfg);updates+=1
   states+=1
 class Fixed:
  def exponential(self,scale):return .01
  def random(self):return .25
 a=np.zeros((7,1),np.int64);a[1,0]=2;a[2,0]=2;cfg=(2.,[(10.,4.)],[],(.2,.3,.4),(.4,.7),(False,False,False,False,True,False));model=Model(1,cfg);assert tick(model,a,.1,Fixed())['events']==1;assert tick(model,a,.01,Fixed())['events']==0
 a=np.zeros((7,2),np.int64);a[0,0]=1;cfg=(2.,[(10.,4.)]*2,[('ActiveLink',0,1,1.),('ActiveLink',1,0,1.)],(1.,1.,1.),(1.,1.),(False,)*6);model=Model(2,cfg);out=tick(model,a,.1,Fixed());assert out['events']==1 and out['counts'][0,1]==1
 BASE.mkdir(parents=True,exist_ok=True);(BASE/'algebra.json').write_text(json.dumps(dict(passed=True,states=states,local_updates_checked=updates,scripted_ticks=3,elapsed_seconds=time.perf_counter()-start,random_draws=0),indent=2)+'\n')

def benchmark():
 start=time.perf_counter();rows=[]
 for p in sorted((ROOT/'results/MG12H3/ring_feasibility').glob('*.npz')):
  with np.load(p) as z:pop=z['counts'];cfg=(float(z['omega']),list(zip(z['a0'],z['b0'])),json.loads(str(z['links'])),tuple(z['diff']),tuple(z['coupling']),(True,)*6)
  clock=time.perf_counter();model=Model(pop.shape[1],cfg);work=model.workspace(pop);setup=time.perf_counter()-clock;clock=time.perf_counter();audit=0
  for step in range(1000):
   if time.perf_counter()-start>60:raise TimeoutError('Benchmark deadline')
   u=((step+1)*.6180339887498949)%1;idx=work.tree.select(u);work.accept(idx)
  elapsed=time.perf_counter()-clock
  fresh=np.array([c.rate(work.available) for c in model.channels]);np.testing.assert_allclose(work.tree.values,fresh,rtol=1e-11,atol=1e-11);np.testing.assert_allclose(work.tree.total(),math.fsum(fresh),rtol=1e-11,atol=1e-6);assert work.output.min()>=0 and work.output.sum()==pop.sum()
  payload=work.tree.values.nbytes+work.tree.tree.nbytes+work.available.nbytes+work.output.nbytes
  rows.append(dict(case=p.stem,channels=len(model.channels),setup_seconds=setup,forced_updates=1000,update_seconds=elapsed,updates_per_second=1000/elapsed,array_payload_bytes=payload,maximum_dependency_channels=max(map(len,model.dependencies))))
 (BASE/'benchmark.json').write_text(json.dumps(dict(complete=True,elapsed_seconds=time.perf_counter()-start,records=rows,random_draws=0,warning='Forced selections, not a simulated trajectory; no actual event throughput or particle runtime inference. Array payload excludes Python objects.'),indent=2)+'\n')
if __name__=='__main__':
 import sys
 check() if '--check' in sys.argv else benchmark()

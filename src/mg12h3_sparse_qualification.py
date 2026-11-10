"""MG137 sparse writer/stress worker and prospective provenance."""
from pathlib import Path
import os
for key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[key]='1'
import numpy as np,json,time,hashlib,math,sys,traceback
from mg12h3_sparse import Model,tick,Tree
from prepare_mg12h3_bridge import fixture,config,ROOT,seed_for
from mg12h3_efficiency_qualification import check,score
from mg10e2c1_pool import atomic_json
from run_mg12h3_sparse_combined import load
BASE=ROOT/'results/MG12H3/sparse_qualification'
FILES=['src/mg12h3_sparse_qualification.py','src/prepare_mg12h3_sparse_qualification.py','src/run_mg12h3_sparse_qualification.py','spec/MG12H3_sparse_qualification.md','src/mg12h3_sparse.py','src/prepare_mg12h3_bridge.py','src/mg12h3_efficiency_qualification.py','src/prepare_mg12h3_efficiency.py','src/mg12h1.py','src/mg12h2.py','src/mg12h1_reference.py','src/mg10e2c1_pool.py','src/run_mg12h3_sparse_combined.py']
def hashes():
 h={p:hashlib.sha256((ROOT/p).read_text(encoding='utf-8').replace('\r\n','\n').encode()).hexdigest() for p in FILES}
 for p in (ROOT/'results/MG12H3/ring_feasibility').glob('*.npz'):h[p.relative_to(ROOT).as_posix()]=hashlib.sha256(p.read_bytes()).hexdigest()
 for c in range(4):
  p=ROOT/f'results/MG12H3/mixed_qualification/references/c{c}_reference.npz';h[p.relative_to(ROOT).as_posix()]=hashlib.sha256(p.read_bytes()).hexdigest()
 return h

def jobs():
 endpoint=[dict(id=f'endpoint_c{c}_b{b}',kind='endpoint',cell=c,chunk=b,reps=4096) for c in range(4) for b in range(8)]
 cases=sorted(p.stem for p in (ROOT/'results/MG12H3/ring_feasibility').glob('*.npz'))
 return endpoint+[dict(id='stress_'+case,kind='stress',case=case,updates=1000000) for case in cases]

def write_endpoint(job,seeds,directory,deadline):
 directory.mkdir(parents=True,exist_ok=True);clock=time.perf_counter();state=fixture(job['cell']);model=Model(len(state.sites),config(state));arrays=[];events=0
 # Save count prefixes on exceptions; incomplete status cannot qualify.
 try:
  with (directory/(job['id']+'.jsonl')).open('x') as f:
   for seed in seeds:
    pop=state.counts();saved=[pop.copy()];rows=[];rng=np.random.default_rng(seed);time0=0.;step=0
    for target in (.05,.1):
     while time0<target:
      end=min(time0+.001,target);out=tick(model,pop,end-time0,rng,deadline=deadline);n=model.n
      for e in out['trace']:
       change=np.zeros(7*n,np.int64);reserve=np.zeros(7*n,np.int64)
       for i,v in e['change']:change[i]+=v
       for i,v in e['reserve']:reserve[i]+=v
       rows.append(dict(time=e['time']+time0,step=step,rule=e['rule'],increment=change.reshape(7,n).tolist(),reserved=reserve.reshape(7,n).tolist()))
      pop=out['counts'];events+=out['events'];time0=end;step+=1
     saved.append(pop.copy())
    arrays.append(saved);f.write(json.dumps(dict(seed=seed,trace=rows))+'\n')
 finally:
  shape=(0,3,7,len(state.sites));data=np.asarray(arrays,dtype=np.int64) if arrays else np.empty(shape,np.int64);np.savez_compressed(directory/(job['id']+'.npz'),counts=data,seeds=seeds[:len(arrays)],job=json.dumps(job))
 path=directory/(job['id']+'.npz');meta=dict(status='COMPLETE',job=job,reps=len(arrays),events=events,elapsed_seconds=time.perf_counter()-clock,sha256=hashlib.sha256(path.read_bytes()).hexdigest());atomic_json(directory/(job['id']+'.json'),meta);check(job,seeds,directory,deadline);return meta

def stress(job,directory,deadline):
 directory.mkdir(parents=True,exist_ok=True);clock=time.perf_counter();pop,model=load(ROOT/f'results/MG12H3/ring_feasibility/{job["case"]}.npz');work=model.workspace(pop);ids=[];rows=[]
 try:
  for step in range(job['updates']):
   if step%1000==0 and time.perf_counter()>deadline:raise TimeoutError('Stress deadline')
   u=((step+1)*.6180339887498949)%1;index=work.tree.select(u);work.accept(index);ids.append(index)
   if (step+1)%1000==0:
    fresh=np.array([c.rate(work.available) for c in model.channels]);np.testing.assert_allclose(work.tree.values,fresh,rtol=1e-10,atol=1e-6);total=math.fsum(fresh);np.testing.assert_allclose(work.tree.total(),total,rtol=1e-10,atol=1e-6);cum=np.cumsum(fresh)
    for probe in (1e-9,.5,1-1e-9):assert work.tree.select(probe)==int(np.searchsorted(cum,probe*work.tree.total(),side='right'))
    assert work.available.min()>=0 and work.output.min()>=0 and work.output.sum()==pop.sum();rows.append(dict(updates=step+1,total=total,tree_total=work.tree.total(),maximum_leaf_error=float(abs(fresh-work.tree.values).max())))
 finally:
  np.savez_compressed(directory/(job['id']+'.npz'),initial=pop,final=work.output.reshape(pop.shape),available=work.available.reshape(pop.shape),channel_ids=np.asarray(ids,np.int32),job=json.dumps(job));atomic_json(directory/(job['id']+'_checkpoints.json'),dict(records=rows))
 meta=dict(status='COMPLETE',job=job,updates=len(ids),elapsed_seconds=time.perf_counter()-clock);atomic_json(directory/(job['id']+'.json'),meta);return meta

def stress_check(job,directory):
 with np.load(directory/(job['id']+'.npz')) as z:initial=z['initial'];ids=z['channel_ids'];final=z['final'];available=z['available'];assert json.loads(str(z['job']))==job
 _,model=load(ROOT/f'results/MG12H3/ring_feasibility/{job["case"]}.npz');assert len(ids)==job['updates'];counts=np.bincount(ids,minlength=len(model.channels));reserve=np.zeros(initial.size,np.int64);change=reserve.copy()
 for c,n in zip(model.channels,counts):
  for i,v in c.reserve:reserve[i]+=int(n)*v
  for i,v in c.change:change[i]+=int(n)*v
 np.testing.assert_array_equal(initial.ravel()-reserve,available.ravel());np.testing.assert_array_equal(initial.ravel()+change,final.ravel());assert available.min()>=0 and final.min()>=0 and final.sum()==initial.sum();return dict(passed=True,updates=len(ids))

def controls():
 tree=Tree([1e9,.1,0.,3.]);tree.update(0,0.);np.testing.assert_allclose(tree.total(),3.1,rtol=1e-10,atol=1e-6)
 for u in (0.,.01,.5,1-1e-9):assert tree.select(u)==int(np.searchsorted(np.cumsum(tree.values),u*tree.total(),side='right'))
 tree=Tree([1.]);tree.update(0,0.);assert tree.total()==0
 try:tree.select(.5);raise RuntimeError('Exhausted selected')
 except ValueError:pass
 return dict(passed=True,scope='Targeted finite dynamic range and exhausted selection, not arbitrary numerical guarantee')

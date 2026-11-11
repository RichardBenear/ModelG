"""MG139 separate saved-array/ID recomputation, no simulation or random draws."""
from pathlib import Path
import json,hashlib,math,numpy as np
from mg12h3_sparse_qualification import ROOT,BASE,hashes,load

def run():
 r=json.loads((BASE/'registration.json').read_text());e=json.loads((BASE/'evaluation.json').read_text());assert hashes()==r['source_hashes'];groups={c:[] for c in range(4)};stress=[];events=0
 for job in r['jobs']:
  p=BASE/'production'/(job['id']+'.npz');meta=json.loads(p.with_suffix('.json').read_text());assert meta['status']=='COMPLETE' and meta['job']==job
  if job['kind']=='endpoint':
   assert hashlib.sha256(p.read_bytes()).hexdigest()==meta['sha256'];events+=meta['events']
   with np.load(p) as z:a=z['counts'];assert z['seeds'].tolist()==r['seeds'][job['id']] and json.loads(str(z['job']))==job
   assert a.shape[0]==4096 and a.min()>=0 and np.all(a.sum(axis=(2,3))==a.sum(axis=(2,3))[:,0,None]);groups[job['cell']].append(a)
  else:
   with np.load(p) as z:initial=z['initial'];final=z['final'];available=z['available'];ids=z['channel_ids'];assert json.loads(str(z['job']))==job
   donor,model=load(ROOT/f'results/MG12H3/ring_feasibility/{job["case"]}.npz');np.testing.assert_array_equal(initial,donor);assert len(ids)==1000000 and ids.min()>=0 and ids.max()<len(model.channels)
   ordered=np.sort(ids);counts=np.searchsorted(ordered,np.arange(len(model.channels)),side='right')-np.searchsorted(ordered,np.arange(len(model.channels)),side='left');reserve=np.zeros(initial.size,np.int64);change=reserve.copy()
   for c,n in zip(model.channels,counts):
    for i,v in c.reserve:reserve[i]+=int(n)*v
    for i,v in c.change:change[i]+=int(n)*v
   np.testing.assert_array_equal(initial.ravel()-reserve,available.ravel());np.testing.assert_array_equal(initial.ravel()+change,final.ravel());assert available.min()>=0 and final.min()>=0 and final.sum()==initial.sum()
   checkpoints=json.loads((BASE/'production'/(job['id']+'_checkpoints.json')).read_text())['records'];assert [x['updates'] for x in checkpoints]==list(range(1000,1000001,1000));assert all(abs(x['tree_total']-x['total'])<=1e-6+1e-10*abs(x['total']) for x in checkpoints);stress.append(dict(case=job['case'],updates=len(ids),checkpoints=len(checkpoints),maximum_leaf_error=max(x['maximum_leaf_error'] for x in checkpoints),maximum_total_error=max(abs(x['tree_total']-x['total']) for x in checkpoints),maximum_relative_total_error=max(abs(x['tree_total']-x['total'])/x['total'] for x in checkpoints),passed=True))
 rows=[]
 for c,parts in groups.items():
  a=np.concatenate(parts);assert len(a)==32768
  with np.load(ROOT/f'results/MG12H3/mixed_qualification/references/c{c}_reference.npz') as z:states=z['states'];prob=z['probabilities']
  lookup={tuple(v):i for i,v in enumerate(states)}
  for t in range(2):
   ids=np.fromiter((lookup[tuple(v)] for v in a[:,t+1].reshape(32768,-1)),dtype=np.int64);freq=np.bincount(ids,minlength=len(states))/32768;distance=math.fsum(abs(float(x)-float(y)) for x,y in zip(freq,prob[t]))/2;upper=distance+math.sqrt(math.log((2**len(states))*8000)/65536);old=next(x for x in e['comparisons'] if x['cell']==c and x['time']==(.05,.1)[t]);assert abs(upper-old['upper'])<1e-12;rows.append(dict(cell=c,time=(.05,.1)[t],observed=distance,upper=upper,passed=upper<=.1))
 result=dict(passed=all(x['passed'] for x in rows+stress),comparisons=rows,stress=stress,endpoint_paths=131072,endpoint_events=events,stress_updates=8000000,source_hashes_unchanged=True,shared='Shared raw arrays/references and compact channel compilation; separate histogram and sorted-ID aggregation. Same author, not external replication. Checkpoint selection/leaf checks rely on registered worker checks.')
 (BASE/'separate_recomputation.json').write_text(json.dumps(result,indent=2)+'\n');files=[dict(path=p.relative_to(ROOT).as_posix(),bytes=p.stat().st_size,sha256=hashlib.sha256(p.read_bytes()).hexdigest()) for p in sorted((BASE/'production').iterdir()) if p.is_file()];(BASE/'raw_manifest.json').write_text(json.dumps(dict(files=files,total_bytes=sum(x['bytes'] for x in files)),indent=2)+'\n');print(json.dumps({k:v for k,v in result.items() if k not in ('comparisons','stress')}));print('max stress relative drift',max(x['maximum_relative_total_error'] for x in stress))
if __name__=='__main__':run()

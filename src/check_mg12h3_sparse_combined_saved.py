"""MG136 separate saved-output aggregation, no simulation or RNG."""
from run_mg12h3_sparse_combined import BASE,ROOT,load,Model,fixture,config,np,json,hashes
from collections import Counter
import hashlib

def aggregate(model,ids,initial):
    frequencies=np.bincount(ids,minlength=len(model.channels));reserve=np.zeros(initial.size,np.int64);change=np.zeros(initial.size,np.int64)
    for c,n in zip(model.channels,frequencies):
        if n:
            for i,v in c.reserve:reserve[i]+=int(n)*v
            for i,v in c.change:change[i]+=int(n)*v
    # Snapshot reservations are monotone: nonnegative final availability also
    # bounds every reservation prefix; products cannot replenish availability.
    available=initial.ravel()-reserve;out=initial.ravel()+change;assert available.min()>=0 and out.min()>=0 and out.sum()==initial.sum()
    return out.reshape(initial.shape),available.reshape(initial.shape)

def run():
 p=json.loads((BASE/'plan.json').read_text());r=json.loads((BASE/'results.json').read_text());assert p['source_hashes']==hashes();actual=[];tiny=[]
 for case,seed in zip(p['cases'],p['actual_seeds']):
  initial,model=load(ROOT/f'results/MG12H3/ring_feasibility/{case}.npz')
  with np.load(BASE/f'actual_{case}.npz') as z:
   assert int(z['seed'])==seed;np.testing.assert_array_equal(z['initial'],initial);ids=z['channel_ids'];times=z['times'];assert len(ids)==len(times) and np.all(np.diff(times)>0) and np.all(times>0) and np.all(times<.001);out,available=aggregate(model,ids,initial);np.testing.assert_array_equal(out,z['final']);np.testing.assert_array_equal(available,z['available']);actual.append(dict(case=case,events=len(ids),passed=True))
 for c in range(4):
  state=fixture(c);model=Model(len(state.sites),config(state))
  with np.load(BASE/f'tiny_c{c}.npz') as z:counts=z['counts'];assert z['seeds'].tolist()==p['tiny_seeds'][c]
  lines=0
  with (BASE/f'tiny_c{c}.jsonl').open() as f:
   for line in f:
    record=json.loads(line);assert record['seed']==p['tiny_seeds'][c][lines];pop=state.counts();np.testing.assert_array_equal(pop,counts[lines,0]);by_step={}
    for event in record['events']:
     assert 0<=event['step']<100 and 0<event['time']<.001;by_step.setdefault(event['step'],[]).append(event['channel'])
    for step in range(100):
     pop,_=aggregate(model,np.asarray(by_step.get(step,[]),np.int64),pop)
     if step in (49,99):np.testing.assert_array_equal(pop,counts[lines,1 if step==49 else 2])
    lines+=1
  assert lines==256
  with np.load(ROOT/f'results/MG12H3/mixed_qualification/references/c{c}_reference.npz') as z:states=z['states'];prob=z['probabilities']
  distances=[]
  for t in range(2):
   d=Counter(map(tuple,counts[:,t+1].reshape(256,-1)));assert set(d)<=set(map(tuple,states));distance=sum(abs(d.get(tuple(v),0)/256-float(pr)) for v,pr in zip(states,prob[t]))/2;assert abs(distance-r['tiny'][c]['observed_TV'][t])<1e-12;distances.append(distance)
  tiny.append(dict(fixture=c,paths=lines,observed_TV=distances,passed=True))
 files=[dict(path=x.relative_to(ROOT).as_posix(),bytes=x.stat().st_size,sha256=hashlib.sha256(x.read_bytes()).hexdigest()) for x in sorted(BASE.iterdir()) if x.is_file()]
 (BASE/'separate_recomputation.json').write_text(json.dumps(dict(passed=True,actual=actual,tiny=tiny,source_hashes_unchanged=True,shared='Same compact model compilation, saved inputs. Independent aggregate reservation/increment arithmetic and categorical frequency calculation; same author, not external replication.'),indent=2)+'\n');(BASE/'raw_manifest.json').write_text(json.dumps(dict(files=files,total_bytes=sum(x['bytes'] for x in files)),indent=2)+'\n')
if __name__=='__main__':run()

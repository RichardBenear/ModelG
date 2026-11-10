"""MG130 actual formal chunk writer/evaluator; registration separate."""
from prepare_mg12h3_bridge import ROOT,fixture,config,seed_for,atomic_json,np,json,time,hashlib
from prepare_mg12h3_efficiency import audit
from mg12h3_exchangeable import simulate
from collections import Counter
import math
BASE=ROOT/'results/MG12H3/efficiency/qualification'
def jobs():return [dict(id=f'c{c}_b{b}',cell=c,chunk=b,reps=4096) for c in range(4) for b in range(8)]
def write(job,seeds,directory,deadline):
    directory.mkdir(parents=True,exist_ok=True);state=fixture(job['cell']);arrays=[];start=time.perf_counter();events=0
    with (directory/(job['id']+'.jsonl')).open('x') as f:
        for seed in seeds:
            out=simulate(state.counts(),config(state),[0.,.05,.1],seed,deadline=deadline);events+=audit(state,out['counts'],out['trace']);arrays.append(out['counts']);f.write(json.dumps(dict(seed=seed,trace=out['trace']))+'\n')
    path=directory/(job['id']+'.npz');np.savez_compressed(path,counts=np.asarray(arrays),seeds=seeds,job=json.dumps(job));atomic_json(directory/(job['id']+'.json'),dict(status='COMPLETE',job=job,reps=len(arrays),events=events,elapsed_seconds=time.perf_counter()-start,sha256=hashlib.sha256(path.read_bytes()).hexdigest()))
def check(job,seeds,directory,deadline):
    path=directory/(job['id']+'.npz');meta=json.loads((directory/(job['id']+'.json')).read_text());assert meta['status']=='COMPLETE' and meta['job']==job and meta['sha256']==hashlib.sha256(path.read_bytes()).hexdigest()
    with np.load(path) as z:a=z['counts'];assert z['seeds'].tolist()==seeds and json.loads(str(z['job']))==job
    state=fixture(job['cell']);assert a.shape==(job['reps'],3,7,len(state.sites)) and a.dtype==np.int64 and np.all(a>=0);lines=0;events=0
    with (directory/(job['id']+'.jsonl')).open() as f:
        for line in f:
            if time.perf_counter()>deadline:raise TimeoutError('Audit deadline')
            row=json.loads(line);assert row['seed']==seeds[lines];events+=audit(state,a[lines],row['trace'])
            for frame,target in enumerate((0.,.05,.1)):
                reconstructed=state.counts().copy()
                for event in row['trace']:
                    if event['time']<=target:reconstructed+=np.asarray(event['increment'])
                np.testing.assert_array_equal(reconstructed,a[lines,frame])
            lines+=1
    assert lines==len(a) and events==meta['events'];return a

def score(arrays,required=32768):
    rows=[]
    for c in range(4):
        a=np.concatenate(arrays[c]);assert len(a)==required
        with np.load(ROOT/f'results/MG12H3/mixed_qualification/references/c{c}_reference.npz') as z:states=z['states'];prob=z['probabilities']
        for t in range(2):
            flat=a[:,t+1].reshape(len(a),-1);d=Counter(map(tuple,flat));assert set(d)<=set(map(tuple,states));frequency=np.array([d.get(tuple(v),0)/required for v in states]);u,n=np.unique(flat,axis=0,return_counts=True);alternative={tuple(v):k/required for v,k in zip(u,n)};np.testing.assert_array_equal(frequency,[alternative.get(tuple(v),0) for v in states]);tv=float(abs(frequency-prob[t]).sum()/2);epsilon=math.sqrt((len(states)*math.log(2)+math.log(8000))/(2*required));rows.append(dict(cell=c,time=(.05,.1)[t],observed=tv,epsilon=epsilon,upper=tv+epsilon,passed=bool(tv+epsilon<=.1)))
    return dict(qualified=all(x['passed'] for x in rows),comparisons=rows)

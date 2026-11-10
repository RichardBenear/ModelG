"""MG132 saved-array recomputation, no executor/evaluator calls or RNG."""
from pathlib import Path
import json,math,hashlib
import numpy as np
ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/'results/MG12H3/efficiency/qualification'
def run():
    registration=json.loads((BASE/'registration.json').read_text());official=json.loads((BASE/'evaluation.json').read_text());groups={c:[] for c in range(4)};manifest=[];events=0
    for job in registration['jobs']:
        p=BASE/'production'/(job['id']+'.npz');meta=json.loads(p.with_suffix('.json').read_text());assert meta['status']=='COMPLETE' and meta['sha256']==hashlib.sha256(p.read_bytes()).hexdigest();events+=meta['events']
        with np.load(p) as z:
            a=z['counts'];assert z['seeds'].tolist()==registration['seeds'][job['id']] and json.loads(str(z['job']))==job
        assert len(a)==job['reps'] and np.all(a>=0);assert np.all(a.sum(axis=(2,3))==a.sum(axis=(2,3))[:,0,None]);groups[job['cell']].append(a)
    rows=[]
    for c,parts in groups.items():
        a=np.concatenate(parts);assert len(a)==32768
        with np.load(ROOT/f'results/MG12H3/mixed_qualification/references/c{c}_reference.npz') as z:states=z['states'];prob=z['probabilities']
        index={tuple(v):i for i,v in enumerate(states)}
        for t in range(2):
            ids=np.fromiter((index[tuple(v)] for v in a[:,t+1].reshape(len(a),-1)),dtype=np.int64);freq=np.bincount(ids,minlength=len(states))/len(a);tv=math.fsum(abs(float(x)-float(y)) for x,y in zip(freq,prob[t]))/2;upper=tv+math.sqrt(math.log((2**len(states))*8000)/(2*len(a)));prior=next(x for x in official['comparisons'] if x['cell']==c and x['time']==(.05,.1)[t]);assert abs(upper-prior['upper'])<1e-12;rows.append(dict(cell=c,time=(.05,.1)[t],observed=tv,upper=upper,passed=upper<=.1))
    for p in sorted((BASE/'production').glob('*')):
        if p.is_file():manifest.append(dict(path=p.relative_to(ROOT).as_posix(),bytes=p.stat().st_size,sha256=hashlib.sha256(p.read_bytes()).hexdigest()))
    result=dict(passed=all(x['passed'] for x in rows),paths=131072,chunks=32,accepted_events_from_metadata=events,comparisons=rows,shared_dependencies='Saved reference probabilities and raw arrays. Separate histogram/bound calculation, same author; not external replication.')
    (BASE/'separate_recomputation.json').write_text(json.dumps(result,indent=2)+'\n');(BASE/'raw_manifest.json').write_text(json.dumps(dict(files=manifest,total_bytes=sum(x['bytes'] for x in manifest)),indent=2)+'\n');print(json.dumps({k:v for k,v in result.items() if k!='comparisons'}))
if __name__=='__main__':run()

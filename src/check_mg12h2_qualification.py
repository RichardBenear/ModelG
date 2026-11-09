"""MG111 saved-data fixed-cell joint-law qualification audit."""
import time,math,json,hashlib
from pathlib import Path
from dataclasses import asdict
import numpy as np
from prepare_mg12h2 import ROOT,read,fixture
from mg12h1 import SPECIES
from mg10d1_pool import atomic_json
BASE=ROOT/'results/MG12H2/qualification'
def check(deadline=math.inf):
    from run_mg12h2_qualification import hashes
    reg=read(BASE/'registration.json');assert reg['source_hashes']==hashes();R=reg['R'];times=np.array(reg['times']);comparisons=[];manifest=[];paths=0;accepted=0
    for job in reg['jobs']:
        if time.perf_counter()>deadline:raise TimeoutError('Audit cap')
        state=fixture(job['cell'])[0];assert reg['fixtures'][job['cell']]==json.loads(json.dumps(asdict(state)));path=BASE/(job['id']+'.npz');meta=read(path.with_suffix('.json'));assert meta['status']=='COMPLETE' and meta['replicas']==R and meta['sha256']==hashlib.sha256(path.read_bytes()).hexdigest()
        with np.load(path) as raw:data=raw['counts'];assert raw['seeds'].tolist()==job['seeds'];np.testing.assert_array_equal(raw['times'],times)
        assert data.dtype==np.int64 and data.shape==(R,3,7,len(state.sites)) and data.min()>=0 and np.all(data.sum(axis=(2,3))==len(state.occupations()));tracepath=BASE/(job['id']+'_traces.jsonl');lines=0;jobaccepted=0;cellset={(rel.kind,frozenset(rel.roles)) for rel in state.relations if rel.kind!='Occupies'}
        with tracepath.open(encoding='utf-8') as f:
            for line in f:
                if time.perf_counter()>deadline:raise TimeoutError('Trace audit cap')
                record=json.loads(line);r=record['replicate'];assert r==lines and record['seed']==job['seeds'][r];current=dict(state.occupations());step=-1;used=set();snapshot=dict(current);last=0.;frame=1;np.testing.assert_array_equal(data[r,0],state.counts())
                def compare(index):
                    expected=np.zeros((7,len(state.sites)),np.int64)
                    for site,sp in current.values():expected[SPECIES.index(sp),site]+=1
                    np.testing.assert_array_equal(expected,data[r,index])
                for row in record['trace']:
                    assert last<=row['time']<=.1;last=row['time']
                    while frame<3 and times[frame]<row['time']:compare(frame);frame+=1
                    if row['step']!=step:assert row['step']>step;step=row['step'];used=set();snapshot=dict(current)
                    tokens=set(row['participants']);before={t:(s,sp) for t,s,sp in row['before']};after={t:(s,sp) for t,s,sp in row['after']};assert len(tokens)==1 and set(before)==set(after)==tokens and all(snapshot.get(t)==v for t,v in before.items());token=next(iter(tokens));source,sp=before[token];dest,newsp=after[token];kind,roles=row['cell'];assert sp==newsp and source!=dest and row['rule']=='H_'+sp and source in roles and dest in roles;assert (kind,frozenset(roles)) in cellset and kind==('ActiveCell' if sp in SPECIES[:3] else 'ACell' if sp=='A' else 'BCell')
                    if row['accepted']:assert not tokens & used;used.update(tokens);current.update(after);jobaccepted+=1
                    else:assert job['method']==1 and tokens & used
                while frame<3:compare(frame);frame+=1
                lines+=1;paths+=1
        assert lines==R and jobaccepted==meta['accepted_events_audited'];accepted+=jobaccepted
        with np.load(BASE/f'c{job["cell"]}_reference.npz') as ref:states=ref['states'];prob=ref['probabilities'];np.testing.assert_array_equal(ref['times'],times[1:])
        lookup={tuple(pop.ravel()):i for i,pop in enumerate(states)};S=len(states)
        for frame in (1,2):
            indices=np.array([lookup[tuple(pop.ravel())] for pop in data[:,frame]]);frequency=np.bincount(indices,minlength=S)/R;observed=float(.5*np.sum(abs(frequency-prob[frame-1])));ordered=np.sort(indices);alternate=(np.searchsorted(ordered,np.arange(S),side='right')-np.searchsorted(ordered,np.arange(S),side='left'))/R;observed2=sum(abs(float(a)-float(b)) for a,b in zip(alternate,prob[frame-1]))/2;assert abs(observed-observed2)<1e-12;epsilon=math.sqrt((S*math.log(2)+math.log(24/.001))/(2*R));eps2=float(np.sqrt(np.log((2.**S)*24/.001)/(2*R)));assert abs(epsilon-eps2)<1e-12;comparisons.append(dict(job_id=job['id'],cell=job['cell'],method=job['method'],time=float(times[frame]),states=S,observed_TV=observed,epsilon=epsilon,upper=observed+epsilon,supported=observed+epsilon<=.10,margin=.10))
        for p in (path,tracepath):manifest.append(dict(file=p.relative_to(ROOT).as_posix(),bytes=p.stat().st_size,sha256=hashlib.sha256(p.read_bytes()).hexdigest()))
    assert len(comparisons)==24 and paths==98304;passed=all(r['supported'] for r in comparisons);result=dict(complete=True,audit_passed=True,qualified=passed,verdict='QUALIFIED WITHIN REGISTERED FIXED-CELL TRANSPORT SCOPE' if passed else 'NOT QUALIFIED UNDER REGISTERED ENDPOINT CRITERIA',paths=paths,accepted_events_replayed=accepted,comparisons=comparisons,source_hashes_unchanged=True,confidence='nominal99.9% simultaneous joint-count endpoint TV upper bounds conditional on sampling/PRNG and computed references',shared_dependencies='same fixtures/raw data/source-hash helper; independent literal cell-role/observation replay and alternate categorical frequency/bound arithmetic; same author not external replication');atomic_json(BASE/'evaluation.json',result);atomic_json(BASE/'raw_manifest.json',dict(records=manifest,files=len(manifest),bytes=sum(r['bytes'] for r in manifest)));return result
if __name__=='__main__':check()

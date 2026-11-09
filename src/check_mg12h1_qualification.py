"""MG107 independent saved relation replay and categorical-law bounds."""
import json,time,math,hashlib
from pathlib import Path
import numpy as np
from prepare_mg12h1 import ROOT,fixture,read
from mg12h1 import SPECIES
from mg10d1_pool import atomic_json
BASE=ROOT/'results/MG12H1/qualification'
def check(deadline=math.inf):
    from run_mg12h1_qualification import hashes
    reg=read(BASE/'registration.json');assert hashes()==reg['source_hashes'];R=reg['R'];times=np.asarray(reg['times']);comparisons=[];diagnostics=[];audited_paths=0;events=0;manifest=[]
    for job in reg['jobs']:
        if time.perf_counter()>deadline:raise TimeoutError('Audit deadline')
        p=BASE/(job['id']+'.npz');meta=read(p.with_suffix('.json'));assert meta['status']=='COMPLETE' and meta['replicas']==R and meta['sha256']==hashlib.sha256(p.read_bytes()).hexdigest();state=fixture(job['cell']);assert read(BASE/'registration.json')['fixtures'][job['cell']]==json.loads(json.dumps(__import__('dataclasses').asdict(state)))
        with np.load(p) as z:counts=z['counts'];assert z['seeds'].tolist()==job['seeds'];np.testing.assert_array_equal(z['times'],times)
        assert counts.shape==(R,3,7,len(state.sites)) and counts.dtype==np.int64 and counts.min()>=0;assert np.all(counts.sum(axis=(2,3))==len(state.occupations()));site_index={s:i for i,s in enumerate(state.sites)};linkset={(r.kind,*r.roles) for r in state.relations if r.kind!='Occupies'}
        tracepath=BASE/(job['id']+'_traces.jsonl');lines=0
        with tracepath.open(encoding='utf-8') as f:
            for line in f:
                if time.perf_counter()>deadline:raise TimeoutError('Trace audit deadline')
                record=json.loads(line);r=record['replicate'];assert r==lines and record['seed']==job['seeds'][r];occ=dict(state.occupations());step=-1;used=set();owners={};parents={};history=[];last_time=0.;snapshot=dict(occ);frame=1
                np.testing.assert_array_equal(counts[r,0],state.counts())
                def compare(index):
                    expected=np.zeros((7,len(state.sites)),np.int64)
                    for site,sp in occ.values():expected[SPECIES.index(sp),site_index[site]]+=1
                    np.testing.assert_array_equal(expected,counts[r,index])
                for index,row in enumerate(record['trace']):
                    assert last_time<=row['time']<=times[-1];last_time=row['time']
                    # Tick boundaries are clipped to observations, so no commit straddles a saved frame.
                    while frame<3 and times[frame]<row['time']:compare(frame);frame+=1
                    if row['step']!=step:assert row['step']>step;step=row['step'];used=set();owners={};snapshot=dict(occ)
                    tokens=set(row['participants']);before={t:(s,sp) for t,s,sp in row['before']};after={t:(s,sp) for t,s,sp in row['after']};assert tokens==set(before) and all(snapshot.get(t)==v for t,v in before.items());assert row['parents']==sorted({parents[t] for t in tokens if t in parents});assert row['blocked_by']==sorted({owners[t] for t in tokens if t in owners})
                    if not row['accepted']:assert job['method']!=0 and tokens & used;history.append(row);continue
                    assert not tokens & used;rule=row['rule']
                    if rule.startswith('H_'):
                        assert len(tokens)==1 and len(after)==1;t=next(iter(tokens));s,sp=before[t];dest,newsp=after[t];kind='ActiveLink' if sp in SPECIES[:3] else 'ALink' if sp=='A' else 'BLink';assert sp==newsp and rule=='H_'+sp and (kind,s,dest) in linkset;expected=after
                    elif rule=='R3':
                        assert sorted(sp for _,sp in before.values())==['B','X'] and len({s for s,_ in before.values()})==1;expected={t:(s,'Z' if sp=='B' else 'Y') for t,(s,sp) in before.items()}
                    elif rule=='R4':
                        assert sorted(sp for _,sp in before.values())==['X','X','Y'] and len({s for s,_ in before.values()})==1;expected={t:(s,'X') for t,(s,sp) in before.items() if sp=='Y'}
                    else:
                        source,dest={'R1':('A','G'),'R2':('G','X'),'R2r':('X','G'),'R5':('X','W')}[rule];assert len(tokens)==1;t=next(iter(tokens));s,sp=before[t];assert sp==source;expected={t:(s,dest)}
                    assert after==expected;occ.update(after);used.update(tokens)
                    for token in tokens:owners[token]=index;parents[token]=index
                    history.append(row);events+=1
                while frame<3:compare(frame);frame+=1
                lines+=1;audited_paths+=1
        assert lines==R
        with np.load(BASE/f'c{job["cell"]}_reference.npz') as ref:states=ref['states'];prob=ref['probabilities'];np.testing.assert_array_equal(ref['times'],times[1:])
        lookup={tuple(x):i for i,x in enumerate(states)};S=len(states)
        for frame in (1,2):
            idx=np.array([lookup[tuple(pop.ravel())] for pop in counts[:,frame]]);frequency=np.bincount(idx,minlength=S)/R;observed=float(.5*np.sum(abs(frequency-prob[frame-1])));sorted_idx=np.sort(idx);frequency2=(np.searchsorted(sorted_idx,np.arange(S),side='right')-np.searchsorted(sorted_idx,np.arange(S),side='left'))/R;tv2=sum(abs(float(a)-float(b)) for a,b in zip(frequency2,prob[frame-1]))/2;assert abs(tv2-observed)<1e-12
            epsilon=math.sqrt((S*math.log(2)+math.log(16/.001))/(2*R));alternative=float(np.sqrt(np.log((2.**S)*16/.001)/(2*R)));assert abs(epsilon-alternative)<1e-12;row=dict(job_id=job['id'],cell=job['cell'],method=job['method'],time=float(times[frame]),states=S,observed_TV=observed,epsilon=epsilon,upper=observed+epsilon,margin=.10,supported=observed+epsilon<=.10)
            (comparisons if job['method'] in (0,3) else diagnostics).append(row)
        for file in (p,tracepath):manifest.append(dict(file=file.relative_to(ROOT).as_posix(),bytes=file.stat().st_size,sha256=hashlib.sha256(file.read_bytes()).hexdigest()))
    assert audited_paths==16*R and len(comparisons)==16;passed=all(row['supported'] for row in comparisons);result=dict(complete=True,audit_passed=True,qualified=passed,verdict='QUALIFIED WITHIN REGISTERED SMALL-SYSTEM ENDPOINT SCOPE' if passed else 'NOT QUALIFIED UNDER REGISTERED ENDPOINT CRITERIA',paths=audited_paths,accepted_events_replayed=events,formal=comparisons,descriptive=diagnostics,source_hashes_unchanged=True,confidence='nominal99.9% simultaneous joint-count endpoint TV upper bounds conditional on sampling/PRNG and computed-reference fidelity; not full temporal law',shared_dependencies='same raw paths/fixtures and source hash helper; literal role replay and separately recomputed categorical frequencies/bounds; no external replication');atomic_json(BASE/'evaluation.json',result);atomic_json(BASE/'raw_manifest.json',dict(records=manifest,files=len(manifest),bytes=sum(r['bytes'] for r in manifest)))
    return result
if __name__=='__main__':check()

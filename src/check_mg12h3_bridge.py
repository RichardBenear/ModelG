"""MG114 literal mixed-rule replay and saved endpoint recomputation."""
import math
from pathlib import Path
import numpy as np
from mg12h1 import SPECIES
from mg10d1_pool import atomic_json
TIMES=np.array([0.,.05,.1])
def audit_trace(initial,counts,trace,method):
    occ=dict(initial.occupations());original=set(occ);step=-1;used=set();snapshot=dict(occ);last=0.;frame=0;accepted=0;cells={(r.kind,frozenset(r.roles)) for r in initial.relations if r.kind!='Occupies'}
    def check(index):
        predicted=np.zeros_like(counts[index])
        for site,sp in occ.values():predicted[SPECIES.index(sp),site]+=1
        np.testing.assert_array_equal(predicted,counts[index]);assert predicted.sum()==len(original)
    for row in trace:
        assert last<=row['time']<=.1;last=row['time']
        while frame<3 and TIMES[frame]<row['time']:check(frame);frame+=1
        if row['step']!=step:assert row['step']>step;step=row['step'];used=set();snapshot=dict(occ)
        participants=set(row['participants']);before={t:(s,sp) for t,s,sp in row['before']};after={t:(s,sp) for t,s,sp in row['after']};assert set(before)==participants and all(snapshot.get(t)==value for t,value in before.items())
        if not row['accepted']:assert method!=0 and participants & used;continue
        assert not participants & used;rule=row['rule']
        if rule.startswith('H_'):
            assert len(participants)==1;token=next(iter(participants));site,sp=before[token];dest,newsp=after[token];kind,roles=row['cell'];assert sp==newsp and site!=dest and site in roles and dest in roles and (kind,frozenset(roles)) in cells;assert kind==('ActiveCell' if sp in SPECIES[:3] else 'ACell' if sp=='A' else 'BCell');assert rule=='H_'+sp;expected=after
        elif rule=='R3':
            assert row['cell'] is None and sorted(sp for _,sp in before.values())==['B','X'] and len({s for s,_ in before.values()})==1;expected={t:(s,'Z' if sp=='B' else 'Y') for t,(s,sp) in before.items()}
        elif rule=='R4':
            assert row['cell'] is None and sorted(sp for _,sp in before.values())==['X','X','Y'] and len({s for s,_ in before.values()})==1;expected={t:(s,'X') for t,(s,sp) in before.items() if sp=='Y'}
        else:
            source,dest={'R1':('A','G'),'R2':('G','X'),'R2r':('X','G'),'R5':('X','W')}[rule];assert row['cell'] is None and len(participants)==1;token=next(iter(participants));site,sp=before[token];assert sp==source;expected={token:(site,dest)}
        assert after==expected;occ.update(after);used.update(participants);assert set(occ)==original;accepted+=1
    while frame<3:check(frame);frame+=1
    return accepted
def main():
    from prepare_mg12h3_bridge import ROOT,BASE,read,hashes,fixture,R
    plan=read(BASE/'plan.json');q=read(BASE/'results.json');assert plan['source_hashes']==hashes()==q['source_hashes'];events=0;paths=0;rows=[]
    for job in plan['jobs']:
        initial=fixture(job['cell']);saved=read(BASE/(job['id']+'_traces.json'));assert saved['job']==job and len(saved['traces'])==R
        with np.load(BASE/(job['id']+'.npz')) as z:counts=z['counts'];assert z['seeds'].tolist()==job['seeds'];np.testing.assert_array_equal(z['times'],TIMES)
        for data,trace in zip(counts,saved['traces']):events+=audit_trace(initial,data,trace,job['method']);paths+=1
        with np.load(BASE/f'c{job["cell"]}_reference.npz') as ref:states=ref['states'];prob=ref['probabilities']
        lookup={tuple(v):i for i,v in enumerate(states)};indices=np.sort([lookup[tuple(pop[-1].ravel())] for pop in counts]);observed=np.searchsorted(indices,np.arange(len(states)),side='right')/R;distance=float(max(abs(observed-np.cumsum(prob))));expected=next(row['CDF_distance'] for row in q['jobs'] if row['id']==job['id']);assert abs(distance-expected)<1e-12
        if job['method']==0:assert distance<=math.sqrt(math.log(800)/(2*R))
        rows.append(dict(id=job['id'],CDF_distance=distance))
    assert paths==1536 and events==q['accepted_events_audited'];atomic_json(ROOT/'results/MG12H3/saved_evidence_check.json',dict(passed=True,paths=paths,accepted_events_replayed=events,endpoint_checks=rows,source_hashes_unchanged=True,shared='same raw data/fixtures and literal-role audit used during sampling; saved reload and sorted CDF arithmetic separate from primary dynamics and count oracle, same author not external replication'))
if __name__=='__main__':main()

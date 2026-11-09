"""MG110 separately reread saved cell transport histories and CDF diagnostics."""
import json,math,hashlib
import numpy as np
from prepare_mg12h2 import ROOT,BASE,read,hashes,fixture,TIMES,R
from mg12h1 import SPECIES
from mg10d1_pool import atomic_json
def main():
    p=read(BASE/'plan.json');q=read(BASE/'results.json');assert p['source_hashes']==hashes()==q['source_hashes'];paths=0;events=0;distances=[]
    for job in p['jobs']:
        initial=fixture(job['cell'])[0];ref=np.load(BASE/f'c{job["cell"]}_reference.npz');lookup={tuple(v.ravel()):i for i,v in enumerate(ref['states'])};records=read(BASE/(job['id']+'_traces.json'));assert records['job']==job and len(records['traces'])==R
        with np.load(BASE/(job['id']+'.npz')) as z:counts=z['counts'];assert z['seeds'].tolist()==job['seeds'];np.testing.assert_array_equal(z['times'],TIMES)
        for replica,trace in enumerate(records['traces']):
            current=dict(initial.occupations());step=-1;used=set();before_tick={};last=0.
            for row in trace:
                assert last<=row['time']<=.1;last=row['time']
                if row['step']!=step:step=row['step'];used=set();before_tick=dict(current)
                parts=set(row['participants']);old={t:(site,species) for t,site,species in row['before']};new={t:(site,species) for t,site,species in row['after']};assert set(old)==set(new)==parts and len(parts)==1 and all(before_tick[t]==v for t,v in old.items());token=next(iter(parts));source,sp=old[token];dest,newsp=new[token];kind,roles=row['cell'];assert sp==newsp and source!=dest and source in roles and dest in roles and row['rule']=='H_'+sp
                assert kind==('ActiveCell' if sp in SPECIES[:3] else 'ACell' if sp=='A' else 'BCell');assert any(rel.kind==kind and set(rel.roles)==set(roles) for rel in initial.relations)
                if row['accepted']:assert not used & parts;used.update(parts);current.update(new);events+=1
                else:assert job['method']==1 and used & parts
            for frame,target in enumerate(TIMES):
                state=dict(initial.occupations())
                for row in trace:
                    if row['accepted'] and row['time']<=target:
                        for t,site,sp in row['after']:state[t]=(site,sp)
                predicted=np.zeros_like(counts[replica,frame])
                for site,sp in state.values():predicted[SPECIES.index(sp),site]+=1
                np.testing.assert_array_equal(predicted,counts[replica,frame]);assert predicted.sum()==len(initial.occupations())
            paths+=1
        indices=np.sort([lookup[tuple(v[-1].ravel())] for v in counts]);empirical=np.searchsorted(indices,np.arange(len(lookup)),side='right')/R;distance=float(max(abs(empirical-np.cumsum(ref['probabilities'][1]))));expected=next(x['CDF_distance'] for x in q['jobs'] if x['id']==job['id']);assert abs(distance-expected)<1e-12
        if job['method']==0:assert distance<=math.sqrt(math.log(1200)/(512))
        distances.append(dict(id=job['id'],CDF_distance=distance))
    assert paths==3072 and events==q['accepted_events_audited'];atomic_json(ROOT/'results/MG12H2/saved_evidence_check.json',dict(passed=True,paths=paths,accepted_events_replayed=events,distances=distances,source_hashes_unchanged=True,shared='same saved data, fixture/hash helpers and reference arrays; separate literal cell-role/observation replay and sorted CDF arithmetic, same author not external replication'))
if __name__=='__main__':main()

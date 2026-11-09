"""MG106 reread saved development evidence; no dynamics/RNG calls."""
from pathlib import Path
import json,hashlib,math
import numpy as np
from prepare_mg12h1 import ROOT,BASE,hashes,fixture,R,TIMES
from mg12h1 import SPECIES
from mg10d1_pool import atomic_json
def read(p):return json.loads(Path(p).read_text(encoding='utf-8'))
def main():
    plan=read(BASE/'plan.json');summary=read(BASE/'results.json');assert plan['source_hashes']==hashes()==summary['source_hashes'];epsilon=math.sqrt(math.log(800)/(512));checks=[];paths=0;events=0
    for result in summary['jobs']:
        j=result['job'];initial=fixture(j['cell']);refs=np.load(BASE/f'c{j["cell"]}_reference.npz');states=refs['states'];prob=refs['probabilities'];lookup={tuple(v):i for i,v in enumerate(states)}
        with np.load(BASE/(j['id']+'.npz')) as raw:counts=raw['counts'];indices=raw['state_index'];assert raw['seeds'].tolist()==j['seeds'];np.testing.assert_array_equal(raw['times'],TIMES)
        traces=read(BASE/(j['id']+'_traces.json'));assert traces['job']==j and len(traces['traces'])==R
        for r,trace in enumerate(traces['traces']):
            occ=dict(initial.occupations());tokens=set(occ);last_step=-1;used=set();pending_times=[];snapshots={};last_time=0.
            for row in trace:
                assert row['time']>=last_time;last_time=row['time']
                if row['step']!=last_step:last_step=row['step'];used=set();pre=dict(occ)
                parts=set(row['participants']);before={t:(site,sp) for t,site,sp in row['before']};assert parts==set(before) and all(pre[t]==v for t,v in before.items())
                for parent in row['parents']:assert parent<len(pending_times) and pending_times[parent]['accepted']
                if row['accepted']:
                    assert not used & parts;used.update(parts)
                    for t,site,sp in row['after']:assert t in parts;occ[t]=(site,sp)
                    events+=1
                else:assert used & parts
                assert set(occ)==tokens;pending_times.append(row)
            # Independently reconstruct counts at observation times from literal saved replacements.
            for frame,target in enumerate(TIMES):
                current=dict(initial.occupations())
                for row in trace:
                    if row['accepted'] and row['time']<=target:
                        for token,site,sp in row['after']:current[token]=(site,sp)
                expected=np.zeros_like(counts[r,frame]);sites={site:i for i,site in enumerate(initial.sites)}
                for site,sp in current.values():expected[SPECIES.index(sp),sites[site]]+=1
                np.testing.assert_array_equal(expected,counts[r,frame])
            assert indices[r]==lookup[tuple(counts[r,-1].ravel())];paths+=1
        # Sorted/searchsorted calculation instead of bincount, same frozen state encoding.
        empirical=np.searchsorted(np.sort(indices),np.arange(len(prob)),side='right')/R;d=float(max(abs(empirical-np.cumsum(prob))));assert abs(d-result['joint_state_index_CDF_distance'])<1e-12
        if j['method']==0:assert d<=epsilon
        checks.append(dict(id=j['id'],paths=R,CDF_distance=d))
    files=[]
    for p in sorted(BASE.iterdir()):
        if p.is_file():files.append(dict(file=p.relative_to(ROOT).as_posix(),sha256=hashlib.sha256(p.read_bytes()).hexdigest(),bytes=p.stat().st_size))
    output=ROOT/'results/MG12H1/saved_evidence_check.json';atomic_json(output,dict(passed=True,paths=paths,accepted_events_replayed=events,jobs=checks,source_hashes_unchanged=True,files=files,shared_dependencies='same saved observations, fixtures/model constants and primary source-hash helper; independent literal replacement replay and sorted CDF arithmetic; same author, no external replication'))
if __name__=='__main__':main()

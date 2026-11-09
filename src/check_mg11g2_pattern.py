"""MG088 separate full-record incidence, dense projection and contrast audit."""
import os
for key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[key]='1'
from pathlib import Path
import json,hashlib,math
import numpy as np
from mg11g2_pattern import factory,jobs,reference,OBS,S0
from mg11g2_pattern_oracle import incidence_check,dense_projection,late
from mg10d1_pool import atomic_json
ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/'results/MG11G2/pattern_response'
def read(p):return json.loads(Path(p).read_text(encoding='utf-8'))
def equal(a,b):np.testing.assert_allclose(a,b,rtol=1e-8,atol=1e-9)
def main(base=BASE,synthetic=False):
    base=Path(base);r=read(base/'registration.json');e=read(base/'evaluation.json');checks=[];complete={};scores={};fields={};matched={}
    assert r['jobs']==jobs()
    if not synthetic:
        from prepare_mg11g2_pattern import hashes
        assert hashes()==r['source_hashes'] and hashlib.sha256((ROOT/'results/MG003/ref_profiles.npz').read_bytes()).hexdigest()==r['reference_sha256']
    for job in r['jobs']:
        path=base/'production'/(job['id']+'.npz');case,g,k,mode=factory(job['config'],job['start'])
        if not path.exists():complete[job['id']]=False;continue
        with np.load(path) as raw:
            for name,value in dict(initial=case['initial'],a0=case['a0'],b0=case['b0'],omega=case['omega'],diff=case['diff'],src=g['src'],dst=g['dst'],weight=g['weight'],kappa=k,requested_times=OBS).items():np.testing.assert_array_equal(raw[name],value)
            assert str(raw['mode'])==mode and json.loads(str(raw['job_json']))==job
            t=raw['times'];state=raw['state'];assert len(t)>=1 and np.array_equal(t,OBS[:len(t)]) and state.shape==(len(t),22*g['n']+2*len(g['src']));np.testing.assert_array_equal(state[0,:7*g['n']],case['initial'].ravel());assert not np.any(state[0,7*g['n']:]);pop=state[:,:7*g['n']].reshape(len(t),7,g['n'])
        v=incidence_check(state,case,g);assert v['passed'];v.update(file=path.relative_to(ROOT).as_posix(),sha256=hashlib.sha256(path.read_bytes()).hexdigest());checks.append(v);row=e['records'][job['id']];assert row['sha256']==v['sha256'] and row['checkpoints']==len(t)
        meta=read(path.with_suffix('.json')) if path.with_suffix('.json').exists() else {};full=np.array_equal(t,OBS) and meta.get('status')=='COMPLETE';complete[job['id']]=full;assert row['complete']==full
        if full:assert meta['sha256']==v['sha256'] and meta['job']==job
        if mode=='FINITE':assert not np.any(state[:,-2*g['n']:])
        if job['start']==0:np.testing.assert_array_equal(pop,np.broadcast_to(pop[:,:,:1],pop.shape))
        score=dense_projection(pop,reference()[0]);equal(score,row['metrics']['projection']);equal(pop[:,:3].mean(axis=-1),row['metrics']['field_mean'])
        stock=pop[:,3:5]/np.array([case['a0'],case['b0']])[None];equal(stock.mean(axis=-1),row['stock_mean']);equal(stock[:,:,(250+np.arange(-25,26))%500].mean(axis=-1),row['stock_local']);equal(stock[:,:,np.arange(-25,26)%500].mean(axis=-1),row['stock_remote']);rx=state[:,7*g['n']:14*g['n']].reshape(len(t),7,g['n']);equal(rx[:,[0,3]].sum(axis=-1),row['total_A_B_uptake'])
        if full:
            scores[job['id']]=late(score,t);equal(scores[job['id']],row['metrics']['late_projection']);fields[job['id']]=pop[:,:3]
            # Classifier shared/frozen; compare raw derived values, not a second
            # physical particle-identity implementation.
            from mg10p1 import continuity
            track=continuity(pop[:,1],t,x0=0.);assert track==row['metrics']['reference_continuity'];matched[job['id']]=track['persistent_reference_signal']
        else:assert 'late_projection' not in row['metrics']
    full=all(complete.values());assert full==e['complete'] and sum(complete.values())==e['jobs_complete'];assert set(e['records'])=={j['id'] for j in r['jobs'] if (base/'production'/(j['id']+'.npz')).exists()}
    refinement=[]
    for c in range(7):
        for s in (0,1):
            a=f'c{c}_s{s}_v0';b=f'c{c}_s{s}_v1'
            if complete[a] and complete[b]:refinement.append(dict(condition=c,start=s,active_error=float(np.max(abs(fields[a]-fields[b])/S0[None,:,None]))))
    assert len(refinement)==len(e['refinement'])
    for a,b in zip(refinement,e['refinement']):assert a['condition']==b['condition'] and a['start']==b['start'];equal(a['active_error'],b['active_error'])
    primary=[];secondary=[]
    for version in (0,1):
        for target,(local,cross) in ((primary,(2,3)),(secondary,(5,6))):
            ids=[f'c{c}_s{s}_v{version}' for c in (local,cross) for s in (0,1)]
            if all(complete[key] for key in ids):target.append((scores[ids[3]]-scores[ids[2]]-scores[ids[1]]+scores[ids[0]]).tolist())
            else:target.append(None)
    for actual,expected in zip(primary,e['primary_GXY_by_setting']):
        if actual is None:assert expected is None
        else:equal(actual,expected)
    for actual,expected in zip(secondary,e['C100_secondary_GXY_by_setting']):
        if actual is None:assert expected is None
        else:equal(actual,expected)
    numerical=full and len(refinement)==14 and max(v['active_error'] for v in refinement)<=1e-4;change=None if any(v is None for v in primary) else abs(primary[1][1]-primary[0][1]);threshold=None if change is None else max(.01,5*change);resolved=bool(numerical and primary[0][1]*primary[1][1]>0 and abs(primary[1][1])>threshold);maintained=all(matched.get(f'c0_s1_v{v}',False) for v in (0,1))
    assert numerical==e['numerical_qualified'] and resolved==e['primary_resolved'] and maintained==e['maintained_reference_control'] and bool(resolved and maintained)==e['pattern_response_supported']
    if threshold is None:assert e['primary_resolution_threshold'] is None
    else:equal(threshold,e['primary_resolution_threshold']);equal(change,e['primary_refinement_difference'])
    expected='INCOMPLETE' if not full else 'NUMERICAL RESOLUTION FAILED' if not numerical else 'RESOLVED CONNECTIVITY CONTRAST' if resolved else 'UNRESOLVED AT REGISTERED RESOLUTION';assert e['verdict']==expected
    charge=None
    if not synthetic:
        ledger=read(base/'production/ledger.json');charge=ledger['charged_worker_seconds'];assert ledger['ceiling']==r['ceiling_worker_seconds'] and ledger['prior_charge']==0 and charge<=r['ceiling_worker_seconds'];equal(charge,math.fsum(v['charged_worker_seconds'] for v in ledger['records']))
        if full:assert not ledger['stop_reason'] and len(ledger['records'])==28 and all(v['returncode']==0 and v['status']=='COMPLETE' for v in ledger['records'])
    atomic_json(base/'independent_check.json',dict(passed=True,verdict=expected,raw_checks=checks,charged_worker_seconds=charge,dependencies='Same author/shared constants/libraries/frozen classifier. Separate sparse directed incidence, real-space dense convolution, explicit trapezoid/interaction/refinement arithmetic; not external replication, statistical confidence or physical identity.'))
    return checks
if __name__=='__main__':main()

"""MG092 separate donor/branch accounting, dense measurement and estimand audit."""
import os
for key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[key]='1'
from pathlib import Path
import json,hashlib,math
import numpy as np
from mg11g3 import factory,donor_path,donor_manifest,jobs,combined,OBS,PARENT_INDEX,SWITCH_TIME,S0
from mg11g3_oracle import accounting,dense_projection,late
from mg11g2_pattern_v2 import reference
from mg10p1 import continuity
from mg10d1_pool import atomic_json
ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/'results/MG11G3'
def read(p):return json.loads(Path(p).read_text(encoding='utf-8'))
def equal(a,b):np.testing.assert_allclose(a,b,rtol=1e-8,atol=1e-9)
def main(base=BASE,synthetic=False):
    base=Path(base);r=read(base/'registration.json');e=read(base/'evaluation.json');assert r['jobs']==jobs();complete={};fields={};scores={};checks=[];shams=[]
    if not synthetic:
        from prepare_mg11g3 import hashes
        assert hashes()==r['source_hashes'] and donor_manifest()==r['donor_hashes']
    for job in jobs():
        p=base/'production'/(job['id']+'.npz');case,ga,gb,k=factory(job)
        if not p.exists():complete[job['id']]=False;continue
        with np.load(p) as raw:
            for name,value in dict(initial=case['initial'],a0=case['a0'],b0=case['b0'],omega=case['omega'],diff=case['diff'],src=ga['src'],dst=ga['dst'],weight=ga['weight'],srcB=gb['src'],dstB=gb['dst'],weightB=gb['weight'],kappa=k,requested_times=OBS,parent_index=PARENT_INDEX,switch_time=SWITCH_TIME).items():np.testing.assert_array_equal(raw[name],value)
            assert json.loads(str(raw['job_json']))==job and str(raw['parent_sha256'])==hashlib.sha256(donor_path(job).read_bytes()).hexdigest();t=raw['times'];state=raw['state'];assert state.shape==(len(t),14000) and len(t)>=1 and np.array_equal(t,OBS[:len(t)]);np.testing.assert_array_equal(state[0,:3500],case['initial'].ravel());assert not np.any(state[0,3500:]);pop=state[:,:3500].reshape(len(t),7,500)
        v=accounting(state,case,ga,gb);assert v['passed'];v.update(file=p.relative_to(ROOT).as_posix(),sha256=hashlib.sha256(p.read_bytes()).hexdigest());checks.append(v);row=e['records'][job['id']];assert row['sha256']==v['sha256'] and row['checkpoints']==len(t)
        if job['arm']==2:assert not np.any(state[:,11500:13000])
        if not job['start']:np.testing.assert_array_equal(pop,np.broadcast_to(pop[:,:,:1],pop.shape))
        meta=read(p.with_suffix('.json')) if p.with_suffix('.json').exists() else {};full=np.array_equal(t,OBS) and meta.get('status')=='COMPLETE';complete[job['id']]=full;assert row['complete']==full
        if full:assert meta['sha256']==v['sha256'] and meta['job']==job
        merged,absolute=combined(job,pop,t);z=dense_projection(merged,reference()[0]);equal(z,row['metrics']['projection']);equal(merged[:,:3].mean(axis=-1),row['metrics']['field_mean'])
        stock=pop[:,3:5]/np.array([case['a0'],case['b0']])[None];equal(stock.mean(axis=-1),row['stock_mean']);equal(stock[:,:,(250+np.arange(-25,26))%500].mean(axis=-1),row['stock_local']);rx=state[:,3500:7000].reshape(len(t),7,500);equal(rx[:,[0,3]].sum(axis=-1),row['total_post_A_B_uptake'])
        if full:
            scores[job['id']]=late(z,absolute);equal(scores[job['id']],row['metrics']['late_projection']);assert continuity(merged[:,1],absolute,x0=0.)==row['metrics']['reference_continuity'];fields[job['id']]=pop[:,:3]
            if job['arm']==job['origin']:
                with np.load(donor_path(job)) as raw:target=raw['state'][40:,:1500].reshape(61,3,500)
                shams.append(dict(job=job['id'],active_error=float(np.max(abs(fields[job['id']]-target)/S0[None,:,None]))))
        else:assert 'late_projection' not in row['metrics']
    full=all(complete.values());assert full==e['complete'] and e['jobs_complete']==sum(complete.values());refinements=[]
    for o in (0,1):
        for s in (0,1):
            for a in (0,1,2):
                ids=[f'o{o}_s{s}_v{v}_a{a}' for v in (0,1)]
                if all(complete[id] for id in ids):refinements.append(dict(origin=o,start=s,arm=a,active_error=float(np.max(abs(fields[ids[1]]-fields[ids[0]])/S0[None,:,None]))))
    for actual,reported in ((refinements,e['refinement']),(shams,e['shams'])):
        assert len(actual)==len(reported)
        for aa,bb in zip(actual,reported):
            for name in aa:
                if name=='active_error':equal(aa[name],bb[name])
                else:assert aa[name]==bb[name]
    origins=[];means=[];off=[]
    for v in (0,1):
        ov=[];of=[]
        for o in (0,1):
            keys=[f'o{o}_s{s}_v{v}_a{a}' for s in (0,1) for a in (0,1)]
            ov.append((scores[keys[3]]-scores[keys[2]]-scores[keys[1]]+scores[keys[0]]).tolist() if all(complete[id] for id in keys) else None)
            keys=[f'o{o}_s{s}_v{v}_a{a}' for s in (0,1) for a in (o,2)]
            of.append((scores[keys[3]]-scores[keys[2]]-scores[keys[1]]+scores[keys[0]]).tolist() if all(complete[id] for id in keys) else None)
        origins.append(ov);off.append(of);means.append(None if any(x is None for x in ov) else ((np.array(ov[0])+np.array(ov[1]))*.5).tolist())
    def compare(actual,reported):
        if actual is None:assert reported is None
        elif isinstance(actual,list) and actual and isinstance(actual[0],(list,type(None))):
            assert len(actual)==len(reported)
            for aa,bb in zip(actual,reported):compare(aa,bb)
        else:equal(actual,reported)
    compare(origins,e['origin_contrasts_GXY']);compare(means,e['primary_GXY_by_setting']);compare(off,e['off_minus_sham_GXY']);numerical=full and len(refinements)==12 and all(x['active_error']<=1e-4 for x in refinements);sham_ok=full and len(shams)==8 and all(x['active_error']<=1e-4 for x in shams);delta=None if any(x is None for x in means) else abs(means[1][1]-means[0][1]);threshold=None if delta is None else max(.01,5*delta);resolved=bool(numerical and sham_ok and means[0][1]*means[1][1]>0 and abs(means[1][1])>threshold)
    assert numerical==e['numerical_qualified'] and sham_ok==e['sham_qualified'] and resolved==e['primary_resolved']
    if delta is None:assert e['primary_resolution_threshold'] is None
    else:equal(delta,e['primary_refinement_difference']);equal(threshold,e['primary_resolution_threshold'])
    verdict='INCOMPLETE' if not full else 'NUMERICAL/SHAM QUALIFICATION FAILED' if not numerical or not sham_ok else 'RESOLVED ONGOING B-TRANSPORT CONTRAST' if resolved else 'UNRESOLVED AT REGISTERED RESOLUTION';assert verdict==e['verdict'];charge=None
    if not synthetic:
        l=read(base/'production/ledger.json');charge=l['charged_worker_seconds'];assert l['ceiling']==3600 and l['prior_charge']==0 and charge<=3600;equal(charge,math.fsum(x['charged_worker_seconds'] for x in l['records']))
        if full:assert not l['stop_reason'] and len(l['records'])==24 and all(x['status']=='COMPLETE' and x['returncode']==0 for x in l['records'])
    atomic_json(base/'independent_check.json',dict(passed=True,raw_checks=checks,verdict=verdict,charged_worker_seconds=charge,dependencies='Same author/shared laws/constants/NumPy/SciPy/frozen classifier. Separate two-graph incidence, dense convolution/explicit trapezoid/conditional contrasts and sham/refinement arithmetic; not external replication.'))
    return checks
if __name__=='__main__':main()

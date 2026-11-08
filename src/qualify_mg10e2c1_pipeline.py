"""MG076 mixed-epoch reuse/replay provenance qualification; algebraic only."""
import os
for key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[key]='1'
from pathlib import Path
import json,time
import numpy as np
from mg10d1_pool import atomic_json
from mg10p2 import fixture
from mg10e2 import activity
from mg10e2c1_replay import PrefixGuard
from run_mg10e2 import save_npz
from run_mg10e2c1 import read,digest,build_manifest,stage_reuse,hashes
from evaluate_mg10e2c1 import main as evaluate
from check_mg10e2c1 import main as check
ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/'results/MG10E2C1/preparation'
def arrays(case,job,frames,parent_sha=''):
    state=case['initial'];data=np.zeros((1,frames,22,500),dtype=np.int64);data[0,:,:7]=state;arm=job.get('arm',0)
    return dict(data=data,initial=state,a0=case['a0'],b0=case['b0'],omega=1000.,clamp_a=True,clamp_b=True,stats=np.zeros((1,3)),times=np.arange(frames)/4,seed=job['seed'],arm=arm,switch=5.,activity_before=activity(arm,0),activity_after=activity(arm,5),parent_sha256=parent_sha)
def main():
    started=time.perf_counter()
    if BASE.exists():raise ValueError('Already attempted; preserve/version failures')
    BASE.mkdir(parents=True);source=BASE/'SYNTHETIC_original';combined=BASE/'SYNTHETIC_combined';source.mkdir();combined.mkdir()
    original=read(ROOT/'results/MG10E2/registration.json');parents=original['parent_jobs'];branches=original['branch_jobs']
    atomic_json(source/'registration.json',dict(parent_jobs=parents,branch_jobs=branches,scope='ALGEBRAIC ONLY; no scientific random draws'))
    first=branches[:8];reused='a1_s1_k0_b0';parent_states={}
    for job in parents:
        case=fixture(0,job['start']);parent_states[job['id']]=case;p=source/'parents'/(job['id']+'.npz');p.parent.mkdir(exist_ok=True);save_npz(p,**arrays(case,job,21))
        atomic_json(p.with_suffix('.json'),dict(status='COMPLETE',job=job,raw_sha256=digest(p),scope='ALGEBRAIC ONLY'))
    for job in first:
        case=parent_states[job['parent']];psha=digest(source/'parents'/(job['parent']+'.npz'));path=source/'branches'/(job['id']+'.npz');path.parent.mkdir(exist_ok=True)
        save_npz(path,**arrays(case,job,101 if job['id']==reused else 96,psha))
        atomic_json(path.with_suffix('.json'),dict(status='COMPLETE' if job['id']==reused else 'BUDGET-SHORT: FORECAST',job=job,raw_sha256=digest(path),parent_sha256=psha,scope='ALGEBRAIC ONLY'))
    atomic_json(source/'parents/ledger.json',dict(charged_worker_seconds=2.5,scope='SYNTHETIC CHARGE'))
    atomic_json(source/'branches/ledger.json',dict(charged_worker_seconds=10.5,scope='SYNTHETIC CHARGE'))
    manifest=build_manifest(source);assert list(manifest['modes'].values()).count('EXACT_REPLAY')==7 and list(manifest['modes'].values()).count('UNSTARTED')==120
    atomic_json(combined/'registration.json',dict(parent_jobs=parents,branch_jobs=branches,**manifest))
    stage_reuse(combined,manifest);records=[]
    for job in branches:
        mode=manifest['modes'][job['id']]
        if mode=='REUSE':continue
        case=parent_states[job['parent']];psha=digest(combined/'parents'/(job['parent']+'.npz'));full=arrays(case,job,101,psha);proof=None
        if mode=='EXACT_REPLAY':
            guard=PrefixGuard(source/'branches'/(job['id']+'.npz'),lambda *a,**kw:None)
            for n in range(1,97):
                prefix=dict(full);prefix['data']=full['data'][:,:n];prefix['times']=full['times'][:n];guard(combined/'not_saved.npz',**prefix)
            guard(combined/'not_saved.npz',**full);proof=guard.finish()
        path=combined/'branches'/(job['id']+'.npz');save_npz(path,**full)
        atomic_json(path.with_suffix('.json'),dict(status='COMPLETE',job=job,raw_sha256=digest(path),parent_sha256=psha,continuation_mode=mode,prefix_proof=proof,scope='ALGEBRAIC ONLY'))
        records.append(dict(id=job['id'],job=job,status='COMPLETE',returncode=0,charged_worker_seconds=.01))
    ledger=dict(records=records,ceiling=518400,prior_charge=10.5,charged_worker_seconds=10.5+sum(v['charged_worker_seconds'] for v in records),stop_reason=None,scope='SYNTHETIC CHARGE')
    atomic_json(combined/'branches/ledger.json',ledger)
    evaluate(combined);check(combined,synthetic=True);v=read(combined/'independent_check.json');assert v['passed'] and v['verified_complete_replays']==7 and v['status']=='COMPLETE RESTORATION ESTIMATION'
    for name in ('evaluation.json','independent_check.json','MG10E2C1_results.md'):(combined/name).rename(combined/('complete_'+name))
    missing=combined/'branches/a0_s0_k1_b0.npz';missing.rename(missing.with_suffix('.npz.withheld'))
    evaluate(combined);check(combined,synthetic=True);v=read(combined/'evaluation.json');assert [v['primary_by_step'][str(k)]['n'] for k in range(2)]==[8,7]
    for name in ('evaluation.json','independent_check.json','MG10E2C1_results.md'):(combined/name).rename(combined/('missing_'+name))
    with np.load(missing.with_suffix('.npz.withheld')) as raw:prefix={key:raw[key] for key in raw.files}
    prefix['data']=prefix['data'][:,:81];prefix['times']=prefix['times'][:81];save_npz(missing,**prefix)
    evaluate(combined);check(combined,synthetic=True);v=read(combined/'evaluation.json');assert v['branches_complete']==127 and [v['primary_by_step'][str(k)]['n'] for k in range(2)]==[8,7]
    for name in ('evaluation.json','independent_check.json','MG10E2C1_results.md'):(combined/name).rename(combined/('prefix_'+name))
    # Corrupt one reused original metadata file; immutable-source check must fail.
    target=source/'parents'/(parents[0]['id']+'.json');saved=target.read_bytes();target.write_bytes(saved+b' ')
    from run_mg10e2c1 import verify_original
    try:verify_original(manifest);raise AssertionError('Changed original accepted')
    except ValueError:pass
    finally:target.write_bytes(saved)
    verify_original(manifest)
    elapsed=time.perf_counter()-started;assert elapsed<300
    atomic_json(BASE/'pipeline_qualification.json',dict(passed=True,source_hashes=hashes(),elapsed_seconds=elapsed,reused_parents=32,reused_complete_branches=1,guarded_prefixes=7,unstarted_jobs=120,synthetic_full_branches=128,synthetic_missing_prefix_counts=[8,7],original_tamper_rejected=True,
      scope='Zero-event algebraic full/missing/prefix mixed-provenance raw-file pipeline; exact seven-prefix guards, copy hashes, original preservation and aggregate prior charges audited. No E2 scientific seed used or original scientific path replayed. Same author/shared kernel/reference/libraries, not external replication.'))
    print(json.dumps({k:v for k,v in read(BASE/'pipeline_qualification.json').items() if k!='source_hashes'},indent=2))
if __name__=='__main__':main()

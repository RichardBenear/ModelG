"""MG074 held-stream production-writer replay and algebraic raw-file pipeline."""
import os
for key in ("OPENBLAS_NUM_THREADS","MKL_NUM_THREADS","OMP_NUM_THREADS","NUMEXPR_NUM_THREADS"):os.environ[key]="1"
from pathlib import Path
import json,time
import numpy as np
from mg10d1_pool import atomic_json
from mg10p2 import fixture,reference
from mg10e2 import activity,fork_case
from prepare_mg10e2 import candidate_jobs
from runner import seed_for
from run_mg10e2 import trajectory,save_npz,digest,hashes
from evaluate_mg10e2 import main as evaluate
from check_mg10e2 import main as check
from check_mg10d1_development import check_raw
ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/'results/MG10E2/preparation'
def main():
    started=time.perf_counter()
    if (BASE/'pipeline_qualification.json').exists() or (BASE/'writer_replays').exists():raise ValueError('Already attempted; preserve/version any correction')
    replays=BASE/'writer_replays';replays.mkdir();count=0;raw_checks=[]
    for s in range(2):
        for k in range(2):
            with np.load(ROOT/f'results/MG10E1/preparation/devparents/parent_s{s}_k{k}.npz') as raw:case=fork_case(fixture(0,s),raw['data'][0,-1,:7])
            for a in range(4):
                source=BASE/f'timing/timing_a{a}_s{s}_k{k}.npz'
                with np.load(source) as raw:expected=raw['data'].copy();st=raw['stats'].copy();seed=int(raw['seed']);times=raw['times'].copy()
                job=dict(id=f'replay_a{a}_s{s}_k{k}',kind='held_development_replay',arm=a,start=s,step=k)
                path=replays/(job['id']+'.npz')
                trajectory(case,(1/24000,1/48000)[k],times,seed,seed_for(10,s,k,7302,a+1,1),a,.125,path,job,90)
                with np.load(path) as raw:np.testing.assert_array_equal(raw['data'],expected);np.testing.assert_array_equal(raw['stats'],st)
                raw_checks.append(check_raw(path));count+=1
    # Real production writer's parent all-one mode reproduces held ungated E1
    # development paths, without using any E2 scientific RNG stream.
    for s in range(2):
        for k in range(2):
            source=ROOT/f'results/MG10E1/preparation/devparents/parent_s{s}_k{k}.npz'
            with np.load(source) as raw:expected=raw['data'].copy();st=raw['stats'].copy();seed=int(raw['seed']);times=raw['times'].copy()
            job=dict(id=f'parent_replay_s{s}_k{k}',kind='held_development_replay',start=s,step=k)
            path=replays/(job['id']+'.npz');trajectory(fixture(0,s),(1/24000,1/48000)[k],times,seed,seed_for(10,s,k,6901,0,1),0,5.,path,job,90)
            with np.load(path) as raw:np.testing.assert_array_equal(raw['data'],expected);np.testing.assert_array_equal(raw['stats'],st)
            raw_checks.append(check_raw(path));count+=1
    synth=BASE/'SYNTHETIC_pipeline_fixture';synth.mkdir();p,b=candidate_jobs()
    atomic_json(synth/'registration.json',dict(parent_jobs=p,branch_jobs=b,scope='ALGEBRAIC ONLY; no E2 scientific random draw'))
    ref=reference()[0]
    for directory,jobs,frames in (('parents',p,21),('branches',b,101)):
        d=synth/directory;d.mkdir()
        for job in jobs:
            case=fixture(0,job['start']);state=case['initial'].copy();arm=job.get('arm',0)
            if directory=='branches':
                amp=job['start']*(.7+.03*arm+.001*job['block']*arm)
                # Some translated fields exercise the new movement diagnostics.
                shift=125 if arm==3 else 0
                state[:3]=np.rint(1000*(case['background'][:,None]+amp*np.roll(ref,shift,axis=-1))).astype(np.int64)
            data=np.zeros((1,frames,22,500),dtype=np.int64);data[0,:,:7]=state;path=d/(job['id']+'.npz')
            save_npz(path,data=data,initial=state,a0=case['a0'],b0=case['b0'],omega=1000.,clamp_a=True,clamp_b=True,stats=np.zeros((1,3)),times=np.arange(frames)/4,
                     arm=arm,switch=5.,activity_before=activity(arm,0),activity_after=activity(arm,5),seed=job['seed'],parent_sha256='ALGEBRAIC fixture: no lineage claim')
            atomic_json(path.with_suffix('.json'),dict(status='COMPLETE',job=job,raw_sha256=digest(path),scope='ALGEBRAIC INPUT ONLY'))
    evaluate(synth);check(synth,synthetic=True)
    value=json.loads((synth/'evaluation.json').read_text());assert value['status']=='COMPLETE RESTORATION ESTIMATION' and all(value['primary_by_step'][str(k)]['n']==8 for k in range(2))
    for name in ('evaluation.json','independent_check.json','MG10E2_results.md'):(synth/name).rename(synth/('complete_'+name))
    missing=synth/'branches/a2_s1_k0_b0.npz';missing.rename(missing.with_suffix('.npz.withheld'))
    evaluate(synth);check(synth,synthetic=True);value=json.loads((synth/'evaluation.json').read_text());assert value['branches_complete']==127 and [value['primary_by_step'][str(k)]['n'] for k in range(2)]==[7,8]
    for name in ('evaluation.json','independent_check.json','MG10E2_results.md'):(synth/name).rename(synth/('missing_'+name))
    with np.load(missing.with_suffix('.npz.withheld')) as raw:arrays={key:raw[key] for key in raw.files}
    arrays['data']=arrays['data'][:,:81];arrays['times']=arrays['times'][:81];save_npz(missing,**arrays)
    evaluate(synth);check(synth,synthetic=True);value=json.loads((synth/'evaluation.json').read_text());assert value['branches_available']==128 and value['branches_complete']==127 and [value['primary_by_step'][str(k)]['n'] for k in range(2)]==[7,8]
    for name in ('evaluation.json','independent_check.json','MG10E2_results.md'):(synth/name).rename(synth/('prefix_'+name))
    # Also exercise an early parent-phase stop: no branches and no primary CI.
    empty=BASE/'SYNTHETIC_parent_stop';empty.mkdir();atomic_json(empty/'registration.json',dict(parent_jobs=p,branch_jobs=b));evaluate(empty);check(empty,synthetic=True)
    value=json.loads((empty/'evaluation.json').read_text());assert value['parents_complete']==0 and value['branches_complete']==0 and all(value['primary_by_step'][str(k)]['n']==0 for k in range(2))
    elapsed=time.perf_counter()-started;prior=json.loads((BASE/'independent_check.json').read_text())['total_preparation_allowance_seconds'];assert prior+elapsed<900
    atomic_json(BASE/'pipeline_qualification.json',dict(passed=True,source_hashes=hashes(),held_development_writer_replays=count,raw_replay_checks=raw_checks,
       synthetic_full_counts=[32,128],synthetic_missing_prefix_primary_counts=[7,8],empty_parent_stop_audited=True,additional_seconds=elapsed,total_preparation_allowance_seconds=prior+elapsed,
       scope='Held development stream replays through actual atomic writer; full/missing/prefix/empty algebraic raw-file pipeline and dense movement/statistical audit. No E2 scientific trajectory or new production seed used.'))
    print(json.dumps({key:value for key,value in json.loads((BASE/'pipeline_qualification.json').read_text()).items() if key not in ('raw_replay_checks','source_hashes')},indent=2))
if __name__=='__main__':main()

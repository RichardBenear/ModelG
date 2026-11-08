"""MG084 held development replay and algebraic full/missing/prefix/law-failure tests."""
import os
for key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[key]='1'
from pathlib import Path
import time,json,traceback
import numpy as np
from scipy.linalg import expm
from mg11g1 import cell,generator
from mg11g1_sample import generate,raw_save
from run_mg11g1 import hashes,digest,read
from prepare_mg11g1 import candidate_jobs
from evaluate_mg11g1 import main as evaluate
from check_mg11g1 import main as check
from mg10d1_pool import atomic_json
ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/'results/MG11G1/preparation/pipeline'

def history(resource,network,reps):
    case,g,k,t,m=cell(resource,network);out=np.zeros((reps,len(t)),dtype=np.int64)
    for obs in range(1,len(t)):
        P=np.maximum(expm(generator(case,g,k,resource)*(t[obs]-t[obs-1])),0);P/=P.sum(axis=1,keepdims=True)
        for a in range(len(P)):
            ids=np.flatnonzero(out[:,obs-1]==a);expected=P[a]*len(ids);counts=np.floor(expected).astype(int);counts[np.argsort(-(expected-counts),kind='stable')[:len(ids)-counts.sum()]]+=1;offset=0
            for b,count in enumerate(counts):out[ids[offset:offset+count],obs]=b;offset+=count
    return out

def manufacture(resource,network,h):
    case,g,k,t,m=cell(resource,network);reps=len(h);data=np.zeros((reps,5,20,6),dtype=np.int64);env=np.zeros((reps,5,2,len(g['src'])),dtype=np.int64);stats=np.zeros((reps,3));rx=np.zeros((reps,7,6),dtype=np.int64);ef=env[:,0].copy();data[:,0,:7]=case['initial'];edges={(int(i),int(j)):e for e,(i,j) in enumerate(zip(g['src'],g['dst']))}
    def route(ids,start,end):
        queue=[(start,[])];seen={start}
        while queue:
            current,path=queue.pop(0)
            if current==end:
                for edge in path:ef[ids,resource,edge]+=1
                return
            for (i,j),edge in edges.items():
                if i==current and j not in seen:seen.add(j);queue.append((j,path+[edge]))
        raise AssertionError('Disconnected fixture')
    for obs in range(1,5):
        for a in range(12 if resource==0 else 7):
            for b in range(12 if resource==0 else 7):
                ids=np.flatnonzero((h[:,obs-1]==a)&(h[:,obs]==b))
                if not len(ids):continue
                if a<6:
                    route(ids,a,(b if b<6 else b-6) if resource==0 else (b if b<6 else 3))
                    if b>=6:rx[ids,0 if resource==0 else 3,b-6 if resource==0 else 3]+=1
                else:assert a==b
        pop=np.zeros((reps,7,6),dtype=np.int64);tag=h[:,obs];alive=tag<6;ids=np.flatnonzero(alive);pop[ids,3+resource,tag[ids]]=1
        if resource==0:ids=np.flatnonzero(~alive);pop[ids,0,tag[ids]-6]=1
        else:pop[:,1,3]=alive;pop[:,2,3]=~alive;pop[:,5,3]=~alive
        data[:,obs,:7]=pop;data[:,obs,7:14]=rx;env[:,obs]=ef
    stats[:,0]=rx.sum(axis=(1,2))+ef.sum(axis=(1,2));return data,env,stats

def archive(base,prefix):
    for name in ('evaluation.json','independent_check.json','MG11G1_results.md'):(base/name).rename(base/(prefix+'_'+name))

def main():
    started=time.perf_counter()
    if BASE.exists():raise ValueError('Already attempted; retain/version failure')
    BASE.mkdir();held=BASE/'held_replays';held.mkdir()
    for resource in range(2):
        for network in range(2):
            old=read(BASE.parent/f'timing/timing_r{resource}_g{network}.json');job=dict(old['job'],id=f'held_r{resource}_g{network}')
            generate(job,held,block=64)
            with np.load(held/(job['id']+'.npz')) as a,np.load(BASE.parent/'timing'/old['raw']) as b:
                for key in a.files:np.testing.assert_array_equal(a[key],b[key])
    synth=BASE/'SYNTHETIC_pipeline';directory=synth/'production';directory.mkdir(parents=True);jobs=candidate_jobs();atomic_json(synth/'registration.json',dict(jobs=jobs,formal_passed=True,scope='ALGEBRAIC ONLY; no scientific random draws'))
    for resource in range(2):
        for network in range(2):
            h=history(resource,network,65536);case,g,k,t,m=cell(resource,network)
            for job in [j for j in jobs if j['resource']==resource and j['network']==network]:
                left=job['chunk']*4096;path=directory/(job['id']+'.npz');arrays=manufacture(resource,network,h[left:left+4096]);raw_save(path,case,g,k,t,m,job['seed'],*arrays,4096)
                raw_save(directory/(job['id']+'_warm.npz'),case,g,k,t,m,job['warmup_seed'],*manufacture(resource,network,h[:1]),1)
                atomic_json(path.with_suffix('.json'),dict(status='COMPLETE',job=job,sha256=digest(path),scope='ALGEBRAIC ONLY'))
    evaluate(synth);check(synth,synthetic=True);v=read(synth/'evaluation.json');assert v['complete'] and v['law_supported'] and v['connectivity_effect_supported'];archive(synth,'complete')
    # Missing B suppresses primary; missing A preserves B CI but suppresses interpretation.
    for name in ('r1_g0_c15','r0_g0_c15'):
        path=directory/(name+'.npz');saved=path.with_suffix('.npz.withheld');path.rename(saved);evaluate(synth);check(synth,synthetic=True);v=read(synth/'evaluation.json');assert not v['complete'] and not v['connectivity_effect_supported'] and all(x['passed'] is None for x in v['analytic']);assert (v['primary']['ci'] is None)==name.startswith('r1');archive(synth,'missing_'+name)
        if name.startswith('r0'):
            with np.load(saved) as raw:prefix={key:raw[key] for key in raw.files}
            for key in ('data','environment','stats'):prefix[key]=prefix[key][:128]
            prefix['reps_completed']=128
            with path.open('wb') as f:np.savez_compressed(f,**prefix)
            evaluate(synth);check(synth,synthetic=True);v=read(synth/'evaluation.json');assert not v['complete'] and v['paths_available']==258176;archive(synth,'prefix');path.rename(path.with_suffix('.npz.prefix_fixture'))
        saved.rename(path)
    # Full, valid-accounting B paths with forced early capture must fail the law family.
    job=next(j for j in jobs if j['id']=='r1_g0_c0');path=directory/(job['id']+'.npz');saved=path.with_suffix('.npz.baseline');path.rename(saved);meta=path.with_suffix('.json');oldmeta=meta.read_bytes();case,g,k,t,m=cell(1,0);h=np.full((4096,5),6);h[:,0]=0
    raw_save(path,case,g,k,t,m,job['seed'],*manufacture(1,0,h),4096);atomic_json(meta,dict(status='COMPLETE',job=job,sha256=digest(path),scope='DELIBERATE VALID-ACCOUNTING LAW FAILURE'))
    evaluate(synth);check(synth,synthetic=True);v=read(synth/'evaluation.json');assert v['complete'] and not v['law_supported'] and not v['connectivity_effect_supported'];archive(synth,'failed_law');path.rename(path.with_suffix('.npz.failed_fixture'));saved.rename(path);meta.write_bytes(oldmeta)
    empty=BASE/'SYNTHETIC_empty';empty.mkdir();atomic_json(empty/'registration.json',dict(jobs=jobs,formal_passed=True));evaluate(empty);check(empty,synthetic=True);assert read(empty/'evaluation.json')['paths_available']==0
    elapsed=time.perf_counter()-started;assert elapsed<300
    atomic_json(BASE.parent/'pipeline_qualification.json',dict(passed=True,source_hashes=hashes(),elapsed_seconds=elapsed,held_stream_replays=4,synthetic_full_paths=262144,synthetic_prefix_paths=258176,complete_missing_B_missing_A_prefix_empty_law_failure=True,scope='Held development replays and deterministic rounded Markov transitions routed through valid graph edges; not scientific samples, not external replication. Chunkwise physical-array processing.'))
    print(json.dumps(dict(passed=True,elapsed_seconds=elapsed)))
if __name__=='__main__':
    try:main()
    except BaseException:atomic_json(BASE.parent/'pipeline_qualification_failure.json',dict(status='FAILED; preserve before correction',error=traceback.format_exc()));raise

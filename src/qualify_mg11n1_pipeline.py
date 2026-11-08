"""MG080 held-stream equivalence and algebraic pass/fail/incomplete pipeline."""
import os
for key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[key]='1'
from pathlib import Path
import json,time,traceback
import numpy as np
from scipy.linalg import expm
from mg11n1_cases import cell
from mg10d1_pool import atomic_json
from run_mg11n1 import generate,raw_save,hashes,digest,read
from prepare_mg11n1 import candidate_jobs
from evaluate_mg11n1 import main as evaluate
from check_mg11n1 import main as check
ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/'results/MG11N1/preparation/pipeline'
def transition(case,g,k,kind,dt):
    n=g['n'];W=np.zeros((n,n))
    for i,j,w in zip(g['src'],g['dst'],g['weight']):W[i,j]=w
    L=W-np.diag(W.sum(axis=1))
    if kind=='transport':return np.kron(expm(dt*k[0]*L),expm(dt*k[1]*L))
    size=2*n if kind=='A_source' else n+1;Q=np.zeros((size,size));s=0 if kind=='A_source' else 1;rate=14*case['omega']/case['a0'] if s==0 else 29*case['initial'][1]/case['b0'];Q[:n,:n]=k[s]*L-np.diag(rate)
    for i in range(n):Q[i,n+i if s==0 else n]=rate[i]
    return expm(dt*Q)
def manufacture(c,reps):
    case,g,k,t,m=cell(c);n=g['n'];data=np.zeros((reps,len(t),20,n),dtype=np.int64);env=np.zeros((reps,len(t),2,len(g['src'])),dtype=np.int64);stats=np.zeros((reps,3));pop=np.broadcast_to(case['initial'],(reps,7,n)).copy();rx=np.zeros((reps,7,n),dtype=np.int64);ef=np.zeros((reps,2,len(g['src'])),dtype=np.int64);data[:,0,:7]=pop
    if c==3:data[:,:,:7]=pop[:,None];return data,env,stats
    labels=np.full(reps,2 if c in (0,2) else 0,dtype=np.int64);edge_index={(int(i),int(j)):edge for edge,(i,j) in enumerate(zip(g['src'],g['dst']))}
    def route(indices,s,start,end):
        for i in range(start,end,1 if end>start else -1):
            j=i+(1 if end>start else -1);ef[indices,s,edge_index[i,j]]+=1
    for obs in range(1,len(t)):
        P=transition(case,g,k,case['name'],float(t[obs]-t[obs-1]));P=np.maximum(P,0);P/=P.sum(axis=1,keepdims=True);new=np.empty_like(labels)
        for a in range(len(P)):
            indices=np.flatnonzero(labels==a);expected=P[a]*len(indices);counts=np.floor(expected).astype(int);remainder=len(indices)-int(counts.sum());order=np.argsort(-(expected-counts),kind='stable');counts[order[:remainder]]+=1;offset=0
            for b,count in enumerate(counts):
                ids=indices[offset:offset+count];offset+=count;new[ids]=b
                if not count:continue
                if c==0:
                    route(ids,0,a//n,b//n);route(ids,1,a%n,b%n)
                elif c==1:
                    if a<n:
                        route(ids,0,a,b if b<n else b-n)
                        if b>=n:rx[ids,0,b-n]+=1
                    else:assert a==b
                else:
                    if a<n:
                        route(ids,1,a,b if b<n else 0)
                        if b==n:rx[ids,3,0]+=1
                    else:assert a==b
        labels=new
        if c==0:
            pop[:,3:5]=0;pop[np.arange(reps),3,labels//n]=1;pop[np.arange(reps),4,labels%n]=1
        elif c==1:
            pop[:,0]=0;pop[:,3]=0;alive=labels<n;ids=np.flatnonzero(alive);pop[ids,3,labels[ids]]=1;ids=np.flatnonzero(~alive);pop[ids,0,labels[ids]-n]=1
        else:
            pop[:]=0;converted=labels==n;pop[:,1,0]=~converted;pop[:,2,0]=converted;pop[:,5,0]=converted;ids=np.flatnonzero(~converted);pop[ids,4,labels[ids]]=1
        data[:,obs,:7]=pop;data[:,obs,7:14]=rx;env[:,obs]=ef
    stats[:,0]=rx.sum(axis=(1,2))+ef.sum(axis=(1,2));return data,env,stats
def archive_reports(base,prefix):
    for name in ('evaluation.json','independent_check.json','MG11N1_results.md'):(base/name).rename(base/(prefix+'_'+name))
def main():
    started=time.perf_counter()
    if BASE.exists():raise ValueError('Already attempted; preserve/version failures')
    BASE.mkdir();held=BASE/'held_replays';held.mkdir();replays=0
    from runner import seed_for
    for c in range(4):
        for engine in range(2 if c==3 else 1):
            meta=read(ROOT/f'results/MG11N1/preparation/timing/timing_c{c}_m{engine}.json');seed=meta['seed'];job=dict(id=f'held_c{c}_m{engine}',cell=c,engine=engine,reps=meta['reps'],seed=seed,warmup_seed=seed_for(11,c,engine,7901,0,1))
            generate(job,held,block=64)
            with np.load(held/(job['id']+'.npz')) as a,np.load(ROOT/f'results/MG11N1/preparation/timing/{meta["raw"]}') as b:
                for key in ('data','environment','stats','initial','a0','b0','times','seed'):np.testing.assert_array_equal(a[key],b[key])
            replays+=1
    synth=BASE/'SYNTHETIC_pipeline';synth.mkdir();jobs=candidate_jobs();atomic_json(synth/'registration.json',dict(jobs=jobs,formal_passed=True,scope='ALGEBRAIC ONLY; no new scientific random draws'))
    arrays={c:manufacture(c,32768 if c<3 else 16384) for c in range(4)}
    for job in jobs:
        c=job['cell'];case,g,k,t,m=cell(c);a,b,stats=arrays[c];left=job['chunk']*4096;right=left+4096;directory=synth/'production';directory.mkdir(exist_ok=True);path=directory/(job['id']+'.npz')
        raw_save(path,case,g,k,t,m,job['seed'],a[left:right],b[left:right],stats[left:right],4096)
        warm=manufacture(c,1);raw_save(directory/(job['id']+'_warm.npz'),case,g,k,t,m,job['warmup_seed'],*warm,1)
        atomic_json(path.with_suffix('.json'),dict(status='COMPLETE',job=job,sha256=digest(path),scope='ALGEBRAIC ONLY'))
    evaluate(synth);check(synth,synthetic=True);assert read(synth/'evaluation.json')['status']=='QUALIFIED WITHIN REGISTERED N1 SCOPE';archive_reports(synth,'complete')
    # Valid-accounting deliberate distribution failure: one quarter of oracle
    # paths undergo an extra A->G reaction at node0 at both positive endpoints.
    path=synth/'production/c3_m1_r0.npz';backup=path.with_suffix('.npz.baseline');path.rename(backup);meta=path.with_suffix('.json');saved_meta=meta.read_bytes()
    with np.load(backup) as raw:modified={key:raw[key] for key in raw.files}
    modified['data'][:,1:,3,0]-=1;modified['data'][:,1:,0,0]+=1;modified['data'][:,1:,7,0]+=1;modified['stats'][:,0]+=1
    with path.open('wb') as f:np.savez_compressed(f,**modified)
    atomic_json(meta,dict(status='COMPLETE',sha256=digest(path),scope='DELIBERATE ALGEBRAIC FAILURE'))
    evaluate(synth);check(synth,synthetic=True);assert read(synth/'evaluation.json')['status']=='NOT QUALIFIED UNDER REGISTERED N1 CRITERIA';archive_reports(synth,'failed_criteria')
    path.rename(path.with_suffix('.npz.failed_fixture'));meta.rename(meta.with_suffix('.json.failed_fixture'));backup.rename(path);meta.write_bytes(saved_meta)
    missing=synth/'production/c0_m0_r7.npz';missing.rename(missing.with_suffix('.npz.withheld'));evaluate(synth);check(synth,synthetic=True)
    v=read(synth/'evaluation.json');assert v['status']=='INCOMPLETE N1 QUALIFICATION' and all(row['passed'] is None for row in v['analytic']);archive_reports(synth,'missing')
    with np.load(missing.with_suffix('.npz.withheld')) as raw:partial={key:raw[key] for key in raw.files}
    for key in ('data','environment','stats'):partial[key]=partial[key][:128]
    partial['reps_completed']=128
    with missing.open('wb') as f:np.savez_compressed(f,**partial)
    evaluate(synth);check(synth,synthetic=True);v=read(synth/'evaluation.json');assert v['status']=='INCOMPLETE N1 QUALIFICATION' and v['jobs_complete']==31 and v['paths_available']==127104;archive_reports(synth,'prefix')
    empty=BASE/'SYNTHETIC_empty';empty.mkdir();atomic_json(empty/'registration.json',dict(jobs=jobs,formal_passed=True));evaluate(empty);check(empty,synthetic=True);assert read(empty/'evaluation.json')['paths_available']==0
    elapsed=time.perf_counter()-started;assert elapsed<300
    atomic_json(BASE.parent/'pipeline_qualification.json',dict(passed=True,source_hashes=hashes(),elapsed_seconds=elapsed,held_stream_writer_replays=replays,synthetic_pass_fail_missing_prefix_empty=True,synthetic_full_paths=131072,synthetic_prefix_paths=127104,
       scope='Held development seeds only; manufactured categorical transitions rounded deterministically and routed through valid edges, full-tiny static states, deliberate valid-accounting law failure. Algebraic pipeline fixtures are not scientific samples. Same author/shared laws/reference/libraries, not external replication.'))
    print(json.dumps({k:v for k,v in read(BASE.parent/'pipeline_qualification.json').items() if k!='source_hashes'},indent=2))
if __name__=='__main__':
    try:main()
    except BaseException:atomic_json(BASE.parent/'pipeline_qualification_failure.json',dict(status='FAILED; preserve before correction',error=traceback.format_exc()));raise

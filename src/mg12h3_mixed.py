"""MG124 bounded actual-writer chunks for frozen mixed incidence engine."""
from prepare_mg12h3_bridge import *
from check_mg12h3_bridge import audit_trace
BASE=ROOT/'results/MG12H3/mixed_qualification';R_FORMAL=32768;CHUNK=4096
FILES=('src/mg12h3_mixed.py','src/check_mg12h3_mixed.py','src/prepare_mg12h3_mixed.py','src/run_mg12h3_mixed.py','spec/MG12H3_mixed_qualification.md','src/mg10e2c1_pool.py')
def mixed_hashes():
    from prepare_mg12h3_bridge import hashes as old_hashes
    result=old_hashes();result.update({name:hashlib.sha256((ROOT/name).read_text(encoding='utf-8').replace('\r\n','\n').encode()).hexdigest() for name in FILES})
    for c in range(4):
        p=ROOT/f'results/MG12H3/bridge/c{c}_reference.npz';result[p.relative_to(ROOT).as_posix()]=hashlib.sha256(p.read_bytes()).hexdigest()
    return result
def canonical_jobs():return [dict(id=f'c{c}_m{m}_b{b}',cell=c,method=m,chunk=b,reps=CHUNK) for c in range(4) for m in range(2) for b in range(8)]
def write_counts(path,data,seeds,job):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True);tmp=path.with_suffix('.npz.tmp')
    with tmp.open('wb') as f:np.savez_compressed(f,counts=data,seeds=np.array(seeds,np.uint32),times=TIMES,job_json=json.dumps(job))
    os.replace(tmp,path)
def generate(job,seeds,directory,deadline):
    directory=Path(directory);directory.mkdir(parents=True,exist_ok=True);state=fixture(job['cell']);tau=None if job['method']==0 else .001;total=job['reps'];assert len(seeds)==total and len(set(seeds))==total;data=np.zeros((total,3,7,len(state.sites)),np.int64);path=directory/(job['id']+'.npz');start=time.perf_counter();write_counts(path,data[:0],[],job);accepted=0
    with (directory/(job['id']+'_traces.jsonl')).open('x',encoding='utf-8') as f:
        for r,seed in enumerate(seeds):
            if time.perf_counter()>deadline:raise TimeoutError('Chunk deadline')
            result=simulate(state,TIMES,seed,tau,deadline);accepted+=audit_trace(state,result['counts'],result['trace'],job['method']);data[r]=result['counts'];f.write(json.dumps(dict(replicate=r,seed=seed,trace=result['trace']),separators=(',',':'))+'\n')
            if (r+1)%512==0 or r+1==total:f.flush();write_counts(path,data[:r+1],seeds[:r+1],job);atomic_json(directory/(job['id']+'_progress.json'),dict(completed=r+1,requested=total))
    meta=dict(status='COMPLETE',job=job,reps=total,elapsed_seconds=time.perf_counter()-start,accepted_events_audited=accepted,sha256=hashlib.sha256(path.read_bytes()).hexdigest());atomic_json(directory/(job['id']+'.json'),meta);return meta

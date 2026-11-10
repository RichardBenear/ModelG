"""MG124 raw-chunk reread and joint-count endpoint-law qualification arithmetic."""
from mg12h3_mixed import *
def audit_chunk(job,seeds,directory,deadline=math.inf):
    directory=Path(directory);p=directory/(job['id']+'.npz');metadata=directory/(job['id']+'.json')
    if not p.exists():return dict(complete=False,missing=True,paths=0),None
    meta=read(metadata) if metadata.exists() else {};state=fixture(job['cell'])
    with np.load(p) as z:
        data=z['counts'];assert json.loads(str(z['job_json']))==job;np.testing.assert_array_equal(z['times'],TIMES);assert z['seeds'].tolist()==seeds[:len(data)]
    assert data.dtype==np.int64 and data.shape==(len(data),3,7,len(state.sites)) and len(data)<=job['reps'] and (len(data)==0 or data.min()>=0)
    if meta.get('status')=='COMPLETE':assert meta['sha256']==hashlib.sha256(p.read_bytes()).hexdigest() and meta['job']==job
    complete=meta.get('status')=='COMPLETE' and len(data)==job['reps'] and meta.get('reps')==job['reps'];accepted=0;lines=0;tracepath=directory/(job['id']+'_traces.jsonl')
    if tracepath.exists():
        with tracepath.open(encoding='utf-8') as f:
            for line in f:
                if time.perf_counter()>deadline:raise TimeoutError('Chunk reread cap')
                record=json.loads(line);r=record['replicate'];assert r==lines and record['seed']==seeds[r]
                if r>=len(data):assert not complete;break
                accepted+=audit_trace(state,data[r],record['trace'],job['method']);lines+=1
    assert lines==len(data)
    if complete:assert accepted==meta['accepted_events_audited']
    return dict(complete=bool(complete),paths=len(data),accepted_events=accepted,job=job),data
def score(tables,references,required=R_FORMAL,all_chunks_complete=True):
    comparisons=[];done=all_chunks_complete
    for c in range(4):
        states,prob=references[c];S=len(states)
        for m in range(2):
            table=tables.get((c,m));n=table['n'] if table else 0
            if n!=required:done=False;continue
            frequency=table['counts']/n
            for t in range(2):
                observed=float(np.sum(abs(frequency[t]-prob[t]))/2);epsilon=math.sqrt((S*math.log(2)+math.log(16/.001))/(2*required));alternative=float(np.sqrt(np.log((2.**S)*16/.001)/(2*required)));assert abs(epsilon-alternative)<1e-12;comparisons.append(dict(cell=c,method=m,time=float(TIMES[t+1]),paths=n,states=S,observed_TV=observed,epsilon=epsilon,upper=observed+epsilon,margin=.10,supported=observed+epsilon<=.10))
    done=bool(done and len(comparisons)==16);return dict(complete=done,qualified=bool(done and all(x['supported'] for x in comparisons)),comparisons=comparisons)
def evaluate(directory,registration,references_dir,output,deadline=math.inf):
    directory=Path(directory);references={};lookups={};tables={};chunks=[];manifest=[]
    for c in range(4):
        with np.load(Path(references_dir)/f'c{c}_reference.npz') as z:states=z['states'];prob=z['probabilities'];np.testing.assert_array_equal(z['times'],TIMES[1:])
        references[c]=(states,prob);lookups[c]={tuple(x):i for i,x in enumerate(states)}
    for job in registration['jobs']:
        result,data=audit_chunk(job,registration['seeds'][job['id']],directory,deadline);chunks.append(result)
        if data is None:continue
        c,m=job['cell'],job['method'];table=tables.setdefault((c,m),dict(n=0,counts=np.zeros((2,len(references[c][0])),np.int64)));table['n']+=len(data)
        for t in (0,1):
            idx=np.array([lookups[c][tuple(pop[t+1].ravel())] for pop in data]);primary=np.bincount(idx,minlength=len(references[c][0]));ordered=np.sort(idx);alternate=np.searchsorted(ordered,np.arange(len(primary)),side='right')-np.searchsorted(ordered,np.arange(len(primary)),side='left');np.testing.assert_array_equal(primary,alternate);table['counts'][t]+=primary
        for suffix in ('.npz','_traces.jsonl'):
            p=directory/(job['id']+suffix)
            if p.exists():manifest.append(dict(file=p.relative_to(ROOT).as_posix(),bytes=p.stat().st_size,sha256=hashlib.sha256(p.read_bytes()).hexdigest()))
    result=score(tables,references,all_chunks_complete=all(c['complete'] for c in chunks));result.update(audit_passed=True,chunks=chunks,paths=sum(x['paths'] for x in chunks),accepted_events=sum(x.get('accepted_events',0) for x in chunks),verdict='QUALIFIED WITHIN REGISTERED TINY MIXED-LAW ENDPOINT SCOPE' if result['qualified'] else 'NOT QUALIFIED OR INCOMPLETE UNDER REGISTERED CRITERIA',shared='same fixtures and literal-role audit as timing/sampling, saved reread and separately assembled count reference/alternate categorical frequencies/bounds; same author not external replication');atomic_json(output,result);atomic_json(Path(output).with_name('raw_manifest.json'),dict(records=manifest,files=len(manifest),bytes=sum(x['bytes'] for x in manifest)));return result

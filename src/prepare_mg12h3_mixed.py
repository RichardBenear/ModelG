"""MG124 prospectively seeded actual-writer timing and saved pipeline preparation."""
from mg12h3_mixed import *
from check_mg12h3_mixed import audit_chunk,score
from scipy.sparse.linalg import expm_multiply
from scipy.linalg import expm
from dataclasses import asdict
import shutil
PREP=BASE/'preparation'
def plan():
    if (PREP/'plan.json').exists():raise ValueError('Already planned')
    from prepare_mg11g4 import collect_seeds
    prior=set();scan=subprocess.run(['rg','-i','-l','-g','*.json','"[^"\\n]*seed[^"\\n]*"','results'],cwd=ROOT,capture_output=True,text=True)
    if scan.returncode not in (0,1):raise ValueError(scan.stderr)
    for name in scan.stdout.splitlines():
        if not name.replace('\\','/').startswith('results/MG12H3/mixed_qualification/'):collect_seeds(read(ROOT/name),prior)
    formal=canonical_jobs();timing=[dict(id=f'timing_c{c}_m{m}',cell=c,method=m,chunk=0,reps=128) for c in range(4) for m in range(2)];seeds={};used=set(prior);collisions=[]
    for group,stage in ((formal,12400),(timing,12401)):
        for job in group:
            values=[]
            for r in range(job['reps']):
                index=job['chunk']*CHUNK+r;nonce=0;value=seed_for(12,job['cell'],job['method'],stage,index,nonce)
                while value in used:collisions.append(dict(id=job['id'],replicate=r,nonce=nonce,rejected=value));nonce+=1;value=seed_for(12,job['cell'],job['method'],stage,index,nonce)
                values.append(value);used.add(value)
            seeds[job['id']]=values
    atomic_json(PREP/'plan.json',dict(stage='MG124',scope='90second preparation only; formal approval pending',source_hashes=mixed_hashes(),jobs=formal,timing_jobs=timing,seeds=seeds,fixtures=[asdict(fixture(c)) for c in range(4)],collision_resolutions=collisions,prior_seed_values=len(prior),ceiling_elapsed_seconds=90,formal_candidate_total_charge=10800,parent_reserve=300))
def build_references(directory,deadline):
    directory=Path(directory);directory.mkdir(parents=True,exist_ok=True);refs={};checks=[]
    for c in range(4):
        if time.perf_counter()>deadline:raise TimeoutError('Reference deadline')
        state=fixture(c);states,_,q=master(state.counts(),*config(state),.1,max_states=8192);p0=np.zeros(len(states));p0[0]=1;prob=np.array([expm_multiply(q.T*t,p0) for t in TIMES[1:]]);dense=np.array([expm(q.toarray().T*t)@p0 for t in TIMES[1:]]);np.testing.assert_allclose(prob,dense,rtol=0,atol=1e-11);assert prob.min()>=-1e-12 and abs(prob.sum(axis=1)-1).max()<=1e-10 and float(abs(q.sum(axis=1)).max())<=1e-10
        with np.load(ROOT/f'results/MG12H3/bridge/c{c}_reference.npz') as old:np.testing.assert_array_equal(np.array(states),old['states']);np.testing.assert_allclose(prob[-1],old['probabilities'],rtol=0,atol=1e-11)
        np.savez_compressed(directory/f'c{c}_reference.npz',states=np.array(states),probabilities=prob,times=TIMES[1:]);refs[c]=(np.array(states),prob);checks.append(dict(cell=c,states=len(states),dense_error=float(abs(prob-dense).max())))
    return refs,checks
def prepare():
    p=read(PREP/'plan.json')
    if p['source_hashes']!=mixed_hashes() or (PREP/'attempt.json').exists():raise ValueError('Changed source or existing attempt')
    subprocess.run(['git','ls-files','--error-unmatch','results/MG12H3/mixed_qualification/preparation/plan.json'],cwd=ROOT,check=True,stdout=subprocess.DEVNULL);start=time.perf_counter();deadline=start+90;atomic_json(PREP/'attempt.json',dict(status='RUNNING'));refs,checks=build_references(PREP/'references',deadline);records=[];forecast={}
    for job in p['timing_jobs']:
        meta=generate(job,p['seeds'][job['id']],PREP/'timing',deadline);checked,data=audit_chunk(job,p['seeds'][job['id']],PREP/'timing',deadline);assert checked['complete'];records.append(meta)
        for formal in [x for x in p['jobs'] if x['cell']==job['cell'] and x['method']==job['method']]:forecast[formal['id']]=2*meta['elapsed_seconds']*(CHUNK/128)+2
    # Algebra-only complete tables; not integer samples or production outcomes.
    tables={(c,m):dict(n=R_FORMAL,counts=refs[c][1]*R_FORMAL) for c in range(4) for m in range(2)};ideal=score(tables,refs);assert ideal['complete'] and ideal['qualified'];atomic_json(PREP/'synthetic_score.json',dict(label='ALGEBRA ONLY: reference probabilities supplied as frequencies, not sampled data',result=ideal));tables[0,0]['n']-=1;incomplete=score(tables,refs);assert not incomplete['complete'] and not incomplete['qualified']
    job=p['timing_jobs'][0];source=PREP/'timing';bad=PREP/'bad_seed';bad.mkdir(exist_ok=True)
    for suffix in ('.npz','.json','_traces.jsonl'):shutil.copyfile(source/(job['id']+suffix),bad/(job['id']+suffix))
    path=bad/(job['id']+'.npz')
    with np.load(path) as z:payload={key:z[key] for key in z.files};payload['seeds']=payload['seeds'].copy();payload['seeds'][0]=np.uint32(int(payload['seeds'][0])^1)
    with path.open('wb') as f:np.savez_compressed(f,**payload)
    meta=read(path.with_suffix('.json'));meta['sha256']=hashlib.sha256(path.read_bytes()).hexdigest();atomic_json(path.with_suffix('.json'),meta)
    try:audit_chunk(job,p['seeds'][job['id']],bad,deadline)
    except AssertionError:pass
    else:raise AssertionError('Wrong raw seed accepted')
    for label in ('bad_method','bad_trace','prefix'):
        target=PREP/label;target.mkdir(exist_ok=True)
        for suffix in ('.npz','.json','_traces.jsonl'):shutil.copyfile(source/(job['id']+suffix),target/(job['id']+suffix))
        rawpath=target/(job['id']+'.npz')
        if label=='bad_method':
            with np.load(rawpath) as z:payload={key:z[key] for key in z.files}
            changed=dict(job,method=1);payload['job_json']=np.array(json.dumps(changed))
            with rawpath.open('wb') as f:np.savez_compressed(f,**payload)
            meta=read(rawpath.with_suffix('.json'));meta['sha256']=hashlib.sha256(rawpath.read_bytes()).hexdigest();atomic_json(rawpath.with_suffix('.json'),meta)
        elif label=='bad_trace':
            tracepath=target/(job['id']+'_traces.jsonl');lines=tracepath.read_text(encoding='utf-8').splitlines();changed=False
            for i,line in enumerate(lines):
                rec=json.loads(line)
                for row in rec['trace']:
                    if row['accepted']:
                        row['after'][0][2]='Z' if row['after'][0][2]!='Z' else 'W';changed=True;break
                if changed:lines[i]=json.dumps(rec);break
            assert changed;tracepath.write_text('\n'.join(lines)+'\n',encoding='utf-8')
        else:
            with np.load(rawpath) as z:subset=z['counts'][:64]
            write_counts(rawpath,subset,p['seeds'][job['id']][:64],job);meta=read(rawpath.with_suffix('.json'));meta['status']='PREFIX';atomic_json(rawpath.with_suffix('.json'),meta);tracepath=target/(job['id']+'_traces.jsonl');tracepath.write_text('\n'.join(tracepath.read_text(encoding='utf-8').splitlines()[:64])+'\n',encoding='utf-8')
        if label=='prefix':prefix,_=audit_chunk(job,p['seeds'][job['id']],target,deadline);assert not prefix['complete'] and prefix['paths']==64
        else:
            try:audit_chunk(job,p['seeds'][job['id']],target,deadline)
            except AssertionError:pass
            else:raise AssertionError('Corrupted method/trace accepted')
    missing,_=audit_chunk(job,p['seeds'][job['id']],PREP/'empty',deadline);assert not missing['complete']
    elapsed=time.perf_counter()-start;assert elapsed<90 and mixed_hashes()==p['source_hashes'];atomic_json(PREP/'results.json',dict(passed=True,source_hashes=mixed_hashes(),plan_sha256=hashlib.sha256((PREP/'plan.json').read_bytes()).hexdigest(),elapsed_seconds=elapsed,paths=1024,reference_checks=checks,timing=records,forecast=forecast,forecast_worker_seconds=sum(forecast.values()),forecast_ideal_8worker_seconds=sum(forecast.values())/8,synthetic_complete_score=True,incomplete_suppressed=True,wrong_seed_rejected=True,wrong_method_rejected=True,corrupt_trace_rejected=True,prefix_not_complete=True,missing_rejected=True));atomic_json(PREP/'attempt.json',dict(status='COMPLETE'))
if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--plan',action='store_true');parser.add_argument('--prepare',action='store_true');a=parser.parse_args()
    if a.plan:plan()
    elif a.prepare:
        try:prepare()
        except BaseException:atomic_json(PREP/'failure.json',dict(error=traceback.format_exc()));raise
    else:parser.error('Choose action')

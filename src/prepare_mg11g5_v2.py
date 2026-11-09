"""MG100 prospective short development; source/seed plan precedes every draw."""
from mg11g5 import *
import sys,json,hashlib,subprocess,argparse,traceback,platform,shutil
from mg11g4 import simulate
from mg11g4_cases import fixture
from mg11g4_io import save,check_raw,digest
from mg11g3_oracle import accounting
from mg11g2_pattern_v2 import reference
from mg11g2_pattern_v2_oracle import dense_projection
from mg10d1_pool import atomic_json
from evaluate_mg11g5 import evaluate
from check_mg11g5 import audit,quantile
BASE=ROOT/'results/MG11G5/preparation_v2'
FILES=('src/mg11g5.py','src/run_mg11g5.py','src/evaluate_mg11g5.py','src/check_mg11g5.py','src/prepare_mg11g5.py','src/prepare_mg11g5_v2.py','spec/MG11G5_pilot.md','src/mg11g4.py','src/mg11g4_io.py','src/mg11g4_cases.py','src/mg11g3.py','src/mg11g3_oracle.py','src/mg11g2.py','src/mg11g2_pattern_v2.py','src/mg11g2_pattern_v2_oracle.py','src/mg10s1.py','src/mg11n1.py','src/mg10p1.py','src/mg10d2_source.py','src/runner.py','src/mg10d1_pool.py','src/mg10e2c1_pool.py','results/MG003/ref_profiles.npz')
def read(p):return json.loads(Path(p).read_text(encoding='utf-8'))
def hashes():
    names=list(FILES)+sorted({donor_path(parent_job(c)).relative_to(ROOT).as_posix() for c in contexts()})
    return {name:digest(ROOT/name) if name.endswith('.npz') else hashlib.sha256((ROOT/name).read_text(encoding='utf-8').replace('\r\n','\n').encode()).hexdigest() for name in names}
def development():
    out=[dict(id=f'identity_c{c}_v{v}',kind='identity',cell=c,setting=v,seeds=[seed_for(11,c,v,10000,r) for r in range(8)]) for c in range(3) for v in range(2)]
    out += [dict(c,id=f'timing_o{c["origin"]}_s{c["start"]}_a{c["arm"]}_v{v}',kind='timing',setting=v,seed=seed_for(11,4*c['origin']+2*c['start']+c['arm'],v,10001,0),warmup_seed=seed_for(11,4*c['origin']+2*c['start']+c['arm'],v,10001,0,1)) for c in contexts() for v in range(2)]
    return out
def plan():
    if (BASE/'plan.json').exists():raise ValueError('Already planned')
    from prepare_mg11g4 import collect_seeds
    values=[j[key] for j in jobs()+development() for key in ('seed','warmup_seed') if key in j]+[s for j in development() for s in j.get('seeds',[])];assert len(values)==len(set(values));prior=set()
    scan=subprocess.run(['rg','-i','-l','-g','*.json','"[^"\\n]*seed[^"\\n]*"','results'],cwd=ROOT,text=True,capture_output=True)
    if scan.returncode not in (0,1):raise ValueError(scan.stderr)
    for name in scan.stdout.splitlines():
        if not name.replace('\\','/').startswith('results/MG11G5/'):collect_seeds(read(ROOT/name),prior)
    assert not set(values)&prior
    atomic_json(BASE/'plan.json',dict(stage='MG100',scope='short development only; full pilot unregistered',source_hashes=hashes(),jobs=jobs(),development=development(),seed_count=len(values),prior_seed_values=len(prior),collisions=0,ceiling_elapsed_seconds=890,threads_per_library=1,workers=1,platform=platform.platform(),python=sys.version))
def algebra():
    # Known interaction, arbitrary background and missing block: no physical draws.
    records={(v,b,o,s,a):100+o+3*a+s*(2+.04*a)+.003*b*s*a for v in range(2) for b in range(4) for o in range(2) for s in range(2) for a in range(2)}
    z=summarize(records);np.testing.assert_allclose(z['settings'][0]['mean'],.0445,atol=1e-12);assert z['direction_consistent']
    from scipy.special import stdtr
    assert abs(stdtr(3,student.ppf(.9875,3))-.9875)<1e-12
    assert abs(stdtr(3,quantile(.9875,3))-.9875)<1e-12
    del records[0,0,0,0,0];z=summarize(records);assert z['settings'][0]['n']==3 and z['settings'][0]['inference'] is None and not z['direction_consistent']
    for bad in ([],[1,2,3],[1,2,3,float('nan')]):
        try:interval(bad)
        except ValueError:pass
        else:raise AssertionError('Invalid primary accepted')
    phi=reference()[0];pop=np.zeros((61,7,500));pop[:,:3]=phi[None];np.testing.assert_allclose(projection(pop,OBS),dense_projection(pop,phi),atol=1e-12)
    return dict(background_cancellation=True,known_interaction=True,missing_block_suppressed=True,dense_projection=True,CDF_quantile=True)
def pipeline(deadline):
    # Constant synthetic raw ledgers test actual evaluator/checker; never scientific data.
    directory=BASE/'synthetic';directory.mkdir(parents=True,exist_ok=True)
    for j in jobs():
        if time.perf_counter()>deadline:raise TimeoutError('Preparation deadline')
        p=directory/(j['id']+'.npz');case,ga,gb,k=factory(j)
        if j['kind']=='stochastic':
            data=np.zeros((1,61,20,500),np.int64);data[0,:,:7]=case['initial'];env=np.zeros((1,61,2,1500),np.int64);stats=np.zeros((1,3));save(p,(case,ga,gb,k,OBS,np.ones(7)),j['seed'],(data,env,stats),1,1,TAUS[j['setting']])
        else:
            state=np.tile(initial(concentration_case(case),ga),(61,1))
            with p.open('wb') as f:np.savez_compressed(f,state=state,times=OBS,job_json=json.dumps(j))
        atomic_json(p.with_suffix('.json'),dict(status='COMPLETE',job=j,sha256=digest(p),SYNTHETIC_NOT_PHYSICAL=True))
    result=evaluate(directory,BASE/'synthetic_complete.json');audit(directory,BASE/'synthetic_complete.json',BASE/'synthetic_complete_check.json');assert result['complete']
    # Preserve original file under separate name, then expose a valid incomplete prefix.
    j=next(j for j in jobs() if j['kind']=='stochastic');p=directory/(j['id']+'.npz');p.rename(p.with_name(p.stem+'_complete_reference.npz'))
    with np.load(p.with_name(p.stem+'_complete_reference.npz')) as raw:save(p,(factory(j)[0],factory(j)[1],factory(j)[2],factory(j)[3],OBS[:5],np.ones(7)),j['seed'],(raw['data'][:,:5],raw['environment'][:,:5],raw['stats']),1,1,TAUS[j['setting']])
    atomic_json(p.with_suffix('.json'),dict(status='PREFIX',job=j,SYNTHETIC_NOT_PHYSICAL=True));result=evaluate(directory,BASE/'synthetic_prefix.json');audit(directory,BASE/'synthetic_prefix.json',BASE/'synthetic_prefix_check.json');assert not result['complete'] and result['settings'][0]['n']==3
    p.rename(p.with_name(p.stem+'_prefix_reference.npz'));result=evaluate(directory,BASE/'synthetic_missing.json');assert result['settings'][0]['n']==3
    empty=BASE/'empty';empty.mkdir(exist_ok=True);result=evaluate(empty,BASE/'synthetic_empty.json');audit(empty,BASE/'synthetic_empty.json',BASE/'synthetic_empty_check.json');assert result['settings'][0]['n']==0
    return dict(complete=True,prefix=True,missing=True,empty=True)
def qualify():
    plan=read(BASE/'plan.json')
    if hashes()!=plan['source_hashes'] or (BASE/'attempt.json').exists():raise ValueError('Changed source or already attempted')
    subprocess.run(['git','ls-files','--error-unmatch','results/MG11G5/preparation_v2/plan.json'],cwd=ROOT,check=True,stdout=subprocess.DEVNULL)
    start=time.perf_counter();deadline=start+890;atomic_json(BASE/'attempt.json',dict(status='RUNNING',source_hashes=hashes()));tests=algebra();checks=[];timings={};detcost={}
    for j in plan['development']:
        if time.perf_counter()>deadline:raise TimeoutError('Preparation deadline')
        v=j['setting'];tau=TAUS[v]
        if j['kind']=='identity':
            case,ga,gb,k,_,mask=fixture(j['cell']);times=np.arange(5)/1000
            expected=[[],[],[]];actual=[[],[],[]]
            for r in range(8):
                # Each path gets a distinct held-replay seed, avoiding public seed=-1.
                seed=j['seeds'][r]
                old=simulate(case,ga,gb,k,times,1,seed,tau);new,_=streamed(case,ga,gb,k,times,seed,tau,deadline,chunk=2)
                for a,b in zip(old,new):np.testing.assert_array_equal(a,b)
                for values,b in zip(expected,old):values.append(b)
                for values,b in zip(actual,new):values.append(b)
            for suffix,z in (('monolithic',tuple(np.concatenate(x) for x in expected)),('segmented',tuple(np.concatenate(x) for x in actual))):
                path=BASE/(j['id']+'_'+suffix+'.npz');save(path,(case,ga,gb,k,times,mask),j['seeds'][0],z,8,1,tau);checks.append(check_raw(path))
        else:
            case,ga,gb,k=factory(j);warmtimes=np.array([0.,.001]);warm,_=streamed(case,ga,gb,k,warmtimes,j['warmup_seed'],tau,deadline);p=BASE/(j['id']+'_warm.npz');save(p,(case,ga,gb,k,warmtimes,np.ones(7)),j['warmup_seed'],warm,1,1,tau);checks.append(check_raw(p));times=np.array([0.,.005,.01]);p=BASE/(j['id']+'.npz')
            def checkpoint(t,z):save(p,(case,ga,gb,k,t,np.ones(7)),j['seed'],z,1,1,tau)
            result,timing=streamed(case,ga,gb,k,times,j['seed'],tau,deadline,checkpoint,chunk=1);clock=time.perf_counter();checks.append(check_raw(p));timing['audit_seconds']=time.perf_counter()-clock
            old=simulate(case,ga,gb,k,times,1,j['seed'],tau)
            for a,b in zip(old,result):np.testing.assert_array_equal(a,b)
            ref=p.with_name(p.stem+'_held_reference.npz');save(ref,(case,ga,gb,k,times,np.ones(7)),j['seed'],old,1,1,tau);checks.append(check_raw(ref));timings[j['id']]=timing
    for j in [j for j in jobs() if j['kind']=='deterministic']:
        clock=time.perf_counter();state,case,ga,gb,k=deterministic(j,np.array([0.,.005,.01]),deadline);detcost[j['id']]=time.perf_counter()-clock;check=accounting(state,case,ga,gb);assert check['passed'];p=BASE/(j['id']+'.npz')
        with p.open('wb') as f:np.savez_compressed(f,state=state,times=np.array([0.,.005,.01]),job_json=json.dumps(j))
    pipeline_check=pipeline(deadline);forecast={}
    for j in jobs():
        if j['kind']=='deterministic':forecast[j['id']]=2*1500*detcost[j['id']]+10
        else:
            q=timings[f'timing_o{j["origin"]}_s{j["start"]}_a{j["arm"]}_v{j["setting"]}'];forecast[j['id']]=1.5*1500*q['compute_seconds']+200*q['io_seconds']+42*q['audit_seconds']+30
    elapsed=time.perf_counter()-start;assert elapsed<890 and hashes()==plan['source_hashes'];result=dict(passed=True,source_hashes=hashes(),tests=tests,pipeline=pipeline_check,integer_checks=checks,timings=timings,deterministic_short_seconds=detcost,forecast=forecast,forecast_worker_seconds=sum(forecast.values()),forecast_ideal_8worker_seconds=sum(forecast.values())/8,elapsed_seconds=elapsed,storage_free_bytes=shutil.disk_usage(ROOT).free,limitations='short-horizon extrapolation; no late outcomes; no full pattern pilot executed');atomic_json(BASE/'results.json',result);atomic_json(BASE/'attempt.json',dict(status='COMPLETE',elapsed_seconds=elapsed))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--plan',action='store_true');p.add_argument('--qualify',action='store_true');a=p.parse_args()
    if a.plan:plan()
    elif a.qualify:
        try:qualify()
        except BaseException:atomic_json(BASE/'failure.json',dict(error=traceback.format_exc(),source_hashes=hashes()));raise
    else:p.error('Choose action')

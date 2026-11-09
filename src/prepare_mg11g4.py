"""MG096 prospective stochastic-coupling development; no formal production CLI."""
import os
for key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[key]='1'
from pathlib import Path
import argparse,time,json,hashlib,subprocess,traceback,math,copy,sys,platform,shutil
import numpy as np,scipy,numba
from scipy.stats import binom
from mg11g4 import simulate,arrays,single_steps,MAX_COUNT
from mg11g4_cases import fixture,one_step,benchmark,TAUS
from mg11g4_io import sample,save,generate,check_raw,jobs,digest
from evaluate_mg11g4 import main as evaluate,cdf_distance
from check_mg11g4 import main as audit
from mg11n1 import graph,simulate as old_network
from mg10s1 import fixture as old_fixture,ssa_chunk as old_exact,leap_chunk as old_leap
from runner import seed_for
from mg10d1_pool import atomic_json
ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/'results/MG11G4/preparation'
FILES=('src/mg11g4.py','src/mg11g4_cases.py','src/mg11g4_oracle.py','src/mg11g4_io.py','src/evaluate_mg11g4.py','src/check_mg11g4.py','src/run_mg11g4.py','src/prepare_mg11g4.py','spec/MG11G4_step1.md','src/mg10s1.py','src/mg10d2_source.py','src/mg11n1.py','src/mg11n1_oracle.py','src/mg11g2_pattern_v2.py','src/mg11g2.py','src/mg10p1.py','src/runner.py','src/mg10d1_pool.py','src/mg10e2c1_pool.py','results/MG003/ref_profiles.npz')
def read(p):return json.loads(Path(p).read_text(encoding='utf-8'))
def hashes():return {name:digest(ROOT/name) if name.endswith('.npz') else hashlib.sha256((ROOT/name).read_text(encoding='utf-8').replace('\r\n','\n').encode()).hexdigest() for name in FILES}
def dev_jobs():
    out=[]
    for c in range(4):
        for e in (0,1):out.append(dict(id=f'legacy_c{c}_e{e}',kind='legacy',cell=c,engine=e,reps=8,seed=seed_for(11,c,e,9600,0),warmup_seed=seed_for(11,c,e,9600,0,1)))
    for c in (1,2):out.append(dict(id=f'common_c{c}',kind='common',cell=c,engine=0,reps=8,seed=seed_for(11,10+c,0,9600,0),warmup_seed=seed_for(11,10+c,0,9600,0,1)))
    for large in (0,1):out.append(dict(id=f'single_large{large}',kind='single',cell=large,engine=1,reps=4096,seed=seed_for(11,20+large,1,9600,0),warmup_seed=seed_for(11,20+large,1,9600,0,1)))
    for c in range(3):
        for e in range(5):out.append(dict(id=f'timing_c{c}_e{e}',kind='timing',cell=c,engine=e,reps=128,seed=seed_for(11,c,e,9601,0),warmup_seed=seed_for(11,c,e,9601,0,1)))
    for cross in (0,1):
        for v,tau in enumerate((1/24000,1/48000)):out.append(dict(id=f'benchmark_g{cross}_v{v}',kind='benchmark',cell=cross,engine=1,reps=2,tau=tau,seed=seed_for(11,100+cross,v,9602,0),warmup_seed=seed_for(11,100+cross,v,9602,0,1)))
    return out
def collect_seeds(value,found,context=False):
    if isinstance(value,dict):
        for key,child in value.items():collect_seeds(child,found,context or 'seed' in key.lower())
    elif isinstance(value,list):
        for child in value:collect_seeds(child,found,context)
    elif context and isinstance(value,int):found.add(value)
def plan():
    if (BASE/'plan.json').exists():raise ValueError('Already planned')
    candidates=jobs();dev=dev_jobs();values=[j[key] for j in candidates+dev for key in ('seed','warmup_seed')];assert len(set(values))==len(values);prior=set()
    scan=subprocess.run(['rg','-i','-l','-g','*.json','"[^"\\n]*seed[^"\\n]*"','results'],cwd=ROOT,text=True,capture_output=True)
    if scan.returncode not in (0,1):raise ValueError(scan.stderr)
    for name in scan.stdout.splitlines():
        if name.replace('\\','/').startswith('results/MG11G4/'):continue
        collect_seeds(read(ROOT/name),prior)
    assert not set(values)&prior
    atomic_json(BASE/'plan.json',dict(stage='MG096',scope='BUILD/DEVELOPMENT/SHORT BENCHMARK ONLY; formal120job qualification approval pending',pi_request='yes, go ahead with step 1',source_hashes=hashes(),jobs=candidates,development_jobs=dev,seed_count=len(values),prior_seed_values=len(prior),collisions=0,held_replays='Reuse timing_c1 engine0/3/4 seeds deliberately for32path/block16 identity; no new streams',ceiling_elapsed_seconds=900,threads_per_library=1,workers=1,formal_candidate_cap_worker_seconds=7200,environment=dict(python=sys.version,platform=platform.platform(),numpy=np.__version__,scipy=scipy.__version__,numba=numba.__version__)))
def invalid_inputs():
    inputs=fixture(1);case,ga,gb,k,t,m=inputs;count=0
    def reject(call):
        nonlocal count
        try:call()
        except (ValueError,TypeError):count+=1;return
        raise AssertionError('Invalid input accepted')
    for bad in (-1,2**32,1.5):reject(lambda bad=bad:simulate(case,ga,gb,k,t,1,bad,.001,m))
    for bad in ([-1.,0.],[float('nan'),0.]):reject(lambda bad=bad:simulate(case,ga,gb,bad,t,1,1,.001,m))
    huge=dict(case,initial=case['initial'].copy());huge['initial'][3,0]=MAX_COUNT+1;reject(lambda:simulate(huge,ga,gb,k,t,1,1,.001,m))
    negative=dict(case,initial=case['initial'].copy());negative['initial'][4,0]=-1;reject(lambda:simulate(negative,ga,gb,k,t,1,1,.001,m));reject(lambda:simulate(case,ga,graph(3,[]),k,t,1,1,.001,m));reject(lambda:simulate(case,ga,gb,k,[0,0],1,1,.001,m));reject(lambda:simulate(case,ga,gb,k,t,1,1,-.01,m))
    bad_graph=dict(gb,weight=gb['weight'].copy());bad_graph['weight'][0]*=1.1;reject(lambda:simulate(case,ga,bad_graph,k,t,1,1,.001,m))
    self_graph=dict(gb,dst=gb['dst'].copy());self_graph['dst'][0]=gb['src'][0];reject(lambda:simulate(case,ga,self_graph,k,t,1,1,.001,m))
    # Conditional first moment identities; no draws or fitted tolerance.
    for x,b,a,source_hazard,hop_hazard in ((4.,1.,8.,.3,.7),(1.,8.,3.,.2,.5),(1000.,100000000.,100000000.,.003,1.)):
        rate=29*x*b/40;dt=min(.001,.1/(1+rate/min(x,b)));mu=rate*dt;assert 0<=mu<min(x,b);np.testing.assert_allclose((x-mu)*(.9*dt)/(1-mu/x),x*.9*dt,rtol=1e-12);np.testing.assert_allclose((b-mu)*(hop_hazard*dt)/(1-mu/b),b*hop_hazard*dt,rtol=1e-12);p=source_hazard*dt;np.testing.assert_allclose((a-a*p)*(hop_hazard*dt)/(1-p),a*hop_hazard*dt,rtol=1e-12)
    return dict(invalid_cases_rejected=count,conditional_mean_algebra=True)
def save_pair(directory,job,inputs,warm,main,tau=None):
    directory.mkdir(parents=True,exist_ok=True);save(directory/(job['id']+'_warm.npz'),inputs,job['warmup_seed'],warm,1,job['engine'],tau);save(directory/(job['id']+'.npz'),inputs,job['seed'],main,job['reps'],job['engine'],tau);return [check_raw(directory/(job['id']+suffix)) for suffix in ('.npz','_warm.npz')]
def formal(dev,deadline):
    directory=BASE/'formal';checks=[];tests=[invalid_inputs()];binomial_result=None
    for job in [j for j in dev if j['kind'] in ('legacy','common','single')]:
        if time.perf_counter()>deadline:raise TimeoutError('Development deadline')
        c,e=job['cell'],job['engine']
        if job['kind']=='legacy':
            old=old_fixture(c);case=dict(initial=old['state'],a0=old['a0'],b0=old['b0'],omega=old['omega'],diff=old['diff']);g=graph(case['initial'].shape[1],[]);inputs=(case,g,g,np.zeros(2),np.array([0.,.01,.05]),np.ones(7));warm=sample(inputs,e,1,job['warmup_seed']);main=sample(inputs,e,8,job['seed']);fn=old_exact if e==0 else old_leap;args=(case['initial'],case['a0'],case['b0'],case['omega'],case['diff'],inputs[-1]);reference=fn(*args,inputs[-2],8,job['seed']) if e==0 else fn(*args,.001,inputs[-2],8,job['seed']);np.testing.assert_array_equal(main[0],reference[0]);np.testing.assert_array_equal(main[2],reference[1]);assert not np.any(main[1]);save(directory/(job['id']+'_reference.npz'),inputs,job['seed'],(reference[0],np.zeros_like(main[1]),reference[1]),8,e)
        elif job['kind']=='common':
            case,ga,gb,k,t,m=fixture(c);inputs=(case,ga,ga,k,t,m);warm=sample(inputs,0,1,job['warmup_seed']);main=sample(inputs,0,8,job['seed']);reference=old_network(case,ga,k,t,8,job['seed'],m)
            for aa,bb in zip(main,reference):np.testing.assert_array_equal(aa,bb)
            save(directory/(job['id']+'_reference.npz'),inputs,job['seed'],reference,8,0)
        else:
            inputs=one_step(bool(c));case,ga,gb,k,t,m=inputs;src,dst,w,degree,offsets,ids=arrays(ga,gb,k)
            args=(case['initial'],case['a0'],case['b0'],case['omega'],case['diff'],m,float(t[-1]));tail=(src,dst,w,k,degree,offsets,ids);warm=single_steps(*args,1,job['warmup_seed'],*tail);main=single_steps(*args,job['reps'],job['seed'],*tail);data,ef,stats=main
            if c==0:
                assert not np.any(data[:,1,7,1]) and not np.any(data[:,1,10,1]);assert not np.any(ef[:,1,:,ga['src']==1]);assert np.all(data[:,1,7,0]+ef[:,1,0].sum(axis=-1)<=1) and np.all(data[:,1,10,0]+ef[:,1,1].sum(axis=-1)<=1)
            else:
                x=data[:,1,3,1];grid=np.unique(np.r_[x,x-1]);empirical=np.searchsorted(np.sort(x),grid,side='right')/len(x);distance=float(abs(empirical-binom.cdf(grid,100000000,1e-5)).max());epsilon=math.sqrt(math.log(4000)/(2*4096));upper=distance+epsilon;binomial_result=dict(distance=distance,epsilon=epsilon,upper=upper,margin=.08,alpha=.0005,passed=upper<=.08);assert upper<=.08
        checks+=save_pair(directory,job,inputs,warm,main,float(inputs[-2][-1]) if job['kind']=='single' else None);tests.append(dict(id=job['id'],passed=True));print('Formal '+job['id'],flush=True)
    return dict(passed=True,tests=tests,raw_checks=checks,large_binomial=binomial_result)
def archive(base,label):
    for name in ('evaluation.json','independent_check.json','MG11G4_results.md'):(base/name).rename(base/(label+'_'+name))
def pipeline(deadline):
    base=BASE/'SYNTHETIC_pipeline';directory=base/'production';directory.mkdir(parents=True);atomic_json(base/'registration.json',dict(jobs=jobs(),formal_passed=True,scope='ALGEBRAIC FIXTURES ONLY'))
    for job in jobs():
        if time.perf_counter()>deadline:raise TimeoutError('Development deadline')
        inputs=fixture(job['cell']);case,ga,gb,k,t,m=inputs;data=np.zeros((job['reps'],len(t),20,ga['n']),np.int64);data[:,:,:7]=case['initial'];env=np.zeros((job['reps'],len(t),2,len(ga['src'])),np.int64);stats=np.zeros((job['reps'],3));p=directory/(job['id']+'.npz');save(p,inputs,job['seed'],(data,env,stats),job['reps'],job['engine']);save(p.with_name(p.stem+'_warm.npz'),inputs,job['warmup_seed'],(data[:1],env[:1],stats[:1]),1,job['engine']);atomic_json(p.with_suffix('.json'),dict(status='COMPLETE',job=job,sha256=digest(p),scope='ALGEBRAIC ONLY'))
    e=evaluate(base);audit(base,synthetic=True);assert e['scientific_pass'] and len(e['formal'])==196;archive(base,'complete')
    # CDF identity at known unequal supports, independent of raw simulation.
    assert cdf_distance(np.array([0,0,1,2]),np.array([0,1,1,2]))==.25
    wrong=copy.deepcopy(e);wrong['formal'][0]['distance']+=.01;atomic_json(base/'evaluation.json',wrong)
    try:audit(base,synthetic=True)
    except AssertionError:pass
    else:raise AssertionError('Wrong CDF accepted')
    (base/'evaluation.json').rename(base/'wrong_statistic_evaluation.json');atomic_json(base/'wrong_statistic_rejected.json',dict(rejected=True))
    p=directory/'c1_e3_r7.npz';backup=p.with_suffix('.npz.withheld');p.rename(backup);evaluate(base);audit(base,synthetic=True);v=read(base/'evaluation.json');assert not v['complete'] and all(x['passed'] is None for x in v['formal']);archive(base,'missing')
    with np.load(backup) as raw:partial={key:raw[key] for key in raw.files}
    for key in ('data','environment','stats'):partial[key]=partial[key][:128]
    partial['reps_completed']=128
    with p.open('wb') as f:np.savez_compressed(f,**partial)
    evaluate(base);audit(base,synthetic=True);v=read(base/'evaluation.json');assert v['paths_available']==243840 and not v['complete'];archive(base,'prefix');p.rename(p.with_suffix('.npz.prefix_fixture'));backup.rename(p)
    with np.load(p) as raw:invalid={key:raw[key] for key in raw.files}
    invalid['data'][:,1:,3,0]-=1;invalid['data'][:,1:,0,0]+=1;bad=base/'invalid_accounting.npz'
    with bad.open('wb') as f:np.savez_compressed(f,**invalid)
    try:check_raw(bad)
    except AssertionError:pass
    else:raise AssertionError('Invalid local ledger accepted')
    atomic_json(base/'invalid_accounting_rejected.json',dict(rejected=True,total_conserved_but_reaction_ledger_missing=True))
    empty=BASE/'SYNTHETIC_empty';empty.mkdir();atomic_json(empty/'registration.json',dict(jobs=jobs(),formal_passed=True));evaluate(empty);audit(empty,synthetic=True);assert read(empty/'evaluation.json')['paths_available']==0
    atomic_json(BASE/'pipeline_qualification.json',dict(passed=True,source_hashes=hashes(),complete_missing_prefix_empty_wrong_statistic_invalid_ledger=True,formal_comparisons=196,scope='Manufactured zero-event balanced fields, not scientific samples'))
def execute():
    start=time.perf_counter();deadline=start+900;p=read(BASE/'plan.json')
    if hashes()!=p['source_hashes'] or jobs()!=p['jobs']:raise ValueError('Source/plan changed')
    subprocess.run(['git','ls-files','--error-unmatch','results/MG11G4/preparation/plan.json'],cwd=ROOT,check=True,stdout=subprocess.DEVNULL)
    if (BASE/'attempt.json').exists():raise ValueError('Already attempted; no retry')
    atomic_json(BASE/'attempt.json',dict(status='RUNNING DEVELOPMENT'));f=formal(p['development_jobs'],deadline);atomic_json(BASE/'formal.json',f);timings={};rawchecks=list(f['raw_checks'])
    for job in [x for x in p['development_jobs'] if x['kind']=='timing']:
        if time.perf_counter()>deadline:raise TimeoutError('Development deadline')
        timings[job['cell'],job['engine']]=generate(job,BASE/'timing');rawchecks+=[check_raw(BASE/'timing'/(job['id']+suffix)) for suffix in ('.npz','_warm.npz')];print('Timed '+job['id'],flush=True)
    held=[]
    for engine in (0,3,4):
        job=dict(timings[1,engine]['job'],id=f'held_c1_e{engine}',reps=32);generate(job,BASE/'held_replays',block=16)
        with np.load(BASE/'held_replays'/(job['id']+'.npz')) as aa,np.load(BASE/'timing'/(timings[1,engine]['job']['id']+'.npz')) as bb:
            for name in ('data','environment','stats'):np.testing.assert_array_equal(aa[name],bb[name][:32])
        held.append(dict(engine=engine,exact_prefix=True,scope='Held development stream, not new science'))
    bench=[]
    for job in [x for x in p['development_jobs'] if x['kind']=='benchmark']:
        if time.perf_counter()>deadline:raise TimeoutError('Development deadline')
        inputs=benchmark(bool(job['cell']));case,ga,gb,k,t,m=inputs;clock=time.perf_counter();warm=simulate(case,ga,gb,k,t,1,job['warmup_seed'],job['tau'],m);warm_seconds=time.perf_counter()-clock;clock=time.perf_counter();main=simulate(case,ga,gb,k,t,2,job['seed'],job['tau'],m);compute=time.perf_counter()-clock;checks=save_pair(BASE/'benchmark',job,inputs,warm,main,job['tau']);rawchecks+=checks;item=dict(job=job,warm_seconds=warm_seconds,main_compute_seconds=compute,seconds_per_path=compute/2,linear_T25_seconds_per_path=compute/2*12500,scope='Microbenchmark only: no long-horizon cost or accuracy guarantee');bench.append(item);atomic_json(BASE/'benchmark'/(job['id']+'.json'),item);print('Bench '+job['id'],flush=True)
    pipeline(deadline);forecast={job['id']:timings[job['cell'],job['engine']]['compute_seconds']*16*2+timings[job['cell'],job['engine']]['io_seconds']*256*2+10 for job in jobs()};total=sum(forecast.values());elapsed=time.perf_counter()-start;assert elapsed<900 and shutil.disk_usage(ROOT).free>=2*1024**3
    result=dict(passed=True,status='IMPLEMENTATION DEVELOPMENT PASS; formal statistical qualification not run',source_hashes=hashes(),elapsed_seconds=elapsed,formal=f,raw_checks=rawchecks,held_replays=held,benchmarks=bench,forecast_seconds_per_job=forecast,forecast_worker_seconds=total,forecast_core_hours=total/3600,eight_worker_wall_minutes=total/480,production_within_reserve=total<=5400,candidate_cap_core_hours=2,scope='Development identities, arrival/binomial probes,128path/method/case timings, and.002 N500benchmark; no R16384formal endpoint qualification or pattern pilot')
    atomic_json(BASE/'preparation_results.json',result);atomic_json(BASE/'attempt.json',dict(status='COMPLETE DEVELOPMENT',elapsed_seconds=elapsed));print(json.dumps({key:result[key] for key in ('passed','elapsed_seconds','forecast_core_hours','eight_worker_wall_minutes','production_within_reserve')},indent=2))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--plan',action='store_true');p.add_argument('--execute',action='store_true');a=p.parse_args()
    if a.plan:plan()
    elif a.execute:
        try:execute()
        except BaseException:atomic_json(BASE/'failure.json',dict(status='FAILED; preserve/version, no automatic retry',error=traceback.format_exc()));raise
    else:p.error('Choose plan/execute')

"""MG080 approved N1 ensemble production, atomic128-path prefixes, no retry."""
import os
for key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[key]='1'
from pathlib import Path
import argparse,json,hashlib,sys,platform,shutil,subprocess,time,traceback
from datetime import datetime,timezone
from mg10d1_pool import atomic_json
from mg10e2c1_pool import run
ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/'results/MG11N1'
EXTRA=('src/run_mg11n1.py','src/mg11n1_stream.py','src/evaluate_mg11n1.py','src/check_mg11n1.py','src/qualify_mg11n1_pipeline.py','src/mg10e2c1_pool.py','spec/MG11N1_prerun.md')
def read(path):return json.loads(Path(path).read_text(encoding='utf-8'))
def digest(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def utc():return datetime.now(timezone.utc).isoformat()
def hashes():
    original=read(BASE/'preparation/preparation_plan.json')['source_hashes'];values={}
    for name in set(original)|set(EXTRA):
        p=ROOT/name;content=p.read_bytes() if name.endswith('.npz') else p.read_text(encoding='utf-8').replace('\r\n','\n').encode();values[name]=hashlib.sha256(content).hexdigest()
        if name in original and values[name]!=original[name]:raise ValueError('Frozen preparation source changed:'+name)
    return values
def raw_save(path,case,g,kappa,times,mask,seed,data,environment,stats,requested):
    import numpy as np
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True);temp=path.with_suffix('.npz.tmp')
    with temp.open('wb') as f:np.savez_compressed(f,data=data,environment=environment,stats=stats,initial=case['initial'],a0=case['a0'],b0=case['b0'],omega=case['omega'],diff=case['diff'],src=g['src'],dst=g['dst'],weight=g['weight'],kappa=kappa,times=times,mask=mask,seed=seed,requested_reps=requested,reps_completed=len(data))
    os.replace(temp,path)
def generate(job,directory,block=128):
    import numpy as np
    from mg11n1_cases import cell
    from mg11n1 import validate
    from mg11n1_stream import stream_chunk,oracle_stream
    from check_mg11n1_preparation import check_raw
    directory=Path(directory);directory.mkdir(parents=True,exist_ok=True);case,g,k,t,m=cell(job['cell']);validate(case,g,k,m);seed=job['seed'];engine=job['engine'];total=job['reps'];rng=np.random.default_rng(seed)
    warm=stream_chunk(case['initial'],case['a0'],case['b0'],case['omega'],case['diff'],m,t,1,job['warmup_seed'],g['src'],g['dst'],g['weight'],k,1000000) if engine==0 else oracle_stream(case,g,k,t,1,np.random.default_rng(job['warmup_seed']),m)
    raw_save(directory/(job['id']+'_warm.npz'),case,g,k,t,m,job['warmup_seed'],*warm,1)
    data=np.zeros((total,len(t),20,g['n']),dtype=np.int64);env=np.zeros((total,len(t),2,len(g['src'])),dtype=np.int64);stats=np.zeros((total,3));path=directory/(job['id']+'.npz');compute=0.
    raw_save(path,case,g,k,t,m,seed,data[:0],env[:0],stats[:0],total)
    for start in range(0,total,block):
        stop=min(total,start+block);clock=time.perf_counter()
        if engine==0:result=stream_chunk(case['initial'],case['a0'],case['b0'],case['omega'],case['diff'],m,t,stop-start,seed if start==0 else -1,g['src'],g['dst'],g['weight'],k,1000000)
        else:result=oracle_stream(case,g,k,t,stop-start,rng,m)
        compute+=time.perf_counter()-clock;data[start:stop],env[start:stop],stats[start:stop]=result
        raw_save(path,case,g,k,t,m,seed,data[:stop],env[:stop],stats[:stop],total)
        atomic_json(directory/(job['id']+'_progress.json'),dict(job=job,replicas_completed=stop,requested=total,updated=utc()))
    checked=check_raw(path)
    atomic_json(path.with_suffix('.json'),dict(status='COMPLETE',job=job,raw=path.name,sha256=digest(path),replicas=total,compute_seconds=compute,raw_check=checked,finished=utc()))
def worker(job,directory):
    r=read(BASE/'registration.json')
    if hashes()!=r['source_hashes'] or job not in r['jobs']:raise ValueError('Unregistered source/job')
    generate(job,directory)
def register():
    import numpy,scipy,numba
    from prepare_mg11n1 import candidate_jobs
    if (BASE/'registration.json').exists():raise ValueError('Already registered')
    q=read(BASE/'preparation/pipeline_qualification.json');prep=read(BASE/'preparation/independent_check.json');plan=read(BASE/'preparation/preparation_plan.json')
    if not q['passed'] or q['source_hashes']!=hashes() or q['elapsed_seconds']>=300 or not prep['passed'] or not prep['within_75_percent_reserve']:raise ValueError('Qualification block')
    jobs=candidate_jobs()
    if jobs!=plan['candidate_jobs'] or len(jobs)!=32 or len({j[k] for j in jobs for k in ('seed','warmup_seed')})!=64:raise ValueError('Seed/denominator matrix changed')
    if sum(prep['forecast_seconds_per_job'].values())>2700 or shutil.disk_usage(ROOT).free<2*1024**3:raise ValueError('Forecast/storage')
    atomic_json(BASE/'registration.json',dict(status='REGISTERED BEFORE N1 SCIENTIFIC QUALIFICATION',stage='MG080',registered=utc(),pi_authorization='yes, i approve of all those items (model,criteria,1core-hour cap and public full-raw/results)',publication_authorized=True,
      jobs=jobs,source_hashes=hashes(),formal_passed=True,formal_evidence={name:digest(BASE/'preparation'/name) for name in ('independent_check.json','qualification/qualification_results.json','qualification/linear_limits.json','qualification/continuum_inputs.npz')},
      forecast_seconds=prep['forecast_seconds_per_job'],ceiling_worker_seconds=3600,job_timeout_seconds=180,workers=8,threads_per_library=1,poll_seconds=.25,block_replicas=128,expected_paths=131072,
      analytical_family=dict(comparisons=6,alpha=.0005,n=32768,margin=.03),two_sample_family=dict(comparisons=28,alpha=.0005,n_per_engine=16384,margin=.08),
      scope='Fixed equal-volume symmetric network, exact tiny endpoint laws and declared linear limits only; no nonlinear particle/full-lattice stochastic/vacuum/energy/entanglement qualification',external_review='Not requested; same-author separate recomputation',environment=dict(python=sys.version,platform=platform.platform(),numpy=numpy.__version__,scipy=scipy.__version__,numba=numba.__version__)))
def production():
    r=read(BASE/'registration.json')
    if hashes()!=r['source_hashes']:raise ValueError('Source changed')
    subprocess.run(['git','ls-files','--error-unmatch','results/MG11N1/registration.json'],cwd=ROOT,check=True,stdout=subprocess.DEVNULL)
    if (BASE/'production/ledger.json').exists():raise ValueError('Already attempted; no automatic retry')
    started=utc();atomic_json(BASE/'run_state.json',dict(status='RUNNING',started=started,parent_pid=os.getpid()))
    l=run(r['jobs'],BASE/'production',Path(__file__).resolve(),ROOT,3600,workers=8,timeout=180,forecast=r['forecast_seconds'],poll_seconds=.25)
    subprocess.run([sys.executable,str(ROOT/'src/evaluate_mg11n1.py')],cwd=ROOT,check=True);subprocess.run([sys.executable,str(ROOT/'src/check_mg11n1.py')],cwd=ROOT,check=True)
    atomic_json(BASE/'run_state.json',dict(status='RESULTS AND SEPARATE CHECK READY; PI closure pending',started=started,finished=utc(),stop_reason=l['stop_reason']))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--register',action='store_true');p.add_argument('--production',action='store_true');p.add_argument('--worker',nargs=2);a=p.parse_args()
    if a.worker:
        j=json.loads(a.worker[0]);d=Path(a.worker[1])
        try:worker(j,d)
        except BaseException:atomic_json(d/(j['id']+'.json'),dict(status='FAILED',job=j,error=traceback.format_exc()));raise
    elif a.register:register()
    elif a.production:
        try:production()
        except BaseException:atomic_json(BASE/'run_state.json',dict(status='DRIVER/ANALYSIS FAILED; data retained',error=traceback.format_exc(),finished=utc()));raise
    else:p.error('Choose register/production/worker')

"""MG088 actual atomic writer and approval-gated bounded production driver."""
import os
for key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[key]='1'
from pathlib import Path
import time,json,hashlib,argparse,traceback,subprocess,sys,shutil
from datetime import datetime,timezone
import numpy as np
from mg11g2_pattern_v2 import factory,solve_pattern,OBS,jobs
from mg11g2_pattern_v2_oracle import incidence_check
from mg10d1_pool import atomic_json
from mg10e2c1_pool import run
ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/'results/MG11G2/pattern_response_v2'
def read(p):return json.loads(Path(p).read_text(encoding='utf-8'))
def digest(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def utc():return datetime.now(timezone.utc).isoformat()

def raw_save(path,job,case,g,k,mode,times,state,requested_times):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True);temp=path.with_suffix('.npz.tmp')
    with temp.open('wb') as f:np.savez_compressed(f,job_json=json.dumps(job,sort_keys=True),initial=case['initial'],a0=case['a0'],b0=case['b0'],omega=case['omega'],diff=case['diff'],src=g['src'],dst=g['dst'],weight=g['weight'],kappa=k,mode=mode,times=times,state=state,requested_times=requested_times)
    os.replace(temp,path)

def generate(job,directory,times=OBS,deadline=None,scope='SCIENTIFIC PRODUCTION'):
    directory=Path(directory);directory.mkdir(parents=True,exist_ok=True);started=time.perf_counter();deadline=deadline or started+590;case,g,k,mode=factory(job['config'],job['start']);path=directory/(job['id']+'.npz');io=0.;writes=0
    def save(t,y):
        nonlocal io,writes
        if len(t)>1 and abs(t[-1]-round(t[-1]))>1e-12 and t[-1]!=times[-1]:return
        before=time.perf_counter();raw_save(path,job,case,g,k,mode,t,y,times);atomic_json(path.with_name(path.stem+'_progress.json'),dict(job=job,checkpoints=len(t),target_checkpoints=len(times),last_time=float(t[-1]),scope=scope));io+=time.perf_counter()-before;writes+=1
    from mg11g2 import initial
    save(np.array([0.]),initial(case,g)[None]);result=solve_pattern(case,g,k,mode,times,job['setting'],deadline,checkpoint=save);check=incidence_check(result['state'],case,g);assert check['passed'];elapsed=time.perf_counter()-started
    meta=dict(status='COMPLETE',job=job,scope=scope,checkpoints=len(result['times']),elapsed_seconds=elapsed,io_seconds=io,non_io_seconds=elapsed-io,atomic_writes=writes,sha256=digest(path),homogeneous_reduction=result['homogeneous_reduction'],nfev=result['nfev'],njev=result['njev'],raw_check=check)
    atomic_json(path.with_suffix('.json'),meta);return meta

def worker(job,directory):
    from prepare_mg11g2_pattern_v2 import hashes
    r=read(BASE/'registration.json')
    if hashes()!=r['source_hashes'] or job not in r['jobs']:raise ValueError('Unregistered source/job')
    generate(job,directory)

def register(approval):
    from prepare_mg11g2_pattern_v2 import hashes
    if not approval:raise ValueError('Record concrete PI design/budget/public-retention approval')
    if (BASE/'registration.json').exists():raise ValueError('Already registered')
    q=read(BASE/'preparation/preparation_results.json');pipeline=read(BASE/'preparation/pipeline_qualification.json')
    if not q['passed'] or not pipeline['passed'] or hashes()!=q['source_hashes'] or hashes()!=pipeline['source_hashes'] or q['forecast_worker_seconds']>5400 or shutil.disk_usage(ROOT).free<2*1024**3:raise ValueError('Qualification/forecast/storage block')
    import scipy,platform
    atomic_json(BASE/'registration.json',dict(status='REGISTERED BEFORE PATTERN SCIENCE',stage='MG089',registered=utc(),pi_authorization=approval,publication_authorized=True,jobs=jobs(),source_hashes=hashes(),reference_sha256=digest(ROOT/'results/MG003/ref_profiles.npz'),preparation_sha256=digest(BASE/'preparation/preparation_results.json'),pipeline_sha256=digest(BASE/'preparation/pipeline_qualification.json'),ceiling_worker_seconds=7200,job_timeout_seconds=600,workers=8,threads_per_library=1,poll_seconds=.25,forecast_seconds=q['forecast_seconds_per_job'],scope='Deterministic full-ring matched-environment field response; no stochastic/physical identity claim',environment=dict(python=sys.version,platform=platform.platform(),numpy=np.__version__,scipy=scipy.__version__)))

def production():
    from prepare_mg11g2_pattern_v2 import hashes
    r=read(BASE/'registration.json')
    if hashes()!=r['source_hashes'] or (BASE/'production/ledger.json').exists():raise ValueError('Source changed or already attempted; no retry')
    subprocess.run(['git','ls-files','--error-unmatch','results/MG11G2/pattern_response_v2/registration.json'],cwd=ROOT,check=True,stdout=subprocess.DEVNULL)
    started=utc();atomic_json(BASE/'run_state.json',dict(status='RUNNING',started=started,parent_pid=os.getpid()))
    ledger=run(r['jobs'],BASE/'production',Path(__file__).resolve(),ROOT,7200,workers=8,timeout=600,forecast=r['forecast_seconds'],poll_seconds=.25)
    subprocess.run([sys.executable,str(ROOT/'src/evaluate_mg11g2_pattern_v2.py')],cwd=ROOT,check=True);subprocess.run([sys.executable,str(ROOT/'src/check_mg11g2_pattern_v2.py')],cwd=ROOT,check=True)
    atomic_json(BASE/'run_state.json',dict(status='RESULTS/AUDIT READY; PI closure pending',started=started,finished=utc(),stop_reason=ledger['stop_reason']))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--register',action='store_true');p.add_argument('--approval');p.add_argument('--production',action='store_true');p.add_argument('--worker',nargs=2);a=p.parse_args()
    if a.worker:
        j=json.loads(a.worker[0]);d=Path(a.worker[1])
        try:worker(j,d)
        except BaseException:atomic_json(d/(j['id']+'.json'),dict(status='FAILED',job=j,error=traceback.format_exc()));raise
    elif a.register:register(a.approval)
    elif a.production:
        try:production()
        except BaseException:atomic_json(BASE/'run_state.json',dict(status='DRIVER/ANALYSIS FAILED; data retained',finished=utc(),error=traceback.format_exc()));raise
    else:p.error('Choose action')

"""MG092 approval-gated B-only branch writer and bounded production controller."""
import os
for key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[key]='1'
from pathlib import Path
import sys,time,json,hashlib,argparse,traceback,subprocess,shutil
from datetime import datetime,timezone
import numpy as np
from mg11g3 import factory,donor_path,donor_manifest,jobs,solve,OBS,SWITCH_TIME,PARENT_INDEX
from mg11g2 import initial
from mg11g3_oracle import accounting
from mg10d1_pool import atomic_json
from mg10e2c1_pool import run
ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/'results/MG11G3'
def read(p):return json.loads(Path(p).read_text(encoding='utf-8'))
def digest(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def raw_save(path,job,case,ga,gb,k,t,state,requested,parent_sha):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True);tmp=path.with_suffix('.npz.tmp')
    with tmp.open('wb') as f:np.savez_compressed(f,job_json=json.dumps(job,sort_keys=True),initial=case['initial'],a0=case['a0'],b0=case['b0'],omega=case['omega'],diff=case['diff'],src=ga['src'],dst=ga['dst'],weight=ga['weight'],srcB=gb['src'],dstB=gb['dst'],weightB=gb['weight'],kappa=k,times=t,state=state,requested_times=requested,parent_sha256=parent_sha,parent_index=PARENT_INDEX,switch_time=SWITCH_TIME)
    os.replace(tmp,path)
def generate(job,directory,times=OBS,deadline=None,scope='SCIENTIFIC PRODUCTION'):
    start=time.perf_counter();directory=Path(directory);directory.mkdir(parents=True,exist_ok=True);case,ga,gb,k=factory(job);p=directory/(job['id']+'.npz');parent_sha=digest(donor_path(job));io=0.;writes=0
    def save(t,state):
        nonlocal io,writes
        if len(t)>1 and abs(t[-1]-round(t[-1]))>1e-12 and t[-1]!=times[-1]:return
        before=time.perf_counter();raw_save(p,job,case,ga,gb,k,t,state,times,parent_sha);atomic_json(p.with_name(p.stem+'_progress.json'),dict(job=job,checkpoints=len(t),target_checkpoints=len(times),last_post_time=float(t[-1]),scope=scope));io+=time.perf_counter()-before;writes+=1
    save(np.array([0.]),initial(case,ga)[None]);result=solve(job,times,deadline or start+290,checkpoint=save);v=accounting(result['state'],case,ga,gb);assert v['passed'];elapsed=time.perf_counter()-start;meta=dict(status='COMPLETE',job=job,scope=scope,checkpoints=len(times),elapsed_seconds=elapsed,io_seconds=io,non_io_seconds=elapsed-io,atomic_writes=writes,sha256=digest(p),parent_sha256=parent_sha,raw_check=v);atomic_json(p.with_suffix('.json'),meta);return meta
def worker(job,directory):
    from prepare_mg11g3 import hashes
    r=read(BASE/'registration.json')
    if hashes()!=r['source_hashes'] or donor_manifest()!=r['donor_hashes'] or job not in r['jobs']:raise ValueError('Unregistered source/donor/job')
    generate(job,directory)
def register(approval):
    from prepare_mg11g3 import hashes
    if not approval or (BASE/'registration.json').exists():raise ValueError('PI approval required or already registered')
    q=read(BASE/'preparation/preparation_results.json')
    if not q['passed'] or hashes()!=q['source_hashes'] or donor_manifest()!=q['donor_hashes'] or q['forecast_worker_seconds']>2700 or shutil.disk_usage(ROOT).free<2*1024**3:raise ValueError('Qualification/forecast/storage block')
    atomic_json(BASE/'registration.json',dict(stage='MG093',registered=datetime.now(timezone.utc).isoformat(),pi_authorization=approval,publication_authorized=True,jobs=jobs(),source_hashes=hashes(),donor_hashes=donor_manifest(),forecast_seconds=q['forecast_seconds_per_job'],ceiling_worker_seconds=3600,job_timeout_seconds=300,workers=8,threads_per_library=1,poll_seconds=.25,scope='Conditional deterministic B-transport intervention; no stochastic/physical identity claim'))
def production():
    from prepare_mg11g3 import hashes
    r=read(BASE/'registration.json')
    if hashes()!=r['source_hashes'] or donor_manifest()!=r['donor_hashes'] or (BASE/'production/ledger.json').exists():raise ValueError('Source/donor changed or attempted; no retry')
    subprocess.run(['git','ls-files','--error-unmatch','results/MG11G3/registration.json'],cwd=ROOT,check=True,stdout=subprocess.DEVNULL);started=datetime.now(timezone.utc).isoformat();atomic_json(BASE/'run_state.json',dict(status='RUNNING',started=started,parent_pid=os.getpid()))
    l=run(r['jobs'],BASE/'production',Path(__file__).resolve(),ROOT,3600,workers=8,timeout=300,forecast=r['forecast_seconds'],poll_seconds=.25)
    for name in ('evaluate_mg11g3.py','check_mg11g3.py'):subprocess.run([sys.executable,str(ROOT/'src'/name)],cwd=ROOT,check=True)
    atomic_json(BASE/'run_state.json',dict(status='RESULTS/AUDIT READY; PI closure pending',started=started,finished=datetime.now(timezone.utc).isoformat(),stop_reason=l['stop_reason']))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--register',action='store_true');p.add_argument('--approval');p.add_argument('--production',action='store_true');p.add_argument('--worker',nargs=2);a=p.parse_args()
    if a.worker:
        j=json.loads(a.worker[0]);d=Path(a.worker[1])
        try:worker(j,d)
        except BaseException:atomic_json(d/(j['id']+'.json'),dict(status='FAILED',job=j,error=traceback.format_exc()));raise
    elif a.register:register(a.approval)
    elif a.production:
        try:production()
        except BaseException:atomic_json(BASE/'run_state.json',dict(status='DRIVER/ANALYSIS FAILED; data retained',error=traceback.format_exc()));raise
    else:p.error('Choose action')

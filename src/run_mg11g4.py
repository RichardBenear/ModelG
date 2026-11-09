"""MG096 approval-gated formal small-system stochastic qualification driver."""
import os
for key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[key]='1'
from pathlib import Path
import sys,json,time,argparse,traceback,subprocess,shutil,platform
from datetime import datetime,timezone
import numpy as np,scipy,numba
from mg11g4_io import generate,jobs,digest
from mg10d1_pool import atomic_json
from mg10e2c1_pool import run
ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/'results/MG11G4'
def read(p):return json.loads(Path(p).read_text(encoding='utf-8'))
def worker(job,directory):
    from prepare_mg11g4 import hashes
    r=read(BASE/'registration.json')
    if hashes()!=r['source_hashes'] or job not in r['jobs']:raise ValueError('Unregistered source/job')
    generate(job,directory)
def register(approval):
    from prepare_mg11g4 import hashes
    if not approval or (BASE/'registration.json').exists():raise ValueError('PI approval required or already registered')
    q=read(BASE/'preparation/preparation_results.json');p=read(BASE/'preparation/pipeline_qualification.json')
    if not q['passed'] or not p['passed'] or hashes()!=q['source_hashes'] or hashes()!=p['source_hashes'] or q['forecast_worker_seconds']>5400 or shutil.disk_usage(ROOT).free<2*1024**3:raise ValueError('Qualification/forecast/storage block')
    atomic_json(BASE/'registration.json',dict(stage='MG097',registered=datetime.now(timezone.utc).isoformat(),pi_authorization=approval,publication_authorized=True,formal_passed=True,jobs=jobs(),source_hashes=hashes(),preparation_sha256=digest(BASE/'preparation/preparation_results.json'),pipeline_sha256=digest(BASE/'preparation/pipeline_qualification.json'),forecast_seconds=q['forecast_seconds_per_job'],ceiling_worker_seconds=7200,job_timeout_seconds=180,workers=8,threads_per_library=1,poll_seconds=.25,formal_comparisons=196,alpha=.0005,R=16384,margin=.08,environment=dict(python=sys.version,platform=platform.platform(),numpy=np.__version__,scipy=scipy.__version__,numba=numba.__version__),scope='Finest tau.00025 marginal endpoint accuracy and exact-oracle validation on three small fixtures; no pattern pilot or full temporal law'))
def production():
    from prepare_mg11g4 import hashes
    r=read(BASE/'registration.json')
    if hashes()!=r['source_hashes'] or (BASE/'production/ledger.json').exists():raise ValueError('Source changed or already attempted; no retry')
    subprocess.run(['git','ls-files','--error-unmatch','results/MG11G4/registration.json'],cwd=ROOT,check=True,stdout=subprocess.DEVNULL);started=datetime.now(timezone.utc).isoformat();atomic_json(BASE/'run_state.json',dict(status='RUNNING',started=started,parent_pid=os.getpid()))
    l=run(r['jobs'],BASE/'production',Path(__file__).resolve(),ROOT,7200,workers=8,timeout=180,forecast=r['forecast_seconds'],poll_seconds=.25)
    for name in ('evaluate_mg11g4.py','check_mg11g4.py'):subprocess.run([sys.executable,str(ROOT/'src'/name)],cwd=ROOT,check=True)
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

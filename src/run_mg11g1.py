"""MG084 approved G1 ensemble production, atomic128-path prefixes, no retry."""
import os
for key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[key]='1'
from pathlib import Path
import argparse,json,hashlib,sys,platform,shutil,subprocess,time,traceback
from datetime import datetime,timezone
from mg10d1_pool import atomic_json
from mg10e2c1_pool import run
ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/'results/MG11G1'
EXTRA=('src/run_mg11g1.py','src/evaluate_mg11g1.py','src/check_mg11g1.py','src/qualify_mg11g1_pipeline.py','src/mg10e2c1_pool.py','spec/MG11G1_prerun.md')

def read(path):return json.loads(Path(path).read_text(encoding='utf-8'))
def digest(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def utc():return datetime.now(timezone.utc).isoformat()
def hashes():
    original=read(BASE/'preparation/preparation_plan.json')['source_hashes'];values={}
    for name in set(original)|set(EXTRA):
        p=ROOT/name;content=p.read_bytes() if name.endswith('.npz') else p.read_text(encoding='utf-8').replace('\r\n','\n').encode();values[name]=hashlib.sha256(content).hexdigest()
        if name in original and values[name]!=original[name]:raise ValueError('Frozen preparation source changed:'+name)
    return values
from mg11g1_sample import generate

def worker(job,directory):
    r=read(BASE/'registration.json')
    if hashes()!=r['source_hashes'] or job not in r['jobs']:raise ValueError('Unregistered source/job')
    generate(job,directory)
def register():
    import numpy,scipy,numba
    from prepare_mg11g1 import candidate_jobs
    if (BASE/'registration.json').exists():raise ValueError('Already registered')
    q=read(BASE/'preparation/pipeline_qualification.json');prep=read(BASE/'preparation/independent_check.json');plan=read(BASE/'preparation/preparation_plan.json')
    if not q['passed'] or q['source_hashes']!=hashes() or q['elapsed_seconds']>=300 or not prep['passed'] or not prep['within_75_percent_reserve']:raise ValueError('Qualification block')
    jobs=candidate_jobs()
    if jobs!=plan['candidate_jobs'] or len(jobs)!=64 or len({j[k] for j in jobs for k in ('seed','warmup_seed')})!=128:raise ValueError('Seed/denominator matrix changed')
    if sum(prep['forecast_seconds_per_job'].values())>2700 or shutil.disk_usage(ROOT).free<4*1024**3:raise ValueError('Forecast/storage')
    atomic_json(BASE/'registration.json',dict(status='REGISTERED BEFORE G1 SCIENTIFIC ESTIMATION',stage='MG084',registered=utc(),pi_authorization='i approve all of them. Go ahead witht he run. (G1 design,criteria,1core-hour cap,public full raw/results)',publication_authorized=True,
      jobs=jobs,source_hashes=hashes(),formal_passed=True,formal_evidence={name:digest(BASE/'preparation'/name) for name in ('independent_check.json','qualification/qualification_results.json','qualification/formal_references.json','qualification/reference_inputs.npz','design_calculation.json')},
      forecast_seconds=prep['forecast_seconds_per_job'],ceiling_worker_seconds=3600,job_timeout_seconds=180,workers=8,threads_per_library=1,poll_seconds=.25,block_replicas=128,expected_paths=262144,
      law_family=dict(comparisons=16,alpha=.001,n=65536,margin=.03),primary=dict(resource='B',time=4.,contrast='bipartite_minus_prism',individual_alpha=.025,joint_coverage=.95),
      scope='Degree/capacity/contact/distance-matched graph delivery probes, not particles/native vacuum/energy/entanglement',external_review='Not requested; same-author separate recomputation',environment=dict(python=sys.version,platform=platform.platform(),numpy=numpy.__version__,scipy=scipy.__version__,numba=numba.__version__)))
def production():
    r=read(BASE/'registration.json')
    if hashes()!=r['source_hashes']:raise ValueError('Source changed')
    subprocess.run(['git','ls-files','--error-unmatch','results/MG11G1/registration.json'],cwd=ROOT,check=True,stdout=subprocess.DEVNULL)
    if (BASE/'production/ledger.json').exists():raise ValueError('Already attempted; no automatic retry')
    started=utc();atomic_json(BASE/'run_state.json',dict(status='RUNNING',started=started,parent_pid=os.getpid()))
    l=run(r['jobs'],BASE/'production',Path(__file__).resolve(),ROOT,3600,workers=8,timeout=180,forecast=r['forecast_seconds'],poll_seconds=.25)
    subprocess.run([sys.executable,str(ROOT/'src/evaluate_mg11g1.py')],cwd=ROOT,check=True);subprocess.run([sys.executable,str(ROOT/'src/check_mg11g1.py')],cwd=ROOT,check=True)
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

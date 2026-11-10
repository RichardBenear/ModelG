"""MG124 compact CLI workers and approval-gated mixed-law qualification."""
from mg12h3_mixed import *
from prepare_mg12h3_mixed import build_references
from check_mg12h3_mixed import evaluate
from mg10e2c1_pool import run as pool_run
from datetime import datetime,timezone
import shutil
PREP=BASE/'preparation'
def register(approval):
    p=read(PREP/'plan.json');q=read(PREP/'results.json')
    if not approval or (BASE/'registration.json').exists() or not q['passed'] or q['source_hashes']!=mixed_hashes() or q['plan_sha256']!=hashlib.sha256((PREP/'plan.json').read_bytes()).hexdigest() or p['jobs']!=canonical_jobs() or q['forecast_worker_seconds']>7875 or shutil.disk_usage(ROOT).free<2*1024**3:raise ValueError('Approval/preparation/source/plan/budget/storage block')
    seeds={j['id']:p['seeds'][j['id']] for j in p['jobs']};values=[v for group in seeds.values() for v in group];assert len(set(values))==len(values)==262144
    atomic_json(BASE/'registration.json',dict(stage='MG125',authorization=approval,registered=datetime.now(timezone.utc).isoformat(),source_hashes=mixed_hashes(),jobs=p['jobs'],seeds=seeds,fixtures=p['fixtures'],R=32768,chunk_reps=4096,formal_comparisons=16,alpha=.001,margin=.10,tick=.001,total_charge_cap=10800,worker_charge_cap=10500,parent_reserve=300,job_timeout=600,workers=8,threads_per_library=1,forecast=q['forecast'],preparation_sha256=hashlib.sha256((PREP/'results.json').read_bytes()).hexdigest(),public_retention=True))
def worker(job,directory):
    r=read(BASE/'registration.json')
    if job not in r['jobs'] or mixed_hashes()!=r['source_hashes']:raise ValueError('Unregistered job or changed source')
    generate(job,r['seeds'][job['id']],directory,time.perf_counter()+590)
def run():
    r=read(BASE/'registration.json')
    if r['source_hashes']!=mixed_hashes() or r['jobs']!=canonical_jobs() or r['preparation_sha256']!=hashlib.sha256((PREP/'results.json').read_bytes()).hexdigest() or (BASE/'run_state.json').exists():raise ValueError('Source/plan/preparation changed or already attempted')
    subprocess.run(['git','ls-files','--error-unmatch','results/MG12H3/mixed_qualification/registration.json'],cwd=ROOT,check=True,stdout=subprocess.DEVNULL);started=datetime.now(timezone.utc).isoformat();clock=time.perf_counter();atomic_json(BASE/'run_state.json',dict(status='RUNNING',started=started,parent_pid=os.getpid()));build_references(BASE/'references',time.perf_counter()+150);reference_seconds=time.perf_counter()-clock;assert reference_seconds<150
    workers=min(8,int(os.environ.get('MODELG_WORKERS','8')));assert 1<=workers<=8;ledger=pool_run(r['jobs'],BASE/'production',Path(__file__).resolve(),ROOT,10500,workers=workers,timeout=600,forecast=r['forecast'],poll_seconds=.25);clock=time.perf_counter();result=evaluate(BASE/'production',r,BASE/'references',BASE/'evaluation.json',time.perf_counter()+300-reference_seconds);parent_seconds=reference_seconds+time.perf_counter()-clock;assert parent_seconds<=300 and ledger['charged_worker_seconds']+parent_seconds<=10800 and mixed_hashes()==r['source_hashes'];atomic_json(BASE/'run_state.json',dict(status='MIXED QUALIFICATION RESULTS/AUDIT READY; PI closure pending',started=started,finished=datetime.now(timezone.utc).isoformat(),qualified=result['qualified'],charged_worker_seconds=ledger['charged_worker_seconds'],parent_numerical_seconds=parent_seconds,total_charge_seconds=ledger['charged_worker_seconds']+parent_seconds,stop_reason=ledger['stop_reason'],workers=workers,threads_per_library=1))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--register',action='store_true');p.add_argument('--approval');p.add_argument('--run',action='store_true');p.add_argument('--worker',nargs=2);a=p.parse_args()
    if a.worker:
        job=json.loads(a.worker[0]);directory=Path(a.worker[1])
        try:worker(job,directory)
        except BaseException:atomic_json(directory/(job['id']+'.json'),dict(status='FAILED',job=job,error=traceback.format_exc()));raise
    elif a.register:register(a.approval)
    elif a.run:
        try:run()
        except BaseException:atomic_json(BASE/'failure.json',dict(error=traceback.format_exc()));atomic_json(BASE/'run_state.json',dict(status='STOPPED/ANALYSIS FAILED; partial data preserved, no retry'));raise
    else:p.error('Choose action')

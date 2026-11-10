"""MG115 full deterministic bridge: concrete PI approval required before launch."""
from prepare_mg12h3_population import *
import shutil
from datetime import datetime,timezone
def register(approval):
    q=read(PREP/'results.json')
    if not approval or (BASE/'registration.json').exists() or not q['passed'] or q['source_hashes']!=hashes() or q['forecast_seconds']>2700 or shutil.disk_usage(ROOT).free<1024**3:raise ValueError('Approval/source/forecast/storage block')
    atomic_json(BASE/'registration.json',dict(stage='MG116',authorization=approval,registered=datetime.now(timezone.utc).isoformat(),source_hashes=hashes(),jobs=jobs(),elapsed_cap=3600,workers=1,threads_per_library=1,public_retention=True,preparation_sha256=hashlib.sha256((PREP/'results.json').read_bytes()).hexdigest()))
def run():
    r=read(BASE/'registration.json')
    if r['source_hashes']!=hashes() or (BASE/'run_state.json').exists():raise ValueError('Source changed or already attempted')
    subprocess.run(['git','ls-files','--error-unmatch','results/MG12H3/population/registration.json'],cwd=ROOT,check=True,stdout=subprocess.DEVNULL);start=time.perf_counter();deadline=start+3600;atomic_json(BASE/'run_state.json',dict(status='RUNNING',started=datetime.now(timezone.utc).isoformat()));records=[]
    for job in r['jobs']:
        clock=time.perf_counter();times=np.arange(int(job['duration']*4)+1)/4
        def checkpoint(t,z):save(BASE/'production'/(job['id']+'.npz'),job,t,z)
        state=solve(job,times,deadline,checkpoint);checked=check(job,state,times);assert checked['passed'];record=dict(job=job,status='COMPLETE',elapsed_seconds=time.perf_counter()-clock,check=checked);records.append(record);atomic_json(BASE/'production'/(job['id']+'.json'),record);atomic_json(BASE/'ledger.json',dict(records=records,elapsed_seconds=time.perf_counter()-start,elapsed_cap=3600))
    from check_mg12h3_population import evaluate
    result=evaluate(BASE/'production',BASE/'evaluation.json');assert result['qualified'];assert hashes()==r['source_hashes'];elapsed=time.perf_counter()-start;assert elapsed<3600;atomic_json(BASE/'run_state.json',dict(status='RESULTS/AUDIT READY; PI closure pending',finished=datetime.now(timezone.utc).isoformat(),elapsed_seconds=elapsed))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--register',action='store_true');p.add_argument('--approval');p.add_argument('--run',action='store_true');a=p.parse_args()
    if a.register:register(a.approval)
    elif a.run:
        try:run()
        except BaseException:atomic_json(BASE/'failure.json',dict(error=traceback.format_exc()));raise
    else:p.error('Choose action')

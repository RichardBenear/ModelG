"""MG120 prospective approval-gated completion, two fresh replays only."""
from mg12h3_completion import *
from check_mg12h3_completion import evaluate
from scipy.integrate import BDF
from datetime import datetime,timezone
import shutil
def register(approval):
    p=read(BASE/'plan.json');q=read(BASE/'preparation.json')
    if not approval or (BASE/'registration.json').exists() or not q['passed'] or completion_hashes()!=q['source_original_hashes'] or p['forecast_seconds']>225 or shutil.disk_usage(ROOT).free<1024**3:raise ValueError('Approval/qualification/source/budget/storage block')
    atomic_json(BASE/'registration.json',dict(stage='MG121',authorization=approval,registered=datetime.now(timezone.utc).isoformat(),source_original_hashes=completion_hashes(),entries=p['entries'],elapsed_cap=300,workers=1,threads_per_library=1,public_retention=True,preparation_sha256=hashlib.sha256((BASE/'preparation.json').read_bytes()).hexdigest()))
def run():
    r=read(BASE/'registration.json')
    if completion_hashes()!=r['source_original_hashes'] or (BASE/'run_state.json').exists():raise ValueError('Changed source/original or existing attempt')
    subprocess.run(['git','ls-files','--error-unmatch','results/MG12H3/completion/registration.json'],cwd=ROOT,check=True,stdout=subprocess.DEVNULL);start=time.perf_counter();deadline=start+300;started=datetime.now(timezone.utc).isoformat();atomic_json(BASE/'run_state.json',dict(status='RUNNING',started=started));verify_reuse(r['entries']);records=[]
    for e in [x for x in r['entries'] if x['origin']=='fresh_corrected_replay']:
        clock=time.perf_counter();job=e['job'];method=e['method'];case,ga,gb,k,mode,reference,parent=factory(job);n=ga['n'];mats=matrices(case,ga,gb);times=np.arange(101)/4;initial=np.r_[case['initial'].ravel(),np.zeros(2*n)]
        def fun(t,y):
            if time.perf_counter()>deadline:raise TimeoutError('Completion cap')
            return rhs(y,case,mats,k,mode)
        solver=BDF(fun,0.,initial,25.,jac_sparsity=sparsity(n,mats),rtol=method['rtol'],atol=method['atol'],max_step=method['max_step']);saved=[initial.copy()];index=1;path=ROOT/e['path'];save_new(path,e,times[:1],np.array(saved))
        while solver.status=='running':
            solver.step()
            if solver.status=='failed':raise RuntimeError('Completion BDF failed')
            while index<len(times) and times[index]<=solver.t:saved.append(solver.dense_output()(times[index]));index+=1;save_new(path,e,times[:index],np.array(saved))
        assert index==101;checked=check(job,np.array(saved),times);record=dict(entry=e,status='COMPLETE' if checked['passed'] else 'NUMERICAL CRITERION FAILED',check=checked,elapsed_seconds=time.perf_counter()-clock,sha256=hashlib.sha256(path.read_bytes()).hexdigest());records.append(record);atomic_json(path.with_suffix('.json'),record);atomic_json(BASE/'ledger.json',dict(records=records,elapsed_seconds=time.perf_counter()-start,elapsed_cap=300));assert checked['passed']
    result=evaluate(r['entries'],BASE/'evaluation.json');assert result['qualified'] and completion_hashes()==r['source_original_hashes'];elapsed=time.perf_counter()-start;assert elapsed<300;atomic_json(BASE/'run_state.json',dict(status='COMPLETION RESULTS/AUDIT READY; PI closure pending',started=started,finished=datetime.now(timezone.utc).isoformat(),elapsed_seconds=elapsed,elapsed_cap=300))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--register',action='store_true');p.add_argument('--approval');p.add_argument('--run',action='store_true');a=p.parse_args()
    if a.register:register(a.approval)
    elif a.run:
        try:run()
        except BaseException:atomic_json(BASE/'failure.json',dict(error=traceback.format_exc()));atomic_json(BASE/'run_state.json',dict(status='COMPLETION STOPPED; partials preserved; no retry'));raise
    else:p.error('Choose action')

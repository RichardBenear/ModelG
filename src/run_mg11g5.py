"""MG100 approval-gated conditional pilot, full raw retention and no retry."""
from mg11g5 import *
import sys,json,argparse,subprocess,traceback,shutil
from datetime import datetime,timezone
from mg11g4_io import save,check_raw,digest
from mg11g3_oracle import accounting
from mg10d1_pool import atomic_json
from mg10e2c1_pool import run
BASE=ROOT/'results/MG11G5'
def read(p):return json.loads(Path(p).read_text(encoding='utf-8'))
def worker(j,directory):
    from prepare_mg11g5_v2 import hashes
    r=read(BASE/'registration.json')
    if hashes()!=r['source_hashes'] or j not in r['jobs']:raise ValueError('Unregistered source/job')
    directory=Path(directory);start=time.perf_counter();path=directory/(j['id']+'.npz');deadline=start+3550
    if j['kind']=='deterministic':
        state,case,ga,gb,k=deterministic(j,OBS,deadline);check=accounting(state,case,ga,gb)
        if not check['passed']:raise ValueError('Deterministic accounting failure')
        tmp=path.with_suffix('.npz.tmp')
        with tmp.open('wb') as handle:np.savez_compressed(handle,state=state,times=OBS,job_json=json.dumps(j))
        os.replace(tmp,path);timing={}
    else:
        case,ga,gb,k=factory(j);mask=np.ones(7);tau=TAUS[j['setting']];warmtimes=np.array([0.,.001]);warm,_=streamed(case,ga,gb,k,warmtimes,j['warmup_seed'],tau,deadline);save(directory/(j['id']+'_warm.npz'),(case,ga,gb,k,warmtimes,mask),j['warmup_seed'],warm,1,1,tau)
        def checkpoint(t,z):
            save(path,(case,ga,gb,k,t,mask),j['seed'],z,1,1,tau)
            atomic_json(path.with_name(path.stem+'_progress.json'),dict(job=j,last_post_time=float(t[-1]),requested_final=15.,requested_frames=61,actual_frames=len(t)))
        _,timing=streamed(case,ga,gb,k,OBS,j['seed'],tau,deadline,checkpoint);check=check_raw(path);check_raw(directory/(j['id']+'_warm.npz'))
    atomic_json(path.with_suffix('.json'),dict(status='COMPLETE',job=j,elapsed_seconds=time.perf_counter()-start,sha256=digest(path),check=check,**timing))
def register(approval):
    from prepare_mg11g5_v2 import hashes
    q=read(BASE/'preparation_v2/results.json')
    if not approval or (BASE/'registration.json').exists() or not q['passed'] or hashes()!=q['source_hashes'] or q['forecast_worker_seconds']>64800 or shutil.disk_usage(ROOT).free<4*1024**3:raise ValueError('Approval/qualification/budget/storage block')
    atomic_json(BASE/'registration.json',dict(stage='MG101',pi_authorization=approval,publication_authorized=True,registered=datetime.now(timezone.utc).isoformat(),source_hashes=hashes(),jobs=jobs(),forecast=q['forecast'],ceiling=86400,timeout=3600,workers=8,preparation_sha256=digest(BASE/'preparation_v2/results.json')))
def production():
    from prepare_mg11g5_v2 import hashes
    r=read(BASE/'registration.json')
    if hashes()!=r['source_hashes'] or (BASE/'production/ledger.json').exists():raise ValueError('Changed source or already attempted')
    subprocess.run(['git','ls-files','--error-unmatch','results/MG11G5/registration.json'],cwd=ROOT,check=True,stdout=subprocess.DEVNULL)
    atomic_json(BASE/'run_state.json',dict(status='RUNNING',started=datetime.now(timezone.utc).isoformat()))
    run(r['jobs'],BASE/'production',Path(__file__).resolve(),ROOT,r['ceiling'],workers=min(8,int(os.environ.get('MODELG_WORKERS','8'))),timeout=r['timeout'],forecast=r['forecast'],poll_seconds=.25)
    for script in ('evaluate_mg11g5.py','check_mg11g5.py'):subprocess.run([sys.executable,str(ROOT/'src'/script)],cwd=ROOT,check=True)
    atomic_json(BASE/'run_state.json',dict(status='RESULTS/AUDIT READY; PI closure pending',finished=datetime.now(timezone.utc).isoformat()))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--register',action='store_true');p.add_argument('--approval');p.add_argument('--production',action='store_true');p.add_argument('--worker',nargs=2);a=p.parse_args()
    if a.worker:
        j=json.loads(a.worker[0]);directory=Path(a.worker[1])
        try:worker(j,directory)
        except BaseException:atomic_json(directory/(j['id']+'.json'),dict(status='FAILED',job=j,error=traceback.format_exc()));raise
    elif a.register:register(a.approval)
    elif a.production:production()
    else:p.error('Choose action')

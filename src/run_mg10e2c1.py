"""MG076 explicitly approved E2-C1; original data immutable, same-seed replay."""
import os
for key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[key]='1'
from pathlib import Path
import argparse,hashlib,json,shutil,subprocess,sys,time,traceback,platform
from datetime import datetime,timezone
from mg10d1_pool import atomic_json
from mg10e2c1_pool import run
ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/'results/MG10E2C1';ORIGINAL=ROOT/'results/MG10E2'
EXTRA=('src/run_mg10e2c1.py','src/evaluate_mg10e2c1.py','src/check_mg10e2c1.py','src/qualify_mg10e2c1_pipeline.py','src/mg10e2c1_pool.py','src/mg10e2c1_replay.py',
       'spec/MG10E2_continuation_amendment_1.md','spec/MG10E2C1_prerun.md')
def read(path):return json.loads(Path(path).read_text(encoding='utf-8'))
def digest(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def utc():return datetime.now(timezone.utc).isoformat()
def hashes():
    old=read(ORIGINAL/'registration.json')['source_hashes'];values={}
    for name in set(old)|set(EXTRA):
        path=ROOT/name;content=path.read_bytes() if name.endswith('.npz') else path.read_text(encoding='utf-8').replace('\r\n','\n').encode();values[name]=hashlib.sha256(content).hexdigest()
        if name in old and values[name]!=old[name]:raise ValueError('Frozen E2 source changed: '+name)
    return values
def build_manifest(original=ORIGINAL):
    original=Path(original);r=read(original/'registration.json');ledger=read(original/'branches/ledger.json');allfiles={};reused=[];modes={}
    def remember(path):
        allfiles[path.relative_to(ROOT).as_posix()]=dict(sha256=digest(path),bytes=path.stat().st_size)
    for job in r['parent_jobs']:
        path=original/'parents'/(job['id']+'.npz');meta=read(path.with_suffix('.json'))
        if meta['status']!='COMPLETE' or digest(path)!=meta['raw_sha256']:raise ValueError('Invalid reusable parent')
        remember(path);remember(path.with_suffix('.json'));reused.append(('parents',job['id']))
    for job in r['branch_jobs']:
        path=original/'branches'/(job['id']+'.npz')
        if path.exists():
            remember(path);remember(path.with_suffix('.json'));meta=read(path.with_suffix('.json'))
            if meta['status']=='COMPLETE':
                if digest(path)!=meta['raw_sha256']:raise ValueError('Invalid reusable branch')
                modes[job['id']]='REUSE';reused.append(('branches',job['id']))
            else:modes[job['id']]='EXACT_REPLAY'
        else:modes[job['id']]='UNSTARTED'
    for name in ('registration.json','parents/ledger.json','branches/ledger.json','evaluation.json','independent_check.json','MG10E2_results.md','run_state.json'):
        if (original/name).exists():remember(original/name)
    return dict(original_relative=original.relative_to(ROOT).as_posix(),original_files=allfiles,reused_files=reused,modes=modes,original_charge_seconds=ledger['charged_worker_seconds'])
def verify_original(r):
    for name,meta in r['original_files'].items():
        if digest(ROOT/name)!=meta['sha256']:raise ValueError('Original provenance changed: '+name)
def stage_reuse(base,r):
    base=Path(base);verify_original(r);source=ROOT/r['original_relative'];resolved={}
    for phase,name in r['reused_files']:
        target=base/phase;target.mkdir(parents=True,exist_ok=True)
        for suffix in ('.npz','.json'):
            dest=target/(name+suffix)
            if dest.exists():raise ValueError('Reuse target exists')
            original=source/phase/(name+suffix);shutil.copy2(original,dest)
            if digest(dest)!=digest(original):raise ValueError('Reuse copy mismatch')
        resolved[name]=dict(mode='REUSE',source=(source/phase/(name+'.npz')).relative_to(ROOT).as_posix(),sha256=digest(target/(name+'.npz')))
    for name,mode in r['modes'].items():
        if mode!='REUSE':resolved[name]=dict(mode=mode,source=(source/'branches'/(name+'.npz')).relative_to(ROOT).as_posix() if mode=='EXACT_REPLAY' else None)
    atomic_json(base/'resolved_manifest.json',resolved)
def worker(job,directory):
    import numpy as np
    import run_mg10e2 as writer
    from mg10p2 import fixture
    from mg10e2 import fork_case
    from mg10e2c1_replay import PrefixGuard
    directory=Path(directory);r=read(BASE/'registration.json')
    if hashes()!=r['source_hashes'] or job not in r['continuation_jobs']:raise ValueError('Unregistered continuation source/job')
    parent=directory.parent/'parents'/(job['parent']+'.npz');meta=read(parent.with_suffix('.json'));sha=digest(parent)
    if sha!=meta['raw_sha256']:raise ValueError('Reusable parent changed')
    source=ROOT/r['original_relative'];original_parent=source/'parents'/parent.name
    if sha!=r['original_files'][original_parent.relative_to(ROOT).as_posix()]['sha256']:raise ValueError('Parent provenance mismatch')
    case=fixture(0,job['start'])
    with np.load(parent) as raw:
        if not np.array_equal(raw['times'],np.arange(21)/4) or not np.array_equal(raw['a0'],case['a0']) or not np.array_equal(raw['b0'],case['b0']):raise ValueError('Parent grid/denominators changed')
        case=fork_case(case,raw['data'][0,-1,:7])
    mode=r['modes'][job['id']];guard=None;save=writer.save_npz;path=directory/(job['id']+'.npz')
    if mode=='EXACT_REPLAY':
        original=source/'branches'/path.name
        if digest(original)!=r['original_files'][original.relative_to(ROOT).as_posix()]['sha256']:raise ValueError('Original prefix changed')
        guard=PrefixGuard(original,save);writer.save_npz=guard
    elif mode!='UNSTARTED':raise ValueError('Completed branch must not run again')
    try:writer.trajectory(case,(1/24000,1/48000)[job['step']],np.arange(101)/4,job['seed'],job['warmup_seed'],job['arm'],5.,path,job,5400,sha)
    finally:writer.save_npz=save
    proof=guard.finish() if guard else None
    m=read(path.with_suffix('.json'));m.update(continuation_mode=mode,prefix_proof=proof,epoch='MG076 E2-C1');atomic_json(path.with_suffix('.json'),m)
def register():
    import numpy,scipy,numba
    if (BASE/'registration.json').exists():raise ValueError('Already registered')
    q=read(BASE/'preparation/pipeline_qualification.json');cost=read(ORIGINAL/'diagnostic/continuation_cost_basis.json');old=read(ORIGINAL/'registration.json')
    if not q['passed'] or q['source_hashes']!=hashes() or q['elapsed_seconds']>=300:raise ValueError('Continuation qualification/source guard')
    m=build_manifest();modes=m['modes']
    if [list(modes.values()).count(k) for k in ('REUSE','EXACT_REPLAY','UNSTARTED')]!=[1,7,120]:raise ValueError('Original denominator changed')
    jobs=[j for j in old['branch_jobs'] if modes[j['id']]!='REUSE'];forecast=cost['per_job_forecast']
    if set(forecast)!={j['id'] for j in jobs} or sum(forecast.values())+m['original_charge_seconds']>.75*518400 or shutil.disk_usage(ROOT).free<4*1024**3:raise ValueError('Forecast/storage')
    atomic_json(BASE/'registration.json',dict(status='REGISTERED CONTINUATION BEFORE SCIENCE',stage='MG076',registered=utc(),pi_authorization='yes, i approve (144core-hour aggregate continuation amendment)',publication_authorized=True,
      source_hashes=hashes(),parent_jobs=old['parent_jobs'],branch_jobs=old['branch_jobs'],continuation_jobs=jobs,forecast_seconds=forecast,ceiling_worker_seconds=518400,job_timeout_seconds=5400,
      workers=8,threads_per_library=1,poll_seconds=.25,original_registration_sha256=digest(ORIGINAL/'registration.json'),**m,
      environment=dict(python=sys.version,platform=platform.platform(),numpy=numpy.__version__,scipy=scipy.__version__,numba=numba.__version__),scope='Original E2 incomplete record preserved; same scientific matrix/endpoints, new combined continuation epoch, no outcome-based replacement',external_review='Not requested; same-author separate recomputation'))
def production():
    r=read(BASE/'registration.json')
    if hashes()!=r['source_hashes']:raise ValueError('Frozen source changed')
    subprocess.run(['git','ls-files','--error-unmatch','results/MG10E2C1/registration.json'],cwd=ROOT,check=True,stdout=subprocess.DEVNULL)
    if (BASE/'branches/ledger.json').exists() or (BASE/'resolved_manifest.json').exists():raise ValueError('Already attempted; no automatic retry')
    verify_original(r);stage_reuse(BASE,r);started=utc()
    atomic_json(BASE/'run_state.json',dict(status='RUNNING CONTINUATION',started=started,parent_pid=os.getpid(),original_charge_seconds=r['original_charge_seconds']))
    ledger=run(r['continuation_jobs'],BASE/'branches',Path(__file__).resolve(),ROOT,518400,workers=8,prior=r['original_charge_seconds'],timeout=5400,forecast=r['forecast_seconds'],poll_seconds=.25)
    atomic_json(BASE/'run_state.json',dict(status='SIMULATION STOPPED; automatic analysis follows',started=started,finished=utc(),stop_reason=ledger['stop_reason']))
    subprocess.run([sys.executable,str(ROOT/'src/evaluate_mg10e2c1.py')],cwd=ROOT,check=True)
    subprocess.run([sys.executable,str(ROOT/'src/check_mg10e2c1.py')],cwd=ROOT,check=True)
    atomic_json(BASE/'run_state.json',dict(status='RESULTS AND SEPARATE CHECK READY; PI closure pending',started=started,finished=utc(),stop_reason=ledger['stop_reason']))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--register',action='store_true');p.add_argument('--production',action='store_true');p.add_argument('--worker',nargs=2);a=p.parse_args()
    if a.worker:
        j=json.loads(a.worker[0]);d=Path(a.worker[1])
        try:worker(j,d)
        except BaseException:atomic_json(d/(j['id']+'.json'),dict(status='FAILED',job=j,error=traceback.format_exc(),finished=utc()));raise
    elif a.register:register()
    elif a.production:
        try:production()
        except BaseException:atomic_json(BASE/'run_state.json',dict(status='CONTINUATION DRIVER/ANALYSIS FAILED; originals retained',error=traceback.format_exc(),finished=utc()));raise
    else:p.error('Choose register/production/worker')

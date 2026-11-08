"""MG075 bounded synthetic admission-drain/hard-limit and prefix guard tests."""
import os
for key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[key]='1'
from pathlib import Path
import json,argparse,time
import numpy as np
from mg10d1_pool import atomic_json
from mg10e2c1_pool import run
from mg10e2c1_replay import PrefixGuard
from run_mg10e2 import save_npz
ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/'results/MG10E2/diagnostic/continuation_controls'
def main():
    started=time.perf_counter()
    if BASE.exists():raise ValueError('Already attempted; preserve/version failures')
    BASE.mkdir();script=Path(__file__).resolve()
    jobs=[dict(id='fast',delay=1.),dict(id='inflight',delay=2.),dict(id='undispatched',delay=.1)]
    l=run(jobs,BASE/'admission',script,ROOT,30,workers=2,forecast={'fast':.01,'inflight':.01,'undispatched':22.4},timeout=10,poll_seconds=.05)
    assert l['stop_reason']=='ADMISSION STOP: FORECAST; INFLIGHT DRAINED' and len(l['records'])==2 and all(v['status']=='COMPLETE' and v['returncode']==0 for v in l['records']) and l['pending_ids']==['undispatched']
    hard=run([dict(id='hard',delay=2.)],BASE/'hard',script,ROOT,.4,workers=1,timeout=5,poll_seconds=.05)
    assert hard['stop_reason']=='BUDGET-SHORT' and hard['records'][0]['returncode']!=0
    timeout=run([dict(id='timeout',delay=2.)],BASE/'timeout',script,ROOT,30,workers=1,timeout=.3,poll_seconds=.05)
    assert timeout['stop_reason']=='TIMEOUT' and timeout['records'][0]['returncode']!=0
    failed=run([dict(id='failure',delay=.1,fail=True)],BASE/'failure',script,ROOT,30,workers=1,timeout=5,poll_seconds=.05)
    assert failed['stop_reason'] and failed['records'][0]['returncode']!=0
    # Replay an already-held E2 schedule path through the actual writer guard.
    from mg10p2 import fixture
    from mg10e2 import fork_case
    from runner import seed_for
    import run_mg10e2 as original
    source=ROOT/'results/MG10E2/preparation/writer_replays/replay_a2_s1_k0.npz'
    guard=PrefixGuard(source,save_npz)
    with np.load(ROOT/'results/MG10E1/preparation/devparents/parent_s1_k0.npz') as raw:case=fork_case(fixture(0,1),raw['data'][0,-1,:7])
    original.save_npz=guard
    try:original.trajectory(case,1/24000,np.array([0.,.125,.25]),seed_for(10,1,0,7302,3),seed_for(10,1,0,7302,3,1),2,.125,BASE/'prefix_replay.npz',dict(id='held_prefix_replay'),60)
    finally:original.save_npz=save_npz
    proof=guard.finish()
    with np.load(source) as raw:arrays={key:raw[key].copy() for key in raw.files}
    corrupted={key:value.copy() for key,value in arrays.items()};corrupted['data'][0,0,0,0]+=1
    reject=PrefixGuard(source,lambda *a,**kw:None)
    try:reject(BASE/'not_written.npz',**corrupted);raise AssertionError('Corrupt prefix accepted')
    except ValueError:pass
    assert not (BASE/'not_written.npz').exists()
    elapsed=time.perf_counter()-started;assert elapsed<30
    atomic_json(BASE/'qualification.json',dict(passed=True,admission_drains_complete_jobs=True,hard_cap_terminates=True,timeout_terminates=True,failure_stops=True,exact_prefix_proof=proof,mismatch_rejected_before_save=True,elapsed_seconds=elapsed,scope='Synthetic controller workers and held-development stream; no scientific restart'))
    print(json.dumps(json.loads((BASE/'qualification.json').read_text()),indent=2))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--worker',nargs=2);a=p.parse_args()
    if a.worker:
        job=json.loads(a.worker[0]);directory=Path(a.worker[1]);time.sleep(job['delay'])
        if job.get('fail'):raise RuntimeError('Deliberate synthetic failure')
        atomic_json(directory/(job['id']+'.json'),dict(status='COMPLETE',job=job))
    else:main()

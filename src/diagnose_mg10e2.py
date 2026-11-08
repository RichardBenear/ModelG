"""MG075 disclosed held-development timing and exact-prefix replay diagnosis."""
import os
for key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[key]='1'
from pathlib import Path
import json,time,argparse,hashlib,traceback
from mg10d1_pool import atomic_json,run
ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/'results/MG10E2/diagnostic'
def digest(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def worker(job,directory):
    import numpy as np
    from mg10p2 import fixture
    from mg10e1 import fork_case,simulate as old
    from mg10e2 import simulate as new
    from runner import seed_for
    directory=Path(directory);s=job['start'];k=job['step'];tau=(1/24000,1/48000)[k]
    with np.load(ROOT/f'results/MG10E1/preparation/devparents/parent_s{s}_k{k}.npz') as raw:case=fork_case(fixture(0,s),raw['data'][0,-1,:7])
    seed=seed_for(10,s,k,7302,1);times=np.array([0.,.125,.25])
    # All-one law/held seeds: no new E2 scientific trajectory or selection.
    old(case,tau,[tau],1,seed);new(case,tau,[0.,tau],seed,0,.125)
    wall=time.perf_counter();cpu=time.process_time();a,sa=old(case,tau,times,1,seed);old_cpu=time.process_time()-cpu;old_wall=time.perf_counter()-wall
    wall=time.perf_counter();cpu=time.process_time();b,sb=new(case,tau,times,seed,0,.125);new_cpu=time.process_time()-cpu;new_wall=time.perf_counter()-wall
    np.testing.assert_array_equal(a,b);np.testing.assert_array_equal(sa,sb)
    path=directory/(job['id']+'.npz');np.savez_compressed(path,old=a,new=b,old_stats=sa,new_stats=sb,times=times,seed=seed)
    atomic_json(path.with_suffix('.json'),dict(status='COMPLETE',job=job,old_wall=old_wall,new_wall=new_wall,old_cpu=old_cpu,new_cpu=new_cpu,ratio_new_old=new_wall/old_wall,byte_identical=True,sha256=digest(path)))
def main():
    if (BASE/'serial/ledger.json').exists():raise ValueError('Already attempted; preserve/version any correction')
    serial=run([dict(id='serial_s1_k0',start=1,step=0)],BASE/'serial',Path(__file__).resolve(),ROOT,300,workers=1,timeout=60,poll_seconds=.25)
    if serial['stop_reason']:raise ValueError('Diagnostic serial failure')
    jobs=[dict(id=f'parallel_s{s}_k{k}_r{r}',start=s,step=k,rep=r) for r in range(2) for k in range(2) for s in range(2)]
    parallel=run(jobs,BASE/'parallel',Path(__file__).resolve(),ROOT,300,workers=8,prior=serial['charged_worker_seconds'],timeout=90,poll_seconds=.25)
    if parallel['stop_reason']:raise ValueError('Diagnostic parallel failure')
    records=[json.loads((BASE/'serial/serial_s1_k0.json').read_text())]+[json.loads((BASE/'parallel'/(j['id']+'.json')).read_text()) for j in jobs]
    atomic_json(BASE/'timing_diagnosis.json',dict(passed=True,records=records,charged_seconds=parallel['charged_worker_seconds'],scope='Held all-one development seed replay, exact old/new data/stats equality; wall vs process CPU measured, no production paths or new scientific RNG streams used'))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--worker',nargs=2);a=p.parse_args()
    if a.worker:
        job=json.loads(a.worker[0]);d=Path(a.worker[1])
        try:worker(job,d)
        except BaseException:atomic_json(d/(job['id']+'.json'),dict(status='FAILED',job=job,error=traceback.format_exc()));raise
    else:main()

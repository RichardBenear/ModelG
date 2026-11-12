"""MG158 prospective entry point; concrete PI budget approval gates registration."""
import os
for k in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):os.environ[k]='1'
import sys,json,hashlib,time,subprocess
from pathlib import Path
import mg12h3_coupling_long as writer
from runner import seed_for
from mg10e2c1_pool import atomic_json,run as pool_run
ROOT=writer.ROOT;BASE=ROOT/'results/MG12H3/tick_accuracy'
FILES=writer.FILES+['src/mg12h3_tick_accuracy.py','src/check_mg12h3_tick_accuracy.py','src/mg12h3_accuracy_readiness.py','spec/MG12H3_tick_accuracy_draft.md']
def read(p):return json.loads(Path(p).read_text(encoding='utf-8'))
def hashes():
    out={p:hashlib.sha256((ROOT/p).read_text(encoding='utf-8').replace('\r\n','\n').encode()).hexdigest() for p in FILES}
    p=ROOT/'results/MG12H3/ring_feasibility/o0_a0_Omega10.npz';out[p.relative_to(ROOT).as_posix()]=hashlib.sha256(p.read_bytes()).hexdigest();return out
def register():
    approval=read(BASE/'approval.json') # Created only following actual PI approval.
    assert approval['approved_worker_ceiling']==300000 and approval['approved_pairs']==[256,64]
    assert not (BASE/'plan.json').exists();prior=set()
    def collect(x):
        if isinstance(x,dict):
            for v in x.values():collect(v)
        elif isinstance(x,list):
            for v in x:collect(v)
        elif type(x) is int and 0<=x<2**32:prior.add(x)
    for p in (ROOT/'results').rglob('*.json'):
        if 'plan' in p.name or 'registration' in p.name:collect(read(p))
    used=set(prior);jobs=[];collisions=[];offsets=[0,0]
    for session,sizes in enumerate(((86,22),(85,21),(85,21)),1):
        for c,size in enumerate(sizes):
            for r in range(offsets[c],offsets[c]+size):
                nonce=0;seed=seed_for(12,c,0,15800,r,nonce)
                while seed in used:collisions.append([c,r,nonce,seed]);nonce+=1;seed=seed_for(12,c,0,15800,r,nonce)
                used.add(seed);jobs.append(dict(id=f'c{c}_paired_r{r:03d}',comparison=c,kind='paired',replicate=r,deltas=writer.DELTAS[c],sentinel=r==0,seeds=[seed],session=session))
            offsets[c]+=size
    atomic_json(BASE/'plan.json',dict(stage='MG158',approval=approval,jobs=jobs,source_hashes=hashes(),pairs=[256,64],workers=8,threads=1,session_worker_ceiling=100000,total_worker_ceiling=300000,job_seconds=1200,collisions_resolved=collisions))
def worker(job,directory):
    # Separate process: redirect only archive location and source binding, not engines.
    old_base,old_hash=writer.BASE,writer.hashes
    try:writer.BASE=BASE;writer.hashes=hashes;writer.worker(job,directory)
    finally:writer.BASE=old_base;writer.hashes=old_hash
def run(session):
    p=read(BASE/'plan.json');assert p['source_hashes']==hashes();assert session in (1,2,3);attempt=BASE/f'session_{session}_attempt.json';assert not attempt.exists()
    if session>1:assert read(BASE/f'session_{session-1}_attempt.json')['status']=='COMPLETE'
    subprocess.run(['git','ls-files','--error-unmatch',str((BASE/'plan.json').relative_to(ROOT))],cwd=ROOT,check=True,stdout=subprocess.DEVNULL)
    jobs=[j for j in p['jobs'] if j['session']==session];forecast={j['id']:1.15*(582.484122 if j['comparison']==0 else 577.835783225) for j in jobs};start=time.perf_counter();atomic_json(attempt,dict(status='RUNNING'));ledger=pool_run(jobs,BASE/f'session_{session}',Path(__file__).resolve(),ROOT,100000,workers=8,timeout=1200,forecast=forecast)
    complete=len(ledger['records'])==len(jobs) and all(r['status']=='COMPLETE' and r['returncode']==0 for r in ledger['records']);atomic_json(attempt,dict(status='COMPLETE' if complete else 'STOPPED',wall_seconds=time.perf_counter()-start,charged_worker_seconds=ledger['charged_worker_seconds'],stop_reason=ledger['stop_reason']))
if __name__=='__main__':
    if sys.argv[1]=='--worker':worker(json.loads(sys.argv[2]),sys.argv[3])
    elif sys.argv[1]=='register':register()
    elif sys.argv[1]=='run':run(int(sys.argv[2]))
    else:raise ValueError('Choose register, run SESSION, or controller worker')

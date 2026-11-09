"""MG107 registered finite joint-law qualification; preserve all paths/traces."""
import os
for key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[key]='1'
import argparse,json,time,math,hashlib,subprocess,shutil,traceback,sys
from pathlib import Path
from datetime import datetime,timezone
from dataclasses import asdict
import numpy as np
from scipy.linalg import expm
from scipy.sparse.linalg import expm_multiply
from prepare_mg12h1 import ROOT,fixture,config,hashes as engine_hashes,trace_audit,read,TIMES,STEPS
from mg12h1 import simulate
from mg12h1_reference import master
from mg10d1_pool import atomic_json
from runner import seed_for
BASE=ROOT/'results/MG12H1/qualification';R=16384;LIMIT=900
FILES=('src/run_mg12h1_qualification.py','src/check_mg12h1_qualification.py','spec/MG12H1_qualification.md','src/check_mg12h1_saved.py')
def hashes():return dict(engine_hashes(),**{name:hashlib.sha256((ROOT/name).read_text(encoding='utf-8').replace('\r\n','\n').encode()).hexdigest() for name in FILES})
def register():
    if (BASE/'registration.json').exists():raise ValueError('Already registered')
    q=read(ROOT/'results/MG12H1/development/results.json');assert q['passed'] and q['source_hashes']==engine_hashes();assert read(ROOT/'results/MG12H1/saved_evidence_check.json')['passed']
    from prepare_mg11g4 import collect_seeds
    prior=set();scan=subprocess.run(['rg','-i','-l','-g','*.json','"[^"\\n]*seed[^"\\n]*"','results'],cwd=ROOT,text=True,capture_output=True)
    if scan.returncode not in (0,1):raise ValueError(scan.stderr)
    for name in scan.stdout.splitlines():
        if not name.replace('\\','/').startswith('results/MG12H1/qualification/'):collect_seeds(read(ROOT/name),prior)
    used=set(prior);jobs=[];collisions=[]
    for c in range(4):
        for m in range(4):
            seeds=[]
            for r in range(R):
                nonce=0;seed=seed_for(12,c,m,10700,r,nonce)
                while seed in used:
                    collisions.append(dict(cell=c,method=m,replicate=r,nonce=nonce,rejected=seed));nonce+=1;seed=seed_for(12,c,m,10700,r,nonce)
                used.add(seed);seeds.append(seed)
            jobs.append(dict(id=f'c{c}_m{m}',cell=c,method=m,seeds=seeds))
    forecast=q['elapsed_seconds']*(16*R/4096)*1.5+30
    if forecast>LIMIT*.75 or shutil.disk_usage(ROOT).free<1024**3:raise ValueError('Forecast/storage block')
    atomic_json(BASE/'registration.json',dict(stage='MG107',authorization='go ahead with the formal qualification',registered=datetime.now(timezone.utc).isoformat(),source_hashes=hashes(),fixtures=[asdict(fixture(c)) for c in range(4)],jobs=jobs,R=R,steps=STEPS,times=TIMES.tolist(),formal_comparisons=16,alpha=.001,margin=.10,elapsed_cap=LIMIT,workers=1,threads_per_library=1,forecast_seconds=forecast,seed_collisions_resolved=collisions,prior_seed_values=len(prior),public_retention=True))
def save_counts(path,counts,seeds):
    tmp=path.with_suffix('.npz.tmp')
    with tmp.open('wb') as f:np.savez_compressed(f,counts=counts,seeds=np.asarray(seeds,np.uint32),times=TIMES)
    os.replace(tmp,path)
def run():
    reg=read(BASE/'registration.json')
    if hashes()!=reg['source_hashes'] or (BASE/'run_state.json').exists():raise ValueError('Source changed or already attempted')
    subprocess.run(['git','ls-files','--error-unmatch','results/MG12H1/qualification/registration.json'],cwd=ROOT,check=True,stdout=subprocess.DEVNULL);start=time.perf_counter();deadline=start+LIMIT;records=[];atomic_json(BASE/'run_state.json',dict(status='RUNNING',started=datetime.now(timezone.utc).isoformat()))
    for c in range(4):
        state=fixture(c);states,_,q=master(state.counts(),*config(state),.1);initial=np.zeros(len(states));initial[0]=1;prob=np.array([expm_multiply(q.T*t,initial) for t in TIMES[1:]]);dense=np.array([expm(q.toarray().T*t)@initial for t in TIMES[1:]]);np.testing.assert_allclose(prob,dense,rtol=0,atol=1e-11);assert prob.min()>=-1e-12 and np.max(abs(prob.sum(axis=1)-1))<=1e-10
        if c==3:
            analytic=np.array([[np.exp(-t),2*(np.exp(-t)-np.exp(-1.5*t)),1-np.exp(-t)-2*(np.exp(-t)-np.exp(-1.5*t))] for t in TIMES[1:]]);np.testing.assert_allclose(prob,analytic,rtol=0,atol=1e-11)
        np.savez_compressed(BASE/f'c{c}_reference.npz',states=np.array(states),probabilities=prob,times=TIMES[1:]);atomic_json(BASE/f'c{c}_reference.json',dict(states=len(states),dense_max_error=float(np.max(abs(prob-dense))),max_row_sum=float(abs(q.sum(axis=1)).max())))
    for j in reg['jobs']:
        state=fixture(j['cell']);tau=STEPS[j['method']];counts=np.zeros((R,len(TIMES),7,len(state.sites)),np.int64);clock=time.perf_counter();path=BASE/(j['id']+'.npz');save_counts(path,counts[:0],[]);accepted=0
        with (BASE/(j['id']+'_traces.jsonl')).open('x',encoding='utf-8') as f:
            for r,seed in enumerate(j['seeds']):
                if time.perf_counter()>deadline:raise TimeoutError('Qualification aggregate cap')
                result=simulate(state,TIMES,seed,tau,deadline=deadline);accepted+=trace_audit(state,result,tau);counts[r]=result['counts'];f.write(json.dumps(dict(replicate=r,seed=seed,trace=result['trace']),separators=(',',':'))+'\n')
                if (r+1)%512==0:
                    f.flush();save_counts(path,counts[:r+1],j['seeds'][:r+1]);atomic_json(BASE/(j['id']+'_progress.json'),dict(job_id=j['id'],completed=r+1,requested=R))
        meta=dict(status='COMPLETE',id=j['id'],replicas=R,elapsed_seconds=time.perf_counter()-clock,accepted_events_audited=accepted,sha256=hashlib.sha256(path.read_bytes()).hexdigest());records.append(meta);atomic_json(BASE/(j['id']+'.json'),meta);atomic_json(BASE/'ledger.json',dict(records=records,elapsed_seconds=time.perf_counter()-start,elapsed_cap=LIMIT))
    from check_mg12h1_qualification import check
    check(deadline);elapsed=time.perf_counter()-start
    if elapsed>=LIMIT:raise TimeoutError('Qualification analysis cap')
    atomic_json(BASE/'run_state.json',dict(status='RESULTS/AUDIT READY; PI closure pending',finished=datetime.now(timezone.utc).isoformat(),elapsed_seconds=elapsed,elapsed_cap=LIMIT,workers=1,threads_per_library=1));atomic_json(BASE/'ledger.json',dict(records=records,elapsed_seconds=elapsed,elapsed_cap=LIMIT))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--register',action='store_true');p.add_argument('--run',action='store_true');a=p.parse_args()
    if a.register:register()
    elif a.run:
        try:run()
        except BaseException:atomic_json(BASE/'failure.json',dict(error=traceback.format_exc()));atomic_json(BASE/'run_state.json',dict(status='STOPPED/FAILED; partial records retained',finished=datetime.now(timezone.utc).isoformat()));raise
    else:p.error('Choose action')

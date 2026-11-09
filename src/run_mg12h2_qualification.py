"""MG111 prospectively registered fixed-cell qualification, no retries."""
import os
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[key]='1'
import argparse,time,math,json,hashlib,subprocess,shutil,traceback
from pathlib import Path
from dataclasses import asdict
from datetime import datetime,timezone
import numpy as np
from scipy.linalg import eigh,expm
from prepare_mg12h2 import ROOT,read,fixture,reference,matrices,audit,hashes as engine_hashes,TIMES
from mg12h2 import simulate
from mg10d1_pool import atomic_json
from runner import seed_for
BASE=ROOT/'results/MG12H2/qualification';R=8192;LIMIT=900
FILES=('src/run_mg12h2_qualification.py','src/check_mg12h2_qualification.py','spec/MG12H2_qualification.md','src/check_mg12h2_saved.py')
def hashes():return dict(engine_hashes(),**{n:hashlib.sha256((ROOT/n).read_text(encoding='utf-8').replace('\r\n','\n').encode()).hexdigest() for n in FILES})
def register():
    if (BASE/'registration.json').exists():raise ValueError('Already registered')
    q=read(ROOT/'results/MG12H2/development/results.json');assert q['passed'] and q['source_hashes']==engine_hashes() and read(ROOT/'results/MG12H2/saved_evidence_check.json')['passed']
    from prepare_mg11g4 import collect_seeds
    prior=set();scan=subprocess.run(['rg','-i','-l','-g','*.json','"[^"\\n]*seed[^"\\n]*"','results'],cwd=ROOT,capture_output=True,text=True)
    if scan.returncode not in (0,1):raise ValueError(scan.stderr)
    for name in scan.stdout.splitlines():
        if not name.replace('\\','/').startswith('results/MG12H2/qualification/'):collect_seeds(read(ROOT/name),prior)
    used=set(prior);jobs=[];collisions=[]
    for c in range(6):
        for m in range(2):
            seeds=[]
            for r in range(R):
                nonce=0;seed=seed_for(12,c,m,11100,r,nonce)
                while seed in used:collisions.append(dict(cell=c,method=m,replicate=r,nonce=nonce,rejected=seed));nonce+=1;seed=seed_for(12,c,m,11100,r,nonce)
                seeds.append(seed);used.add(seed)
            jobs.append(dict(id=f'c{c}_m{m}',cell=c,method=m,seeds=seeds))
    forecast=q['elapsed_seconds']*32*1.5+30
    if forecast>.75*LIMIT or shutil.disk_usage(ROOT).free<1024**3:raise ValueError('Forecast/storage block')
    atomic_json(BASE/'registration.json',dict(stage='MG111',authorization='go ahead with H2',registered=datetime.now(timezone.utc).isoformat(),source_hashes=hashes(),jobs=jobs,fixtures=[asdict(fixture(c)[0]) for c in range(6)],R=R,times=TIMES.tolist(),tick=.005,comparisons=24,margin=.10,alpha=.001,elapsed_cap=LIMIT,workers=1,threads_per_library=1,forecast_seconds=forecast,prior_seed_values=len(prior),collision_resolutions=collisions,public_retention=True))
def save(path,data,seeds):
    tmp=path.with_suffix('.npz.tmp')
    with tmp.open('wb') as f:np.savez_compressed(f,counts=data,seeds=np.asarray(seeds,np.uint32),times=TIMES)
    os.replace(tmp,path)
def run():
    p=read(BASE/'registration.json')
    if hashes()!=p['source_hashes'] or (BASE/'run_state.json').exists():raise ValueError('Source changed or existing attempt')
    subprocess.run(['git','ls-files','--error-unmatch','results/MG12H2/qualification/registration.json'],cwd=ROOT,check=True,stdout=subprocess.DEVNULL);start=time.perf_counter();deadline=start+LIMIT;started=datetime.now(timezone.utc).isoformat();atomic_json(BASE/'run_state.json',dict(status='RUNNING',started=started));records=[]
    for c in range(6):
        state,coords,d=fixture(c);states,prob,q=reference(state);assert abs(q.sum(axis=1)).max()<1e-12 and abs(q.sum(axis=0)).max()<1e-12 and abs(q-q.T).max()<1e-12;initial=next(i for i,v in enumerate(states) if np.array_equal(v,state.counts()));w,v=eigh(q);alternate=np.array([((v*np.exp(w*t))@v.T)[initial] for t in TIMES[1:]]);np.testing.assert_allclose(prob,alternate,atol=1e-11,rtol=0)
        analytic=None
        if c==0:analytic=np.array([[.25*(1+2*np.exp(-2*t)+np.exp(-4*t)),.25*(1-np.exp(-4*t)),.25*(1-2*np.exp(-2*t)+np.exp(-4*t)),.25*(1-np.exp(-4*t))] for t in TIMES[1:]])
        elif c==1:analytic=np.array([[1/3+2/3*np.exp(-1.5*t),(1-np.exp(-1.5*t))/3,(1-np.exp(-1.5*t))/3] for t in TIMES[1:]])
        elif c in (3,4):analytic=np.array([[math.prod(1/3+((1 if x==0 else 0)-1/3)*np.exp(-3*t) for x in coord) for coord in coords] for t in TIMES[1:]])
        elif c==5:
            mats=matrices(state);analytic=np.array([np.kron(expm(mats['A']*t)[0],expm(mats['B']*t)[3]) for t in TIMES[1:]])
        if analytic is not None:np.testing.assert_allclose(prob,analytic,atol=1e-11,rtol=0)
        np.savez_compressed(BASE/f'c{c}_reference.npz',states=states,probabilities=prob,generator=q,times=TIMES[1:]);atomic_json(BASE/f'c{c}_reference.json',dict(states=len(states),eigen_max_error=float(abs(prob-alternate).max()),analytic_max_error=float(abs(prob-analytic).max()) if analytic is not None else None))
    for job in p['jobs']:
        state=fixture(job['cell'])[0];tau=None if job['method']==0 else .005;data=np.zeros((R,3,7,len(state.sites)),np.int64);clock=time.perf_counter();path=BASE/(job['id']+'.npz');save(path,data[:0],[]);accepted=0
        with (BASE/(job['id']+'_traces.jsonl')).open('x',encoding='utf-8') as f:
            for r,seed in enumerate(job['seeds']):
                if time.perf_counter()>deadline:raise TimeoutError('Qualification cap')
                result=simulate(state,TIMES,seed,tau,deadline);accepted+=audit(state,result,job['method']);data[r]=result['counts'];f.write(json.dumps(dict(replicate=r,seed=seed,trace=result['trace']),separators=(',',':'))+'\n')
                if (r+1)%512==0:f.flush();save(path,data[:r+1],job['seeds'][:r+1]);atomic_json(BASE/(job['id']+'_progress.json'),dict(completed=r+1,requested=R))
        meta=dict(id=job['id'],status='COMPLETE',replicas=R,elapsed_seconds=time.perf_counter()-clock,accepted_events_audited=accepted,sha256=hashlib.sha256(path.read_bytes()).hexdigest());records.append(meta);atomic_json(path.with_suffix('.json'),meta);atomic_json(BASE/'ledger.json',dict(records=records,elapsed_seconds=time.perf_counter()-start,elapsed_cap=LIMIT))
    from check_mg12h2_qualification import check
    check(deadline);elapsed=time.perf_counter()-start
    if elapsed>=LIMIT:raise TimeoutError('Audit cap')
    atomic_json(BASE/'run_state.json',dict(status='RESULTS/AUDIT READY; PI closure pending',started=started,finished=datetime.now(timezone.utc).isoformat(),elapsed_seconds=elapsed,elapsed_cap=LIMIT));atomic_json(BASE/'ledger.json',dict(records=records,elapsed_seconds=elapsed,elapsed_cap=LIMIT))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--register',action='store_true');p.add_argument('--run',action='store_true');a=p.parse_args()
    if a.register:register()
    elif a.run:
        try:run()
        except BaseException:atomic_json(BASE/'failure.json',dict(error=traceback.format_exc()));atomic_json(BASE/'run_state.json',dict(status='STOPPED/FAILED; partial records retained'));raise
    else:p.error('Choose action')

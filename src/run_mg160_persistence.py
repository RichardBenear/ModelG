"""MG160 prospective P1 writer/controller; no run without explicit approval."""
import os
for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[k]='1'
import json,sys,time,hashlib,traceback,subprocess,platform
from pathlib import Path
import numpy as np
from mg12h3_pattern_persistence import Coherence,homogenize
from mg12h3_finer_precision import setup
from mg12h3_tick_coupling import pair
from mg12h3_coupling_long import save_trace
from runner import seed_for
from mg10e2c1_pool import atomic_json,run as pool_run
ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/'results/MG12H3/pattern_persistence'
FILES=['src/run_mg160_persistence.py','src/mg12h3_pattern_persistence.py','src/prepare_mg160_persistence.py','src/mg12h3_sparse.py','src/mg12h3_tick_coupling.py','src/mg12h3_coupling_long.py','src/mg12h3_finer_precision.py','src/mg12h3_contrast_estimator.py','src/measure.py','src/runner.py','spec/MG12H3_pattern_persistence_plan.md']
def read(p):return json.loads(Path(p).read_text(encoding='utf-8'))
def hashes():
    out={p:hashlib.sha256((ROOT/p).read_text(encoding='utf-8').replace('\r\n','\n').encode()).hexdigest() for p in FILES}
    for p in (BASE/'prepared_inputs.npz',ROOT/'results/MG12H3/ring_feasibility/o0_a0_Omega10.npz'):out[p.relative_to(ROOT).as_posix()]=hashlib.sha256(p.read_bytes()).hexdigest()
    return out
def register():
    approval=read(BASE/'P1_approval.json');assert approval['worker_cap']==100000 and approval['pairs']==8 and not (BASE/'P1_plan.json').exists();used=set()
    def collect(x):
        if isinstance(x,dict):
            for v in x.values():collect(v)
        elif isinstance(x,list):
            for v in x:collect(v)
        elif type(x) is int and 0<=x<2**32:used.add(x)
    for p in (ROOT/'results').rglob('*.json'):
        if 'plan' in p.name or 'registration' in p.name:collect(read(p))
    jobs=[]
    for condition in ('prepared','homogeneous'):
        for r in range(4):
            nonce=0;c=int(condition=='homogeneous');seed=seed_for(12,c,0,16000,r,nonce)
            while seed in used:nonce+=1;seed=seed_for(12,c,0,16000,r,nonce)
            used.add(seed);jobs.append(dict(id=f'{condition}_r{r}',condition=condition,replicate=r,seed=seed,seed_nonce=nonce,sentinel=r==0))
    atomic_json(BASE/'P1_plan.json',dict(approval=approval,jobs=jobs,source_hashes=hashes(),duration=.1,segment=.001,deltas=[.0001,.00002],workers=8,worker_cap=100000,parent_cap=1800,job_cap=12000,environment=dict(python=sys.version,numpy=np.__version__,platform=platform.platform(),numerical_library_threads=1)))
def worker(job,directory):
    directory=Path(directory);plan=read(BASE/'P1_plan.json');assert plan['source_hashes']==hashes() and job in plan['jobs'];path=directory/job['id'];path.mkdir(exist_ok=False);started=time.perf_counter();deadline=started+11980
    try:
        canonical,contrast=setup();model=contrast.model
        with np.load(BASE/'prepared_inputs.npz') as z:initial=z['initial' if job['condition']=='prepared' else 'homogeneous'];detector=Coherence(z['reference'],z['initial'],float(z['omega']))
        pop=np.stack([initial,initial]);rng=np.random.default_rng(job['seed']);frames=[pop.copy()];metrics=[[detector.measure(x) for x in pop]];events=0;atomic_json(path/'metadata.json',dict(job=job,source_hashes=hashes()))
        for segment in range(100):
            previous=pop.copy();before=json.dumps(rng.bit_generator.state);out=pair(model,*pop,[.0001,.00002],[0,.001],job['seed'],deadline,rng=rng,sentinel=job['sentinel']);pop=out['counts'][:,-1];hist=out['histogram'];events+=out['events']
            for side in range(2):
                flat=previous[side].ravel().copy()
                for i,ch in enumerate(model.channels):
                    for coordinate,v in ch.change:flat[coordinate]+=int(hist[0,i]+hist[side+1,i])*v
                np.testing.assert_array_equal(flat,pop[side].ravel());assert sum(map(int,flat))==sum(map(int,initial.ravel()))
            np.savez_compressed(path/f'checkpoint_{segment:03d}.npz',initial=previous,final=pop,available=out['final_available_before_reset'],histogram=hist,rng_before=before,rng_after=out['rng_state'],events=out['events'])
            if job['sentinel']:save_trace(path/f'trace_{segment:03d}.npz',out)
            del out;frames.append(pop.copy());metrics.append([detector.measure(x) for x in pop]);atomic_json(path/'progress.json',dict(segments_complete=segment+1,required=100,events=events))
        np.savez_compressed(path/'observations.npz',counts=frames,times=np.arange(101)*.001,seed=job['seed']);atomic_json(path/'metrics.json',dict(values=metrics));atomic_json(directory/(job['id']+'.json'),dict(status='COMPLETE',job=job,seconds=time.perf_counter()-started,events=events,source_hashes=hashes()))
    except BaseException:
        atomic_json(directory/(job['id']+'.json'),dict(status='FAILED',job=job,seconds=time.perf_counter()-started,error=traceback.format_exc()));raise
def production():
    plan=read(BASE/'P1_plan.json');assert hashes()==plan['source_hashes'] and not (BASE/'P1_attempt.json').exists();subprocess.run(['git','ls-files','--error-unmatch',str((BASE/'P1_plan.json').relative_to(ROOT))],cwd=ROOT,check=True,stdout=subprocess.DEVNULL);atomic_json(BASE/'P1_attempt.json',dict(status='RUNNING'));start=time.perf_counter();forecast=read(BASE/'preparation_results.json')['candidate_P1']['forecast_worker_seconds']/8
    ledger=pool_run(plan['jobs'],BASE/'P1_production',Path(__file__).resolve(),ROOT,100000,workers=8,timeout=12000,forecast={j['id']:forecast for j in plan['jobs']});complete=len(ledger['records'])==8 and all(r['status']=='COMPLETE' and r['returncode']==0 for r in ledger['records']);atomic_json(BASE/'P1_attempt.json',dict(status='COMPLETE' if complete else 'STOPPED',wall_seconds=time.perf_counter()-start,charged_worker_seconds=ledger['charged_worker_seconds'],stop_reason=ledger['stop_reason']))
if __name__=='__main__':
    if sys.argv[1]=='--worker':worker(json.loads(sys.argv[2]),sys.argv[3])
    elif sys.argv[1]=='register':register()
    elif sys.argv[1]=='run':production()
    else:raise ValueError('Choose register or run')

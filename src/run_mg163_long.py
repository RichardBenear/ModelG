"""MG163 parameterized long-horizon driver; prior frozen engines remain unchanged."""
import os
for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[k]='1'
import json,sys,time,hashlib,traceback,shutil,platform,subprocess
from pathlib import Path
import numpy as np
from run_mg160_persistence import FILES as OLD_FILES,Coherence,setup,pair,save_trace,seed_for,read,atomic_json,pool_run
ROOT=Path(__file__).resolve().parents[1];PUBLIC=ROOT/'results/MG12H3/long_horizon_persistence';BASE=ROOT/'local_raw/MG163_long_horizon'
FILES=OLD_FILES+['src/run_mg163_long.py','src/check_mg163_long.py','src/qualify_mg163_long.py','src/score_mg163_long.py','src/check_mg163_shape.py','src/mg161_local_shape.py','spec/MG163_long_horizon_persistence.md']
def hashes():
    h={p:hashlib.sha256((ROOT/p).read_text(encoding='utf-8').replace('\r\n','\n').encode()).hexdigest() for p in FILES}
    for p in ['results/MG12H3/pattern_persistence/prepared_inputs.npz','results/MG12H3/ring_feasibility/o0_a0_Omega10.npz','results/MG12H3/noise_aware_shape/calibration.json']:h[p]=hashlib.sha256((ROOT/p).read_bytes()).hexdigest()
    return h

def trace_segment(job,segment):return job['sentinel'] and segment%10==0

def storage_guard(directory,cap=4000000000):
    if sum(p.stat().st_size for p in Path(directory).rglob('*') if p.is_file())>=cap:raise RuntimeError('LOCAL STORAGE CAP')

def resource_guard(free):
    if free<8000000000:raise RuntimeError('Require 8 GB free local storage')

def register():
    approval=read(PUBLIC/'approval.json');assert approval['worker_cap']==800000 and approval['pairs']==8
    assert read(PUBLIC/'qualification.json')['passed'] and not (PUBLIC/'plan.json').exists()
    used=set()
    def collect(x):
        if isinstance(x,dict):
            for v in x.values():collect(v)
        elif isinstance(x,list):
            for v in x:collect(v)
        elif type(x) is int and 0<=x<2**32:used.add(x)
    for p in (ROOT/'results').rglob('*.json'):
        if 'plan' in p.name or 'registration' in p.name:collect(read(p))
    jobs=[]
    for c,condition in enumerate(('prepared','homogeneous')):
        for r in range(4):
            nonce=0;seed=seed_for(12,c,0,16300,r,nonce)
            while seed in used:nonce+=1;seed=seed_for(12,c,0,16300,r,nonce)
            used.add(seed);jobs.append(dict(id=f'{condition}_r{r}',condition=condition,replicate=r,seed=seed,seed_nonce=nonce,sentinel=r==0))
    atomic_json(PUBLIC/'plan.json',dict(stage='MG163',approval=approval,jobs=jobs,source_hashes=hashes(),segments=1000,segment=.001,duration=1.,deltas=[.0001,.00002],workers=8,worker_cap=800000,job_cap=100000,parent_cap=3600,trace_stride=10,storage_cap=4000000000,forecast_worker_seconds=563540.728079,environment=dict(python=sys.version,numpy=np.__version__,platform=platform.platform(),threads=1)))

def run():
    p=read(PUBLIC/'plan.json');assert p['source_hashes']==hashes() and not (PUBLIC/'attempt.json').exists() and not BASE.exists();resource_guard(shutil.disk_usage(ROOT).free)
    subprocess.run(['git','ls-files','--error-unmatch',str((PUBLIC/'plan.json').relative_to(ROOT))],cwd=ROOT,check=True,stdout=subprocess.DEVNULL)
    BASE.mkdir();shutil.copy2(ROOT/'results/MG12H3/pattern_persistence/prepared_inputs.npz',BASE/'prepared_inputs.npz');atomic_json(BASE/'plan.json',p);atomic_json(PUBLIC/'attempt.json',dict(status='RUNNING'));start=time.perf_counter()
    ledger=pool_run(p['jobs'],BASE/'production',Path(__file__).resolve(),ROOT,p['worker_cap'],workers=8,timeout=p['job_cap'],forecast={j['id']:p['forecast_worker_seconds']/8 for j in p['jobs']})
    done=len(ledger['records'])==8 and all(r['status']=='COMPLETE' and r['returncode']==0 for r in ledger['records']);atomic_json(PUBLIC/'attempt.json',dict(status='COMPLETE' if done else 'STOPPED',wall_seconds=time.perf_counter()-start,charged_worker_seconds=ledger['charged_worker_seconds'],stop_reason=ledger['stop_reason']))
def worker(job,directory):
    directory=Path(directory);plan=read(BASE/'plan.json');assert plan['source_hashes']==hashes() and job in plan['jobs'];path=directory/job['id'];path.mkdir(exist_ok=False);started=time.perf_counter();deadline=started+plan['job_cap']-20
    try:
        canonical,contrast=setup();model=contrast.model
        with np.load(BASE/'prepared_inputs.npz') as z:initial=z['initial' if job['condition']=='prepared' else 'homogeneous'];detector=Coherence(z['reference'],z['initial'],float(z['omega']))
        pop=np.stack([initial,initial]);rng=np.random.default_rng(job['seed']);frames=[pop.copy()];metrics=[[detector.measure(x) for x in pop]];events=0;atomic_json(path/'metadata.json',dict(job=job,source_hashes=hashes()))
        for segment in range(plan['segments']):
            storage_guard(BASE)
            previous=pop.copy();before=json.dumps(rng.bit_generator.state);out=pair(model,*pop,[.0001,.00002],[0,.001],job['seed'],deadline,rng=rng,sentinel=trace_segment(job,segment));pop=out['counts'][:,-1];hist=out['histogram'];events+=out['events']
            for side in range(2):
                flat=previous[side].ravel().copy()
                for i,ch in enumerate(model.channels):
                    for coordinate,v in ch.change:flat[coordinate]+=int(hist[0,i]+hist[side+1,i])*v
                np.testing.assert_array_equal(flat,pop[side].ravel());assert sum(map(int,flat))==sum(map(int,initial.ravel()))
            np.savez_compressed(path/f'checkpoint_{segment:03d}.npz',initial=previous,final=pop,available=out['final_available_before_reset'],histogram=hist,rng_before=before,rng_after=out['rng_state'],events=out['events'])
            if trace_segment(job,segment):save_trace(path/f'trace_{segment:03d}.npz',out)
            del out;frames.append(pop.copy());metrics.append([detector.measure(x) for x in pop]);atomic_json(path/'progress.json',dict(segments_complete=segment+1,required=plan['segments'],events=events))
        np.savez_compressed(path/'observations.npz',counts=frames,times=np.arange(plan['segments']+1)*.001,seed=job['seed']);atomic_json(path/'metrics.json',dict(values=metrics));atomic_json(directory/(job['id']+'.json'),dict(status='COMPLETE',job=job,seconds=time.perf_counter()-started,events=events,source_hashes=hashes()))
    except BaseException:
        atomic_json(directory/(job['id']+'.json'),dict(status='FAILED',job=job,seconds=time.perf_counter()-started,error=traceback.format_exc()));raise

if __name__=='__main__':
    if sys.argv[1]=='--worker':worker(json.loads(sys.argv[2]),sys.argv[3])
    elif sys.argv[1]=='register':register()
    elif sys.argv[1]=='run':run()
    else:raise ValueError('Choose register/run')

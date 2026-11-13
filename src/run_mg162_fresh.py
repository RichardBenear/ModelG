"""MG162B prospective independent-data pilot; explicit concrete approval required."""
import os
for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):os.environ[k]='1'
from pathlib import Path
import json,sys,time,hashlib,shutil,subprocess,traceback,platform
import numpy as np
import run_mg160_persistence as writer
from mg161_local_shape import LocalShape
from runner import seed_for
from mg10e2c1_pool import atomic_json,run as pool_run
ROOT=writer.ROOT;BASE=ROOT/'results/MG12H3/shape_real_noise_challenge';RAW=ROOT/'local_raw/MG162_fresh_model'
FILES=writer.FILES+['src/run_mg162_fresh.py','src/mg161_local_shape.py','src/check_mg160_pilot.py','spec/MG162_shape_validation.md','spec/MG162B_fresh_pilot.md']
def read(p):return json.loads(Path(p).read_text(encoding='utf-8'))
def hashes():
    out={p:hashlib.sha256((ROOT/p).read_text(encoding='utf-8').replace('\r\n','\n').encode()).hexdigest() for p in FILES}
    for p in (ROOT/'results/MG12H3/pattern_persistence/prepared_inputs.npz',ROOT/'results/MG12H3/ring_feasibility/o0_a0_Omega10.npz',ROOT/'results/MG12H3/noise_aware_shape/calibration.json'):out[p.relative_to(ROOT).as_posix()]=hashlib.sha256(p.read_bytes()).hexdigest()
    return out
def register():
    approval=read(BASE/'fresh_approval.json');assert approval['pairs']==8 and approval['worker_cap']==100000;assert not (BASE/'fresh_plan.json').exists();assert read(BASE/'challenge_results.json')['supported_within_finite_challenge'];used=set()
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
            nonce=0;seed=seed_for(12,c,0,16200,r,nonce)
            while seed in used:nonce+=1;seed=seed_for(12,c,0,16200,r,nonce)
            used.add(seed);jobs.append(dict(id=f'fresh_{condition}_r{r}',condition=condition,replicate=r,seed=seed,seed_nonce=nonce,sentinel=r==0))
    atomic_json(BASE/'fresh_plan.json',dict(stage='MG162B',approval=approval,jobs=jobs,source_hashes=hashes(),worker_cap=100000,parent_cap=1800,job_cap=12000,workers=8,threads=1,duration=.1,observation=.001,deltas=[.0001,.00002],raw_location='local_raw/MG162_fresh_model',environment=dict(python=sys.version,numpy=np.__version__,platform=platform.platform(),threads=1),primary_rule='All4preparedpairs core-compatible atall101frames/bothmarginals and no boundaryadvisoryframes in preparedpaths and nohomogeneouspath compatible atanyframe; otherwise notestablished. Finitepilot only, not population qualification.'))
def worker(job,directory):
    assert hashlib.sha256((RAW/'prepared_inputs.npz').read_bytes()).hexdigest()==hashes()['results/MG12H3/pattern_persistence/prepared_inputs.npz']
    previous_base,previous_hash=writer.BASE,writer.hashes
    try:writer.BASE=RAW;writer.hashes=hashes;writer.worker(job,directory)
    finally:writer.BASE=previous_base;writer.hashes=previous_hash
def run():
    plan=read(BASE/'fresh_plan.json');assert plan['source_hashes']==hashes();assert not (BASE/'fresh_attempt.json').exists();subprocess.run(['git','ls-files','--error-unmatch',str((BASE/'fresh_plan.json').relative_to(ROOT))],cwd=ROOT,check=True,stdout=subprocess.DEVNULL);RAW.mkdir(parents=True,exist_ok=False);shutil.copy2(ROOT/'results/MG12H3/pattern_persistence/prepared_inputs.npz',RAW/'prepared_inputs.npz');atomic_json(RAW/'P1_plan.json',plan);atomic_json(BASE/'fresh_attempt.json',dict(status='RUNNING'));start=time.perf_counter()
    ledger=pool_run(plan['jobs'],RAW/'P1_production',Path(__file__).resolve(),ROOT,100000,workers=8,timeout=12000,forecast={j['id']:36912.6878852/8*1.25 for j in plan['jobs']});complete=len(ledger['records'])==8 and all(r['status']=='COMPLETE' and r['returncode']==0 for r in ledger['records']);atomic_json(BASE/'fresh_attempt.json',dict(status='COMPLETE' if complete else 'STOPPED',wall_seconds=time.perf_counter()-start,charged_worker_seconds=ledger['charged_worker_seconds'],stop_reason=ledger['stop_reason']));atomic_json(RAW/'P1_attempt.json',read(BASE/'fresh_attempt.json'))
def audit():
    import check_mg160_pilot as checker
    prior_base,prior_hash=checker.BASE,checker.hashes
    try:checker.BASE=RAW;checker.hashes=hashes;checker.check()
    finally:checker.BASE=prior_base;checker.hashes=prior_hash
def score():
    started=time.perf_counter();attempt=BASE/'fresh_score_attempt.json';assert not attempt.exists();atomic_json(attempt,dict(status='RUNNING'))
    try:
        plan=read(BASE/'fresh_plan.json');assert hashes()==plan['source_hashes'] and read(BASE/'fresh_attempt.json')['status']=='COMPLETE';assert read(RAW/'P1_analysis_results.json')['passed'];cal=read(ROOT/'results/MG12H3/noise_aware_shape/calibration.json')
        with np.load(RAW/'prepared_inputs.npz') as z:model=LocalShape(z['reference'][1],z['initial'][1],float(z['omega']))
        rows=[]
        for j in plan['jobs']:
            with np.load(RAW/'P1_production'/j['id']/'observations.npz') as z:counts=z['counts'];times=z['times'];assert counts.shape==(101,2,7,500) and int(z['seed'])==j['seed']
            sides=[]
            for s in (0,1):
                fits=[model.fit(frame[s,1]/model.omega) for frame in counts];labels=[model.classify(f,cal['detection_threshold'],cal['residual_threshold']) for f in fits];bad=[i for i,l in enumerate(labels) if l!='core_compatible'];edge=[i for i,f in enumerate(fits) if labels[i]=='core_compatible' and (f['amplitude']/model.initial_amplitude<=.6 or min(abs(f['scale']-.5),abs(f['scale']-2))<=.125)];sides.append(dict(side=s,compatible_frames=labels.count('core_compatible'),all_compatible=not bad,any_compatible='core_compatible' in labels,first_noncompatible_bracket=[float(times[max(0,bad[0]-1)]),float(times[bad[0]])] if bad else None,boundary_advisory_frames=edge,labels=labels,fits=fits))
            rows.append(dict(id=j['id'],condition=j['condition'],outcomes=sides))
        supported=all(s['all_compatible'] for r in rows if r['condition']=='prepared' for s in r['outcomes']) and all(not s['any_compatible'] for r in rows if r['condition']=='homogeneous' for s in r['outcomes']) and all(not s['boundary_advisory_frames'] for r in rows if r['condition']=='prepared' for s in r['outcomes']);atomic_json(BASE/'fresh_candidate_scores.json',dict(pilot_separation_supported=supported,jobs=rows,source_hashes=hashes(),qualification_status='Candidate scores only; causal/raw audit and separate metric recomputation required before final report.',boundary_advisory='Engineering sensitivity flags only; no fitted-parameter confidence intervals and no hidden score changes.',elapsed_seconds=time.perf_counter()-started));atomic_json(attempt,dict(status='COMPLETE',elapsed_seconds=time.perf_counter()-started))
    except BaseException:atomic_json(attempt,dict(status='FAILED',elapsed_seconds=time.perf_counter()-started,error=traceback.format_exc()));raise
if __name__=='__main__':
    if sys.argv[1]=='--worker':worker(json.loads(sys.argv[2]),sys.argv[3])
    elif sys.argv[1]=='register':register()
    elif sys.argv[1]=='run':run()
    elif sys.argv[1]=='audit':audit()
    elif sys.argv[1]=='score':score()
    else:raise ValueError('Choose register/run/score')

"""MG073 E2 prospective qualification and short timing only, no production."""
import os
for key in ("OPENBLAS_NUM_THREADS","MKL_NUM_THREADS","OMP_NUM_THREADS","NUMEXPR_NUM_THREADS"):os.environ[key]="1"
from pathlib import Path
import argparse,hashlib,json,platform,sys,time,traceback
from mg10d1_pool import run,atomic_json
ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/"results/MG10E2/preparation"
FILES=("src/mg10e2.py","src/prepare_mg10e2.py","src/check_mg10e2_preparation.py","spec/MG10E2_proposal.md","src/mg10e1.py","src/mg10e1_stream.py",
 "src/mg10d2.py","src/mg10d2_source.py","src/mg10s1.py","src/mg10s2.py","src/mg10p1.py","src/mg10p2.py","src/runner.py","src/mg10d1_pool.py","src/check_mg10d1_development.py")
def digest(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def hashes():
    values={p:hashlib.sha256((ROOT/p).read_text(encoding="utf-8").replace("\r\n","\n").encode()).hexdigest() for p in FILES}
    for s in range(2):
        for k in range(2):
            name=f"results/MG10E1/preparation/devparents/parent_s{s}_k{k}.npz";values[name]=digest(ROOT/name)
            for b in range(8):
                name=f"results/MG10E1/parents/parent_s{s}_k{k}_b{b}.json";values[name]=hashlib.sha256((ROOT/name).read_text(encoding="utf-8").replace("\r\n","\n").encode()).hexdigest()
    values["results/MG003/ref_profiles.npz"]=digest(ROOT/"results/MG003/ref_profiles.npz");return values
def candidate_jobs():
    from runner import seed_for
    p=[dict(id=f"parent_s{s}_k{k}_b{b}",kind="parent",start=s,step=k,block=b,seed=seed_for(10,s,k,7303,b),warmup_seed=seed_for(10,s,k,7303,b,1)) for b in range(8) for k in range(2) for s in range(2)]
    branches=[dict(id=f"a{a}_s{s}_k{k}_b{b}",kind="branch",parent=f"parent_s{s}_k{k}_b{b}",arm=a,start=s,step=k,block=b,seed=seed_for(10,s,k,7304,b,a+1),warmup_seed=seed_for(10,s,k,7304,b,101+a)) for b in range(8) for k in range(2) for s in range(2) for a in range(4)]
    return p,branches
def plan():
    import numpy,scipy,numba
    from runner import seed_for
    from mg10e2 import activity
    if (BASE/"preparation_plan.json").exists():raise ValueError("Already planned; preserve/version changes")
    p,b=candidate_jobs();dev=[seed_for(10,a,k,7300,0) for a in range(4) for k in range(2)]
    dev += [seed_for(10,s,k,7302,a+1,seg) for s in range(2) for k in range(2) for a in range(4) for seg in (0,1)]
    seeds=[j[key] for j in p+b for key in ("seed","warmup_seed")]
    assert len(set(seeds))==320 and len(set(dev))==40 and not set(seeds)&set(dev)
    prior=set();count=0
    def collect(v):
        if isinstance(v,dict):
            for key,value in v.items():
                if "seed" in key.lower() and isinstance(value,int):prior.add(value)
                else:collect(value)
        elif isinstance(v,list):
            for value in v:collect(value)
    for directory in (ROOT/"results").glob("MG10*"):
        if directory.name=="MG10E2":continue
        for path in directory.rglob("*.json"):
            collect(json.loads(path.read_text(encoding="utf-8")));count+=1
    assert not (set(seeds)|set(dev))&prior
    atomic_json(BASE/"preparation_plan.json",dict(status="DRAFT E2; PREPARATION ONLY",source_hashes=hashes(),candidate_parents=p,candidate_branches=b,development_seeds=dev,
      prior_recorded_seed_count=len(prior),prior_json_files=count,production_seed_values=320,development_seed_values=40,collisions=0,
      activities_before={str(a):activity(a,0).tolist() for a in range(4)},activities_after={str(a):activity(a,5).tolist() for a in range(4)},
      parent_time=5,post_time=25,switch_time=5,late_window=[20,25],development_ceiling_seconds=900,candidate_production_ceiling_seconds=144000,
      timing_time=.25,timing_switch=.125,environment=dict(python=sys.version,platform=platform.platform(),numpy=numpy.__version__,scipy=scipy.__version__,numba=numba.__version__)))
def save_raw(path,data,case,stats,times,seed,arm,switch):
    import numpy as np
    np.savez_compressed(path,data=data,initial=case["initial"],a0=case["a0"],b0=case["b0"],omega=case["omega"],clamp_a=case["clamp_a"],clamp_b=case["clamp_b"],stats=stats,times=times,seed=seed,arm=arm,switch=switch)
def qualify(directory):
    import numpy as np
    import numba as nb
    import mg10e2 as e2
    from mg10e1 import simulate as constant
    from mg10e1 import _gated_leap
    from mg10p2 import fixture,reference
    from runner import seed_for
    from check_mg10d1_development import check_raw
    # Reference checks schedule selection every leap, independent of wrapper
    # segment construction. Shared frozen event kernel explicitly disclosed.
    @nb.njit
    def oracle(initial,a0,b0,diff,tau,times,seed,before,after,switch):
        np.random.seed(seed);state=initial.copy();f=np.zeros((7,500),np.int64);h=np.zeros((6,500),np.int64);x=np.zeros((2,500),np.int64)
        out=np.zeros((1,len(times),22,500),np.int64);out[0,0,:7]=state;stats=np.zeros((1,3));t=0.
        for i in range(1,len(times)):
            target=times[i]
            while t<target:
                stop=min(target,switch) if t<switch else target
                requested=min(tau,stop-t)
                dt,prob=_gated_leap(state,f,h,a0,b0,1000.,diff,np.ones(7),requested,True,True,x,before if t<switch else after)
                stats[0,0]+=1;stats[0,1]+=int(dt<requested*(1-1e-12));stats[0,2]=max(stats[0,2],prob)
                if dt==stop-t:t=stop
                else:t+=dt
            out[0,i,:7]=state;out[0,i,7:14]=f;out[0,i,14:20]=h;out[0,i,20:22]=x
        return out,stats
    tests=[];checks=[]
    for a in range(4):
        for k in range(2):
            case=fixture(0,1);case["initial"][5]=123;case["initial"][6]=456;tau=(1/24000,1/48000)[k];times=np.array([0.,.001,.003,.005]);switch=.002;seed=seed_for(10,a,k,7300,0)
            data,stats=e2.simulate(case,tau,times,seed,a,switch)
            expected,st=oracle(case["initial"],case["a0"],case["b0"],case["diff"],tau,times,seed,e2.activity(a,0,switch),e2.activity(a,switch,switch),switch if a in (2,3) else -1.)
            np.testing.assert_array_equal(data,expected);np.testing.assert_array_equal(stats,st)
            if a in (0,1):
                old,st=constant(case,tau,times,1,seed,a_activity=e2.activity(a,0,switch));np.testing.assert_array_equal(data,old);np.testing.assert_array_equal(stats,st)
            if a==2:
                noop,ns=e2.simulate(case,tau,times,seed,a,0);old,st=constant(case,tau,times,1,seed);np.testing.assert_array_equal(noop,old);np.testing.assert_array_equal(ns,st)
            path=directory/f"schedule_a{a}_k{k}.npz";save_raw(path,data,case,stats,times,seed,a,switch);checks.append(check_raw(path));tests.append(dict(name=f"per_leap_schedule_and_identity_a{a}_k{k}",passed=True))
    for args in ((-1,0,5),(4,0,5),(0,-1,5),(0,0,-1),(0,float('nan'),5)):
        try:e2.activity(*args);raise AssertionError("Invalid schedule accepted")
        except ValueError:pass
    tests.append(dict(name="invalid_schedule_rejection",passed=True))
    ref=reference()[0];background=np.array([140/9,140/9,29/(140/9)])[:,None];fields=[];labels=[]
    for label,f in (("flat",background+np.zeros((3,500))),("reference",background+ref),("translated",background+np.roll(ref,125,axis=-1)),("scaled_amplitude",background+1.25*ref),("mean_offset",background+ref+3),("deformed",background+ref+np.roll(ref,25,axis=-1)*.5),("two_peaks",background+.5*ref+.5*np.roll(ref,200,axis=-1))):labels.append(label);fields.append(f)
    fields=np.array(fields);metrics=e2.field_metrics(fields)
    assert metrics["x_phase_grid"][0]==0 and abs(metrics["projection"][0,1])<1e-9
    np.testing.assert_allclose(metrics["projection"][1],1,atol=1e-10);np.testing.assert_allclose(metrics["x_aligned_projection"][2],1,atol=1e-10);assert metrics["x_phase_grid"][2]==25
    np.testing.assert_allclose(metrics["projection"][3],1.25,atol=1e-10);np.testing.assert_allclose(metrics["projection"][4],metrics["projection"][1],atol=1e-10)
    assert metrics["x_aligned_residual_rms"][5,1]>0 and metrics["x_aligned_residual_rms"][6,1]>0
    np.savez_compressed(directory/"metric_inputs.npz",fields=fields)
    atomic_json(directory/"metric_results.json",dict(labels=labels,metrics={k:v.tolist() for k,v in metrics.items()}));tests.append(dict(name="translation_amplitude_background_deformation_and_tie",passed=True))
    times=e2.OBS;linear=np.stack([2+3*times,-1+2*times,0*times],axis=-1);np.testing.assert_allclose(e2.window_average(linear,times),[69.5,44,0])
    late=np.arange(24,dtype=float).reshape(2,4,3)/10;late[1,2]+=np.array([.5,1.,1.5]);shift=np.arange(12).reshape(4,3)/10
    effects=e2.paired_effects(late);shifted=e2.paired_effects(late+shift[None])
    for key in effects:np.testing.assert_allclose(effects[key],shifted[key],atol=1e-12)
    records={}
    for k in range(2):
        for b in range(8):
            for s in range(2):
                for a in range(4):records[f"a{a}_s{s}_k{k}_b{b}"]=dict(complete=True,late_projection=(late[s,a]+s*(a==2)*b*.01).tolist())
    ef,full=e2.balanced_analysis(records);assert full['0']['n']==full['1']['n']==8
    missing=dict(records);del missing['a2_s1_k0_b0'];_,partial=e2.balanced_analysis(missing);assert partial['0']['n']==7 and partial['1']['n']==8
    prefix=dict(records);prefix['a2_s1_k0_b0']=dict(complete=False);_,truncated=e2.balanced_analysis(prefix);assert truncated==partial
    np.savez_compressed(directory/"contrast_inputs.npz",times=times,linear=linear,late=late,shift=shift)
    atomic_json(directory/"contrast_results.json",dict(effects={k:v.tolist() for k,v in effects.items()},late_average=e2.window_average(linear,times).tolist(),records=records,full=full,missing=partial,prefix=truncated))
    tests.append(dict(name="integral_background_cancellation_balanced_missing_and_prefix",passed=True))
    atomic_json(directory/"qualification_results.json",dict(passed=True,tests=tests,raw_checks=checks,status="DEVELOPMENT ONLY; NO E2 SCIENCE"))
def worker(job,directory):
    directory=Path(directory)
    if hashes()!=json.loads((BASE/"preparation_plan.json").read_text())["source_hashes"]:raise ValueError("Frozen development source changed")
    if job['kind']=='qualification':qualify(directory);atomic_json(directory/(job['id']+'.json'),dict(status='COMPLETE',job=job));return
    import numpy as np
    import mg10e2 as e2
    from mg10p2 import fixture
    from runner import seed_for
    from check_mg10d1_development import check_raw
    with np.load(ROOT/f"results/MG10E1/preparation/devparents/parent_s{job['start']}_k{job['step']}.npz") as raw:case=e2.fork_case(fixture(0,job['start']),raw['data'][0,-1,:7])
    tau=(1/24000,1/48000)[job['step']];times=np.array([0.,.125,.25]);seed=seed_for(10,job['start'],job['step'],7302,job['arm']+1)
    e2.simulate(case,tau,[0.,tau],seed_for(10,job['start'],job['step'],7302,job['arm']+1,1),job['arm'],.125)
    clock=time.perf_counter();data,stats=e2.simulate(case,tau,times,seed,job['arm'],.125);elapsed=time.perf_counter()-clock
    path=directory/(job['id']+'.npz');save_raw(path,data,case,stats,times,seed,job['arm'],.125);checked=check_raw(path)
    atomic_json(directory/(job['id']+'.json'),dict(status='COMPLETE',job=job,seed=seed,raw=path.name,sha256=digest(path),compute_seconds=elapsed,raw_check=checked))
def execute():
    if (BASE/'qualification/ledger.json').exists():raise ValueError('Already attempted; no automatic retry')
    q=run([dict(id='qualification',kind='qualification')],BASE/'qualification',Path(__file__).resolve(),ROOT,900,workers=1,timeout=300,poll_seconds=.25)
    if q['stop_reason']:raise ValueError('Qualification failed')
    jobs=[dict(id=f'timing_a{a}_s{s}_k{k}',kind='timing',arm=a,start=s,step=k) for a in range(4) for s in range(2) for k in range(2)]
    timing=run(jobs,BASE/'timing',Path(__file__).resolve(),ROOT,900,workers=8,prior=q['charged_worker_seconds'],timeout=90,poll_seconds=.25)
    atomic_json(BASE/'execution.json',dict(qualification=q,timing=timing))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--plan',action='store_true');p.add_argument('--execute',action='store_true');p.add_argument('--worker',nargs=2);a=p.parse_args()
    if a.plan:plan()
    elif a.execute:execute()
    elif a.worker:
        j=json.loads(a.worker[0]);d=Path(a.worker[1])
        try:worker(j,d)
        except BaseException:atomic_json(d/(j['id']+'.json'),dict(status='FAILED',job=j,error=traceback.format_exc()));raise
    else:p.error('Preparation only: plan/execute/worker')

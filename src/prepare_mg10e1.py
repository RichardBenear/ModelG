"""MG069 E1 source-gate development/qualification/timing, not production."""
import os
for key in ("OPENBLAS_NUM_THREADS","MKL_NUM_THREADS","OMP_NUM_THREADS","NUMEXPR_NUM_THREADS"):
    os.environ[key]="1"
from pathlib import Path
import argparse,hashlib,json,platform,sys,time,traceback
from mg10d1_pool import run,atomic_json

ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/"results/MG10E1/preparation"
FILES=("src/mg10e1.py","src/prepare_mg10e1.py","src/check_mg10e1_preparation.py","spec/MG10E1_proposal.md",
       "src/mg10d2.py","src/mg10d2_source.py","src/mg10p2.py","src/mg10p1.py","src/mg10s1.py","src/mg10s2.py","src/mg10d1_pool.py","src/runner.py")


def digest(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def hashes():
    d={p:hashlib.sha256((ROOT/p).read_text(encoding="utf-8").replace("\r\n","\n").encode()).hexdigest() for p in FILES}
    d["results/MG003/ref_profiles.npz"]=digest(ROOT/"results/MG003/ref_profiles.npz")
    for s in range(2):
        for k in range(2):
            for rep in range(2):
                p=f"results/MG10P3/preparation/timing/timing_s{s}_k{k}_r{rep}.json"
                d[p]=hashlib.sha256((ROOT/p).read_text(encoding="utf-8").replace("\r\n","\n").encode()).hexdigest()
    return d


def candidate_jobs():
    from runner import seed_for
    parents=[dict(id=f"parent_s{s}_k{k}_b{b}",start=s,step=k,block=b,seed=seed_for(10,s,k,6903,b),warmup_seed=seed_for(10,s,k,6903,b,1)) for b in range(8) for k in range(2) for s in range(2)]
    branches=[dict(id=f"a{a}_s{s}_k{k}_b{b}",parent=f"parent_s{s}_k{k}_b{b}",arm=a,start=s,step=k,block=b,
                  seed=seed_for(10,s,k,6904,b,a+1),warmup_seed=seed_for(10,s,k,6904,b,101+a)) for b in range(8) for k in range(2) for s in range(2) for a in range(4)]
    return parents,branches


def plan():
    import numpy,scipy,numba
    from mg10e1 import activity
    if (BASE/"preparation_plan.json").exists():raise ValueError("Already planned; preserve before revision")
    p,b=candidate_jobs();seeds=[j[key] for j in p+b for key in ("seed","warmup_seed")];assert len(set(seeds))==320
    atomic_json(BASE/"preparation_plan.json",dict(status="DRAFT E1; PREPARATION ONLY",source_hashes=hashes(),candidate_parents=p,candidate_branches=b,
        masks={str(a):activity(a).tolist() for a in range(4)},parent_time=5,post_time=20,late_window=[15,20],blocks_per_step=8,
        development_ceiling_seconds=900,production_candidate_ceiling_seconds=129600,workers_max=8,threads_per_library=1,
        preparation_scope="Short.25 development parents and.25 timing forks only; no production5/20 trajectories or matcher-based selection",
        environment=dict(python=sys.version,platform=platform.platform(),numpy=numpy.__version__,scipy=scipy.__version__,numba=numba.__version__)))


def save_raw(path,data,case,stats,**extra):
    import numpy as np
    np.savez_compressed(path,data=data,initial=case["initial"],a0=case["a0"],b0=case["b0"],omega=case["omega"],clamp_a=case["clamp_a"],clamp_b=case["clamp_b"],stats=stats,**extra)


def qualify(directory):
    import numpy as np
    from scipy.stats import poisson
    import mg10e1 as new
    import mg10d2 as old
    import mg10s1 as small
    from runner import seed_for
    from check_mg10d1_development import check_raw
    tests=[];raw_checks=[]
    cases=[]
    for i in range(4):
        c=small.fixture(i);cases.append((f"small{i}",dict(initial=c["state"],a0=c["a0"],b0=c["b0"],omega=c["omega"],diff=c["diff"],clamp_a=False,clamp_b=False)))
    for start in range(2):
        cases.append((f"full_s{start}",new.baseline(0,start)))
        c=new.baseline(0,start);c["initial"][5]=123;c["initial"][6]=456;cases.append((f"prepared_like_s{start}",c))
    for i,(name,case) in enumerate(cases):
        seed=seed_for(10,i,0,6900,0);obs=np.array([0.,.002,.005]);tau=1/24000
        a,sa=old.simulate(case,tau,obs,1,seed);b,sb=new.simulate(case,tau,obs,1,seed)
        assert np.array_equal(a,b) and np.array_equal(sa,sb),name
        path=directory/(name+".npz");save_raw(path,b,case,sb,seed=seed,times=obs,a_activity=np.ones(case["initial"].shape[1]));raw_checks.append(check_raw(path))
        tests.append(dict(name="all_one_untouched_identity_"+name,passed=True,seed=seed))
    for arm in range(4):
        activity=new.activity(arm);assert activity.shape==(500,) and np.all((activity>=0)&(activity<=1))
        assert abs(activity.sum()-(500 if arm==0 else 449))<1e-10
        if arm in (1,2):assert np.sum(activity==0)==51
        case=new.baseline(0,1);seed=seed_for(10,arm,1,6900,1)
        data,stats=new.simulate(case,1/24000,np.array([0.,.005]),1,seed,a_activity=activity)
        assert new.ledger(data,case["initial"])["pass_check"]
        assert not np.any(data[:,:,7,activity==0]) and not np.any(data[:,:,20,activity==0])
        path=directory/f"mask_a{arm}.npz";save_raw(path,data,case,stats,seed=seed,times=np.array([0.,.005]),a_activity=activity);raw_checks.append(check_raw(path))
        tests.append(dict(name=f"spatial_mask_dose_off_events_a{arm}",passed=True,seed=seed))
    assert not np.any((new.activity(1)==0)&(new.activity(2)==0))
    case=new.baseline(0,1);state=case["initial"].copy();state[5]=123;state[6]=456
    fork=new.fork_case(case,state);assert np.array_equal(fork["initial"],state) and np.array_equal(fork["a0"],case["a0"])
    fork["initial"][0,0]+=1;fork["a0"][0]+=1;assert not np.array_equal(fork["initial"],state) and case["a0"][0]==4200000
    tests.append(dict(name="fork_preserves_all_species_denominators_without_alias",passed=True))
    for bad in (np.ones(499),np.full(500,-1.),np.full(500,1.01),np.full(500,np.nan)):
        try:new.simulate(case,1/24000,[.001],1,1,a_activity=bad);raise AssertionError("Bad activity accepted")
        except ValueError:pass
    tests.append(dict(name="activity_domain_rejection",passed=True))
    source=[];tau=1/24000;R=32768;epsilon=float(np.sqrt(np.log(6000)/(2*R)))
    for i,alpha in enumerate((0.,.898,1.)):
        c=dict(initial=np.array([[0],[0],[0],[4200000],[135333333],[0],[0]],dtype=np.int64),a0=np.array([4200000.]),b0=np.array([135333333.]),omega=1000.,diff=np.zeros(3),clamp_a=True,clamp_b=True)
        mask=np.zeros(7);mask[0]=1;seed=seed_for(10,i,2,6900,2)
        data,stats=new.simulate(c,tau,np.array([tau]),R,seed,mask=mask,a_activity=np.array([alpha]))
        counts=data[:,0,7,0];grid=np.arange(int(counts.max())+1);cdf=np.cumsum(np.bincount(counts,minlength=len(grid)))/R
        lam=14000*tau*alpha;distance=float(np.max(np.abs(cdf-poisson.cdf(grid,lam))));upper=distance+epsilon
        assert upper<=.03 and (alpha!=0 or not np.any(counts))
        path=directory/f"source_alpha{i}.npz";save_raw(path,data,c,stats,seed=seed,a_activity=np.array([alpha]),times=np.array([tau]));raw_checks.append(check_raw(path))
        source.append(dict(alpha=alpha,seed=seed,expected_mean=lam,mean=float(counts.mean()),distance=distance,epsilon=epsilon,upper=upper,margin=.03,raw=path.name,sha256=digest(path)))
    atomic_json(directory/"source_tests.json",source)
    times=new.POST_OBS;linear=np.stack([2+3*times,-1+2*times,times*0],axis=1)
    np.testing.assert_allclose(new.late_average(linear,times),[54.5,34.,0.])
    late=np.zeros((2,4,3));late[0]=np.array([0.,3.,1.,.5])[:,None];late[1]=np.array([1.,3.2,2.,1.2])[:,None]
    effect=new.paired_effects(late);np.testing.assert_allclose(effect["primary_local_remote"],[-.8]*3)
    shift=np.arange(12).reshape(4,3)/10;shifted=new.paired_effects(late+shift[None])
    for key in effect:np.testing.assert_allclose(effect[key],shifted[key],atol=1e-12)
    bumped=late.copy();bumped[1,1]+=2;np.testing.assert_allclose(new.paired_effects(bumped)["primary_local_remote"],effect["primary_local_remote"]+2)
    sample=np.arange(8)-3.5;ci=new.primary_interval(sample);assert ci["mean"]==0 and ci["ci_nominal"][0]<0<ci["ci_nominal"][1]
    assert new.primary_interval([])["mean"] is None and new.primary_interval([1])["ci_nominal"] is None
    np.savez_compressed(directory/"contrast_inputs.npz",times=times,linear=linear,late=late,shift=shift,sample=sample)
    atomic_json(directory/"contrast_results.json",dict(effects={k:v.tolist() for k,v in effect.items()},late_average=new.late_average(linear,times).tolist(),primary_interval=ci))
    tests.append(dict(name="late_integral_paired_background_cancellation_nominal_interval",passed=True))
    atomic_json(directory/"qualification_results.json",dict(passed=True,status="PREPARATION QUALIFICATION PASS; E1 SCIENCE NOT RUN",tests=tests,raw_checks=raw_checks))


def worker(job,directory):
    directory=Path(directory)
    if hashes()!=json.loads((BASE/"preparation_plan.json").read_text(encoding="utf-8"))["source_hashes"]:raise ValueError("Preparation source changed")
    if job["kind"]=="qualification":qualify(directory);atomic_json(directory/(job["id"]+".json"),dict(status="COMPLETE",job=job));return
    import numpy as np
    from runner import seed_for
    import mg10e1 as new
    import mg10d2 as old
    case=new.baseline(0,job["start"]);tau=(1/24000,1/48000)[job["step"]]
    if job["kind"]=="parent":
        seed=seed_for(10,job["start"],job["step"],6901,0)
        old.simulate(case,tau,np.array([tau]),1,seed_for(10,job["start"],job["step"],6901,0,1))
        clock=time.perf_counter();data,stats=old.simulate(case,tau,np.array([0.,.25]),1,seed);elapsed=time.perf_counter()-clock
        extra={}
    else:
        with np.load(BASE/f"devparents/parent_s{job['start']}_k{job['step']}.npz") as parent:case=new.fork_case(case,parent["data"][0,-1,:7])
        alpha=new.activity(job["arm"]);seed=seed_for(10,job["start"],job["step"],6902,job["arm"]+1)
        new.simulate(case,tau,np.array([tau]),1,seed_for(10,job["start"],job["step"],6902,job["arm"]+1,1),a_activity=alpha)
        clock=time.perf_counter();data,stats=new.simulate(case,tau,np.array([0.,.25]),1,seed,a_activity=alpha);elapsed=time.perf_counter()-clock
        assert not np.any(data[:,:,7,alpha==0]) and not np.any(data[:,:,20,alpha==0]);extra=dict(a_activity=alpha)
    if not new.ledger(data,case["initial"])["pass_check"] or np.max(stats[:,2])>.25+1e-14:raise ValueError("Timing ledger/safety")
    path=directory/(job["id"]+".npz");save_raw(path,data,case,stats,seed=seed,times=np.array([0.,.25]),**extra)
    atomic_json(directory/(job["id"]+".json"),dict(status="COMPLETE",job=job,seed=seed,compute_seconds=elapsed,raw=path.name,sha256=digest(path)))


def execute():
    if (BASE/"qualification/ledger.json").exists():raise ValueError("Already attempted; no automatic retry")
    q=run([dict(id="qualification",kind="qualification")],BASE/"qualification",Path(__file__).resolve(),ROOT,900,workers=1,timeout=300,poll_seconds=.25)
    if q["stop_reason"]:raise ValueError("Blocking E1 qualification failure")
    parents=[dict(id=f"parent_s{s}_k{k}",kind="parent",start=s,step=k) for s in range(2) for k in range(2)]
    p=run(parents,BASE/"devparents",Path(__file__).resolve(),ROOT,900,workers=4,prior=q["charged_worker_seconds"],timeout=60,poll_seconds=.25)
    if p["stop_reason"]:raise ValueError("Development parent failure")
    timing=[dict(id=f"timing_a{a}_s{s}_k{k}",kind="timing",arm=a,start=s,step=k) for a in range(4) for s in range(2) for k in range(2)]
    t=run(timing,BASE/"timing",Path(__file__).resolve(),ROOT,900,workers=8,prior=p["charged_worker_seconds"],timeout=90,poll_seconds=.25)
    atomic_json(BASE/"execution.json",dict(qualification=q,parents=p,timing=t))


if __name__=="__main__":
    p=argparse.ArgumentParser();p.add_argument("--plan",action="store_true");p.add_argument("--execute",action="store_true");p.add_argument("--worker",nargs=2);args=p.parse_args()
    if args.plan:plan()
    elif args.execute:execute()
    elif args.worker:
        j=json.loads(args.worker[0]);d=Path(args.worker[1])
        try:worker(j,d)
        except BaseException:atomic_json(d/(j["id"]+".json"),dict(status="FAILED",job=j,error=traceback.format_exc()));raise
    else:p.error("Preparation only: choose plan/execute/worker")

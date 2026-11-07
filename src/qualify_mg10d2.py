"""MG055 preregistered scalar-law diagnostic and corrected-code qualification."""
import os
for key in ("OPENBLAS_NUM_THREADS","MKL_NUM_THREADS","OMP_NUM_THREADS","NUMEXPR_NUM_THREADS"):
    os.environ[key]="1"
from pathlib import Path
import argparse,hashlib,json,math,platform,subprocess,sys,time,traceback
from mg10d1_pool import run,atomic_json

ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/"results/MG10D2"
FILES=("src/mg10d2_source.py","src/mg10d2.py","src/mg10d2_stream.py","src/qualify_mg10d2.py",
       "src/check_mg10d2_source.py","spec/MG10D2_proposal.md","src/mg10d1.py","src/mg10s1.py",
       "src/mg10s2.py","src/mg10d1_pool.py","src/runner.py")


def digest(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def hashes():return {p:hashlib.sha256((ROOT/p).read_text(encoding="utf-8").replace("\r\n","\n").encode()).hexdigest() for p in FILES}


def primary_jobs():
    from runner import seed_for
    return [dict(id=f"source_c{c}_m{m}",case=c,method=m,n=4200000000000,p=14000*(1/24000,1/48000)[c]/4200000000000,reps=2**24,
                 seed=seed_for(10,c,m,5501,0)) for c in range(2) for m in range(3)]


def register():
    import numpy,scipy,numba
    file=BASE/"diagnostic_registration.json"
    if file.exists():raise ValueError("Already registered")
    vendor=Path(numba.__file__).parent/"np/random/legacy.py"
    atomic_json(file,dict(status="REGISTERED BEFORE SOURCE DIAGNOSTIC",stage="MG055",source_hashes=hashes(),jobs=primary_jobs(),
        dkw_family_alpha=.01,primary_comparisons=6,revised_cdf_upper_margin=.0008,supplement_alpha=.001,supplement_margin=.02,
        diagnostic_ceiling_seconds=900,production_candidate_ceiling_seconds=86400,
        prior_D1_charge_seconds=120000.33734830027,combined_worst_seconds=120000.33734830027+900+86400,
        pi_authorization="looks ok, proceed with the next steps; allocation remains within prior approved216000worker-seconds",
        environment=dict(python=sys.version,platform=platform.platform(),numpy=numpy.__version__,scipy=scipy.__version__,numba=numba.__version__),
        vendor_source=dict(path=str(vendor),sha256=digest(vendor),relevant_expression="q=1-p; qn=q**n")))


def histogram(n,p,reps,seed,method):
    import numpy as np
    import numba as nb
    from mg10d2_source import finite_source,geometric_binomial
    @nb.njit
    def sample(n,p,reps,seed,method):
        np.random.seed(seed);h=np.zeros(65,dtype=np.int64);overflow=[];total=0;squares=0
        for _ in range(reps):
            k=np.random.binomial(n,p) if method==0 else finite_source(n,p) if method==1 else geometric_binomial(n,p)
            h[min(k,64)]+=1;total+=k;squares+=k*k
            if k>=64:overflow.append(k)
        return h,np.asarray(overflow,dtype=np.int64),total,squares
    if method in (0,1,3):return sample(n,p,reps,seed,method)
    rng=np.random.default_rng(seed);h=np.zeros(65,dtype=np.int64);overflow=[];total=0;squares=0
    for offset in range(0,reps,2**20):
        v=rng.binomial(n,p,size=min(2**20,reps-offset));h+=np.bincount(np.minimum(v,64),minlength=65)
        overflow.extend(v[v>=64]);total+=int(v.sum());squares+=int((v*v).sum())
    return h,np.asarray(overflow,dtype=np.int64),total,squares


def theory(n,p,maximum):
    import numpy as np
    prob=math.exp(n*math.log1p(-p));pmf=[prob]
    for k in range(1,maximum+1):pmf.append(pmf[-1]*(n-k+1)*p/(k*(1-p)))
    return np.cumsum(pmf)


def test_law(job,directory,alpha,comparisons,margin):
    import numpy as np
    from scipy.stats import binom
    clock=time.perf_counter();h,overflow,total,squares=histogram(job["n"],job["p"],job["reps"],job["seed"],job["method"])
    path=directory/(job["id"]+".npz");np.savez_compressed(path,hist=h,overflow=overflow,total=total,squares=squares)
    full=h[:64].copy()
    if len(overflow):full=np.pad(full,(0,max(0,int(overflow.max())+1-len(full))));full+=np.bincount(overflow,minlength=len(full))
    empirical=np.cumsum(full)/job["reps"];target=theory(job["n"],job["p"],len(full)-1)
    distance=float(max(np.max(np.abs(empirical-target)),1-float(target[-1])))
    epsilon=math.sqrt(math.log(2*comparisons/alpha)/(2*job["reps"]))
    result=dict(job=job,raw=path.name,sha256=digest(path),cdf_distance=distance,epsilon=epsilon,upper=distance+epsilon,lower=max(0,distance-epsilon),
        margin=margin,passed=bool(distance+epsilon<=margin),mean=total/job["reps"],expected_mean=job["n"]*job["p"],
        variance=squares/job["reps"]-(total/job["reps"])**2,expected_variance=job["n"]*job["p"]*(1-job["p"]),
        scipy_max_cdf_difference=float(np.max(np.abs(target-binom.cdf(np.arange(len(full)),job["n"],job["p"])))),elapsed_seconds=time.perf_counter()-clock)
    atomic_json(directory/(job["id"]+".json"),result)
    return result


def qualification(directory):
    import numpy as np
    import numba as nb
    import mg10d1 as old
    import mg10d2 as new
    import mg10s1 as small
    from mg10d2_source import finite_source,geometric_binomial
    from mg10d2_stream import seed_rng,advance
    from runner import seed_for
    from check_mg10d1_development import check_raw
    @nb.njit
    def powers(n,p):return (1-p)**n,np.exp(n*np.log(1-p)),np.exp(n*np.log1p(-p))
    results=[test_law(j,directory,.01,6,.0008) for j in primary_jobs()]
    rounding=[dict(n=j["n"],p=j["p"],numba_power_values=list(powers(j["n"],j["p"]))) for j in primary_jobs() if j["method"]==0]
    supplements=[]
    for i,(n,p) in enumerate(((3,.2),(40,.3),(4200000000000,30/4200000000000),(2100000000000,14000/48000/4200000000000))):
        supplements.append(test_law(dict(id=f"supplement_{i}",n=n,p=p,reps=2**18,method=3,seed=seed_for(10,i,3,5501,1)),directory,.001,1,.02))
    tests=[];raw=[]
    cases=[]
    for i in range(4):
        c=small.fixture(i);cases.append((f"small{i}",dict(initial=c["state"],a0=c["a0"],b0=c["b0"],omega=c["omega"],diff=c["diff"],clamp_a=False,clamp_b=False)))
    for arm in (0,2,3,4):
        for start in range(2):cases.append((f"unchanged_a{arm}_s{start}",new.fixture(arm,start)))
    for i,(name,case) in enumerate(cases):
        seed=seed_for(10,i,0,5501,2);obs=np.array([0.,.005,.01]);tau=1/24000
        a,astat=old.simulate(case,tau,obs,1,seed);b,bstat=new.simulate(case,tau,obs,1,seed)
        assert np.array_equal(a,b) and np.array_equal(astat,bstat),name
        tests.append(dict(name=name,passed=True,seed=seed))
    for arm in range(2):
        for start in range(2):
            case=new.fixture(arm,start);seed=seed_for(10,2*arm+start,0,5501,3);obs=np.array([0.,.005,.01]);tau=1/24000
            expected,stats=new.simulate(case,tau,obs,1,seed)
            state=case["initial"].copy();flux=np.zeros((7,500),dtype=np.int64);hops=np.zeros((6,500),dtype=np.int64);external=np.zeros((2,500),dtype=np.int64)
            seed_rng(seed);actual=np.zeros_like(expected);s=np.zeros((1,3));previous=0.
            for k,target in enumerate(obs):
                v=advance(state,flux,hops,external,case["a0"],case["b0"],1000.,case["diff"],np.ones(7),tau,previous,target,case["clamp_a"],case["clamp_b"])
                previous=target;s[0,:2]+=v[:2];s[0,2]=max(s[0,2],v[2]);actual[0,k]=np.concatenate([state,flux,hops,external])
            assert np.array_equal(actual,expected) and np.array_equal(s,stats)
            path=directory/f"stream_a{arm}_s{start}.npz"
            np.savez_compressed(path,data=actual,initial=case["initial"],a0=case["a0"],b0=case["b0"],clamp_a=case["clamp_a"],clamp_b=case["clamp_b"],omega=1000.,stats=stats,times=obs,seed=seed)
            raw.append(check_raw(path));tests.append(dict(name=f"stream_a{arm}_s{start}",passed=True,seed=seed))
    assert geometric_binomial(0,.1)==0 and geometric_binomial(10,0.)==0 and geometric_binomial(10,1.)==10
    for n,p in ((-1,.1),(4200000000000,.01),(9007199254740993,1e-15)):
        try:finite_source(n,p);raise AssertionError("Domain not rejected")
        except ValueError:pass
    passed=all(r["passed"] for r in results if r["job"]["method"]==1) and all(r["passed"] for r in supplements)
    atomic_json(directory/"qualification.json",dict(passed=passed,primary=results,supplements=supplements,tests=tests,raw=raw,rounding=rounding,
        status="CORRECTED SOURCE QUALIFICATION PASS" if passed else "QUALIFICATION FAILED",legacy_diagnostic="Original failures preserved, not replaced or pooled with revised outcomes"))
    if not passed:raise ValueError("Blocking corrected-source qualification failure")


def worker(job,directory):
    directory=Path(directory);registration=json.loads((BASE/"diagnostic_registration.json").read_text(encoding="utf-8"))
    if hashes()!=registration["source_hashes"]:raise ValueError("Frozen diagnostic source changed")
    if job["kind"]=="diagnostic":qualification(directory)
    else:
        import numpy as np
        import mg10d2 as model
        from runner import seed_for
        case=model.fixture(job["arm"],job["start"]);tau=(1/24000,1/48000)[job["step"]]
        seed=seed_for(10,2*job["arm"]+job["start"],job["step"],5502,0)
        model.simulate(case,tau,np.array([tau]),1,seed_for(10,2*job["arm"]+job["start"],job["step"],5502,0,1))
        clock=time.perf_counter();data,stats=model.simulate(case,tau,np.array([0.,.1]),1,seed);elapsed=time.perf_counter()-clock
        if not model.ledger(data,case["initial"])["pass_check"]:raise ValueError("Timing ledger")
        path=directory/(job["id"]+".npz");np.savez_compressed(path,data=data,initial=case["initial"],a0=case["a0"],b0=case["b0"],omega=1000.,clamp_a=case["clamp_a"],clamp_b=case["clamp_b"],stats=stats,seed=seed,times=np.array([0.,.1]))
        atomic_json(directory/(job["id"]+".json"),dict(status="COMPLETE",job=job,seed=seed,compute_seconds=elapsed,raw=path.name,sha256=digest(path)))


def execute():
    if (BASE/"development/diagnostic/ledger.json").exists():raise ValueError("Already attempted; no automatic diagnostic retry")
    d=run([dict(id="diagnostic",kind="diagnostic")],BASE/"development/diagnostic",Path(__file__).resolve(),ROOT,900,workers=1,timeout=900,poll_seconds=.25)
    if d["stop_reason"] or any(r["returncode"]!=0 for r in d["records"]):raise ValueError("Diagnostic stopped")
    matrix=[dict(id=f"timing_a{a}_s{s}_k{k}",kind="timing",arm=a,start=s,step=k) for a in range(2) for s in range(2) for k in range(2)]
    t=run(matrix,BASE/"development/timing",Path(__file__).resolve(),ROOT,900,workers=8,prior=d["charged_worker_seconds"],timeout=60,poll_seconds=.25)
    if t["stop_reason"] or any(r["returncode"]!=0 for r in t["records"]):raise ValueError("Timing stopped")
    subprocess.run([sys.executable,str(ROOT/"src/check_mg10d2_source.py")],cwd=ROOT,check=True)


if __name__=="__main__":
    parser=argparse.ArgumentParser();parser.add_argument("--register",action="store_true");parser.add_argument("--execute",action="store_true");parser.add_argument("--worker",nargs=2);args=parser.parse_args()
    if args.register:register()
    elif args.execute:execute()
    elif args.worker:
        j=json.loads(args.worker[0]);directory=Path(args.worker[1])
        try:
            worker(j,directory)
            if j["kind"]=="diagnostic":atomic_json(directory/"diagnostic.json",dict(status="COMPLETE",job=j))
        except BaseException:
            atomic_json(directory/(j["id"]+".json"),dict(status="FAILED",job=j,error=traceback.format_exc()));raise
    else:parser.error("Choose register/execute/worker")

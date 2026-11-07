"""MG061 preparation/qualification/timing only; no scientific-production command."""
import os
for key in ("OPENBLAS_NUM_THREADS","MKL_NUM_THREADS","OMP_NUM_THREADS","NUMEXPR_NUM_THREADS"):
    os.environ[key]="1"
from pathlib import Path
import argparse,hashlib,json,platform,sys,time,traceback
from mg10d1_pool import run,atomic_json

ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/"results/MG10P2/preparation"
FILES=("src/mg10p2.py","src/prepare_mg10p2.py","src/check_mg10p2_preparation.py","spec/MG10P2_proposal.md",
       "src/mg10p1.py","src/mg10d2.py","src/mg10d2_source.py","src/mg10s1.py","src/mg10s2.py","src/mg10d1_pool.py","src/runner.py")


def digest(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def hashes():
    d={name:hashlib.sha256((ROOT/name).read_text(encoding="utf-8").replace("\r\n","\n").encode()).hexdigest() for name in FILES}
    d["results/MG003/ref_profiles.npz"]=digest(ROOT/"results/MG003/ref_profiles.npz")
    return d


def candidate_jobs():
    from runner import seed_for
    return [dict(id=f"c{c}_s{s}_k{k}_b{b}",condition=c,start=s,step=k,replicates=list(range(8*b,8*b+8)),
                 seeds=[seed_for(10,2*c+s,k,6103,r) for r in range(8*b,8*b+8)])
            for b in range(16) for k in range(2) for s in range(2) for c in range(2)]


def plan():
    import numpy as np
    import scipy,numba
    from mg10p2 import reference,upper_bound
    if (BASE/"preparation_plan.json").exists():raise ValueError("Already planned/attempted; preserve before any revision")
    matrix=candidate_jobs();seeds=[v for j in matrix for v in j["seeds"]]
    if len(matrix)!=128 or len(set(seeds))!=1024:raise ValueError("Candidate stream denominator")
    _,template=reference()
    atomic_json(BASE/"preparation_plan.json",dict(status="DEVELOPMENT ONLY; DRAFT PRODUCTION NOT REGISTERED",source_hashes=hashes(),
        candidate_jobs=matrix,candidate_paths=1024,candidate_groups=8,candidate_R=128,production_candidate_ceiling_seconds=14400,
        primary_H_groups=4,family_alpha=.05,match_frequency_cap=.05,zero_match_upper=upper_bound(0,128),
        reference_smooth_rms=np.sqrt(np.mean(template**2,axis=1)).tolist(),
        development_ceiling_seconds=900,workers_max=8,threads_per_library=1,
        timing="T.025 in8 arm/start/step classes; cost only; one development path/class plus separate warmup",
        environment=dict(python=sys.version,platform=platform.platform(),numpy=np.__version__,scipy=scipy.__version__,numba=numba.__version__)))


def qualify(directory):
    import numpy as np
    import mg10p2 as model
    from mg10p1 import assess
    tests=[]
    for c in range(2):
        for s in range(2):
            case=model.fixture(c,s)
            assert case["initial"].dtype==np.int64 and np.all(case["initial"]>=0)
            assert np.all(case["a0"]==4200000) and np.all(case["b0"]==135333333)
            assert case["clamp_a"]==(c==0) and case["clamp_b"]==(c==0)
            factor=1. if c==0 else .82
            assert np.all(case["initial"][3]==round(factor*4200000)) and np.all(case["initial"][4]==round(factor*135333333))
            assessment=assess(case["initial"][1]/1000)
            assert assessment["status"]==("REFERENCE_MATCH" if s else "FLAT")
            tests.append(dict(name=f"prepared_state_c{c}_s{s}",passed=True,initial_X_status=assessment["status"]))
    profile,template=model.reference()
    scales=np.array([.75,1.1,1.3]);background=np.array([25.,30.,5.])
    synthetic=np.zeros((1,1,22,500))
    synthetic[0,0,:3]=1000*(background[:,None]+scales[:,None]*profile)
    metrics,_=model.measure(synthetic)
    np.testing.assert_allclose(metrics["projection"][0,0],scales,atol=1e-10)
    np.testing.assert_allclose(metrics["smooth_rms"][0,0],scales*np.sqrt(np.mean(template**2,axis=-1)),atol=1e-10)
    assert np.max(metrics["residual_rms"])<1e-10
    np.savez_compressed(directory/"synthetic_metric_inputs.npz",data=synthetic)
    atomic_json(directory/"synthetic_metric_results.json",dict(metrics={k:v.tolist() for k,v in metrics.items()},expected_projection=scales.tolist()))
    tests.append(dict(name="known_three_field_template_metrics",passed=True))
    bounds=[]
    for n,k in ((128,0),(128,1),(128,2),(128,64),(128,127),(128,128)):
        bounds.append(dict(n=n,k=k,upper=model.upper_bound(k,n),interval=model.binomial_interval(k,n)))
    assert model.upper_bound(0,128)<.05 and model.upper_bound(128,128)==1
    atomic_json(directory/"binomial_bounds.json",bounds)
    for c,s in ((2,0),(0,2)):
        try:model.fixture(c,s);raise AssertionError("Invalid state accepted")
        except ValueError:pass
    tests.append(dict(name="state_domain_rejection",passed=True))
    atomic_json(directory/"qualification_results.json",dict(passed=True,status="PREPARATION CHECKS PASS; no scientific production",tests=tests))


def worker(job,directory):
    directory=Path(directory)
    if hashes()!=json.loads((BASE/"preparation_plan.json").read_text(encoding="utf-8"))["source_hashes"]:raise ValueError("Preparation source changed")
    if job["kind"]=="qualification":qualify(directory)
    else:
        import numpy as np
        from mg10p2 import fixture
        from mg10d2 import simulate,ledger
        from runner import seed_for
        case=fixture(job["condition"],job["start"]);tau=(1/24000,1/48000)[job["step"]]
        cell=2*job["condition"]+job["start"]
        seed=seed_for(10,cell,job["step"],6102,0)
        simulate(case,tau,np.array([tau]),1,seed_for(10,cell,job["step"],6102,0,1))
        clock=time.perf_counter();data,stats=simulate(case,tau,np.array([0.,.025]),1,seed);elapsed=time.perf_counter()-clock
        check=ledger(data,case["initial"])
        if not check["pass_check"] or stats[0,2]>.25+1e-14:raise ValueError("Timing ledger/safety")
        path=directory/(job["id"]+".npz")
        np.savez_compressed(path,data=data,initial=case["initial"],a0=case["a0"],b0=case["b0"],clamp_a=case["clamp_a"],clamp_b=case["clamp_b"],omega=1000.,stats=stats,times=np.array([0.,.025]),seed=seed)
        atomic_json(directory/(job["id"]+".json"),dict(status="COMPLETE",job=job,seed=seed,compute_seconds=elapsed,raw=path.name,sha256=digest(path),ledger=check))
    if job["kind"]=="qualification":atomic_json(directory/(job["id"]+".json"),dict(status="COMPLETE",job=job))


def execute():
    if (BASE/"qualification/ledger.json").exists():raise ValueError("Already attempted; no automatic retry")
    previous=BASE.parent/"preparation_v1/timing/ledger.json"
    prior=json.loads(previous.read_text(encoding="utf-8"))["charged_worker_seconds"]+7 if previous.exists() else 0
    q=run([dict(id="qualification",kind="qualification")],BASE/"qualification",Path(__file__).resolve(),ROOT,900,workers=1,prior=prior,timeout=300,poll_seconds=.25)
    if q["stop_reason"]:raise ValueError("Blocking preparation failure")
    matrix=[dict(id=f"timing_c{c}_s{s}_k{k}",kind="timing",condition=c,start=s,step=k) for c in range(2) for s in range(2) for k in range(2)]
    t=run(matrix,BASE/"timing",Path(__file__).resolve(),ROOT,900,workers=8,prior=q["charged_worker_seconds"],timeout=60,poll_seconds=.25)
    atomic_json(BASE/"execution.json",dict(qualification=q,timing=t))


if __name__=="__main__":
    parser=argparse.ArgumentParser();parser.add_argument("--plan",action="store_true");parser.add_argument("--execute",action="store_true");parser.add_argument("--worker",nargs=2);args=parser.parse_args()
    if args.plan:plan()
    elif args.execute:execute()
    elif args.worker:
        job=json.loads(args.worker[0]);directory=Path(args.worker[1])
        try:worker(job,directory)
        except BaseException:atomic_json(directory/(job["id"]+".json"),dict(status="FAILED",job=job,error=traceback.format_exc()));raise
    else:parser.error("Preparation only: choose plan/execute/worker")

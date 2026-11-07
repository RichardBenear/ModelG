"""MG065 development only; no P3 scientific production command."""
import os
for key in ("OPENBLAS_NUM_THREADS","MKL_NUM_THREADS","OMP_NUM_THREADS","NUMEXPR_NUM_THREADS"):
    os.environ[key]="1"
from pathlib import Path
import argparse,hashlib,json,platform,sys,time,traceback
from mg10d1_pool import run,atomic_json

ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/"results/MG10P3/preparation"
FILES=("src/mg10p3.py","src/prepare_mg10p3.py","src/check_mg10p3_preparation.py","spec/MG10P3_proposal.md",
       "src/mg10p1.py","src/mg10p2.py","src/mg10d2.py","src/mg10d2_stream.py","src/mg10d2_source.py",
       "src/mg10s1.py","src/mg10s2.py","src/mg10d1_pool.py","src/runner.py")


def digest(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def hashes():
    out={p:hashlib.sha256((ROOT/p).read_text(encoding="utf-8").replace("\r\n","\n").encode()).hexdigest() for p in FILES}
    out["results/MG003/ref_profiles.npz"]=digest(ROOT/"results/MG003/ref_profiles.npz")
    return out


def candidate_jobs():
    from runner import seed_for
    return [dict(id=f"s{s}_k{k}_r{r}",start=s,step=k,replicate=r,seed=seed_for(10,s,k,6503,r),warmup_seed=seed_for(10,s,k,6503,r,1)) for r in range(16) for k in range(2) for s in range(2)]


def plan():
    import numpy,scipy,numba
    if (BASE/"preparation_plan.json").exists():raise ValueError("Already planned/attempted; preserve before revision")
    matrix=candidate_jobs();seeds=[j[key] for j in matrix for key in ("seed","warmup_seed")]
    assert len(matrix)==64 and len(set(seeds))==128
    atomic_json(BASE/"preparation_plan.json",dict(status="DRAFT PRODUCTION; DEVELOPMENT ONLY",source_hashes=hashes(),candidate_jobs=matrix,
        paths=64,R=16,groups=4,horizon=25,checkpoint_interval=.25,window=[0,20],production_candidate_ceiling_seconds=86400,
        per_path_timeout_seconds=1800,development_ceiling_seconds=900,workers_max=8,threads_per_library=1,
        timing="Two T1 paths per H/S/step class, full maintained conditions, stages6502; no endpoint matcher selection/inference",
        environment=dict(python=sys.version,platform=platform.platform(),numpy=numpy.__version__,scipy=scipy.__version__,numba=numba.__version__)))


def qualify(directory):
    import numpy as np
    from mg10p2 import reference
    from mg10p3 import episodes,average,summarize,OBS
    from mg10d2 import simulate
    from mg10d2_stream import seed_rng,advance
    from mg10p2 import fixture
    from runner import seed_for
    from check_mg10d1_development import check_raw
    tests=[];episode_inputs=[]
    times=OBS.copy();match=lambda p:dict(status="REFERENCE_MATCH",position=float((p+50)%100-50))
    for kind in ("continuous","wrap","temporary_gap","late20","late_short","motion_split","ambiguity","sampling_gap"):
        t=times.copy();items=[match(49+x if kind=="wrap" else 0.) for x in t]
        if kind=="sampling_gap":keep=(t<10)|(t>=11.5);t=t[keep];items=[match(0.) for _ in t]
        if kind=="temporary_gap":items[40]=dict(status="FLAT",position=None)
        if kind=="late20":items=[v if x>=5 else dict(status="FLAT",position=None) for x,v in zip(t,items)]
        if kind=="late_short":items=[v if x>=5.25 else dict(status="FLAT",position=None) for x,v in zip(t,items)]
        if kind=="motion_split":items=[match(0. if x<10 else 10.) for x in t]
        if kind=="ambiguity":items[40]=dict(status="AMBIGUOUS_MATCH",position=None)
        runs=episodes(items,t);expected=kind in ("continuous","wrap","late20")
        assert any(r["persistent"] for r in runs)==expected
        tests.append(dict(name="episodes_"+kind,passed=True));episode_inputs.append(dict(kind=kind,times=t.tolist(),observations=items,result=runs,expected_persistent=expected))
    atomic_json(directory/"episode_inputs.json",episode_inputs)
    profile,_=reference();times20=np.arange(81)/4;averages=[]
    def shift(f,pos):return np.fft.ifft(np.fft.fft(f)*np.exp(-2j*np.pi*np.fft.fftfreq(500)*pos/.2)).real
    wave=np.cos(2*np.pi*7*(np.arange(500)-250)/500)
    for kind in ("clean","fractional_motion","constant_disturbance","alternating_disturbance"):
        phases=.1*times20 if kind=="fractional_motion" else np.zeros(81)
        fields=np.array([[30+shift(f,pos) for f in profile] for pos in phases])
        if kind in ("constant_disturbance","alternating_disturbance"):
            for i in range(81):fields[i]+=np.array([.02,.03,.01])[:,None]*wave[None]*(1 if kind=="constant_disturbance" else (-1)**i)
        fixed=average(fields);aligned=average(fields,phases)
        path=directory/(kind+".npz");np.savez_compressed(path,fields=fields,phases=phases)
        def jsonable(value):return {k:v.tolist() if isinstance(v,np.ndarray) else v for k,v in value.items()}
        averages.append(dict(kind=kind,raw=path.name,sha256=digest(path),fixed=jsonable(fixed),aligned=jsonable(aligned)))
        if kind in ("clean","fractional_motion"):
            np.testing.assert_allclose(aligned["projection"],[1,1,1],atol=1e-9)
            assert np.max(aligned["residual_rms"])<1e-9 and all(v is None for v in aligned["gain"])
        else:np.testing.assert_allclose(aligned["gain"],[1 if kind=="constant_disturbance" else 81]*3,rtol=1e-7,atol=1e-7)
        if kind=="fractional_motion":assert np.linalg.norm(fixed["profile"]-aligned["profile"])>1e-3
        tests.append(dict(name="average_"+kind,passed=True))
    atomic_json(directory/"average_results.json",averages)
    # A temporary early gap must not be mistaken for the initial seed even
    # when a later newly observed episode itself spans20units.
    fields=np.array([30+profile for _ in times]);fields[4]=30.
    data=np.zeros((101,22,500));data[:,:3]=fields*1000
    summary=summarize(data,times,1)
    assert not summary["track"]["persistent_reference_signal"] and summary["any_observed_episode20"] and summary["aligned"] is None
    np.savez_compressed(directory/"censored_seed_input.npz",data=data,times=times)
    tests.append(dict(name="initial_seed_not_reacquired_despite_later_episode20",passed=True))
    try:summarize(data[:-1],times[:-1],1);raise AssertionError("Partial path accepted")
    except ValueError:pass
    tests.append(dict(name="partial_path_rejected_as_full_estimate",passed=True))
    raw_checks=[]
    for start in range(2):
        for step in range(2):
            case=fixture(0,start);tau=(1/24000,1/48000)[step];obs=np.array([0.,.001,.002]);seed=seed_for(10,start,step,6501,0)
            expected,stats=simulate(case,tau,obs,1,seed)
            state=case["initial"].copy();flux=np.zeros((7,500),dtype=np.int64);hops=np.zeros((6,500),dtype=np.int64);external=np.zeros((2,500),dtype=np.int64)
            actual=np.zeros_like(expected);ss=np.zeros((1,3));seed_rng(seed);previous=0.
            for i,target in enumerate(obs):
                d=advance(state,flux,hops,external,case["a0"],case["b0"],1000.,case["diff"],np.ones(7),tau,previous,target,True,True);previous=target
                ss[0,:2]+=d[:2];ss[0,2]=max(ss[0,2],d[2]);actual[0,i]=np.concatenate([state,flux,hops,external])
            assert np.array_equal(actual,expected) and np.array_equal(ss,stats)
            path=directory/f"stream_s{start}_k{step}.npz";np.savez_compressed(path,data=actual,initial=case["initial"],a0=case["a0"],b0=case["b0"],clamp_a=True,clamp_b=True,omega=1000.,stats=stats,times=obs,seed=seed)
            raw_checks.append(check_raw(path));tests.append(dict(name=f"stream_s{start}_k{step}",passed=True,seed=seed))
    atomic_json(directory/"qualification_results.json",dict(passed=True,status="DEVELOPMENT QUALIFICATION PASS; NO SCIENTIFIC PRODUCTION",tests=tests,raw_checks=raw_checks))


def worker(job,directory):
    directory=Path(directory)
    if hashes()!=json.loads((BASE/"preparation_plan.json").read_text(encoding="utf-8"))["source_hashes"]:raise ValueError("Preparation source changed")
    if job["kind"]=="qualification":qualify(directory);atomic_json(directory/(job["id"]+".json"),dict(status="COMPLETE",job=job));return
    import numpy as np
    from mg10p2 import fixture
    from mg10d2 import simulate,ledger
    from runner import seed_for
    case=fixture(0,job["start"]);tau=(1/24000,1/48000)[job["step"]];seed=seed_for(10,job["start"],job["step"],6502,job["replicate"])
    simulate(case,tau,np.array([tau]),1,seed_for(10,job["start"],job["step"],6502,job["replicate"],1))
    clock=time.perf_counter();data,stats=simulate(case,tau,np.arange(5)/4,1,seed);elapsed=time.perf_counter()-clock
    check=ledger(data,case["initial"])
    if not check["pass_check"] or stats[0,2]>.25+1e-14:raise ValueError("Timing ledger/safety")
    path=directory/(job["id"]+".npz");np.savez_compressed(path,data=data,initial=case["initial"],a0=case["a0"],b0=case["b0"],clamp_a=True,clamp_b=True,omega=1000.,stats=stats,times=np.arange(5)/4,seed=seed)
    atomic_json(directory/(job["id"]+".json"),dict(status="COMPLETE",job=job,seed=seed,compute_seconds=elapsed,raw=path.name,sha256=digest(path)))


def execute():
    if (BASE/"qualification/ledger.json").exists():raise ValueError("Already attempted; no automatic retry")
    q=run([dict(id="qualification",kind="qualification")],BASE/"qualification",Path(__file__).resolve(),ROOT,900,workers=1,timeout=300,poll_seconds=.25)
    if q["stop_reason"]:raise ValueError("Blocking qualification failure")
    matrix=[dict(id=f"timing_s{s}_k{k}_r{r}",kind="timing",start=s,step=k,replicate=r) for r in range(2) for k in range(2) for s in range(2)]
    t=run(matrix,BASE/"timing",Path(__file__).resolve(),ROOT,900,workers=8,prior=q["charged_worker_seconds"],timeout=90,poll_seconds=.25)
    atomic_json(BASE/"execution.json",dict(qualification=q,timing=t))


if __name__=="__main__":
    parser=argparse.ArgumentParser();parser.add_argument("--plan",action="store_true");parser.add_argument("--execute",action="store_true");parser.add_argument("--worker",nargs=2);args=parser.parse_args()
    if args.plan:plan()
    elif args.execute:execute()
    elif args.worker:
        j=json.loads(args.worker[0]);d=Path(args.worker[1])
        try:worker(j,d)
        except BaseException:atomic_json(d/(j["id"]+".json"),dict(status="FAILED",job=j,error=traceback.format_exc()));raise
    else:parser.error("Preparation only: choose plan/execute/worker")

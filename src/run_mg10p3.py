"""MG066 approved P3 production: one uninterrupted stream, quarter checkpoints."""
import os
for key in ("OPENBLAS_NUM_THREADS","MKL_NUM_THREADS","OMP_NUM_THREADS","NUMEXPR_NUM_THREADS"):
    os.environ[key]="1"
from pathlib import Path
import argparse,hashlib,json,platform,shutil,subprocess,sys,time,traceback
from datetime import datetime,timezone
from mg10d1_pool import run,atomic_json

ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/"results/MG10P3"
SOURCES=("src/run_mg10p3.py","src/mg10p3.py","src/evaluate_mg10p3.py","src/check_mg10p3.py","src/mg10p3_oracle.py",
         "src/qualify_mg10p3_pipeline.py","src/prepare_mg10p3.py","src/check_mg10p3_preparation.py",
         "src/mg10p2.py","src/mg10p1.py","src/mg10p2_oracle.py","src/check_mg10p2.py","src/check_mg10d1_development.py",
         "src/mg10d2.py","src/mg10d2_stream.py","src/mg10d2_source.py","src/mg10s1.py","src/mg10s2.py","src/runner.py","src/mg10d1_pool.py",
         "spec/MG10P3_proposal.md","spec/MG10P3_prerun.md","results/MG10P1/registration.json")


def utc():return datetime.now(timezone.utc).isoformat()
def digest(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def hashes():
    d={p:hashlib.sha256((ROOT/p).read_text(encoding="utf-8").replace("\r\n","\n").encode()).hexdigest() for p in SOURCES}
    d["results/MG003/ref_profiles.npz"]=digest(ROOT/"results/MG003/ref_profiles.npz")
    return d


def jobs():
    from prepare_mg10p3 import candidate_jobs
    return candidate_jobs()


def save_npz(path,**arrays):
    import numpy as np
    temp=Path(path).with_suffix(".npz.tmp")
    with temp.open("wb") as handle:np.savez_compressed(handle,**arrays)
    os.replace(temp,path)


def worker(job,directory):
    import numpy as np
    from mg10p2 import fixture
    from mg10p3 import OBS
    from mg10d2 import simulate,ledger
    from mg10d2_stream import seed_rng,advance
    directory=Path(directory);started=time.perf_counter()
    registration=json.loads((BASE/"registration.json").read_text(encoding="utf-8"))
    if hashes()!=registration["source_hashes"] or job not in registration["jobs"]:raise ValueError("Unregistered source/job")
    case=fixture(0,job["start"]);tau=(1/24000,1/48000)[job["step"]]
    simulate(case,tau,np.array([tau]),1,job["warmup_seed"])
    state=case["initial"].copy();flux=np.zeros((7,500),dtype=np.int64);hops=np.zeros((6,500),dtype=np.int64);external=np.zeros((2,500),dtype=np.int64)
    advance(state,flux,hops,external,case["a0"],case["b0"],1000.,case["diff"],np.ones(7),tau,0.,0.,True,True)
    seed_rng(job["seed"])
    data=np.zeros((1,101,22,500),dtype=np.int64);data[0,0,:7]=state;stats=np.zeros((1,3));compute=0.
    path=directory/(job["id"]+".npz")
    for k,target in enumerate(OBS):
        if time.perf_counter()-started>=1800:raise TimeoutError("Path elapsed limit; atomic checkpoints retained")
        if k:
            clock=time.perf_counter();d=advance(state,flux,hops,external,case["a0"],case["b0"],1000.,case["diff"],np.ones(7),tau,OBS[k-1],target,True,True);compute+=time.perf_counter()-clock
            stats[0,:2]+=d[:2];stats[0,2]=max(stats[0,2],d[2]);data[0,k]=np.concatenate([state,flux,hops,external])
        check=ledger(data[:,k:k+1],case["initial"])
        if not check["pass_check"] or stats[0,2]>.25+1e-14:raise ValueError("Checkpoint accounting/probability safety")
        if not np.array_equal(state[3:5],case["initial"][3:5]):raise ValueError("Maintained stocks changed")
        save_npz(path,data=data[:,:k+1],initial=case["initial"],a0=case["a0"],b0=case["b0"],omega=1000.,clamp_a=True,clamp_b=True,
                 stats=stats,times=OBS[:k+1],seed=job["seed"])
        atomic_json(directory/(job["id"]+"_progress.json"),dict(job=job,last_time=float(target),checkpoints=k+1,elapsed_seconds=time.perf_counter()-started,updated=utc()))
    atomic_json(directory/(job["id"]+".json"),dict(status="COMPLETE",job=job,raw=path.name,raw_sha256=digest(path),checkpoints=101,
                 simulation_seconds=compute,elapsed_seconds=time.perf_counter()-started,finished=utc(),ledger=check))


def register():
    import numpy,scipy,numba
    from runner import seed_for
    if (BASE/"registration.json").exists():raise ValueError("Already registered")
    prep=json.loads((BASE/"preparation/independent_check.json").read_text(encoding="utf-8"))
    pipeline=json.loads((BASE/"preparation/pipeline_qualification.json").read_text(encoding="utf-8"))
    if not prep["passed"] or not prep["within_75_percent_reserve"] or not pipeline["passed"] or pipeline["total_preparation_allowance_seconds"]>=900:raise ValueError("Qualification/budget blocks")
    matrix=jobs();allseeds=[j[k] for j in matrix for k in ("seed","warmup_seed")]
    if len(matrix)!=64 or len(set(allseeds))!=128:raise ValueError("Denominator/seed collisions")
    excluded=set()
    def flatten(v):
        if isinstance(v,int):excluded.add(v)
        elif isinstance(v,list):
            for x in v:flatten(x)
        elif isinstance(v,dict):
            for x in v.values():flatten(x)
    def collect(v):
        if isinstance(v,dict):
            for k,x in v.items():
                if k in ("seed","seeds","warmup_seed","warmup_seeds"):flatten(x)
                else:collect(x)
        elif isinstance(v,list):
            for x in v:collect(x)
    for old in ("MG10S1","MG10S2","MG10D1","MG10D2","MG10P1","MG10P2"):collect(json.loads((ROOT/f"results/{old}/registration.json").read_text(encoding="utf-8")))
    for p in (BASE/"preparation").rglob("*.json"):
        if "SYNTHETIC" in str(p) or p.name in ("preparation_plan.json","stream_and_factory_audit.json"):continue
        collect(json.loads(p.read_text(encoding="utf-8")))
    for s in range(2):
        for k in range(2):
            for r in range(2):excluded.add(seed_for(10,s,k,6502,r,1))
    if set(allseeds)&excluded:raise ValueError("Prior/development stream collision")
    forecast={j["id"]:prep["forecast_seconds_per_path"][f"s{j['start']}_k{j['step']}"] for j in matrix}
    if sum(forecast.values())>64800 or shutil.disk_usage(ROOT).free<3_000_000_000:raise ValueError("Forecast/retention")
    atomic_json(BASE/"registration.json",dict(status="REGISTERED BEFORE P3 SCIENTIFIC PRODUCTION",stage="MG066",registered=utc(),
        pi_authorization="i approve (P3 protocol64paths/T25,20-unit sampled persistence and new24-core-hour ceiling)",jobs=matrix,source_hashes=hashes(),
        excluded_seed_count=len(excluded),forecast_seconds=forecast,ceiling_worker_seconds=86400,job_timeout_seconds=1800,workers=8,threads_per_library=1,poll_seconds=.25,
        observations=[i/4 for i in range(101)],expected_paths=64,expected_cell_n=16,average_window=[0,20],average_frames=81,
        inference="Descriptive exact per-cell fractions only, no rate threshold/5% cap/family equivalence; initial S continuity versus any H episode distinct",
        retention="All complete/partial atomic state/reaction/hop/input checkpoints, averages/inputs/config/streams/stats/progress/failures/ledger/evaluation/audit retained and uploaded after inspection",
        external_review="Not requested with PI approval; optional, same-author checks only",available_disk_bytes=shutil.disk_usage(ROOT).free,
        environment=dict(python=sys.version,platform=platform.platform(),numpy=numpy.__version__,scipy=scipy.__version__,numba=numba.__version__)))


def production():
    r=json.loads((BASE/"registration.json").read_text(encoding="utf-8"))
    if hashes()!=r["source_hashes"]:raise ValueError("Frozen source changed")
    subprocess.run(["git","ls-files","--error-unmatch","results/MG10P3/registration.json"],cwd=ROOT,check=True,stdout=subprocess.DEVNULL)
    if (BASE/"production/ledger.json").exists():raise ValueError("Already attempted; no automatic restart")
    atomic_json(BASE/"run_state.json",dict(status="RUNNING",started=utc(),parent_pid=os.getpid(),registration_sha256=digest(BASE/"registration.json")))
    l=run(r["jobs"],BASE/"production",Path(__file__).resolve(),ROOT,86400,workers=8,timeout=1800,forecast=r["forecast_seconds"],poll_seconds=.25)
    atomic_json(BASE/"run_state.json",dict(status="SIMULATION FINISHED; automatic analysis follows",finished=utc(),stop_reason=l["stop_reason"]))
    subprocess.run([sys.executable,str(ROOT/"src/evaluate_mg10p3.py")],cwd=ROOT,check=True)
    subprocess.run([sys.executable,str(ROOT/"src/check_mg10p3.py")],cwd=ROOT,check=True)
    atomic_json(BASE/"run_state.json",dict(status="RESULTS AND INDEPENDENT CHECK READY; PI closure pending",finished=utc(),stop_reason=l["stop_reason"]))


if __name__=="__main__":
    p=argparse.ArgumentParser();p.add_argument("--register",action="store_true");p.add_argument("--production",action="store_true");p.add_argument("--worker",nargs=2);args=p.parse_args()
    if args.worker:
        j=json.loads(args.worker[0]);d=Path(args.worker[1])
        try:worker(j,d)
        except BaseException:atomic_json(d/(j["id"]+".json"),dict(status="FAILED",job=j,error=traceback.format_exc(),finished=utc()));raise
    elif args.register:register()
    elif args.production:
        try:production()
        except BaseException:atomic_json(BASE/"run_state.json",dict(status="DRIVER/ANALYSIS FAILED; data not replaced",error=traceback.format_exc(),finished=utc()));raise
    else:p.error("Choose register/production/worker")

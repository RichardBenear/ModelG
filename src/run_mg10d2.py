"""MG055 registered D2 driver: checkpointed paths, no reseeding or replacement."""
import os
for key in ("OPENBLAS_NUM_THREADS","MKL_NUM_THREADS","OMP_NUM_THREADS","NUMEXPR_NUM_THREADS"):
    os.environ[key]="1"
from pathlib import Path
import argparse, hashlib, json, platform, shutil, subprocess, sys, time, traceback
from datetime import datetime, timezone
from mg10d1_pool import run, atomic_json

ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/"results/MG10D2"
SOURCES=("src/mg10d2.py","src/mg10d2_source.py","src/mg10d2_stream.py","src/run_mg10d2.py",
         "src/evaluate_mg10d2.py","src/check_mg10d2.py","src/qualify_mg10d2.py","src/check_mg10d2_source.py",
         "src/qualify_mg10d2_pipeline.py","src/mg10d1_pool.py","src/bounded_pool.py",
         "src/check_mg10d1_development.py","src/mg10d1.py","src/mg10s1.py","src/mg10s2.py","src/runner.py",
         "spec/MG10D2_proposal.md","spec/MG10D2_prerun.md")


def utc():return datetime.now(timezone.utc).isoformat()
def digest(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def hashes():
    result={name:hashlib.sha256((ROOT/name).read_text(encoding="utf-8").replace("\r\n","\n").encode()).hexdigest() for name in SOURCES}
    result["results/MG003/ref_profiles.npz"]=digest(ROOT/"results/MG003/ref_profiles.npz")
    return result


def jobs():
    from runner import seed_for
    # Interleave arms/starts in each replicate/step block; order is frozen.
    return [dict(id=f"a{a}_s{s}_k{k}_r{r}",arm=a,start=s,step=k,replicate=r,
                 seed=seed_for(10,2*a+s,k,5503,r),horizon=60,tau=(1/24000,1/48000)[k])
            for r in range(4) for k in range(2) for s in range(2) for a in range(2)]


def save_npz(path,**arrays):
    import numpy as np
    path=Path(path);temp=path.with_suffix(".npz.tmp")
    with temp.open("wb") as handle:np.savez_compressed(handle,**arrays)
    os.replace(temp,path)


def worker(job,directory):
    import numpy as np
    import mg10d2 as model
    from mg10d2_stream import seed_rng, advance
    directory=Path(directory);started=time.perf_counter()
    registration=json.loads((BASE/"registration.json").read_text(encoding="utf-8"))
    if hashes()!=registration["source_hashes"]:raise ValueError("Registered source changed")
    if job not in registration["jobs"]:raise ValueError("Unregistered job")
    case=model.fixture(job["arm"],job["start"]);n=case["initial"].shape[1]
    state=case["initial"].copy();flux=np.zeros((7,n),dtype=np.int64)
    hops=np.zeros((6,n),dtype=np.int64);external=np.zeros((2,n),dtype=np.int64)
    # Load compiled specializations without consuming scientific events.
    advance(state,flux,hops,external,case["a0"],case["b0"],case["omega"],case["diff"],np.ones(7),job["tau"],0.,0.,case["clamp_a"],case["clamp_b"])
    seed_rng(job["seed"])
    data=np.zeros((1,job["horizon"]+1,22,n),dtype=np.int64)
    data[0,0,:7]=state;stats=np.zeros((1,3));simulation_seconds=0.
    path=directory/(job["id"]+".npz")
    for checkpoint in range(job["horizon"]+1):
        if time.perf_counter()-started>=3600:raise TimeoutError("Worker elapsed limit; checkpoints retained")
        if checkpoint:
            clock=time.perf_counter()
            delta=advance(state,flux,hops,external,case["a0"],case["b0"],case["omega"],case["diff"],np.ones(7),job["tau"],float(checkpoint-1),float(checkpoint),case["clamp_a"],case["clamp_b"])
            simulation_seconds+=time.perf_counter()-clock
            stats[0,:2]+=delta[:2];stats[0,2]=max(stats[0,2],delta[2])
            data[0,checkpoint,:7]=state;data[0,checkpoint,7:14]=flux
            data[0,checkpoint,14:20]=hops;data[0,checkpoint,20:22]=external
        check=model.ledger(data[:,checkpoint:checkpoint+1],case["initial"])
        if not check["pass_check"] or stats[0,2]>.25+1e-14:raise ValueError("Checkpoint accounting/safety")
        for species,flag in ((3,case["clamp_a"]),(4,case["clamp_b"])):
            if flag and not np.array_equal(state[species],case["initial"][species]):raise ValueError("Maintained stock changed")
        save_npz(path,data=data[:,:checkpoint+1],initial=case["initial"],a0=case["a0"],b0=case["b0"],
                 omega=case["omega"],clamp_a=case["clamp_a"],clamp_b=case["clamp_b"],stats=stats,
                 times=np.arange(checkpoint+1,dtype=float),seed=job["seed"])
        atomic_json(directory/(job["id"]+"_progress.json"),dict(job=job,checkpoint=checkpoint,
                    elapsed_seconds=time.perf_counter()-started,simulation_seconds=simulation_seconds,updated=utc()))
    atomic_json(directory/(job["id"]+".json"),dict(status="COMPLETE",job=job,raw=path.name,raw_sha256=digest(path),
                checkpoints=job["horizon"]+1,ledger=check,stats=stats.tolist(),simulation_seconds=simulation_seconds,
                elapsed_seconds=time.perf_counter()-started,finished=utc(),source_hashes=hashes()))


def register():
    import numpy, scipy, numba
    from runner import seed_for
    if (BASE/"registration.json").exists():raise ValueError("Already registered")
    q=json.loads((BASE/"development/independent_source_check.json").read_text(encoding="utf-8"))
    pq=json.loads((BASE/"development/pipeline_qualification.json").read_text(encoding="utf-8"))
    if not q["passed"] or not q["within_24_hour_75_percent_reserve"] or not pq["passed"]:raise ValueError("Qualification/forecast blocks")
    if pq["combined_development_seconds_with_startup_allowance"]>=900:raise ValueError("Qualification budget")
    if shutil.disk_usage(ROOT).free<2_000_000_000:raise ValueError("Retention space")
    matrix=jobs();seeds=[j["seed"] for j in matrix]
    if len(seeds)!=32 or len(set(seeds))!=32:raise ValueError("Seed denominator")
    excluded=set()
    def flatten(value):
        if isinstance(value,int):excluded.add(value)
        elif isinstance(value,list):
            for item in value:flatten(item)
        elif isinstance(value,dict):
            for item in value.values():flatten(item)
    def collect(value):
        if isinstance(value,dict):
            for key,item in value.items():
                if key in ("seed","seeds"):flatten(item)
                else:collect(item)
        elif isinstance(value,list):
            for item in value:collect(item)
    for previous in ("MG10S1","MG10S2","MG10D1"):
        collect(json.loads((ROOT/f"results/{previous}/registration.json").read_text(encoding="utf-8")))
    collect(json.loads((BASE/"diagnostic_registration.json").read_text(encoding="utf-8")))
    for path in (BASE/"development").rglob("*.json"):
        if "SYNTHETIC" in str(path):continue
        collect(json.loads(path.read_text(encoding="utf-8")))
    for arm in range(2):
        for start in range(2):
            for step in range(2):excluded.add(seed_for(10,2*arm+start,step,5502,0,1))
    if set(seeds)&excluded:raise ValueError("Seed collision")
    forecast={j["id"]:q["forecast_by_class_seconds"][f"a{j['arm']}_s{j['start']}_k{j['step']}"] for j in matrix}
    prior=json.loads((ROOT/"results/MG10D1/production/ledger.json").read_text(encoding="utf-8"))["charged_worker_seconds"]
    if sum(forecast.values())>.75*86400 or prior+86400+900>216000:raise ValueError("Approved aggregate budget exceeded")
    atomic_json(BASE/"registration.json",dict(status="REGISTERED BEFORE CORRECTED TRAJECTORIES",stage="MG055",registered=utc(),
        pi_authorization="looks ok, proceed with the next steps; corrected allocation stays within prior approved60core-hours",
        prior_D1_charge_seconds=prior,new_qualification_ceiling_seconds=900,combined_worst_seconds=prior+86400+900,
        external_review="Not requested with proceed instruction; optional, same-author checks only",jobs=matrix,
        source_hashes=hashes(),excluded_seed_count=len(excluded),production_seed_collisions=0,forecast_seconds=forecast,
        ceiling_worker_seconds=86400,workers=8,threads_per_library=1,job_timeout_seconds=3600,poll_seconds=.25,checkpoints=list(range(61)),
        diagnostic_registration_sha256=digest(BASE/"diagnostic_registration.json"),qualified_diagnostic_source_hashes=json.loads((BASE/"diagnostic_registration.json").read_text(encoding="utf-8"))["source_hashes"],
        qualification="independent_source_check.json plus pipeline_qualification.json; original D1 controls preserved separately",
        retention="All32 full spatial/event/external-input checkpoint archives and progress retained in results/MG10D2/production; repository upload after inspection; D1 raw files unchanged",
        available_disk_bytes=shutil.disk_usage(ROOT).free,
        environment=dict(python=sys.version,platform=platform.platform(),numpy=numpy.__version__,scipy=scipy.__version__,numba=numba.__version__)))


def production():
    registration=json.loads((BASE/"registration.json").read_text(encoding="utf-8"))
    if hashes()!=registration["source_hashes"]:raise ValueError("Frozen source changed")
    # Require registration tracked in the local commit; operator verifies push
    # before invoking. Never silently resume a completed/failed batch.
    subprocess.run(["git","ls-files","--error-unmatch","results/MG10D2/registration.json"],cwd=ROOT,check=True,stdout=subprocess.DEVNULL)
    directory=BASE/"production"
    if (directory/"ledger.json").exists():raise ValueError("Batch already attempted; no automatic resume")
    atomic_json(BASE/"run_state.json",dict(status="RUNNING",started=utc(),parent_pid=os.getpid(),registration_sha256=digest(BASE/"registration.json")))
    ledger=run(registration["jobs"],directory,Path(__file__).resolve(),ROOT,86400,workers=8,timeout=3600,
               forecast=registration["forecast_seconds"],poll_seconds=.25)
    atomic_json(BASE/"run_state.json",dict(status="SIMULATION FINISHED; analysis follows",finished=utc(),stop_reason=ledger["stop_reason"],completed=sum(r.get("status")=="COMPLETE" and r.get("returncode")==0 for r in ledger["records"])))
    # Produce the report and independent audit unattended, including incomplete
    # batches. Commit/upload is done separately, after inspection.
    subprocess.run([sys.executable,str(ROOT/"src/evaluate_mg10d2.py")],cwd=ROOT,check=True)
    subprocess.run([sys.executable,str(ROOT/"src/check_mg10d2.py")],cwd=ROOT,check=True)
    atomic_json(BASE/"run_state.json",dict(status="RESULTS AND INDEPENDENT CHECK READY; PI closure pending",finished=utc(),stop_reason=ledger["stop_reason"]))


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--worker",nargs=2)
    parser.add_argument("--register",action="store_true")
    parser.add_argument("--production",action="store_true")
    args=parser.parse_args()
    if args.worker:
        job=json.loads(args.worker[0]);directory=Path(args.worker[1])
        try:worker(job,directory)
        except BaseException:
            atomic_json(directory/(job["id"]+".json"),dict(status="FAILED",job=job,error=traceback.format_exc(),finished=utc()))
            raise
    elif args.register:register()
    elif args.production:
        try:production()
        except BaseException:
            atomic_json(BASE/"run_state.json",dict(status="DRIVER/ANALYSIS FAILED; inspect retained ledger and raw",error=traceback.format_exc(),finished=utc()))
            raise
    else:parser.error("Choose register, production or worker")

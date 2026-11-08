"""MG070 approved E1; parent then gated-fork phases, aggregate owned budget."""
import os
for key in ("OPENBLAS_NUM_THREADS","MKL_NUM_THREADS","OMP_NUM_THREADS","NUMEXPR_NUM_THREADS"):os.environ[key]="1"
from pathlib import Path
import argparse,hashlib,json,platform,shutil,subprocess,sys,time,traceback
from datetime import datetime,timezone
from mg10d1_pool import run,atomic_json
ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/"results/MG10E1"
SOURCES=("src/run_mg10e1.py","src/mg10e1_stream.py","src/mg10e1.py","src/evaluate_mg10e1.py","src/check_mg10e1.py","src/qualify_mg10e1_pipeline.py",
 "src/prepare_mg10e1.py","src/check_mg10e1_preparation.py","src/mg10d2.py","src/mg10d2_stream.py","src/mg10d2_source.py","src/mg10s1.py","src/mg10s2.py",
 "src/mg10p2.py","src/mg10p1.py","src/check_mg10d1_development.py","src/runner.py","src/mg10d1_pool.py","spec/MG10E1_proposal.md","spec/MG10E1_prerun.md")
def utc():return datetime.now(timezone.utc).isoformat()
def digest(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def hashes():
    values={p:hashlib.sha256((ROOT/p).read_text(encoding="utf-8").replace("\r\n","\n").encode()).hexdigest() for p in SOURCES}
    values["results/MG003/ref_profiles.npz"]=digest(ROOT/"results/MG003/ref_profiles.npz");return values
def save_npz(path,**arrays):
    import numpy as np
    temp=Path(path).with_suffix(".npz.tmp")
    with temp.open("wb") as f:np.savez_compressed(f,**arrays)
    os.replace(temp,path)
def worker(job,directory):
    import numpy as np
    from mg10p2 import fixture
    from mg10e1 import activity,fork_case,ledger
    from mg10e1_stream import seed_rng,advance
    directory=Path(directory);started=time.perf_counter();r=json.loads((BASE/"registration.json").read_text())
    if hashes()!=r["source_hashes"] or job not in r[job["kind"]+"_jobs"]:raise ValueError("Unregistered source/job")
    case=fixture(0,job["start"]);parent_sha=None
    if job["kind"]=="branch":
        parent_path=BASE/"parents"/(job["parent"]+".npz");meta=json.loads(parent_path.with_suffix(".json").read_text())
        if meta["status"]!="COMPLETE" or digest(parent_path)!=meta["raw_sha256"]:raise ValueError("Invalid parent")
        with np.load(parent_path) as raw:
            if not np.array_equal(raw["times"],np.arange(21)/4):raise ValueError("Incomplete parent")
            case=fork_case(case,raw["data"][0,-1,:7])
        parent_sha=meta["raw_sha256"]
    alpha=activity(job["arm"]) if job["kind"]=="branch" else np.ones(500)
    times=np.arange(81 if job["kind"]=="branch" else 21)/4;tau=(1/24000,1/48000)[job["step"]]
    limit=1800 if job["kind"]=="branch" else 600
    state=case["initial"].copy();flux=np.zeros((7,500),dtype=np.int64);hops=np.zeros((6,500),dtype=np.int64);external=np.zeros((2,500),dtype=np.int64)
    seed_rng(job["warmup_seed"]);advance(state.copy(),flux.copy(),hops.copy(),external.copy(),case["a0"],case["b0"],1000.,case["diff"],np.ones(7),tau,0.,tau,True,True,alpha)
    seed_rng(job["seed"]);data=np.zeros((1,len(times),22,500),dtype=np.int64);stats=np.zeros((1,3));compute=0.;path=directory/(job["id"]+".npz")
    for k,target in enumerate(times):
        if time.perf_counter()-started>=limit:raise TimeoutError("Path elapsed limit; prefix retained")
        if k:
            clock=time.perf_counter();d=advance(state,flux,hops,external,case["a0"],case["b0"],1000.,case["diff"],np.ones(7),tau,times[k-1],target,True,True,alpha);compute+=time.perf_counter()-clock
            stats[0,:2]+=d[:2];stats[0,2]=max(stats[0,2],d[2])
        data[0,k]=np.concatenate([state,flux,hops,external]);check=ledger(data[:,k:k+1],case["initial"])
        if not check["pass_check"] or stats[0,2]>.25+1e-14:raise ValueError("Checkpoint safety/accounting")
        if not np.array_equal(state[3:5],case["initial"][3:5]) or np.any(flux[0,alpha==0]) or np.any(external[0,alpha==0]):raise ValueError("Stock/source intervention error")
        save_npz(path,data=data[:,:k+1],initial=case["initial"],a0=case["a0"],b0=case["b0"],omega=1000.,clamp_a=True,clamp_b=True,stats=stats,times=times[:k+1],seed=job["seed"],a_activity=alpha)
        atomic_json(directory/(job["id"]+"_progress.json"),dict(job=job,last_time=float(target),checkpoints=k+1,elapsed_seconds=time.perf_counter()-started,updated=utc()))
    atomic_json(directory/(job["id"]+".json"),dict(status="COMPLETE",job=job,raw=path.name,raw_sha256=digest(path),parent_sha256=parent_sha,checkpoints=len(times),simulation_seconds=compute,elapsed_seconds=time.perf_counter()-started,finished=utc()))
def register():
    import numpy,scipy,numba
    from prepare_mg10e1 import candidate_jobs
    if (BASE/"registration.json").exists():raise ValueError("Already registered")
    prep=json.loads((BASE/"preparation/independent_check.json").read_text());q=json.loads((BASE/"preparation/pipeline_qualification.json").read_text())
    if not prep["passed"] or not prep["within_75_percent_reserve"] or not q["passed"] or q["total_preparation_allowance_seconds"]>=900:raise ValueError("Qualification/budget")
    p,b=candidate_jobs();p=[dict(j,kind="parent") for j in p];b=[dict(j,kind="branch") for j in b]
    if len(p)!=32 or len(b)!=128 or len({j[k] for j in p+b for k in ("seed","warmup_seed")})!=320:raise ValueError("Denominator/seeds")
    forecast={j["id"]:prep["parent_forecast_seconds"][f"s{j['start']}_k{j['step']}"] for j in p}
    forecast.update({j["id"]:prep["branch_forecast_seconds"][f"a{j['arm']}_s{j['start']}_k{j['step']}"] for j in b})
    if sum(forecast.values())>97200 or shutil.disk_usage(ROOT).free<4*1024**3:raise ValueError("Forecast/storage")
    atomic_json(BASE/"registration.json",dict(status="REGISTERED BEFORE E1 SCIENCE",stage="MG070",registered=utc(),pi_authorization="I approve the E1 study. Proceed running it.",
      parent_jobs=p,branch_jobs=b,source_hashes=hashes(),forecast_seconds=forecast,ceiling_worker_seconds=129600,workers=8,threads_per_library=1,poll_seconds=.25,
      parent_timeout_seconds=600,branch_timeout_seconds=1800,blocks_per_step=8,parent_time=5,post_time=20,primary_window=[15,20],
      external_review="Not requested; same-author separate recomputation, not external replication",retention="All complete/partial checkpoint populations/events/hops/inputs, parent/fork hashes, seed/config/masks/progress/errors, ledgers and full profile/contrast analysis retained and uploaded",environment=dict(python=sys.version,platform=platform.platform(),numpy=numpy.__version__,scipy=scipy.__version__,numba=numba.__version__)))
def production():
    r=json.loads((BASE/"registration.json").read_text())
    if hashes()!=r["source_hashes"]:raise ValueError("Frozen source changed")
    subprocess.run(["git","ls-files","--error-unmatch","results/MG10E1/registration.json"],cwd=ROOT,check=True,stdout=subprocess.DEVNULL)
    if (BASE/"parents/ledger.json").exists() or (BASE/"branches/ledger.json").exists():raise ValueError("Already attempted; no automatic restart")
    atomic_json(BASE/"run_state.json",dict(status="RUNNING PARENTS",started=utc(),parent_pid=os.getpid(),registration_sha256=digest(BASE/"registration.json")))
    p=run(r["parent_jobs"],BASE/"parents",Path(__file__).resolve(),ROOT,129600,workers=8,timeout=600,forecast=r["forecast_seconds"],poll_seconds=.25)
    if not p["stop_reason"] and len(p["records"])==32:
        atomic_json(BASE/"run_state.json",dict(status="RUNNING BRANCHES",updated=utc(),parent_pid=os.getpid()))
        b=run(r["branch_jobs"],BASE/"branches",Path(__file__).resolve(),ROOT,129600,workers=8,prior=p["charged_worker_seconds"],timeout=1800,forecast=r["forecast_seconds"],poll_seconds=.25)
    atomic_json(BASE/"run_state.json",dict(status="SIMULATION FINISHED; automatic analysis follows",finished=utc()))
    subprocess.run([sys.executable,str(ROOT/"src/evaluate_mg10e1.py")],cwd=ROOT,check=True)
    subprocess.run([sys.executable,str(ROOT/"src/check_mg10e1.py")],cwd=ROOT,check=True)
    atomic_json(BASE/"run_state.json",dict(status="RESULTS AND SEPARATE CHECK READY; PI closure pending",finished=utc()))
if __name__=="__main__":
    p=argparse.ArgumentParser();p.add_argument("--register",action="store_true");p.add_argument("--production",action="store_true");p.add_argument("--worker",nargs=2);a=p.parse_args()
    if a.worker:
        j=json.loads(a.worker[0]);d=Path(a.worker[1])
        try:worker(j,d)
        except BaseException:atomic_json(d/(j["id"]+".json"),dict(status="FAILED",job=j,error=traceback.format_exc(),finished=utc()));raise
    elif a.register:register()
    elif a.production:
        try:production()
        except BaseException:atomic_json(BASE/"run_state.json",dict(status="DRIVER/ANALYSIS FAILED; data retained",error=traceback.format_exc(),finished=utc()));raise
    else:p.error("Choose register/production/worker")

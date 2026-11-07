"""MG062 approved P2 run, individually seeded paths in bounded atomic chunks."""
import os
for key in ("OPENBLAS_NUM_THREADS","MKL_NUM_THREADS","OMP_NUM_THREADS","NUMEXPR_NUM_THREADS"):
    os.environ[key]="1"
from pathlib import Path
import argparse,hashlib,json,platform,shutil,subprocess,sys,time,traceback
from datetime import datetime,timezone
from mg10d1_pool import run,atomic_json

ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/"results/MG10P2"
SOURCES=("src/mg10p2.py","src/run_mg10p2.py","src/evaluate_mg10p2.py","src/check_mg10p2.py","src/mg10p2_oracle.py",
         "src/qualify_mg10p2_pipeline.py","src/prepare_mg10p2.py","src/check_mg10p2_preparation.py",
         "src/mg10p1.py","src/mg10d2.py","src/mg10d2_source.py","src/mg10d1_pool.py","src/bounded_pool.py",
         "src/check_mg10d1_development.py","src/mg10s1.py","src/mg10s2.py","src/runner.py",
         "spec/MG10P2_proposal.md","spec/MG10P2_prerun.md","results/MG10P1/registration.json")


def utc():return datetime.now(timezone.utc).isoformat()
def digest(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def hashes():
    result={p:hashlib.sha256((ROOT/p).read_text(encoding="utf-8").replace("\r\n","\n").encode()).hexdigest() for p in SOURCES}
    result["results/MG003/ref_profiles.npz"]=digest(ROOT/"results/MG003/ref_profiles.npz")
    return result


def jobs():
    from prepare_mg10p2 import candidate_jobs
    return candidate_jobs()


def save_npz(path,**arrays):
    import numpy as np
    temp=Path(path).with_suffix(".npz.tmp")
    with temp.open("wb") as handle:np.savez_compressed(handle,**arrays)
    os.replace(temp,path)


def worker(job,directory):
    import numpy as np
    from mg10p2 import fixture,OBS
    from mg10d2 import simulate,ledger
    from runner import seed_for
    directory=Path(directory);started=time.perf_counter()
    registration=json.loads((BASE/"registration.json").read_text(encoding="utf-8"))
    if hashes()!=registration["source_hashes"] or job not in registration["jobs"]:raise ValueError("Unregistered source/job")
    case=fixture(job["condition"],job["start"]);tau=(1/24000,1/48000)[job["step"]]
    simulate(case,tau,np.array([tau]),1,seed_for(10,2*job["condition"]+job["start"],job["step"],6103,job["replicates"][0],1))
    data=np.zeros((8,4,22,500),dtype=np.int64);stats=np.zeros((8,3));compute=0.
    path=directory/(job["id"]+".npz")
    for index,seed in enumerate(job["seeds"]):
        if time.perf_counter()-started>=180:raise TimeoutError("Worker elapsed limit; complete prefix retained")
        clock=time.perf_counter();out,s=simulate(case,tau,OBS,1,seed);compute+=time.perf_counter()-clock
        check=ledger(out,case["initial"])
        if not check["pass_check"] or s[0,2]>.25+1e-14:raise ValueError("Integer/probability safety check")
        for species,flag in ((3,case["clamp_a"]),(4,case["clamp_b"])):
            if flag and not np.all(out[:,:,species]==case["initial"][species]):raise ValueError("Maintained stock changed")
        data[index]=out[0];stats[index]=s[0]
        save_npz(path,data=data[:index+1],initial=case["initial"],a0=case["a0"],b0=case["b0"],omega=1000.,
                 clamp_a=case["clamp_a"],clamp_b=case["clamp_b"],stats=stats[:index+1],times=OBS,
                 seeds=np.array(job["seeds"][:index+1],dtype=np.uint32),replicates=np.array(job["replicates"][:index+1]))
        atomic_json(directory/(job["id"]+"_progress.json"),dict(job=job,completed_paths=index+1,elapsed_seconds=time.perf_counter()-started,updated=utc()))
    atomic_json(directory/(job["id"]+".json"),dict(status="COMPLETE",job=job,raw=path.name,raw_sha256=digest(path),
        paths=8,compute_seconds=compute,elapsed_seconds=time.perf_counter()-started,finished=utc(),ledger=check))


def register():
    import numpy,scipy,numba
    if (BASE/"registration.json").exists():raise ValueError("Already registered")
    prep=json.loads((BASE/"preparation/independent_check.json").read_text(encoding="utf-8"))
    pipeline=json.loads((BASE/"preparation/pipeline_qualification.json").read_text(encoding="utf-8"))
    if not prep["passed"] or not prep["forecast_within_75_percent_ceiling"] or not pipeline["passed"]:raise ValueError("Qualification/forecast blocks")
    if pipeline["total_preparation_allowance_seconds"]>=900:raise ValueError("Preparation budget")
    matrix=jobs();seeds=[s for j in matrix for s in j["seeds"]]
    if len(matrix)!=128 or len(seeds)!=1024 or len(set(seeds))!=1024:raise ValueError("Denominator/streams")
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
                if k in ("seed","seeds"):flatten(x)
                else:collect(x)
        elif isinstance(v,list):
            for x in v:collect(x)
    for name in ("MG10S1","MG10S2","MG10D1","MG10D2","MG10P1"):collect(json.loads((ROOT/f"results/{name}/registration.json").read_text(encoding="utf-8")))
    for p in (BASE/"preparation").rglob("*.json"):
        if "SYNTHETIC" in str(p) or p.name in ("preparation_plan.json","stream_and_replay_audit.json"):continue
        collect(json.loads(p.read_text(encoding="utf-8")))
    from runner import seed_for
    for c in range(2):
        for start in range(2):
            for step in range(2):excluded.add(seed_for(10,2*c+start,step,6102,0,1))
    warm={j["id"]:seed_for(10,2*j["condition"]+j["start"],j["step"],6103,j["replicates"][0],1) for j in matrix}
    if set(seeds)&excluded or len(set(warm.values()))!=128 or set(warm.values())&(set(seeds)|excluded):raise ValueError("Prior/development/warmup seed collision")
    forecast={j["id"]:prep["forecast_by_cell_seconds"][f"c{j['condition']}_s{j['start']}_k{j['step']}"]/16 for j in matrix}
    if sum(forecast.values())>10800 or shutil.disk_usage(ROOT).free<2_000_000_000:raise ValueError("Forecast/retention")
    atomic_json(BASE/"registration.json",dict(status="REGISTERED BEFORE P2 SCIENTIFIC PRODUCTION",stage="MG062",registered=utc(),
        pi_authorization="I approve the plan. Proceed. (P2 approved1024-path protocol,5% frequency criterion and4-core-hour budget)",
        jobs=matrix,source_hashes=hashes(),unique_seeds=1024,warmup_seeds=warm,excluded_seed_count=len(excluded),forecast_seconds=forecast,
        ceiling_worker_seconds=14400,job_timeout_seconds=180,workers=8,threads_per_library=1,poll_seconds=.25,
        observations=[0.,.05,.1,.25],expected_paths=1024,expected_jobs=128,expected_cell_n=128,
        family_alpha=.05,primary_H_groups=4,one_sided_group_alpha=.0125,match_frequency_cap=.05,
        primary_scope="Endpoint .25 H reference-match frequency only; S retention/earlier times/GXY metrics descriptive, not physical false-particle rate",
        partial_data="Complete individual T.25 paths saved as atomic chunk prefixes retained/included; any incomplete batch cannot score the primary family",
        external_review="Not requested with PI approval; optional, same-author qualification/checks only",
        retention="All full population/reaction/hop/input checkpoints for complete paths, partial prefixes/progress/failures, seeds/config/ledger/evaluation/audit retained and uploaded after inspection",
        available_disk_bytes=shutil.disk_usage(ROOT).free,
        environment=dict(python=sys.version,platform=platform.platform(),numpy=numpy.__version__,scipy=scipy.__version__,numba=numba.__version__)))


def production():
    registration=json.loads((BASE/"registration.json").read_text(encoding="utf-8"))
    if hashes()!=registration["source_hashes"]:raise ValueError("Frozen source changed")
    subprocess.run(["git","ls-files","--error-unmatch","results/MG10P2/registration.json"],cwd=ROOT,check=True,stdout=subprocess.DEVNULL)
    if (BASE/"production/ledger.json").exists():raise ValueError("Already attempted; no automatic restart")
    atomic_json(BASE/"run_state.json",dict(status="RUNNING",started=utc(),parent_pid=os.getpid(),registration_sha256=digest(BASE/"registration.json")))
    l=run(registration["jobs"],BASE/"production",Path(__file__).resolve(),ROOT,14400,workers=8,timeout=180,forecast=registration["forecast_seconds"],poll_seconds=.25)
    atomic_json(BASE/"run_state.json",dict(status="SIMULATION FINISHED; automatic analysis follows",finished=utc(),stop_reason=l["stop_reason"]))
    subprocess.run([sys.executable,str(ROOT/"src/evaluate_mg10p2.py")],cwd=ROOT,check=True)
    subprocess.run([sys.executable,str(ROOT/"src/check_mg10p2.py")],cwd=ROOT,check=True)
    atomic_json(BASE/"run_state.json",dict(status="RESULTS AND INDEPENDENT CHECK READY; PI closure pending",finished=utc(),stop_reason=l["stop_reason"]))


if __name__=="__main__":
    parser=argparse.ArgumentParser();parser.add_argument("--register",action="store_true");parser.add_argument("--production",action="store_true");parser.add_argument("--worker",nargs=2);args=parser.parse_args()
    if args.worker:
        j=json.loads(args.worker[0]);d=Path(args.worker[1])
        try:worker(j,d)
        except BaseException:atomic_json(d/(j["id"]+".json"),dict(status="FAILED",job=j,error=traceback.format_exc(),finished=utc()));raise
    elif args.register:register()
    elif args.production:
        try:production()
        except BaseException:atomic_json(BASE/"run_state.json",dict(status="DRIVER/ANALYSIS FAILED; retained data not replaced",error=traceback.format_exc(),finished=utc()));raise
    else:parser.error("Choose register/production/worker")

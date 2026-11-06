"""MG044 local owned-process controller and registered prerequisite sequence."""
import os
for key in ("OPENBLAS_NUM_THREADS","MKL_NUM_THREADS","OMP_NUM_THREADS","NUMEXPR_NUM_THREADS"): os.environ[key]="1"
from pathlib import Path
import argparse,hashlib,json,platform,subprocess,sys,time,traceback

ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/"results/MG10C1"


def write(path,obj):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    temp=path.with_suffix(path.suffix+".tmp")
    temp.write_text(json.dumps(obj,indent=2)+"\n",encoding="utf-8");os.replace(temp,path)


def digest(path): return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def hashes():
    names=("src/mg10c1.py","src/run_mg10c1.py","src/evaluate_mg10c1.py","src/check_mg10c1.py",
           "src/mg10_reservoir.py","src/check_mg10.py","src/measure.py","src/runner.py")
    return {name:digest(ROOT/name) for name in names}


def manifest():
    jobs=[dict(id=kind,kind=kind,T=300) for kind in ("reduced_primary","reduced_reference","linear_primary","linear_reference")]
    jobs += [dict(id=f"full_{label}_{method}",kind="full",method=method,epsilon=epsilon,
                  mode="SELECT_FROM_PRIMARY_LINEAR",T=300) for label,epsilon in (("zero",0.),("eps6",1e-6),("eps8",1e-8)) for method in ("E2","Oplus")]
    return jobs


def worker(job,directory):
    began=time.perf_counter()
    import numpy as np
    import scipy
    import mg10c1 as m
    directory=Path(directory)
    partial=directory/(job["id"]+".checkpoint.npz")
    def save_partial(data):
        temp=partial.with_suffix(".tmp")
        with temp.open("wb") as f: np.savez_compressed(f,**data)
        os.replace(temp,partial)
    start_compute=time.perf_counter()
    try:
        if job["kind"].startswith("reduced_"): data=m.reduced(job)
        elif job["kind"].startswith("linear_"): data=m.linear(job,directory)
        elif job["kind"]=="full": data=m.full(job,save_partial)
        else: raise ValueError(job["kind"])
        compute=time.perf_counter()-start_compute
        save_partial(data)
        final=directory/(job["id"]+".npz");os.replace(partial,final)
        audit={"pass":True}
        if job["kind"]=="full":
            audit=m.old.accounting(data["times"],data["state"],dict(capacity=1,mode="FINITE",method="E2" if job["method"]=="E2" else "O"))
        elif job["kind"].startswith("reduced_"):
            error=float(np.max(abs(1-data["state"][:,3]-np.exp(-data["times"]/300))))
            audit={"A_oracle_error":error,"pass":error<=1e-10}
        record=dict(status="COMPLETE" if audit["pass"] else "QUALIFICATION FAILED",job=job,
                    compute_seconds=compute,audit=audit,sha256=digest(final))
    except Exception:
        record=dict(status="FAILED",job=job,compute_seconds=time.perf_counter()-start_compute,traceback=traceback.format_exc())
    record.update(worker_elapsed_seconds=time.perf_counter()-began,python=sys.version,numpy=np.__version__,scipy=scipy.__version__,
                  platform=platform.platform(),threads_per_library=1)
    write(directory/(job["id"]+".json"),record)
    return 0 if record["status"]=="COMPLETE" else 2


def terminate(process):
    process.terminate()
    try: process.wait(timeout=1)
    except subprocess.TimeoutExpired: process.kill();process.wait()


def pool(jobs,directory,ceiling,workers=1,prior=0,timeout=600,forecast=None,future_cost=0):
    directory=Path(directory);directory.mkdir(parents=True,exist_ok=True)
    lock=directory/"driver.lock"
    with lock.open("x") as f:f.write(str(os.getpid()))
    file=directory/"ledger.json"
    ledger=json.loads(file.read_text(encoding="utf-8")) if file.exists() else dict(records=[],ceiling=ceiling,prior_charge=prior,poll_seconds=.5)
    active={};stop=None
    pending=[j for j in jobs if j["id"] not in {r["id"] for r in ledger["records"]}]
    if ledger["ceiling"]!=ceiling or ledger["prior_charge"]!=prior:
        lock.unlink();raise ValueError("Budget changed")
    def charged():return prior+sum(r["charged_worker_seconds"] for r in ledger["records"])
    def save():
        ledger.update(charged_worker_seconds=charged(),inflight_worker_seconds=sum(time.perf_counter()-a["start"] for a in active.values()),
                      pending_ids=[j["id"] for j in pending],active_ids=list(active),stop_reason=stop,workers=workers,threads_per_library=1)
        write(file,ledger)
    try:
        while pending or active:
            now=time.perf_counter();elapsed=sum(now-a["start"] for a in active.values())
            if charged()+elapsed>=ceiling:stop="BUDGET-SHORT"
            if forecast and pending and len(active)<workers and not stop:
                prediction=sum(forecast[j["id"]] for j in pending)
                prediction+=sum(max(0,forecast[k]-(now-a["start"])) for k,a in active.items())
                if prediction+future_cost>.75*(ceiling-charged()-elapsed):stop="BUDGET-SHORT: FORECAST"
            while pending and len(active)<workers and not stop:
                job=pending.pop(0)
                log=(directory/(job["id"]+".log")).open("w",encoding="utf-8")
                process=subprocess.Popen([sys.executable,str(Path(__file__)),"--worker",json.dumps(job),str(directory)],
                                         cwd=ROOT,env=os.environ.copy(),stdout=log,stderr=subprocess.STDOUT)
                active[job["id"]]=dict(job=job,process=process,log=log,start=time.perf_counter())
            for key,a in list(active.items()):
                process=a["process"]
                reason=stop or ("TIMEOUT" if time.perf_counter()-a["start"]>=timeout else None)
                if reason and process.poll() is None:terminate(process)
                if process.poll() is not None:
                    meta=directory/(key+".json")
                    record=json.loads(meta.read_text(encoding="utf-8")) if meta.exists() else dict(status=reason or "FAILED",job=a["job"])
                    record.update(id=key,charged_worker_seconds=time.perf_counter()-a["start"],returncode=process.returncode)
                    if reason:record["termination_reason"]=reason
                    if not meta.exists():write(meta,record)
                    ledger["records"].append(record)
                    a["log"].close();del active[key]
                    if record["status"] in ("FAILED","QUALIFICATION FAILED"):stop=record["status"]
            save()
            if stop and not active:break
            if active:time.sleep(.5)
    except BaseException:
        stop="INTERRUPTED"
        for key,a in active.items():
            if a["process"].poll() is None:terminate(a["process"])
            ledger["records"].append(dict(id=key,status="INTERRUPTED",job=a["job"],charged_worker_seconds=time.perf_counter()-a["start"]))
            a["log"].close()
        active.clear();save();raise
    finally:lock.unlink(missing_ok=True)
    save();return ledger


def benchmark():
    tests=json.loads((BASE/"qualification/tests.json").read_text())
    if tests["returncode"]:raise RuntimeError("Qualification failed")
    prior=6.1446315+1.761948+tests["elapsed_worker_seconds"]
    directory=BASE/"timing"
    jobs=[dict(id=k,kind=k,T=5) for k in ("reduced_primary","reduced_reference","linear_primary","linear_reference")]
    jobs += [dict(id="development_"+method,kind="full",method=method,epsilon=1e-5,mode=7,T=5) for method in ("E2","Oplus")]
    for phase in (jobs[:2],jobs[2:4],jobs[4:]):
        ledger=pool(phase,directory,360,2,prior=prior)
        if any(j["id"] not in {r["id"] for r in ledger["records"] if r["status"]=="COMPLETE"} for j in phase):raise RuntimeError("Development timing incomplete")
    measured={r["id"]:r for r in ledger["records"]}
    overhead=max(5.,2*max(r["charged_worker_seconds"]-r["compute_seconds"] for r in measured.values()))
    cost={}
    for job in manifest():
        key=job["id"] if job["kind"]!="full" else "development_"+job["method"]
        cost[job["id"]]=2*60*measured[key]["compute_seconds"]+overhead
    result=dict(projected_worker_seconds=sum(cost.values()),per_job=cost,overhead_allowance=overhead,
                ceiling=1800,reserve=.25,fits=sum(cost.values())<=1350 and max(cost.values())<=600,
                qualification_timing_seconds=ledger["charged_worker_seconds"],
                formula="2x compute-only T5->T300 linear extrapolation per solver class + max(5s,2x measured startup/serialization/poll overhead) per job")
    write(BASE/"timing/projection.json",result);print(json.dumps(result,indent=2))


def production():
    import numpy as np
    registration=json.loads((BASE/"registration.json").read_text())
    if registration["status"]!="REGISTERED" or registration["source_hashes"]!=hashes():raise RuntimeError("Registration/code mismatch")
    if digest(ROOT/"spec/MG10C1_proposal.md")!=registration["spec_sha256"]:raise RuntimeError("Specification changed")
    projection=json.loads((BASE/"timing/projection.json").read_text())
    if not projection["fits"]:raise RuntimeError("Forecast failed")
    jobs=manifest();directory=BASE/"production"
    workers=min(8,max(1,int(os.environ.get("MODELG_WORKERS","2"))))
    for phase,future in ((jobs[:2],jobs[2:]),(jobs[2:4],jobs[4:])):
        future_cost=sum(projection["per_job"][j["id"]] for j in future)
        ledger=pool(phase,directory,1800,workers,forecast=projection["per_job"],future_cost=future_cost)
        if any(j["id"] not in {r["id"] for r in ledger["records"] if r["status"]=="COMPLETE"} for j in phase):return
    with np.load(directory/"linear_primary.npz") as raw:
        peaks=np.max(abs(raw["response"][:,1]),axis=0)
        mode=int(raw["modes"][np.argmax(peaks)])
    write(directory/"selection.json",dict(mode=mode,rule="Smallest mode maximizing primary integer-checkpoint peak X response",peak_gain=float(np.max(peaks))))
    for job in jobs[4:]:job["mode"]=mode
    write(directory/"resolved_manifest.json",dict(jobs=jobs))
    ledger=pool(jobs[4:],directory,1800,workers,forecast=projection["per_job"])
    print(dict(records=len(ledger["records"]),charged_worker_seconds=ledger["charged_worker_seconds"],stop=ledger["stop_reason"]))


if __name__=="__main__":
    p=argparse.ArgumentParser();p.add_argument("--worker",nargs=2);p.add_argument("--benchmark",action="store_true");p.add_argument("--production",action="store_true")
    a=p.parse_args()
    if a.worker:sys.exit(worker(json.loads(a.worker[0]),a.worker[1]))
    elif a.benchmark:benchmark()
    elif a.production:production()
    else:p.error("Choose an action")

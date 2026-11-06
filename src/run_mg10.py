"""MG042 bounded local process pool, qualification timings and fixed production.

All processes are owned subprocesses; no affinity or global settings change.
No automatic retry, replacement, capacity extension, or seed selection.
"""
import os
for _key in ("OPENBLAS_NUM_THREADS","MKL_NUM_THREADS","OMP_NUM_THREADS","NUMEXPR_NUM_THREADS"):
    os.environ[_key]="1"
import argparse
import hashlib
import json
from pathlib import Path
import platform
import subprocess
import sys
import time
import traceback

ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/"results/MG10"


def write_json(path, obj):
    path=Path(path)
    temporary=path.with_suffix(path.suffix+".tmp")
    temporary.write_text(json.dumps(obj,indent=2)+"\n",encoding="utf-8")
    os.replace(temporary,path)


def sha(path): return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def source_hashes():
    names=("src/mg10_reservoir.py","src/run_mg10.py","src/evaluate_mg10.py",
           "src/measure.py","src/runner.py","src/lattice.py")
    return {name:sha(ROOT/name) for name in names}


def jobs():
    from mg10_reservoir import CAPACITIES
    from evaluate_mg10 import run_id
    return [dict(id=run_id(mode,start,method,c),mode=mode,start=start,method=method,capacity=c,T=300,n=500)
            for method in ("E0","E1","E2","O") for start in ("H","S")
            for mode,c in [("CHEMOSTAT",1)]+[("FINITE",c) for c in CAPACITIES]]


def budget_state(charged, inflight, ceiling, elapsed, timeout):
    if charged+sum(inflight)>=ceiling: return "BUDGET-SHORT"
    if elapsed>=timeout: return "TIMEOUT"
    return "RUNNING"


def worker(job, directory):
    began=time.perf_counter()
    import numpy as np
    from mg10_reservoir import trajectory, accounting
    directory=Path(directory)
    file=directory/(job["id"]+".npz")
    partial=directory/(job["id"]+".checkpoint.npz")
    def checkpoint(times,data):
        temporary=partial.with_suffix(".tmp")
        with temporary.open("wb") as handle:
            np.savez_compressed(handle,times=times,state=data)
        os.replace(temporary,partial)
    try:
        times,data=trajectory(job,checkpoint)
        checkpoint(times,data)
        os.replace(partial,file)
        audit=accounting(times,data,job)
        record={"status":"COMPLETE" if audit["pass"] else "QUALIFICATION FAILED",
                "job":job,"audit":audit,"sha256":sha(file),"last_checkpoint":float(times[-1])}
    except Exception:
        record={"status":"FAILED","job":job,"traceback":traceback.format_exc()}
    record.update(elapsed_worker_seconds=time.perf_counter()-began,platform=platform.platform(),
                  python=sys.version,threads_per_library=1)
    write_json(directory/(job["id"]+".json"),record)
    return 0 if record["status"]=="COMPLETE" else 2


def stop_process(process):
    process.terminate()
    try: process.wait(timeout=1)
    except subprocess.TimeoutExpired:
        process.kill();process.wait()


def execute(matrix,directory,ceiling,timeout,workers,previous_charge=0,projection=None):
    directory=Path(directory);directory.mkdir(parents=True,exist_ok=True)
    lock=directory/"driver.lock"
    with lock.open("x",encoding="utf-8") as f: f.write(str(os.getpid()))
    ledger_file=directory/"ledger.json"
    ledger={"ceiling_worker_seconds":ceiling,"previous_charge":previous_charge,
            "records":[],"workers":workers,"poll_seconds":.5,"projection":projection,
            "threads_per_library":1,"started_utc":time.strftime("%Y-%m-%dT%H:%M:%SZ",time.gmtime())}
    if ledger_file.exists():
        ledger=json.loads(ledger_file.read_text(encoding="utf-8"))
        if ledger["ceiling_worker_seconds"]!=ceiling: raise ValueError("Ceiling changed")
    done={r["id"] for r in ledger["records"]}
    pending=[j for j in matrix if j["id"] not in done]
    active={}
    stop_all=None
    def charged(): return previous_charge+sum(r["charged_worker_seconds"] for r in ledger["records"])
    def save():
        now=time.perf_counter()
        ledger["charged_worker_seconds"]=charged()
        ledger["inflight_worker_seconds"]=sum(now-a["start"] for a in active.values())
        ledger["pending_ids"]=[j["id"] for j in pending]
        ledger["active_ids"]=list(active)
        ledger["stop_reason"]=stop_all
        write_json(ledger_file,ledger)
    try:
        while pending or active:
            now=time.perf_counter()
            if charged()+sum(now-a["start"] for a in active.values())>=ceiling:
                stop_all="BUDGET-SHORT"
            # Recheck the measured forecast before each dispatch wave, including
            # elapsed work on active jobs. Reserve is 25% of remaining budget.
            if projection and pending and len(active)<workers and not stop_all:
                inflight=sum(now-a["start"] for a in active.values())
                remaining=sum(projection[j["method"]] for j in pending)
                remaining+=sum(max(0,projection[a["job"]["method"]]-(now-a["start"])) for a in active.values())
                if remaining > .75*(ceiling-charged()-inflight): stop_all="BUDGET-SHORT: FORECAST"
            while pending and len(active)<workers and not stop_all:
                job=pending.pop(0)
                log=(directory/(job["id"]+".log")).open("w",encoding="utf-8")
                process=subprocess.Popen([sys.executable,str(Path(__file__)),"--worker",json.dumps(job),str(directory)],
                                         cwd=ROOT,stdout=log,stderr=subprocess.STDOUT,env=os.environ.copy())
                active[job["id"]]={"job":job,"process":process,"start":time.perf_counter(),"log":log}
            for key,a in list(active.items()):
                elapsed=time.perf_counter()-a["start"]
                process=a["process"]
                state=budget_state(charged(),[time.perf_counter()-v["start"] for v in active.values()],ceiling,elapsed,timeout)
                termination=stop_all or (state if state!="RUNNING" else None)
                if termination and process.poll() is None:
                    stop_process(process)
                    meta=directory/(key+".json")
                    if not meta.exists(): write_json(meta,{"status":termination,"job":a["job"]})
                if process.poll() is not None:
                    meta=directory/(key+".json")
                    record=json.loads(meta.read_text(encoding="utf-8")) if meta.exists() else {"status":"FAILED","job":a["job"]}
                    record.update(id=key,charged_worker_seconds=time.perf_counter()-a["start"],returncode=process.returncode)
                    if termination: record["termination_reason"]=termination
                    ledger["records"].append(record)
                    if record["status"]=="QUALIFICATION FAILED": stop_all="QUALIFICATION FAILED"
                    a["log"].close();del active[key]
            save()
            if stop_all and not active: break
            if active: time.sleep(.5)
    except BaseException:
        stop_all="INTERRUPTED"
        for key,a in list(active.items()):
            if a["process"].poll() is None: stop_process(a["process"])
            ledger["records"].append({"id":key,"status":"INTERRUPTED","job":a["job"],
                                      "charged_worker_seconds":time.perf_counter()-a["start"]})
            a["log"].close()
        active.clear();save()
        raise
    finally:
        lock.unlink(missing_ok=True)
    save()
    return ledger


def benchmark():
    test=json.loads((BASE/"qualification/full_tests.json").read_text(encoding="utf-8"))
    if test["returncode"]!=0: raise RuntimeError("Qualification tests failed")
    prior=5.7361898+14.19872620003298+2.2223602+test["elapsed_worker_seconds"]
    matrix=[dict(id="development_"+method,mode="FINITE",start="S",method=method,capacity=1,T=10,n=500)
            for method in ("E2","O")]
    directory=BASE/"timing"
    write_json(directory.parent/"qualification/timing_manifest.json",{"jobs":matrix})
    ledger=execute(matrix,directory,900,600,1,previous_charge=prior)
    if len(ledger["records"])!=2 or any(r["status"]!="COMPLETE" for r in ledger["records"]):
        raise RuntimeError("Timing/qualification incomplete")
    measured={r["job"]["method"]:r["charged_worker_seconds"] for r in ledger["records"]}
    # Conservative 2x extrapolation, plus startup/serialization allowance.
    forecast={"E0":2*measured["E2"]*30/4+5,"E1":2*measured["E2"]*30/2+5,
              "E2":2*measured["E2"]*30+5,"O":2*measured["O"]*30+5}
    total=16*sum(forecast.values())
    projection={"measured_T10_worker_seconds":measured,"per_T300_job_seconds":forecast,
                "total_projected_worker_seconds":total,"production_ceiling":14400,
                "reserve_fraction":.25,"fits_with_reserve":total<=10800 and max(forecast.values())<=600,
                "qualification_charged_worker_seconds":ledger["charged_worker_seconds"],
                "model":"2x linear-in-horizon extrapolation from finite C=1 S; Euler step-count scaling; +5 s/job. Same upper estimate for every start/capacity, capacity does not multiply work.",
                "source_hashes":source_hashes()}
    write_json(BASE/"timing/projection.json",projection)
    print(json.dumps(projection,indent=2))
    return projection


def production():
    directory=BASE/"production"
    reg=json.loads((BASE/"registration.json").read_text(encoding="utf-8"))
    if reg["status"]!="REGISTERED" or reg["source_hashes"]!=source_hashes():
        raise RuntimeError("Missing registration or changed implementation")
    projection=json.loads((BASE/"timing/projection.json").read_text(encoding="utf-8"))
    if not projection["fits_with_reserve"]: raise RuntimeError("Forecast does not fit")
    manifest=json.loads((directory/"manifest.json").read_text(encoding="utf-8"))
    if manifest["jobs"]!=jobs(): raise RuntimeError("Matrix changed")
    workers=min(8,max(1,int(os.environ.get("MODELG_WORKERS","2"))))
    ledger=execute(manifest["jobs"],directory,14400,600,workers,
                   projection=projection["per_T300_job_seconds"])
    print(json.dumps({"records":len(ledger["records"]),"charged_worker_seconds":ledger["charged_worker_seconds"],"stop_reason":ledger["stop_reason"]}))


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--worker",nargs=2)
    parser.add_argument("--benchmark",action="store_true")
    parser.add_argument("--production",action="store_true")
    args=parser.parse_args()
    if args.worker: sys.exit(worker(json.loads(args.worker[0]),args.worker[1]))
    elif args.benchmark: benchmark()
    elif args.production: production()
    else: parser.error("Choose --benchmark or --production")

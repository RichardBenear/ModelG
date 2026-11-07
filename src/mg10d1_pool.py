"""MG053 D1 owned-process controller, derived from bounded_pool at2397a0c.

Historical MG10/MG10-C1 driver sources are unchanged. Workers accept --worker
JSON DIRECTORY. No retry/replacement; elapsed time includes failed/terminated jobs.
"""
from pathlib import Path
import json,os,subprocess,sys,time


def atomic_json(path,obj):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    temp=path.with_suffix(path.suffix+".tmp")
    temp.write_text(json.dumps(obj,indent=2)+"\n",encoding="utf-8");os.replace(temp,path)


def terminate(process):
    process.terminate()
    try:process.wait(timeout=1)
    except subprocess.TimeoutExpired:process.kill();process.wait()


def run(jobs,directory,worker_script,cwd,ceiling,workers=1,prior=0,timeout=600,forecast=None,poll_seconds=.5):
    if not 1<=workers<=8:raise ValueError("Worker limit")
    if not 0<poll_seconds<=.5:raise ValueError("Poll interval")
    directory=Path(directory);directory.mkdir(parents=True,exist_ok=True)
    lock=directory/"driver.lock"
    with lock.open("x") as handle:handle.write(str(os.getpid()))
    file=directory/"ledger.json";active={};stop=None
    try:
        ledger=json.loads(file.read_text(encoding="utf-8")) if file.exists() else dict(records=[],ceiling=ceiling,prior_charge=prior,poll_seconds=poll_seconds)
        if ledger["ceiling"]!=ceiling or ledger["prior_charge"]!=prior:raise ValueError("Budget changed")
        ids=[j["id"] for j in jobs]
        if len(ids)!=len(set(ids)):raise ValueError("Duplicate job id")
        pending=[j for j in jobs if j["id"] not in {r["id"] for r in ledger["records"]}]
        def charged():return prior+sum(r["charged_worker_seconds"] for r in ledger["records"])
        def save():
            ledger.update(charged_worker_seconds=charged(),inflight_worker_seconds=sum(time.perf_counter()-a["start"] for a in active.values()),
                          pending_ids=[j["id"] for j in pending],active_ids=list(active),stop_reason=stop,workers=workers,threads_per_library=1)
            atomic_json(file,ledger)
        try:
            previous_cycle=None
            while pending or active:
                cycle=time.perf_counter()
                if previous_cycle is not None:ledger["max_observed_poll_gap_seconds"]=max(ledger.get("max_observed_poll_gap_seconds",0.),cycle-previous_cycle)
                previous_cycle=cycle
                now=cycle;inflight=sum(now-a["start"] for a in active.values())
                if charged()+inflight>=ceiling:stop="BUDGET-SHORT"
                if forecast and pending and len(active)<workers and not stop:
                    prediction=sum(forecast[j["id"]] for j in pending)
                    prediction+=sum(max(0,forecast[k]-(now-a["start"])) for k,a in active.items())
                    if prediction>.75*(ceiling-charged()-inflight):stop="BUDGET-SHORT: FORECAST"
                while pending and len(active)<workers and not stop:
                    job=pending.pop(0)
                    log=(directory/(job["id"]+".log")).open("w",encoding="utf-8")
                    process=subprocess.Popen([sys.executable,str(worker_script),"--worker",json.dumps(job),str(directory)],
                                             cwd=cwd,env=os.environ.copy(),stdout=log,stderr=subprocess.STDOUT)
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
                        if not meta.exists():atomic_json(meta,record)
                        ledger["records"].append(record);a["log"].close();del active[key]
                        if reason or process.returncode!=0 or record["status"] in ("FAILED","QUALIFICATION FAILED"):
                            stop=reason if reason else (record["status"] if record["status"]!="COMPLETE" else "FAILED: NONZERO EXIT")
                save()
                if stop and not active:break
                if active:time.sleep(max(0.,poll_seconds-(time.perf_counter()-cycle)))
        except BaseException:
            stop="INTERRUPTED"
            for key,a in active.items():
                if a["process"].poll() is None:terminate(a["process"])
                ledger["records"].append(dict(id=key,status="INTERRUPTED",job=a["job"],charged_worker_seconds=time.perf_counter()-a["start"]))
                a["log"].close()
            active.clear();save();raise
        save();return ledger
    finally:lock.unlink(missing_ok=True)

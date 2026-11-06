"""MG049 bounded actual-lattice qualification, reconstructible raw shards."""
import os
for key in ("OPENBLAS_NUM_THREADS","MKL_NUM_THREADS","OMP_NUM_THREADS","NUMEXPR_NUM_THREADS"):os.environ[key]="1"
from pathlib import Path
import argparse,hashlib,json,platform,sys,time,traceback
from bounded_pool import run,atomic_json

ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/"results/MG10S2"


def digest(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def text_digest(path):
    return hashlib.sha256(Path(path).read_text(encoding="utf-8").replace("\r\n","\n").encode("utf-8")).hexdigest()


def hashes():
    files=("src/mg10s2.py","src/mg10s2_measure.py","src/run_mg10s2.py","src/evaluate_mg10s2.py","src/check_mg10s2.py",
           "src/mg10s1.py","src/bounded_pool.py","src/runner.py")
    return {f:text_digest(ROOT/f) for f in files}


def classes():
    return [(case,method) for case in range(6) for method in range(6 if case<4 else 4)]


def jobs(timing=False):
    out=[]
    for case,method in classes():
        for block in range(1 if case<4 or timing else 8):
            kind="snapshot" if case<4 else "short"
            reps=(64 if case<4 else 8) if timing else (16384 if case<4 else 2048)
            chunk=reps if timing else (512 if case<4 else 256)
            out.append(dict(id=f"{kind}_c{case}_m{method}_b{block}",kind=kind,case=case,method=method,block=block,
                            reps=reps,chunk=chunk,offset=0 if timing or case<4 else block*8,
                            stage=4902 if timing else 4903,warmup=timing))
    return out


def seed_allocation(matrix):
    from runner import seed_for
    return {j["id"]:[seed_for(10,j["case"],j["method"],j["stage"],j["offset"]+k) for k in range((j["reps"]+j["chunk"]-1)//j["chunk"])] for j in matrix}


def worker(job,directory):
    started=time.perf_counter()
    import numpy as np
    import scipy,numba
    import mg10s2 as m
    from mg10s2_measure import measure
    from runner import seed_for
    directory=Path(directory);case=m.fixture(job["case"])
    shards=[];generated=0;compute=0.
    def sample(reps,seed):
        if job["kind"]=="snapshot":return m.snapshot(case,job["method"]//2,job["method"]%2,reps,seed)
        return m.short(case,job["method"],reps,seed)
    try:
        if job.get("warmup"):
            sample(1,seed_for(10,job["case"],job["method"],job["stage"],job["offset"],1))
        if job["kind"]=="snapshot":
            tau=m.TAUS[job["method"]//2];mean=m.conditional_mean(case,tau)
            q=m.covariance(case,tau,case["modes"]);q0=m.covariance(case,tau,case["modes"],True)
        while generated<job["reps"]:
            size=min(job["chunk"],job["reps"]-generated)
            global_chunk=job["offset"]+len(shards)
            seed=seed_for(10,job["case"],job["method"],job["stage"],global_chunk)
            begin=time.perf_counter();raw=sample(size,seed)
            path=directory/f"{job['id']}_part{global_chunk:03d}.npz"
            temporary=path.with_suffix(".tmp")
            common=dict(initial=case["initial"],a0=case["a0"],b0=case["b0"],omega=np.array(case["omega"]),
                        modes=case["modes"],seed=np.array(seed,dtype=np.uint32))
            if job["kind"]=="snapshot":
                events,populations,stats=raw
                values,_=measure(populations[:,:3],mean,case["omega"],tau,case["modes"])
                payload=dict(events=events,stats=stats,mean=mean,time=np.array(tau),measurements=values,q_snapshot=q,q_instant=q0,**common)
                population_sha=m.population_hash(populations)
            else:
                data,stats=raw
                measured=[]
                for k,t in enumerate(m.OBS):
                    values,_=measure(data[:,k,:3],case["initial"][:3],case["omega"],t,case["modes"])
                    measured.append(values)
                payload=dict(data=data,stats=stats,observations=m.OBS,measurements=np.stack(measured,axis=1),**common)
                population_sha=None
            with temporary.open("wb") as handle:np.savez_compressed(handle,**payload)
            os.replace(temporary,path)
            shard=dict(file=path.name,replicates=size,chunk=global_chunk,seed=seed,sha256=digest(path),population_sha256=population_sha)
            compute+=time.perf_counter()-begin;shards.append(shard);generated+=size
            # Checkpoint is a list of atomically completed raw shards.
            atomic_json(directory/(job["id"]+".progress.json"),dict(job=job,replicates=generated,shards=shards))
        record=dict(status="COMPLETE",job=job,replicates=generated,shards=shards,compute_seconds=compute)
    except Exception:
        record=dict(status="QUALIFICATION FAILED",job=job,replicates=generated,shards=shards,compute_seconds=compute,traceback=traceback.format_exc())
    record.update(worker_elapsed_seconds=time.perf_counter()-started,platform=platform.platform(),python=sys.version,
                  numpy=np.__version__,scipy=scipy.__version__,numba=numba.__version__,threads_per_library=1)
    atomic_json(directory/(job["id"]+".json"),record)
    return 0 if record["status"]=="COMPLETE" else 2


def benchmark():
    test=json.loads((BASE/"qualification/tests.json").read_text())
    if test["returncode"]:raise RuntimeError("Unit qualification failed")
    ledger=run(jobs(True),BASE/"timing",Path(__file__),ROOT,900,workers=min(8,int(os.environ.get("MODELG_WORKERS","2"))),prior=test["elapsed_worker_seconds"])
    if len(ledger["records"])!=32 or any(r["status"]!="COMPLETE" for r in ledger["records"]):raise RuntimeError("Timing incomplete")
    return project(ledger)


def project(ledger=None):
    # Pure cost arithmetic. Reuse class timings if only the unregistered sample
    # allocation changes; never generate extra physical trajectories here.
    if ledger is None:ledger=json.loads((BASE/"timing/ledger.json").read_text())
    overhead=max(5.,2*max(r["charged_worker_seconds"]-r["compute_seconds"] for r in ledger["records"]))
    records={(r["job"]["case"],r["job"]["method"]):r for r in ledger["records"]}
    forecast={}
    for job in jobs():
        r=records[job["case"],job["method"]]
        forecast[job["id"]]=2*r["compute_seconds"]*job["reps"]/r["replicates"]+overhead
    projection=dict(projected_worker_seconds=sum(forecast.values()),per_job=forecast,overhead=overhead,ceiling=7200,reserve=.25,
                    fits=sum(forecast.values())<=5400 and max(forecast.values())<=600,qualification_timing_seconds=ledger["charged_worker_seconds"],
                    formula="2x measured warm compute/measurement/serialization scaling + max(5s,2x observed startup/warm-up overhead) per fixed job")
    atomic_json(BASE/"timing/projection.json",projection)
    print({k:v for k,v in projection.items() if k!="per_job"})


def production():
    registration=json.loads((BASE/"registration.json").read_text())
    if registration["status"]!="REGISTERED" or registration["source_hashes"]!=hashes():raise RuntimeError("Registration/source mismatch")
    if registration["spec_sha256"]!=text_digest(ROOT/"spec/MG10S2_proposal.md"):raise RuntimeError("Specification changed")
    if registration["seeds"]!=seed_allocation(jobs()):raise RuntimeError("Seed allocation changed")
    projection=json.loads((BASE/"timing/projection.json").read_text())
    if not projection["fits"]:raise RuntimeError("Forecast does not fit")
    ledger=run(jobs(),BASE/"production",Path(__file__),ROOT,7200,workers=min(8,int(os.environ.get("MODELG_WORKERS","2"))),forecast=projection["per_job"],poll_seconds=.25)
    print(dict(jobs=len(ledger["records"]),charged_worker_seconds=ledger["charged_worker_seconds"],stop=ledger["stop_reason"]))


if __name__=="__main__":
    p=argparse.ArgumentParser();p.add_argument("--worker",nargs=2);p.add_argument("--benchmark",action="store_true");p.add_argument("--project-only",action="store_true");p.add_argument("--production",action="store_true");a=p.parse_args()
    if a.worker:sys.exit(worker(json.loads(a.worker[0]),a.worker[1]))
    elif a.benchmark:benchmark()
    elif a.project_only:project()
    elif a.production:production()
    else:p.error("Choose action")

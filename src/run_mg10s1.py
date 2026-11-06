"""MG046 finite stochastic qualification ensembles and bounded execution."""
import os
for key in ("OPENBLAS_NUM_THREADS","MKL_NUM_THREADS","OMP_NUM_THREADS","NUMEXPR_NUM_THREADS"):os.environ[key]="1"
from pathlib import Path
import argparse,hashlib,json,platform,sys,time,traceback
from bounded_pool import run,atomic_json

ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/"results/MG10S1"


def digest(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def hashes():
    names=("src/mg10s1.py","src/run_mg10s1.py","src/evaluate_mg10s1.py","src/check_mg10s1.py","src/bounded_pool.py","src/runner.py","src/lattice.py")
    return {name:digest(ROOT/name) for name in names}


def jobs(stage=4603,reps=32768):
    return [dict(id=f"case{case}_method{method}",case=case,method=method,stage=stage,reps=reps,chunk=512) for case in range(4) for method in range(4)]


def seed_allocation(matrix):
    from runner import seed_for
    return {j["id"]:[seed_for(10,j["case"],j["method"],j["stage"],c) for c in range((j["reps"]+j["chunk"]-1)//j["chunk"])] for j in matrix}


def worker(job,directory):
    began=time.perf_counter()
    import numpy as np
    import scipy,numba
    import mg10s1 as m
    from runner import seed_for
    directory=Path(directory)
    case=m.fixture(job["case"]);data=[];stats=[];seeds=[];compute=0.;generated=0
    partial=directory/(job["id"]+".checkpoint.npz")
    def checkpoint():
        temp=partial.with_suffix(".tmp")
        with temp.open("wb") as handle:
            np.savez_compressed(handle,data=np.concatenate(data),stats=np.concatenate(stats),initial=case["state"],
                                observations=m.OBS,seeds=np.array(seeds,dtype=np.uint32),a0=case["a0"],b0=case["b0"])
        os.replace(temp,partial)
    try:
        total_audit=dict(max_integer_residual=0,max_total_count_error=0,nonnegative=True,pass_check=True)
        while generated<job["reps"]:
            chunk=len(seeds);size=min(job["chunk"],job["reps"]-generated)
            seed=seed_for(10,job["case"],job["method"],job["stage"],chunk)
            start=time.perf_counter();d,s=m.simulate(case,job["method"],size,seed);compute+=time.perf_counter()-start
            audit=m.ledger(d,case["state"])
            if not audit["pass_check"]:raise ValueError(f"Ledger {audit}")
            data.append(d);stats.append(s);seeds.append(seed);generated+=size
            if job["method"]>0 and np.max(s[:,2])>.25+1e-14:raise ValueError("Probability bound")
            if len(seeds)%8==0 or generated==job["reps"]:checkpoint()
        final=directory/(job["id"]+".npz");os.replace(partial,final)
        record=dict(status="COMPLETE",job=job,replicates=generated,compute_seconds=compute,audit=total_audit,sha256=digest(final),
                    total_events_or_steps=float(sum(np.sum(s[:,0]) for s in stats)),adapted_steps=float(sum(np.sum(s[:,1]) for s in stats)),
                    max_probability=float(max(s[:,2].max() for s in stats)))
    except Exception:
        record=dict(status="QUALIFICATION FAILED",job=job,replicates=generated,compute_seconds=compute,traceback=traceback.format_exc())
    record.update(worker_elapsed_seconds=time.perf_counter()-began,platform=platform.platform(),python=sys.version,numpy=np.__version__,scipy=scipy.__version__,numba=numba.__version__,threads_per_library=1)
    atomic_json(directory/(job["id"]+".json"),record)
    return 0 if record["status"]=="COMPLETE" else 2


def benchmark():
    test=json.loads((BASE/"qualification/tests.json").read_text())
    if test["returncode"]:raise RuntimeError("Unit qualification failed")
    matrix=jobs(4602,512)
    ledger=run(matrix,BASE/"timing",Path(__file__),ROOT,360,workers=min(8,int(os.environ.get("MODELG_WORKERS","2"))),prior=test["elapsed_worker_seconds"])
    if len(ledger["records"])!=16 or any(r["status"]!="COMPLETE" for r in ledger["records"]):raise RuntimeError("Timing/qualification incomplete")
    overhead=max(5.,2*max(r["charged_worker_seconds"]-r["compute_seconds"] for r in ledger["records"]))
    estimates={r["id"]:2*64*r["compute_seconds"]+overhead for r in ledger["records"]}
    result=dict(projected_worker_seconds=sum(estimates.values()),per_job=estimates,overhead=overhead,
                ceiling=3600,reserve=.25,fits=sum(estimates.values())<=2700 and max(estimates.values())<=600,
                qualification_timing_seconds=ledger["charged_worker_seconds"],
                formula="2x compute scaling from512 to32768 per case/method + max(5s,2x observed startup/serialization/poll overhead)")
    atomic_json(BASE/"timing/projection.json",result);print(json.dumps(result,indent=2))


def production():
    reg=json.loads((BASE/"registration.json").read_text())
    if reg["status"]!="REGISTERED" or reg["source_hashes"]!=hashes() or reg["spec_sha256"]!=digest(ROOT/"spec/MG10S1_proposal.md"):raise RuntimeError("Registration/code/spec mismatch")
    projection=json.loads((BASE/"timing/projection.json").read_text())
    if not projection["fits"]:raise RuntimeError("Forecast does not fit")
    if reg["seeds"]!=seed_allocation(jobs()):raise RuntimeError("Seed allocation changed")
    ledger=run(jobs(),BASE/"production",Path(__file__),ROOT,3600,workers=min(8,int(os.environ.get("MODELG_WORKERS","2"))),forecast=projection["per_job"])
    print(dict(groups=len(ledger["records"]),charged_worker_seconds=ledger["charged_worker_seconds"],stop=ledger["stop_reason"]))


if __name__=="__main__":
    p=argparse.ArgumentParser();p.add_argument("--worker",nargs=2);p.add_argument("--benchmark",action="store_true");p.add_argument("--production",action="store_true");args=p.parse_args()
    if args.worker:sys.exit(worker(json.loads(args.worker[0]),args.worker[1]))
    elif args.benchmark:benchmark()
    elif args.production:production()
    else:p.error("Choose --benchmark or --production")

"""MG052 qualification and cost evidence only; no production entry point."""
import os
for key in ("OPENBLAS_NUM_THREADS","MKL_NUM_THREADS","OMP_NUM_THREADS","NUMEXPR_NUM_THREADS"):
    os.environ[key]="1"
from pathlib import Path
import argparse, hashlib, json, platform, sys, time
from bounded_pool import run, atomic_json

ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/"results/MG10D1/development"


def hashes():
    files=("src/mg10d1.py","src/prepare_mg10d1.py","src/check_mg10d1_development.py",
           "src/mg10s1.py","src/mg10s2.py","src/bounded_pool.py","src/runner.py",
           "spec/MG10D1_proposal.md")
    out={f:hashlib.sha256((ROOT/f).read_text(encoding="utf-8").replace("\r\n","\n").encode()).hexdigest() for f in files}
    out["results/MG003/ref_profiles.npz"]=hashlib.sha256((ROOT/"results/MG003/ref_profiles.npz").read_bytes()).hexdigest()
    return out


def store(directory,name,data,case,stats):
    import numpy as np
    path=directory/(name+".npz")
    np.savez_compressed(path,data=data,initial=case["initial"],a0=case["a0"],b0=case["b0"],
                        omega=case["omega"],clamp_a=case["clamp_a"],clamp_b=case["clamp_b"],stats=stats)
    return dict(path=path.name,sha256=hashlib.sha256(path.read_bytes()).hexdigest())


def qualify(directory):
    import numpy as np
    from scipy.stats import poisson, binom
    import mg10s1 as old
    import mg10d1 as model
    from runner import seed_for
    tests=[];raw=[]
    for index in range(4):
        oldcase=old.fixture(index)
        case=dict(initial=oldcase["state"],a0=oldcase["a0"],b0=oldcase["b0"],omega=oldcase["omega"],diff=oldcase["diff"],clamp_a=False,clamp_b=False)
        seed=seed_for(10,index,0,5201,0)
        obs=np.array([.005,.01]);mask=np.ones(7)
        a,astat=old.leap_chunk(case["initial"],case["a0"],case["b0"],case["omega"],case["diff"],mask,.00025,obs,4,seed)
        b,bstat=model.simulate(case,.00025,obs,4,seed)
        ok=np.array_equal(a,b[:,:,:20]) and np.array_equal(astat,bstat) and not np.any(b[:,:,20:22]) and model.ledger(b,case["initial"])["pass_check"]
        tests.append(dict(name=f"finite_bit_regression_{index}",passed=bool(ok),seed=seed))
        raw.append(store(directory,f"finite_{index}",b,case,bstat))
    for arm in range(5):
        for start in range(2):
            case=model.fixture(arm,start)
            seed=seed_for(10,2*arm+start,0,5201,1)
            data,stats=model.simulate(case,1/24000,np.array([0.,.005,.01]),1,seed)
            checks=model.ledger(data,case["initial"])
            ok=checks["pass_check"] and stats[:,2].max()<=.25+1e-14
            for species,flag in ((3,case["clamp_a"]),(4,case["clamp_b"])):
                if flag:ok=ok and np.all(data[:,:,species]==case["initial"][species])
            tests.append(dict(name=f"full_arm{arm}_start{start}",passed=bool(ok),seed=seed,ledger=checks))
            raw.append(store(directory,f"full_a{arm}_s{start}",data,case,stats))
    for channel in (0,3):
        c=old.fixture(0)
        case=dict(initial=c["state"],a0=c["a0"],b0=c["b0"],omega=c["omega"],diff=np.zeros(3),clamp_a=True,clamp_b=True)
        # Scarce nominal B tests that replenishment, rather than finite stock,
        # determines conversion. X initially16, B initially1.
        case["initial"][4]=1;case["b0"][:]=1
        mask=np.zeros(7);mask[channel]=1
        seed=seed_for(10,channel,0,5201,2)
        data,stats=model.simulate(case,.001,np.array([.001]),32768,seed,mask)
        values=data[:,0,7+channel,0]
        grid=np.arange(int(values.max())+1)
        empirical=np.cumsum(np.bincount(values,minlength=len(grid)))/len(values)
        target=poisson.cdf(grid,.014) if channel==0 else binom.cdf(grid,16,.029)
        distance=float(np.max(np.abs(empirical-target)))
        upper=distance+float(np.sqrt(np.log(2000)/(2*len(values))))
        tests.append(dict(name="maintained_A_poisson" if channel==0 else "maintained_B_binomial",passed=bool(upper<=.05 and model.ledger(data,case["initial"])["pass_check"]),seed=seed,cdf_distance=distance,dkw_upper=upper,margin=.05,alpha=.001))
        raw.append(store(directory,f"isolated_channel{channel}",data,case,stats))
    # Meaningful domain rejection checks protect long jobs from silent bad input.
    bad=[lambda:model.fixture(5,0),lambda:model.fixture(0,2),lambda:model.simulate(model.fixture(0,0),0,[.1],1,1),lambda:model.simulate(model.fixture(0,0),.001,[.1,.05],1,1)]
    for index,call in enumerate(bad):
        try:call();ok=False
        except ValueError:ok=True
        tests.append(dict(name=f"invalid_input_{index}",passed=ok))
    result=dict(status="DEVELOPMENT ONLY",tests=tests,raw=raw,passed=all(t["passed"] for t in tests),source_hashes=hashes())
    atomic_json(directory/"qualification.json",result)
    if not result["passed"]:raise ValueError("Blocking D1 qualification failure")


def worker(job,directory):
    directory=Path(directory);start=time.perf_counter()
    if job["kind"]=="qualification":
        qualify(directory)
    else:
        import numpy as np
        import mg10d1 as model
        from runner import seed_for
        case=model.fixture(job["arm"],job["start"])
        tau=(1/24000,1/48000)[job["step"]]
        seed=seed_for(10,2*job["arm"]+job["start"],job["step"],5202,0)
        # Compile/load before measuring steady-state sampling cost. The parent
        # still charges all startup, compilation, analysis and serialization.
        model.simulate(case,tau,np.array([tau]),1,seed_for(10,2*job["arm"]+job["start"],job["step"],5202,0,1))
        clock=time.perf_counter()
        data,stats=model.simulate(case,tau,np.array([0.,.1]),1,seed)
        compute=time.perf_counter()-clock
        check=model.ledger(data,case["initial"])
        if not check["pass_check"]:raise ValueError("Timing ledger")
        raw=store(directory,job["id"],data,case,stats)
        atomic_json(directory/(job["id"]+".json"),dict(status="COMPLETE",job=job,tau=tau,seed=seed,horizon=.1,compute_seconds=compute,ledger=check,raw=raw,stats=stats.tolist(),source_hashes=hashes()))
    atomic_json(directory/(job["id"]+"_complete.json"),dict(job=job,elapsed_seconds=time.perf_counter()-start))


def develop():
    import numpy, scipy, numba
    BASE.mkdir(parents=True,exist_ok=True)
    if (BASE/"control_plan.json").exists():raise ValueError("Development already attempted; preserve it and explicitly version any retry")
    # Earlier unchanged-kernel throughput work is retained and charged here.
    prior=sum(r["compute_seconds"] for r in json.loads((BASE/"timing.json").read_text(encoding="utf-8")))+15.
    matrix=[dict(id=f"timing_a{arm}_s{start}_k{step}",kind="timing",arm=arm,start=start,step=step) for arm in range(5) for start in range(2) for step in range(2)]
    atomic_json(BASE/"control_plan.json",dict(status="DEVELOPMENT ONLY; no scientific production",jobs=matrix,qualification=dict(id="qualification",kind="qualification"),source_hashes=hashes(),ceiling_seconds=900,prior_charge_seconds=prior,prior_charge_note="Recorded unchanged-kernel compute plus conservative15seconds for startup/JIT/serialization; no exact parent ledger existed for that initial inline timing.",workers=8,threads=1,environment=dict(python=sys.version,platform=platform.platform(),numpy=numpy.__version__,scipy=scipy.__version__,numba=numba.__version__)))
    q=run([dict(id="qualification",kind="qualification")],BASE/"qualification_run",Path(__file__).resolve(),ROOT,900,workers=1,prior=prior,timeout=300,poll_seconds=.25)
    if q.get("stop_reason") or any(r.get("returncode")!=0 for r in q["records"]):raise ValueError("Qualification controller stopped")
    # Qualification output is inside its controlled directory.
    if not json.loads((BASE/"qualification_run/qualification.json").read_text(encoding="utf-8"))["passed"]:raise ValueError("Qualification did not pass")
    t=run(matrix,BASE/"control_timing",Path(__file__).resolve(),ROOT,900,workers=8,prior=q["charged_worker_seconds"],timeout=60,poll_seconds=.25)
    atomic_json(BASE/"development_summary.json",dict(qualification_ledger=q,timing_ledger=t))


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--worker",nargs=2,metavar=("JSON","DIRECTORY"))
    parser.add_argument("--develop",action="store_true")
    args=parser.parse_args()
    if args.worker:worker(json.loads(args.worker[0]),args.worker[1])
    elif args.develop:develop()
    else:parser.error("Choose development only; production is not implemented/registered")

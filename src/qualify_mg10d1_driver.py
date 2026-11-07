"""MG053 prospective checkpoint/controller checks; development streams only."""
import os
for key in ("OPENBLAS_NUM_THREADS","MKL_NUM_THREADS","OMP_NUM_THREADS","NUMEXPR_NUM_THREADS"):
    os.environ[key]="1"
from pathlib import Path
import argparse, json, sys, time
from mg10d1_pool import run, atomic_json

ROOT=Path(__file__).resolve().parents[1];DIRECTORY=ROOT/"results/MG10D1/development/driver_checks"


def fake(job,directory):
    if job["kind"]=="sleep":time.sleep(3)
    status="FAILED" if job["kind"]=="fail" else "COMPLETE"
    atomic_json(Path(directory)/(job["id"]+".json"),dict(status=status,job=job))
    if job["kind"] in ("fail","bad_exit"):sys.exit(1)


def qualify():
    import numpy as np
    import mg10d1 as model
    import mg10s1 as old
    from mg10d1_stream import seed_rng, advance
    from run_mg10d1 import save_npz, digest
    from runner import seed_for
    from check_mg10d1_development import check_raw
    from check_mg10d1 import independent_metrics
    from evaluate_mg10d1 import describe
    started=time.perf_counter()
    DIRECTORY.mkdir(parents=True,exist_ok=True)
    output=DIRECTORY.parent/"driver_qualification.json"
    if output.exists():raise ValueError("Qualification already attempted; preserve before any revision")
    cases=[]
    for index in range(4):
        c=old.fixture(index)
        cases.append((f"small{index}",dict(initial=c["state"],a0=c["a0"],b0=c["b0"],omega=c["omega"],diff=c["diff"],clamp_a=False,clamp_b=False),.00025,np.array([0.,.005,.01])))
    for arm in range(5):
        for start in range(2):cases.append((f"full_a{arm}_s{start}",model.fixture(arm,start),1/24000,np.array([0.,.005,.01])))
    stiff=old.fixture(0)
    stiff["state"][1]=1000
    cases.append(("adaptive",dict(initial=stiff["state"],a0=stiff["a0"],b0=stiff["b0"],omega=1.,diff=np.zeros(3),clamp_a=False,clamp_b=False),.00025,np.array([0.,.000001,.000002])))
    tests=[];raw=[]
    try:
        for index,(name,case,tau,obs) in enumerate(cases):
            seed=seed_for(10,index,0,5301,0)
            expected,stats=model.simulate(case,tau,obs,1,seed)
            state=case["initial"].copy();n=state.shape[1]
            flux=np.zeros((7,n),dtype=np.int64);hops=np.zeros((6,n),dtype=np.int64);external=np.zeros((2,n),dtype=np.int64)
            seed_rng(seed);actual=np.zeros_like(expected);observed_stats=np.zeros((1,3));previous=0.
            for k,target in enumerate(obs):
                d=advance(state,flux,hops,external,case["a0"],case["b0"],case["omega"],case["diff"],np.ones(7),tau,previous,target,case["clamp_a"],case["clamp_b"])
                observed_stats[0,:2]+=d[:2];observed_stats[0,2]=max(observed_stats[0,2],d[2]);previous=target
                actual[0,k,:7]=state;actual[0,k,7:14]=flux;actual[0,k,14:20]=hops;actual[0,k,20:22]=external
            assert np.array_equal(actual,expected) and np.array_equal(observed_stats,stats),name
            if name=="adaptive":assert stats[0,1]>0
            path=DIRECTORY/(name+".npz")
            save_npz(path,data=actual,initial=case["initial"],a0=case["a0"],b0=case["b0"],omega=case["omega"],clamp_a=case["clamp_a"],clamp_b=case["clamp_b"],stats=observed_stats,times=obs,seed=seed)
            checked=check_raw(path);assert checked["sha256"]==digest(path)
            if name.startswith("full"):
                independent=independent_metrics(actual,case["a0"],case["b0"])
                for metric,v in model.measurements(actual,case).items():np.testing.assert_allclose(v[0],independent[metric],rtol=2e-12,atol=2e-11)
            tests.append(dict(name="streaming_bit_identity_"+name,passed=True,seed=seed,adaptive_steps=float(stats[0,1])))
            raw.append(checked)
        for kind,kwargs,expected_reason in (("fail",dict(ceiling=10,timeout=5),"FAILED"),
                                           ("bad_exit",dict(ceiling=10,timeout=5),"FAILED: NONZERO EXIT"),
                                           ("timeout",dict(ceiling=10,timeout=.5),"TIMEOUT"),
                                           ("budget",dict(ceiling=.5,timeout=5),"BUDGET-SHORT")):
            folder=DIRECTORY/("controller_"+kind)
            jobkind=kind if kind in ("fail","bad_exit") else "sleep"
            ledger=run([dict(id="first",kind=jobkind),dict(id="never_dispatch",kind="complete")],folder,Path(__file__).resolve(),ROOT,workers=1,poll_seconds=.1,**kwargs)
            assert ledger["stop_reason"]==expected_reason,(kind,ledger)
            assert len(ledger["records"])==1 and ledger["pending_ids"]==["never_dispatch"]
            assert not ledger["active_ids"] and not (folder/"driver.lock").exists()
            tests.append(dict(name="controller_"+kind,passed=True,charge=ledger["charged_worker_seconds"],stop_reason=ledger["stop_reason"]))
        description=describe([[1.,10.],[2.,20.],[3.,30.],[4.,40.]])
        np.testing.assert_allclose(description["mean"],[2.5,25.])
        half=3.182446305284263*np.sqrt(5/3)/2
        np.testing.assert_allclose(description["ci95_approx"],[[2.5-half,25-10*half],[2.5+half,25+10*half]])
        assert describe([])["mean"] is None and describe([[1.]])["ci95_approx"] is None
        tests.append(dict(name="descriptive_statistics_known_values_missing_denominators",passed=True))
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig,ax=plt.subplots();ax.plot([0,1],[0,1]);fig.savefig(DIRECTORY/"plot_backend.png");plt.close(fig)
        tests.append(dict(name="plot_backend_available",passed=True))
        result=dict(passed=True,status="DRIVER QUALIFICATION PASS; DEVELOPMENT ONLY",tests=tests,raw=raw,
                    elapsed_seconds=time.perf_counter()-started,prior_development_charge=136.54094180004904,
                    total_development_charge_with_5s_startup_allowance=136.54094180004904+time.perf_counter()-started+5)
        assert result["total_development_charge_with_5s_startup_allowance"]<900
    except BaseException:
        import traceback
        atomic_json(output,dict(passed=False,status="QUALIFICATION FAILED; STOP",tests=tests,raw=raw,error=traceback.format_exc(),elapsed_seconds=time.perf_counter()-started));raise
    atomic_json(output,result)
    print(f"Driver qualification PASS: {len(tests)} checks; {result['elapsed_seconds']:.3f}s",flush=True)


if __name__=="__main__":
    parser=argparse.ArgumentParser();parser.add_argument("--worker",nargs=2);args=parser.parse_args()
    if args.worker:fake(json.loads(args.worker[0]),args.worker[1])
    else:qualify()

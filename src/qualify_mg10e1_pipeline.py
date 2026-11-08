"""MG070 disclosed replay and algebraic complete/incomplete pipeline fixtures."""
import os
for key in ("OPENBLAS_NUM_THREADS","MKL_NUM_THREADS","OMP_NUM_THREADS","NUMEXPR_NUM_THREADS"):os.environ[key]="1"
from pathlib import Path
import json,time
import numpy as np
from mg10d1_pool import atomic_json
from mg10p2 import fixture,reference
from mg10e1 import fork_case,activity
from mg10e1_stream import seed_rng,advance
from prepare_mg10e1 import candidate_jobs
from run_mg10e1 import save_npz
from evaluate_mg10e1 import main as evaluate
from check_mg10e1 import main as check
ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/"results/MG10E1/preparation"
def main():
    started=time.perf_counter()
    if (BASE/"pipeline_qualification.json").exists():raise ValueError("Already attempted; preserve and version corrections")
    matches=0
    for s in range(2):
        for k in range(2):
            for a in range(4):
                with np.load(BASE/f"devparents/parent_s{s}_k{k}.npz") as raw:case=fork_case(fixture(0,s),raw["data"][0,-1,:7])
                with np.load(BASE/f"timing/timing_a{a}_s{s}_k{k}.npz") as raw:expected=raw["data"].copy();expected_stats=raw["stats"].copy();seed=int(raw["seed"])
                state=case["initial"].copy();flux=np.zeros((7,500),dtype=np.int64);hops=np.zeros((6,500),dtype=np.int64);external=np.zeros((2,500),dtype=np.int64);tau=(1/24000,1/48000)[k]
                # Compile before resetting the held development stream.
                advance(state.copy(),flux.copy(),hops.copy(),external.copy(),case["a0"],case["b0"],1000.,case["diff"],np.ones(7),tau,0.,0.,True,True,activity(a))
                seed_rng(seed);stats=advance(state,flux,hops,external,case["a0"],case["b0"],1000.,case["diff"],np.ones(7),tau,0.,.25,True,True,activity(a))
                np.testing.assert_array_equal(np.concatenate([state,flux,hops,external]),expected[0,-1]);np.testing.assert_array_equal(stats,expected_stats[0]);matches+=1
    # Parent stream all-one regression against four held ungated .25 records.
    for s in range(2):
        for k in range(2):
            case=fixture(0,s)
            with np.load(BASE/f"devparents/parent_s{s}_k{k}.npz") as raw:expected=raw["data"].copy();stats_expected=raw["stats"].copy();seed=int(raw["seed"])
            state=case["initial"].copy();flux=np.zeros((7,500),dtype=np.int64);hops=np.zeros((6,500),dtype=np.int64);external=np.zeros((2,500),dtype=np.int64)
            seed_rng(seed);stats=advance(state,flux,hops,external,case["a0"],case["b0"],1000.,case["diff"],np.ones(7),(1/24000,1/48000)[k],0.,.25,True,True,np.ones(500))
            np.testing.assert_array_equal(np.concatenate([state,flux,hops,external]),expected[0,-1]);np.testing.assert_array_equal(stats,stats_expected[0]);matches+=1
    synth=BASE/"SYNTHETIC_pipeline_fixture";synth.mkdir(exist_ok=False);p,b=candidate_jobs()
    atomic_json(synth/"registration.json",dict(parent_jobs=p,branch_jobs=b))
    for directory,jobs,frames in (("parents",p,21),("branches",b,81)):
        d=synth/directory;d.mkdir();ref=reference()[0]
        for job in jobs:
            case=fixture(0,job["start"]);state=case["initial"].copy()
            if directory=="branches":
                # Algebraic response fixture, not a trajectory or scientific run.
                amp=job["start"]*(1+.05*job["arm"]+.002*job["block"]*job["arm"])
                state[:3]=np.rint(1000*(case["background"][:,None]+amp*ref)).astype(np.int64)
            data=np.zeros((1,frames,22,500),dtype=np.int64);data[0,:,:7]=state
            save_npz(d/(job["id"]+".npz"),data=data,initial=state,a0=case["a0"],b0=case["b0"],omega=1000.,clamp_a=True,clamp_b=True,stats=np.zeros((1,3)),times=np.arange(frames)/4,a_activity=activity(job["arm"]) if directory=="branches" else np.ones(500))
            atomic_json(d/(job["id"]+".json"),dict(status="COMPLETE",scope="ALGEBRAIC SYNTHETIC ONLY"))
    evaluate(synth);check(synth,synthetic=True)
    v=json.loads((synth/"evaluation.json").read_text());assert v["status"]=="COMPLETE INTERVENTION ESTIMATION" and all(v["primary_by_step"][str(k)]["n"]==8 for k in range(2))
    # Preserve complete evidence before imposing an explicit absent branch.
    for name in ("evaluation.json","independent_check.json","MG10E1_results.md"):(synth/name).rename(synth/("complete_"+name))
    missing=synth/"branches/a1_s1_k0_b0.npz";missing.rename(missing.with_suffix(".npz.withheld"))
    evaluate(synth);check(synth,synthetic=True);v=json.loads((synth/"evaluation.json").read_text());assert v["primary_by_step"]["0"]["n"]==7 and v["primary_by_step"]["1"]["n"]==8 and v["status"].startswith("INCOMPLETE")
    for name in ("evaluation.json","independent_check.json","MG10E1_results.md"):(synth/name).rename(synth/("incomplete_"+name))
    # A short available prefix is also excluded, not counted as physical loss.
    with np.load(missing.with_suffix(".npz.withheld")) as raw:arrays={name:raw[name] for name in raw.files}
    arrays["data"]=arrays["data"][:,:61];arrays["times"]=arrays["times"][:61];save_npz(missing,**arrays)
    evaluate(synth);check(synth,synthetic=True);v=json.loads((synth/"evaluation.json").read_text());assert v["branches_available"]==128 and v["branches_complete"]==127 and v["primary_by_step"]["0"]["n"]==7
    elapsed=time.perf_counter()-started;prior=json.loads((BASE/"independent_check.json").read_text())["total_preparation_allowance_seconds"];assert prior+elapsed<900
    atomic_json(BASE/"pipeline_qualification.json",dict(passed=True,held_development_stream_replays=matches,synthetic_complete_branches=128,synthetic_missing_branches=127,synthetic_prefix_excluded=True,primary_counts_incomplete=[7,8],additional_seconds=elapsed,total_preparation_allowance_seconds=prior+elapsed,scope="Replayed held development streams, no fresh scientific trajectory; complete/incomplete/prefix algebraic fixtures and separate dense-projection/statistical audit. Not physical simulations or external replication."))
    print(json.dumps(json.loads((BASE/"pipeline_qualification.json").read_text()),indent=2))
if __name__=="__main__":main()

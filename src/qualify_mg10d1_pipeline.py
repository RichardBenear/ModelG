"""MG053 end-to-end SYNTHETIC analysis fixture; not simulated evidence."""
import os
for key in ("OPENBLAS_NUM_THREADS","MKL_NUM_THREADS","OMP_NUM_THREADS","NUMEXPR_NUM_THREADS"):
    os.environ[key]="1"
from pathlib import Path
import json,time,traceback
import numpy as np
import mg10d1 as model
import evaluate_mg10d1 as evaluator
import check_mg10d1 as checker
from run_mg10d1 import jobs, save_npz, digest
from bounded_pool import atomic_json

ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/"results/MG10D1/development"


def main():
    started=time.perf_counter();output=BASE/"pipeline_qualification.json"
    if output.exists():raise ValueError("Already attempted; preserve any failed qualification")
    folder=BASE/"SYNTHETIC_pipeline_fixture"
    (folder/"production").mkdir(parents=True,exist_ok=True)
    matrix=jobs();records=[]
    atomic_json(folder/"registration.json",dict(status="SYNTHETIC TEST ONLY; NOT REGISTERED SCIENTIFIC DATA",source_hashes={},jobs=matrix,
                ceiling_worker_seconds=216000,workers=8,forecast_seconds={j["id"]:0 for j in matrix}))
    try:
        for job in matrix:
            case=model.fixture(job["arm"],job["start"])
            data=np.zeros((1,61,22,500),dtype=np.int64)
            data[0,:,:7]=case["initial"]
            # An algebraic source-only fixture, no stochastic simulation.
            # Replicate-dependent slopes exercise nonzero summary intervals.
            for t in range(61):
                count=t*job["replicate"]
                data[0,t,0]+=count;data[0,t,7]=count
                if case["clamp_a"]:data[0,t,20]=count
                else:data[0,t,3]-=count
            path=folder/"production"/(job["id"]+".npz")
            save_npz(path,data=data,initial=case["initial"],a0=case["a0"],b0=case["b0"],omega=1000.,
                     clamp_a=case["clamp_a"],clamp_b=case["clamp_b"],stats=np.array([[int(round(60/job["tau"])),0,0]]),times=np.arange(61),seed=job["seed"])
            records.append(dict(id=job["id"],status="COMPLETE",returncode=0,raw=path.name,raw_sha256=digest(path),job=job,charged_worker_seconds=0))
        atomic_json(folder/"production/ledger.json",dict(status="SYNTHETIC TEST ONLY",records=records,charged_worker_seconds=0,stop_reason=None))
        evaluator.BASE=folder;checker.BASE=folder
        evaluator.evaluate();checker.main()
        result=json.loads((folder/"production/evaluation.json").read_text(encoding="utf-8"))
        assert result["completed_paths"]==80 and result["verdict"]=="COMPLETE DESCRIPTIVE PILOT"
        for group in result["groups"].values():
            assert group["n"]==4
            initial=model.fixture(0,1 if group["start"]=="S" else 0)["initial"][0].mean()/1000
            assert abs(group["summary"]["field_mean"]["endpoint"]["mean"][0]-initial-.09)<1e-11
            assert abs(group["summary"]["field_mean"]["time_average"]["mean"][0]-initial-.045)<1e-11
        assert len(result["contrasts"])==26
        for contrast in result["contrasts"].values():
            np.testing.assert_allclose(contrast["mean_difference"]["field_mean"]["endpoint"],[0,0,0],atol=1e-11)
        # Exercise incomplete-denominator reporting/audit with one recorded
        # timeout. Its raw file remains as partial/failed evidence.
        records[-1]["status"]="TIMEOUT";records[-1]["returncode"]=1
        atomic_json(folder/"production/ledger.json",dict(status="SYNTHETIC TEST ONLY",records=records,charged_worker_seconds=0,stop_reason="TIMEOUT"))
        evaluator.evaluate();checker.main()
        partial=json.loads((folder/"production/evaluation.json").read_text(encoding="utf-8"))
        assert partial["completed_paths"]==79 and partial["verdict"]=="INCOMPLETE" and len(partial["excluded"])==1
        tests=["80-path raw reconstruction and complete report","Known endpoint/time-average means with nonzero intervals","26 supply/step contrasts","79/80 incomplete report and independent check","Plot generation, source/reference/seed and raw hash checks"]
        elapsed=time.perf_counter()-started
        previous=json.loads((BASE/"driver_qualification.json").read_text(encoding="utf-8"))["total_development_charge_with_5s_startup_allowance"]
        result=dict(passed=True,status="SYNTHETIC PIPELINE QUALIFICATION PASS; not scientific evidence",tests=tests,elapsed_seconds=elapsed,
                    combined_development_seconds_with_startup_allowance=previous+elapsed+5,
                    final_fixture_state="Deliberately incomplete79/80 to test missing-denominator handling; no scientific trajectories generated")
        assert result["combined_development_seconds_with_startup_allowance"]<900
        atomic_json(output,result)
    except BaseException:
        atomic_json(output,dict(passed=False,status="QUALIFICATION FAILED; no scientific production",error=traceback.format_exc(),elapsed_seconds=time.perf_counter()-started));raise
    print(json.dumps(result,indent=2),flush=True)


if __name__=="__main__":main()

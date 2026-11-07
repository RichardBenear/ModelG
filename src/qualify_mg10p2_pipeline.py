"""MG062 full/partial synthetic P2 pipeline fixtures, no stochastic draws."""
import os
for key in ("OPENBLAS_NUM_THREADS","MKL_NUM_THREADS","OMP_NUM_THREADS","NUMEXPR_NUM_THREADS"):
    os.environ[key]="1"
from pathlib import Path
import json,shutil,time,traceback
import numpy as np
from mg10d1_pool import atomic_json
import mg10p2 as model
import evaluate_mg10p2 as evaluator
import check_mg10p2 as checker
from run_mg10p2 import jobs,save_npz,digest

ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/"results/MG10P2"


def main():
    started=time.perf_counter();output=BASE/"preparation/pipeline_qualification.json"
    if output.exists():raise ValueError("Already attempted; preserve before correction")
    folder=BASE/"preparation/SYNTHETIC_pipeline_fixture";(folder/"production").mkdir(parents=True,exist_ok=True)
    matrix=jobs();records=[]
    atomic_json(folder/"registration.json",dict(status="SYNTHETIC ONLY; NO SCIENTIFIC TRAJECTORIES",source_hashes={},jobs=matrix,ceiling_worker_seconds=14400,expected_paths=1024))
    try:
        for job in matrix:
            case=model.fixture(job["condition"],job["start"]);data=np.zeros((8,4,22,500),dtype=np.int64);data[:,:,:7]=case["initial"]
            planted=model.fixture(job["condition"],1)["initial"][1]
            for index,rep in enumerate(job["replicates"]):
                for t in range(4):
                    source=rep*t
                    data[index,t,0]+=source;data[index,t,7]=source
                    if case["clamp_a"]:data[index,t,20]=source
                    else:data[index,t,3]-=source
                delta=np.zeros(500,dtype=np.int64)
                if not job["start"] and rep<2:delta=planted-case["initial"][1]
                if job["start"] and rep<3:delta=round(1000*case["background"][1])-case["initial"][1]
                data[index,-1,0]-=delta;data[index,-1,1]+=delta
                data[index,-1,8]=np.maximum(delta,0);data[index,-1,9]=np.maximum(-delta,0)
            assert np.all(data>=0)
            path=folder/"production"/(job["id"]+".npz")
            stats=np.tile([int(round(.25/(1/24000,1/48000)[job["step"]])),0,0],(8,1))
            save_npz(path,data=data,initial=case["initial"],a0=case["a0"],b0=case["b0"],omega=1000.,clamp_a=case["clamp_a"],clamp_b=case["clamp_b"],
                     stats=stats,times=model.OBS,seeds=np.array(job["seeds"],dtype=np.uint32),replicates=np.array(job["replicates"]))
            records.append(dict(id=job["id"],status="COMPLETE",returncode=0,job=job,raw=path.name,raw_sha256=digest(path),charged_worker_seconds=0))
        complete_ledger=dict(status="SYNTHETIC ONLY",records=records,charged_worker_seconds=0,stop_reason=None)
        atomic_json(folder/"production/ledger.json",complete_ledger);atomic_json(folder/"complete_ledger.json",complete_ledger)
        evaluator.BASE=folder;checker.BASE=folder
        evaluator.evaluate();checker.main()
        full=json.loads((folder/"production/evaluation.json").read_text(encoding="utf-8"))
        assert full["complete"] and full["completed_paths"]==1024 and full["verdict"]=="MATCH CAP NOT QUALIFIED"
        assert all(v["matches"]==2 and v["n"]==128 and v["upper"]>.05 and v["scored"] for v in full["primary"])
        for g in full["groups"].values():
            if g["start"]:assert g["match_counts"][-1]==125 and g["n"]==128
        for name in ("evaluation.json","independent_check.json"):
            shutil.copy2(folder/"production"/name,folder/("complete_"+name))
        shutil.copy2(folder/"MG10P2_results.md",folder/"complete_report.md")
        last=matrix[-1];path=folder/"production"/(last["id"]+".npz");shutil.copy2(path,folder/"complete_last_chunk.npz")
        with np.load(path) as a:arrays={k:a[k].copy() for k in a.files}
        for key in ("data","stats","seeds","replicates"):arrays[key]=arrays[key][:-1]
        save_npz(path,**arrays);records[-1]["status"]="TIMEOUT";records[-1]["returncode"]=1
        atomic_json(folder/"production/ledger.json",dict(status="SYNTHETIC ONLY",records=records,charged_worker_seconds=0,stop_reason="TIMEOUT"))
        evaluator.evaluate();checker.main()
        partial=json.loads((folder/"production/evaluation.json").read_text(encoding="utf-8"))
        assert partial["completed_paths"]==1023 and partial["verdict"]=="INCOMPLETE" and not partial["complete"]
        assert all(not v["scored"] and v["cap_met"] is None for v in partial["primary"])
        assert partial["groups"][f"c{last['condition']}_s{last['start']}_k{last['step']}"]["n"]==127
        elapsed=time.perf_counter()-started
        prep=json.loads((BASE/"preparation/independent_check.json").read_text(encoding="utf-8"))["development_charged_seconds"]
        result=dict(passed=True,status="FULL SYNTHETIC P2 PIPELINE QUALIFICATION PASS; NO STOCHASTIC PRODUCTION",elapsed_seconds=elapsed,
            total_preparation_allowance_seconds=prep+elapsed+5,
            checks=["1024 algebraic paths/128jobs reconstructed;4096 direct instrument frames","Exact primary cap failure with2/128 matches per H group","Seeded retention125/128 and binomial-CDF-inverted intervals","Signed GXY projections/residuals/means/ranges/t summaries and reference/noise proxies","1023-path partial prefix retained; incomplete family not scored","Both descriptive plots/reports and complete-case backups retained"],
            data_scope="All fixture paths are algebraic arrays with fixed seed labels only; no RNG draws or scientific trajectories. Complete and deliberately incomplete phases retained.")
        assert result["total_preparation_allowance_seconds"]<900
        atomic_json(output,result);print(json.dumps(result,indent=2),flush=True)
    except BaseException:atomic_json(output,dict(passed=False,status="QUALIFICATION FAILED; PRODUCTION BLOCKED",error=traceback.format_exc(),elapsed_seconds=time.perf_counter()-started));raise


if __name__=="__main__":main()

"""MG066 quarter-stream regression and full/partial algebraic pipeline checks."""
import os
for key in ("OPENBLAS_NUM_THREADS","MKL_NUM_THREADS","OMP_NUM_THREADS","NUMEXPR_NUM_THREADS"):
    os.environ[key]="1"
from pathlib import Path
import argparse,json,shutil,time,traceback
import numpy as np
from mg10d1_pool import run,atomic_json
from mg10p2 import fixture
from mg10p3 import OBS
from run_mg10p3 import jobs,save_npz,digest
import evaluate_mg10p3 as evaluator
import check_mg10p3 as checker

ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/"results/MG10P3"


def stream_check(job,directory):
    from mg10d2_stream import seed_rng,advance
    from check_mg10d1_development import check_raw
    directory=Path(directory);s,k=job["start"],job["step"]
    source=BASE/f"preparation/timing/timing_s{s}_k{k}_r0.npz"
    with np.load(source) as old:expected={key:old[key].copy() for key in old.files}
    case=fixture(0,s);state=case["initial"].copy();flux=np.zeros((7,500),dtype=np.int64);hops=np.zeros((6,500),dtype=np.int64);external=np.zeros((2,500),dtype=np.int64)
    times=np.arange(5)/4;data=np.zeros((1,5,22,500),dtype=np.int64);stats=np.zeros((1,3));seed_rng(int(expected["seed"]));previous=0.
    for i,target in enumerate(times):
        d=advance(state,flux,hops,external,case["a0"],case["b0"],1000.,case["diff"],np.ones(7),(1/24000,1/48000)[k],previous,target,True,True)
        previous=target;stats[0,:2]+=d[:2];stats[0,2]=max(stats[0,2],d[2]);data[0,i]=np.concatenate([state,flux,hops,external])
    assert np.array_equal(data,expected["data"]) and np.array_equal(stats,expected["stats"])
    path=directory/(job["id"]+".npz");save_npz(path,data=data,initial=case["initial"],a0=case["a0"],b0=case["b0"],omega=1000.,clamp_a=True,clamp_b=True,stats=stats,times=times,seed=expected["seed"])
    checked=check_raw(path)
    atomic_json(directory/(job["id"]+".json"),dict(status="COMPLETE",job=job,bit_identical_to=str(source.relative_to(ROOT)),reference_sha256=digest(source),raw=checked))


def qualify():
    started=time.perf_counter();output=BASE/"preparation/pipeline_qualification.json"
    if output.exists():raise ValueError("Already attempted; preserve any correction")
    prep=json.loads((BASE/"preparation/independent_check.json").read_text(encoding="utf-8"))
    try:
        matrix=[dict(id=f"stream_s{s}_k{k}",start=s,step=k) for s in range(2) for k in range(2)]
        stream_ledger=run(matrix,BASE/"preparation/quarter_stream_checks",Path(__file__).resolve(),ROOT,900,workers=4,prior=prep["total_preparation_allowance_seconds"],timeout=150,poll_seconds=.25)
        if stream_ledger["stop_reason"] or any(v["returncode"]!=0 for v in stream_ledger["records"]):raise ValueError("Quarter-stream qualification stopped")
        analysis_started=time.perf_counter()
        folder=BASE/"preparation/SYNTHETIC_pipeline_fixture";(folder/"production").mkdir(parents=True,exist_ok=True)
        matrix=jobs();records=[]
        atomic_json(folder/"registration.json",dict(status="SYNTHETIC ONLY; NO SCIENTIFIC DRAWS",source_hashes={},jobs=matrix,expected_paths=64,ceiling_worker_seconds=86400))
        with np.load(ROOT/"results/MG003/ref_profiles.npz") as raw:profile=raw["h0.2_b29.0"][1:]
        def shifted(f,pos):return np.fft.ifft(np.fft.fft(f)*np.exp(-2j*np.pi*np.fft.fftfreq(500)*pos/.2)).real
        for job in matrix:
            case=fixture(0,job["start"]);planted=fixture(0,1)["initial"][1]
            data=np.zeros((1,101,22,500),dtype=np.int64);positive=np.zeros(500,dtype=np.int64);negative=np.zeros(500,dtype=np.int64);last=case["initial"][1].copy()
            for i,t in enumerate(OBS):
                target=case["initial"][1].copy();rep=job["replicate"]
                if not job["start"] and (rep==0 and t>=5 or rep==1 and t>=5.25):target=planted
                if job["start"] and (rep==0 and t==1 or rep==1 and t>=20.25):target=np.full(500,round(1000*case["background"][1]),dtype=np.int64)
                if job["start"] and rep==2:target=np.rint(1000*(case["background"][1]+shifted(profile[1],.05*t))).astype(np.int64)
                delta=target-last;positive+=np.maximum(delta,0);negative+=np.maximum(-delta,0);last=target
                source=rep*i
                data[0,i,:7]=case["initial"];data[0,i,0]+=source-(target-case["initial"][1]);data[0,i,1]=target
                data[0,i,7]=source;data[0,i,8]=positive;data[0,i,9]=negative;data[0,i,20]=source
            assert np.all(data>=0)
            path=folder/"production"/(job["id"]+".npz");stats=np.array([[int(round(25/(1/24000,1/48000)[job["step"]])),0,0]])
            save_npz(path,data=data,initial=case["initial"],a0=case["a0"],b0=case["b0"],omega=1000.,clamp_a=True,clamp_b=True,stats=stats,times=OBS,seed=job["seed"])
            records.append(dict(id=job["id"],status="COMPLETE",returncode=0,job=job,raw=path.name,raw_sha256=digest(path),charged_worker_seconds=0))
        complete_ledger=dict(status="SYNTHETIC ONLY",records=records,charged_worker_seconds=0,stop_reason=None)
        atomic_json(folder/"production/ledger.json",complete_ledger);atomic_json(folder/"complete_ledger.json",complete_ledger)
        evaluator.BASE=folder;checker.BASE=folder;evaluator.evaluate();checker.main()
        complete=json.loads((folder/"production/evaluation.json").read_text(encoding="utf-8"))
        assert complete["complete"] and complete["completed_paths"]==64
        for key,g in complete["groups"].items():
            assert g["n"]==16 and g["events20"]==(15 if g["start"] else 1)
            assert g["fixed"]["n"]==16 and g["aligned"]["n"]==(15 if g["start"] else 0)
        for k in range(2):
            censored=complete["paths"][f"s1_k{k}_r0"]
            assert not censored["event20"] and censored["any_observed_episode20"] and censored["aligned"] is None
            late=complete["paths"][f"s1_k{k}_r1"];assert late["event20"] and late["track"]["matched_through"]==20 and late["aligned"] is not None
            moving=complete["paths"][f"s1_k{k}_r2"];assert np.linalg.norm(np.array(moving["fixed"]["profile"])[1]-np.array(moving["aligned"]["profile"])[1])>1e-3
        for name in ("evaluation.json","independent_check.json"):shutil.copy2(folder/"production"/name,folder/("complete_"+name))
        shutil.copy2(folder/"MG10P3_results.md",folder/"complete_report.md")
        last=matrix[-1];path=folder/"production"/(last["id"]+".npz");shutil.copy2(path,folder/"complete_last_path.npz")
        with np.load(path) as raw:arrays={k:raw[k].copy() for k in raw.files}
        arrays["data"]=arrays["data"][:,:-1];arrays["times"]=arrays["times"][:-1];arrays["stats"][0,0]=round(24.75/(1/24000,1/48000)[last["step"]])
        save_npz(path,**arrays);records[-1]["status"]="TIMEOUT";records[-1]["returncode"]=1
        atomic_json(folder/"production/ledger.json",dict(status="SYNTHETIC ONLY",records=records,charged_worker_seconds=0,stop_reason="TIMEOUT"))
        evaluator.evaluate();checker.main()
        partial=json.loads((folder/"production/evaluation.json").read_text(encoding="utf-8"))
        assert partial["completed_paths"]==63 and partial["verdict"]=="INCOMPLETE" and last["id"] not in partial["paths"]
        assert partial["groups"]["s1_k1"]["n"]==15 and partial["groups"]["s1_k1"]["expected_n"]==16
        total=stream_ledger["charged_worker_seconds"]+time.perf_counter()-analysis_started+5
        assert total<900
        result=dict(passed=True,status="P3 STREAM/COMPLETE/PARTIAL SYNTHETIC PIPELINE PASS; NO SCIENTIFIC PRODUCTION",elapsed_seconds=time.perf_counter()-started,
            total_preparation_allowance_seconds=total,quarter_stream_ledger=stream_ledger,
            checks=["Four T1 quarter-grid streams bit-identical to held development chunks","64 full paths and6464 direct instrument frames","S temporary gap does not regain initial event; late20.25censor still eligible","H onset5 qualifies20; onset5.25 does not","Conditional15/16 aligned versus16/16 fixed and moving X-only alignment","Direct Dirichlet profiles/gains/signed metrics, exactn16 intervals and masked missing gains","63-path case excludes24.75prefix despite available20window; full backups retained"],
            synthetic_scope="Only algebraic arrays with production seed labels; no production RNG draws. All original expected/development inputs retained.")
        atomic_json(output,result);print(json.dumps({k:v for k,v in result.items() if k!="quarter_stream_ledger"},indent=2),flush=True)
    except BaseException:atomic_json(output,dict(passed=False,status="QUALIFICATION FAILED; SCIENTIFIC PRODUCTION BLOCKED",error=traceback.format_exc(),elapsed_seconds=time.perf_counter()-started));raise


if __name__=="__main__":
    p=argparse.ArgumentParser();p.add_argument("--worker",nargs=2);args=p.parse_args()
    if args.worker:
        j=json.loads(args.worker[0]);d=Path(args.worker[1])
        try:stream_check(j,d)
        except BaseException:atomic_json(d/(j["id"]+".json"),dict(status="FAILED",job=j,error=traceback.format_exc()));raise
    else:qualify()

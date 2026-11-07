"""MG066 full raw, temporal, direct-average and per-trajectory stats audit."""
import os
for key in ("OPENBLAS_NUM_THREADS","MKL_NUM_THREADS","OMP_NUM_THREADS","NUMEXPR_NUM_THREADS"):
    os.environ[key]="1"
from pathlib import Path
import json,hashlib
import numpy as np
from mg10p3_oracle import summarize
from mg10p2_oracle import circular
from check_mg10p2 import interval,summary
from check_mg10d1_development import check_raw

ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/"results/MG10P3"


def compare_average(actual,saved):
    if actual is None:assert saved is None;return
    assert saved is not None and actual["alignment"]==saved["alignment"] and actual["frames"]==saved["frames"]==81
    for key in ("profile","projection","residual_rms","contrast_rms","mean_individual_residual_rms","individual_projection","individual_residual_rms"):
        np.testing.assert_allclose(actual[key],saved[key],rtol=1e-8,atol=1e-9)
    for a,b in zip(actual["gain"],saved["gain"]):
        if a is None:assert b is None
        else:np.testing.assert_allclose(a,b,rtol=1e-7,atol=1e-7)


def compare_stats(values,saved):
    assert saved["n"]==len(values)
    if not values:
        assert all(saved[k] is None for k in ("mean","min","max","ci95_approx"));return
    mean,mini,maxi,ci=summary(values)
    for a,b in ((mean,saved["mean"]),(mini,saved["min"]),(maxi,saved["max"])):np.testing.assert_allclose(a,np.array(b).ravel(),rtol=1e-7,atol=1e-8)
    if ci is None:assert saved["ci95_approx"] is None
    else:np.testing.assert_allclose(ci,np.array(saved["ci95_approx"]).reshape(2,-1),rtol=1e-7,atol=1e-8)


def main():
    r=json.loads((BASE/"registration.json").read_text(encoding="utf-8"));e=json.loads((BASE/"production/evaluation.json").read_text(encoding="utf-8"));l=json.loads((BASE/"production/ledger.json").read_text(encoding="utf-8"))
    assert len(r["jobs"])==64 and r["ceiling_worker_seconds"]==86400 and r["expected_paths"]==64
    for name,expected in r["source_hashes"].items():
        path=ROOT/name;raw=path.read_bytes() if name.endswith(".npz") else path.read_text(encoding="utf-8").replace("\r\n","\n").encode()
        assert hashlib.sha256(raw).hexdigest()==expected,name
    with np.load(ROOT/"results/MG003/ref_profiles.npz") as raw:profile=raw["h0.2_b29.0"][1:]
    records={v["id"]:v for v in l["records"]};raw_checks=[];paths={};observations=0;average_checks=0;seed_values=[]
    for job in r["jobs"]:
        for segment,key in ((0,"seed"),(1,"warmup_seed")):
            seed=int(np.random.SeedSequence(entropy=20261027,spawn_key=(10,job["start"],job["step"],6503,job["replicate"],segment)).generate_state(1,dtype=np.uint32)[0]);assert seed==job[key];seed_values.append(seed)
        file=BASE/"production"/(job["id"]+".npz")
        if not file.exists():assert job["id"] not in e["paths"];continue
        checked=check_raw(file);raw_checks.append(checked)
        with np.load(file) as raw:
            data=raw["data"];n=data.shape[1];times=raw["times"]
            assert data.shape==(1,n,22,500) and 1<=n<=101 and np.array_equal(times,np.arange(n)/4)
            assert int(raw["seed"])==job["seed"] and raw["clamp_a"] and raw["clamp_b"]
            background=np.array([14/.9,14/.9,29/(14/.9)])
            fields=np.repeat(background[:,None],500,axis=1)+(profile if job["start"] else 0)
            initial=np.zeros((7,500),dtype=np.int64);initial[:3]=np.rint(fields*1000).astype(np.int64);initial[3]=4200000;initial[4]=135333333
            assert np.array_equal(initial,raw["initial"]) and np.array_equal(data[0,0,:7],initial) and not np.any(data[:,0,7:])
            assert np.all(raw["a0"]==4200000) and np.all(raw["b0"]==135333333)
            assert raw["stats"][0,0]>=int(times[-1]/(1/24000,1/48000)[job["step"]])-1
            record=records.get(job["id"],{})
            if record.get("status")!="COMPLETE" or record.get("returncode")!=0 or n!=101:
                assert job["id"] not in e["paths"];continue
            assert checked["sha256"]==record["raw_sha256"]
            actual=summarize(data[0],times,job["start"]);saved=e["paths"][job["id"]]
            assert actual["event20"]==saved["event20"] and actual["episodes"]==saved["episodes"] and actual["any_observed_episode20"]==saved["any_observed_episode20"]
            if actual["track"] is None:assert saved["track"] is None
            else:
                for key in ("status","matched_through","visible_duration","persistent_reference_signal","censor_reason"):assert actual["track"][key]==saved["track"][key],(job["id"],key)
            for a,b in zip(actual["observations"],saved["observations"]):
                assert a["status"]==b["status"]
                for key in ("amplitude","error","local_energy"):
                    if a[key] is None:assert b[key] is None
                    else:np.testing.assert_allclose(a[key],b[key],rtol=1e-9,atol=1e-9)
                if a["position"] is None:assert b["position"] is None
                else:assert circular(a["position"],b["position"])<=1e-7
                observations+=1
            compare_average(actual["fixed"],saved["fixed"]);average_checks+=1
            compare_average(actual["aligned"],saved["aligned"]);average_checks+=actual["aligned"] is not None
            paths[job["id"]]=dict(job=job,**actual)
    assert len(set(seed_values))==128 and len(paths)==e["completed_paths"]
    complete=len(paths)==64 and len(records)==64 and not l["stop_reason"] and all(v.get("status")=="COMPLETE" and v.get("returncode")==0 for v in records.values())
    assert e["complete"]==complete and e["verdict"]==("COMPLETE TEMPORAL ESTIMATION" if complete else "INCOMPLETE")
    bundles=0
    for start in range(2):
        for step in range(2):
            items=[p for p in paths.values() if (p["job"]["start"],p["job"]["step"])==(start,step)];g=e["groups"][f"s{start}_k{step}"];n=len(items);k=sum(p["event20"] for p in items)
            assert (g["n"],g["expected_n"],g["events20"])==(n,16,k)
            if n:
                assert abs(g["fraction"]-k/n)<1e-12;np.testing.assert_allclose(interval(k,n),g["ci95_descriptive"],atol=1e-10)
            else:assert g["fraction"] is None and g["ci95_descriptive"] is None
            if start:
                assert [p["track"]["matched_through"] for p in items]==g["matched_through"]
                reasons={}
                for p in items:
                    reason=p["track"]["censor_reason"] or p["track"]["status"];reasons[reason]=reasons.get(reason,0)+1
                assert reasons==g["censor_counts"]
            for key in ("fixed","aligned"):
                selected=[p[key] for p in items if p[key] is not None];bundle=g[key];assert bundle["n"]==len(selected)
                if not selected:
                    assert bundle["profile_mean"] is None and not bundle["stats"]
                else:
                    np.testing.assert_allclose(np.mean([v["profile"] for v in selected],axis=0),bundle["profile_mean"],rtol=1e-8,atol=1e-9)
                    for metric in ("projection","residual_rms","contrast_rms","mean_individual_residual_rms"):compare_stats([v[metric] for v in selected],bundle["stats"][metric])
                for j in range(3):compare_stats([v["gain"][j] for v in selected if v["gain"][j] is not None],bundle["gain"][j])
                bundles+=1
            assert g["fixed"]["n"]==n and g["aligned"]["n"]==(k if start else 0)
    result=dict(status="INDEPENDENT RAW/TEMPORAL/DIRECT-AVERAGE/STATISTICAL AUDIT PASS",verdict=e["verdict"],complete_paths=len(paths),
        raw_archives=raw_checks,instrument_frame_comparisons=observations,temporal_average_comparisons=average_checks,group_average_bundles=bundles,
        charged_worker_seconds=l["charged_worker_seconds"],ceiling_overshoot_seconds=max(0,l["charged_worker_seconds"]-86400),
        dependencies="Same author; shared reference/constants/NumPy/SciPy and earlier explicit integer reconstruction. Separate dense/direct classification, candidate-start episodes, real Dirichlet interpolation/averaging, binomial-CDF inversion and trajectory statistics; not physical lineage or external replication.")
    (BASE/"production/independent_check.json").write_text(json.dumps(result,indent=2)+"\n",encoding="utf-8")
    print(f"Independent P3 audit PASS: {len(paths)} paths, {observations} frames, {average_checks} temporal averages",flush=True)


if __name__=="__main__":main()

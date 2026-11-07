"""MG053 independent D1 raw arithmetic and descriptive-statistics audit.

Shares prior explicit integer reconstruction and NumPy/SciPy libraries;
does not import the D1 measurement/evaluation functions or sampler.
"""
import os
for key in ("OPENBLAS_NUM_THREADS","MKL_NUM_THREADS","OMP_NUM_THREADS","NUMEXPR_NUM_THREADS"):
    os.environ[key]="1"
from pathlib import Path
import hashlib, json, math
import numpy as np
from scipy.stats import t as student_t
from check_mg10d1_development import check_raw

ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/"results/MG10D1"


def statistics(values):
    values=np.asarray(values,dtype=float);n=len(values)
    mean=np.array([math.fsum(values[:,i])/n for i in range(values.shape[1])])
    interval=None
    if n>=2:
        variance=np.array([math.fsum((float(v)-mean[i])**2 for v in values[:,i])/(n-1) for i in range(values.shape[1])])
        half=float(student_t.ppf(.975,n-1))*np.sqrt(variance/n)
        interval=np.stack([mean-half,mean+half])
    return mean,interval


def independent_metrics(data,a0,b0):
    counts=data[0,:,:3].astype(float)
    means=np.sum(counts,axis=-1)/500/1000
    # Center integer-valued counts before scaling, independent of the
    # evaluator's reference-site centering in concentration units.
    centered=counts-counts.mean(axis=-1,keepdims=True)
    rms=np.sqrt(np.sum(centered**2,axis=-1)/500)/1000
    ratios=data[0,:,3:5]/np.stack([a0,b0])[None]
    return dict(field_mean=means,field_rms=rms,availability_mean=ratios.sum(axis=-1)/500,
        availability_min=np.amin(ratios,axis=-1),availability_max=np.amax(ratios,axis=-1),
        availability_sd=np.sqrt(np.mean((ratios-ratios.mean(axis=-1,keepdims=True))**2,axis=-1)),
        reaction_totals=data[0,:,7:14].sum(axis=-1),external_totals=data[0,:,20:22].sum(axis=-1))


def main():
    registration=json.loads((BASE/"registration.json").read_text(encoding="utf-8"))
    evaluation=json.loads((BASE/"production/evaluation.json").read_text(encoding="utf-8"))
    ledger=json.loads((BASE/"production/ledger.json").read_text(encoding="utf-8"))
    for name,expected in registration["source_hashes"].items():
        path=ROOT/name
        actual=hashlib.sha256(path.read_bytes() if name.endswith(".npz") else path.read_text(encoding="utf-8").replace("\r\n","\n").encode()).hexdigest()
        assert actual==expected,("source changed",name)
    assert registration["ceiling_worker_seconds"]==216000 and registration["workers"]==8
    assert len(registration["jobs"])==80 and len({j["id"] for j in registration["jobs"]})==80
    assert sum(registration["forecast_seconds"].values())<=.75*216000
    for j in registration["jobs"]:
        seed=int(np.random.SeedSequence(entropy=20261027,spawn_key=(10,2*j["arm"]+j["start"],j["step"],5203,j["replicate"],0)).generate_state(1,dtype=np.uint32)[0])
        assert seed==j["seed"]
    raw_records=[];metrics={};records={r["id"]:r for r in ledger["records"]}
    with np.load(ROOT/"results/MG003/ref_profiles.npz") as ref:phi=ref["h0.2_b29.0"][1:]
    for job in registration["jobs"]:
        path=BASE/"production"/(job["id"]+".npz")
        if not path.exists():continue
        checked=check_raw(path);raw_records.append(checked)
        with np.load(path) as raw:
            expected=np.zeros((7,500),dtype=np.int64)
            values=np.repeat(np.array([140/9,140/9,29/(140/9)])[:,None],500,axis=1)
            if job["start"]:values+=phi
            expected[:3]=np.rint(1000*values).astype(np.int64)
            factor=1000000 if job["arm"]==1 else 1
            expected[3]=round(1000*14*300)*factor
            expected[4]=round(1000*29*(140/9)*300)*factor
            assert np.array_equal(raw["initial"],expected)
            assert np.array_equal(raw["a0"],expected[3]) and np.array_equal(raw["b0"],expected[4])
            assert raw["clamp_a"]==(job["arm"] in (0,4)) and raw["clamp_b"]==(job["arm"] in (0,3))
            assert int(raw["seed"])==job["seed"]
            assert np.array_equal(raw["times"],np.arange(raw["data"].shape[1]))
            rec=records.get(job["id"])
            if not rec or rec.get("status")!="COMPLETE" or rec.get("returncode")!=0:
                assert job["id"] not in evaluation["paths"];continue
            assert checked["sha256"]==rec["raw_sha256"]
            assert raw["data"].shape==(1,61,22,500)
            assert raw["stats"][0,0]>=int(60/job["tau"])-1 and 0<=raw["stats"][0,1]<=raw["stats"][0,0]
            metric=independent_metrics(raw["data"],raw["a0"],raw["b0"])
            for name,value in metric.items():
                np.testing.assert_allclose(value,evaluation["paths"][job["id"]]["metrics"][name],rtol=2e-12,atol=2e-11)
            metrics[job["id"]]=metric
    assert len(metrics)==evaluation["completed_paths"]
    group_checks=0
    for arm in range(5):
        for start in range(2):
            for step in range(2):
                selected=[metrics[j["id"]] for j in registration["jobs"] if (j["arm"],j["start"],j["step"])==(arm,start,step) and j["id"] in metrics]
                group=evaluation["groups"][f"a{arm}_s{start}_k{step}"]
                assert group["n"]==len(selected)
                if not selected:continue
                for name in selected[0]:
                    series=np.array([m[name] for m in selected])
                    averages=np.array([[math.fsum((path[t,i]+path[t+1,i])/2 for t in range(60))/60 for i in range(path.shape[1])] for path in series])
                    for label,values in (("time_series",series.reshape(len(series),-1)),("endpoint",series[:,-1]),("time_average",averages)):
                        mean,interval=statistics(values)
                        target=group["summary"][name][label]
                        np.testing.assert_allclose(mean,np.array(target["mean"]).ravel(),rtol=2e-12,atol=2e-11)
                        np.testing.assert_allclose(values.min(axis=0),np.array(target["min"]).ravel(),rtol=2e-12,atol=2e-11)
                        np.testing.assert_allclose(values.max(axis=0),np.array(target["max"]).ravel(),rtol=2e-12,atol=2e-11)
                        if interval is not None:np.testing.assert_allclose(interval,np.array(target["ci95_approx"]).reshape(2,-1),rtol=2e-12,atol=2e-11)
                        else:assert target["ci95_approx"] is None
                        group_checks+=1
    for contrast in evaluation["contrasts"].values():
        left=evaluation["groups"][contrast["left"]];right=evaluation["groups"][contrast["right"]]
        assert contrast["n_left"]==left["n"] and contrast["n_right"]==right["n"]
        if not left["n"] or not right["n"]:assert not contrast["mean_difference"];continue
        for metric,differences in contrast["mean_difference"].items():
            for label,value in differences.items():
                expected=np.array(left["summary"][metric][label]["mean"])-np.array(right["summary"][metric][label]["mean"])
                np.testing.assert_allclose(value,expected,rtol=2e-12,atol=2e-11)
    complete=len(metrics)==80 and not ledger["stop_reason"]
    expected_verdict="COMPLETE DESCRIPTIVE PILOT" if complete else "BUDGET-SHORT / INCOMPLETE" if "BUDGET" in str(ledger["stop_reason"]) else "INCOMPLETE"
    assert evaluation["verdict"]==expected_verdict
    # Polling can overshoot the cap by at most active workers times the observed
    # polling/termination interval; report actual charge rather than hide it.
    result=dict(status="INDEPENDENT RAW AND SUMMARY CHECK PASS",verdict=expected_verdict,completed_paths=len(metrics),
                raw_records=raw_records,summary_checks=group_checks,charged_worker_seconds=ledger["charged_worker_seconds"],
                ceiling_overshoot_seconds=max(0,ledger["charged_worker_seconds"]-216000),
                dependencies="Same author; shares earlier explicit integer reconstruction and NumPy/SciPy/theoretical definitions, not D1 sampler/measurement/evaluation functions. Not external replication.")
    (BASE/"production/independent_check.json").write_text(json.dumps(result,indent=2)+"\n",encoding="utf-8")
    print(f"Independent checker PASS: {len(metrics)} complete paths, {group_checks} summaries, {len(raw_records)} raw archives",flush=True)


if __name__=="__main__":main()

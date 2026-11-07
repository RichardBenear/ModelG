"""MG062 separate raw/state/metric/classification and exact-bound recomputation."""
import os
for key in ("OPENBLAS_NUM_THREADS","MKL_NUM_THREADS","OMP_NUM_THREADS","NUMEXPR_NUM_THREADS"):
    os.environ[key]="1"
from pathlib import Path
import json,hashlib,math
import numpy as np
from scipy.stats import binom,t as student_t
from scipy.optimize import brentq
from check_mg10d1_development import check_raw
from mg10p2_oracle import independent_measure,circular

ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/"results/MG10P2"


def upper(k,n,alpha):
    if k==n:return 1.
    if k==0:return -math.expm1(math.log(alpha)/n)
    return brentq(lambda p:float(binom.cdf(k,n,p))-alpha,0.,1.,xtol=1e-14)


def interval(k,n):
    lo=0. if k==0 else brentq(lambda p:float(binom.sf(k-1,n,p))-.025,0.,1.,xtol=1e-14)
    return [lo,upper(k,n,.025)]


def summary(values):
    a=np.asarray(values,dtype=float);n=len(a);flat=a.reshape(n,-1)
    mean=np.array([math.fsum(flat[:,i])/n for i in range(flat.shape[1])])
    ci=None
    if n>1:
        var=np.array([math.fsum((v-mean[i])**2 for v in flat[:,i])/(n-1) for i in range(flat.shape[1])])
        half=student_t.ppf(.975,n-1)*np.sqrt(var/n);ci=np.stack([mean-half,mean+half])
    return mean,flat.min(axis=0),flat.max(axis=0),ci


def main():
    registration=json.loads((BASE/"registration.json").read_text(encoding="utf-8"));e=json.loads((BASE/"production/evaluation.json").read_text(encoding="utf-8"))
    l=json.loads((BASE/"production/ledger.json").read_text(encoding="utf-8"));records={r["id"]:r for r in l["records"]}
    assert registration["ceiling_worker_seconds"]==14400 and registration["expected_paths"]==1024 and len(registration["jobs"])==128
    for name,expected in registration["source_hashes"].items():
        p=ROOT/name;raw=p.read_bytes() if name.endswith(".npz") else p.read_text(encoding="utf-8").replace("\r\n","\n").encode()
        assert hashlib.sha256(raw).hexdigest()==expected,name
    with np.load(ROOT/"results/MG003/ref_profiles.npz") as raw:profile=raw["h0.2_b29.0"][1:]
    paths={};raw_checks=[];frame_checks=0;all_keys=set()
    for job in registration["jobs"]:
        c,s,k=job["condition"],job["start"],job["step"]
        for rep,seed in zip(job["replicates"],job["seeds"]):
            expected=int(np.random.SeedSequence(entropy=20261027,spawn_key=(10,2*c+s,k,6103,rep,0)).generate_state(1,dtype=np.uint32)[0]);assert seed==expected
            assert (c,s,k,rep) not in all_keys;all_keys.add((c,s,k,rep))
        p=BASE/"production"/(job["id"]+".npz")
        if not p.exists():continue
        checked=check_raw(p);raw_checks.append(checked)
        rec=records.get(job["id"],{})
        if rec.get("status")=="COMPLETE" and rec.get("returncode")==0:assert checked["sha256"]==rec["raw_sha256"]
        with np.load(p) as raw:
            data=raw["data"];n=len(data);assert data.shape==(n,4,22,500) and 1<=n<=8
            assert list(raw["seeds"])==job["seeds"][:n] and list(raw["replicates"])==job["replicates"][:n]
            assert np.array_equal(raw["times"],[0,.05,.1,.25]) and not np.any(data[:,0,7:])
            factor=1 if c==0 else .82;A=round(4200000*factor);B=round(135333333*factor)
            ae=14*A/4200000;be=29*B/135333333;bg=np.array([ae/.9,ae/.9,be/(ae/.9)])
            initial=np.zeros((7,500),dtype=np.int64);fields=np.repeat(bg[:,None],500,axis=1)+(profile if s else 0)
            initial[:3]=np.rint(1000*fields).astype(np.int64);initial[3]=A;initial[4]=B
            assert np.array_equal(raw["initial"],initial) and np.all(data[:,0,:7]==initial)
            assert np.all(raw["a0"]==4200000) and np.all(raw["b0"]==135333333)
            assert raw["clamp_a"]==(c==0) and raw["clamp_b"]==(c==0)
            assert np.all(raw["stats"][:,0]>=int(.25/(1/24000,1/48000)[k])-1)
            metrics,matches=independent_measure(data)
            for index,rep in enumerate(job["replicates"][:n]):
                key=f"c{c}_s{s}_k{k}_r{rep}";item=e["paths"][key]
                for name,v in metrics.items():np.testing.assert_allclose(v[index],item["metrics"][name],rtol=1e-9,atol=1e-9)
                for actual,saved in zip(matches[index],item["matches"]):
                    assert actual["status"]==saved["status"],(key,actual,saved)
                    for name in ("amplitude","error","local_energy"):
                        if actual[name] is None:assert saved[name] is None
                        else:np.testing.assert_allclose(actual[name],saved[name],rtol=1e-9,atol=1e-9)
                    if actual["position"] is None:assert saved["position"] is None
                    else:assert circular(actual["position"],saved["position"])<=1e-7
                    frame_checks+=1
                paths[key]=dict(condition=c,start=s,step=k,metrics={name:v[index] for name,v in metrics.items()},matches=matches[index])
    assert len(all_keys)==1024 and len(paths)==e["completed_paths"]
    complete=len(paths)==1024 and len(records)==128 and not l["stop_reason"] and all(v.get("status")=="COMPLETE" and v.get("returncode")==0 for v in records.values())
    assert e["complete"]==complete
    summary_checks=0;primary=[]
    for c in range(2):
        for s in range(2):
            for k in range(2):
                key=f"c{c}_s{s}_k{k}";group=e["groups"][key];items=[v for v in paths.values() if (v["condition"],v["start"],v["step"])==(c,s,k)];n=len(items)
                assert group["n"]==n and group["expected_n"]==128
                counts=[sum(v["matches"][t]["status"]=="REFERENCE_MATCH" for v in items) for t in range(4)];assert counts==group["match_counts"]
                if n:
                    np.testing.assert_allclose(group["match_rates"],np.array(counts)/n,atol=1e-12)
                    for count,saved in zip(counts,group["match_ci95_descriptive"]):np.testing.assert_allclose(interval(count,n),saved,atol=1e-10)
                    for metric in items[0]["metrics"]:
                        mean,mini,maxi,ci=summary([v["metrics"][metric] for v in items]);saved=group["summaries"][metric]
                        for a,b in ((mean,saved["mean"]),(mini,saved["min"]),(maxi,saved["max"])):np.testing.assert_allclose(a,np.array(b).ravel(),rtol=1e-9,atol=1e-9)
                        if ci is None:assert saved["ci95_approx"] is None
                        else:np.testing.assert_allclose(ci,np.array(saved["ci95_approx"]).reshape(2,-1),rtol=1e-9,atol=1e-9)
                        summary_checks+=1
                    if not s:
                        # Recompute the descriptive reference/noise proxy from
                        # independent smooth RMS means rather than trusting it.
                        _,smooth=__import__('mg10p2_oracle').tools()
                        reference=np.array([smooth(f) for f in profile]);ref_rms=np.sqrt(np.mean(reference**2,axis=1))
                        avg=np.mean([v["metrics"]["smooth_rms"] for v in items],axis=0)
                        for t in range(4):
                            for species in range(3):
                                saved=group["reference_over_mean_H_smooth_rms"][t][species]
                                if avg[t,species]==0:assert saved is None
                                else:np.testing.assert_allclose(saved,ref_rms[species]/avg[t,species],rtol=1e-8,atol=1e-8)
                else:assert group["match_rates"] is None and group["match_ci95_descriptive"] is None
                if not s:
                    bound=upper(counts[-1],n,.0125) if n else None;p=group["primary"]
                    assert p["scored"]==complete and p["matches"]==counts[-1] and p["n"]==n
                    if bound is None:assert p["upper"] is None
                    else:assert abs(bound-p["upper"])<1e-10
                    assert p["cap_met"]==(bound<=.05 if complete else None)
                    primary.append(dict(group=key,n=n,matches=counts[-1],upper=bound,cap_met=bound<=.05 if complete else None))
    verdict="INCOMPLETE" if not complete else "SHORT-HORIZON UNSEEDED MATCH CAP SUPPORTED ON TESTED CONDITIONS" if all(v["cap_met"] for v in primary) else "MATCH CAP NOT QUALIFIED"
    assert verdict==e["verdict"]
    result=dict(status="INDEPENDENT RAW, INSTRUMENT AND STATISTICAL CHECK PASS",verdict=verdict,complete_paths=len(paths),raw_archives=raw_checks,
        instrument_frame_comparisons=frame_checks,scalar_metric_summaries=summary_checks,primary=primary,charged_worker_seconds=l["charged_worker_seconds"],
        ceiling_overshoot_seconds=max(0,l["charged_worker_seconds"]-14400),
        dependencies="Same author; shares reference/constants/NumPy/SciPy and earlier explicit integer reconstruction. Dense convolution/direct dots and separate classification/metric/statistical/CP-CDF-inversion paths, not P1/P2 measurement functions. Not external replication.")
    (BASE/"production/independent_check.json").write_text(json.dumps(result,indent=2)+"\n",encoding="utf-8")
    print(f"Independent P2 audit PASS: {len(paths)} paths, {frame_checks} instrument frames, {summary_checks} scalar summaries",flush=True)


if __name__=="__main__":main()

"""Separate MG052 raw arithmetic checker; no model reconstruction helper."""
import os
for key in ("OPENBLAS_NUM_THREADS","MKL_NUM_THREADS","OMP_NUM_THREADS","NUMEXPR_NUM_THREADS"):
    os.environ[key]="1"
from pathlib import Path
import hashlib, json, math
import numpy as np
from scipy.stats import poisson, binom

ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/"results/MG10D1/development"


def check_raw(path):
    with np.load(path) as raw:
        data=raw["data"];initial=raw["initial"];stats=raw["stats"]
        pop=data[:,:,:7];events=data[:,:,7:14];external=data[:,:,20:22]
        expected=np.broadcast_to(initial,pop.shape).copy()
        # Explicit channel mutations, independent of model STOICH/einsum.
        for source,destination,channel in ((3,0,0),(0,1,1),(1,0,2),(4,5,3),(2,1,4),(1,6,5),(0,3,6)):
            expected[:,:,source]-=events[:,:,channel]
            expected[:,:,destination]+=events[:,:,channel]
        expected[:,:,1]-=events[:,:,3];expected[:,:,2]+=events[:,:,3]
        expected[:,:,3]+=external[:,:,0];expected[:,:,4]+=external[:,:,1]
        n=initial.shape[1]
        for node in range(n):
            for species in range(3):
                left=data[:,:,14+2*species,node];right=data[:,:,15+2*species,node]
                expected[:,:,species,node]-=left+right
                expected[:,:,species,(node-1)%n]+=left
                expected[:,:,species,(node+1)%n]+=right
        assert np.array_equal(pop,expected),path
        assert np.all(data>=0) and np.all(stats[:,2]<=.25+1e-14),path
        assert np.array_equal(pop.sum(axis=(2,3))-initial.sum(),external.sum(axis=(2,3))),path
        for species,key,ext in ((3,"clamp_a",0),(4,"clamp_b",1)):
            if raw[key]:assert np.all(pop[:,:,species]==initial[species]),path
            else:assert not np.any(external[:,:,ext]),path
        return dict(file=str(path.relative_to(ROOT)),sha256=hashlib.sha256(path.read_bytes()).hexdigest(),replicas=len(data),checkpoints=data.shape[1],integer_reconstruction=True)


def main():
    q=json.loads((BASE/"qualification_run/qualification.json").read_text(encoding="utf-8"))
    assert q["passed"] and all(t["passed"] for t in q["tests"])
    records=[]
    for rec in q["raw"]:
        checked=check_raw(BASE/"qualification_run"/rec["path"])
        assert checked["sha256"]==rec["sha256"]
        records.append(checked)
    marginals=[]
    for channel in (0,3):
        with np.load(BASE/f"qualification_run/isolated_channel{channel}.npz") as raw:
            sample=sorted(int(x) for x in raw["data"][:,0,7+channel,0])
        # Independent merging CDF, compared at jumps and left limits.
        target=poisson(.014) if channel==0 else binom(16,.029)
        count=0;distance=0.
        for value in sorted(set(sample)):
            before=count
            while count<len(sample) and sample[count]==value:count+=1
            distance=max(distance,abs(before/len(sample)-float(target.cdf(value-1))),abs(count/len(sample)-float(target.cdf(value))))
        upper=distance+math.sqrt(math.log(2000)/(2*len(sample)))
        original=next(t for t in q["tests"] if t["name"]==("maintained_A_poisson" if channel==0 else "maintained_B_binomial"))
        assert abs(upper-original["dkw_upper"])<1e-14 and upper<=.05
        marginals.append(dict(channel=channel,distance=distance,upper=upper))
    timing=[]
    for arm in range(5):
        for start in range(2):
            for step in range(2):
                path=BASE/f"control_timing/timing_a{arm}_s{start}_k{step}.json"
                item=json.loads(path.read_text(encoding="utf-8"))
                checked=check_raw(path.parent/item["raw"]["path"])
                assert checked["sha256"]==item["raw"]["sha256"]
                records.append(checked);timing.append(item)
    ledger=json.loads((BASE/"control_timing/ledger.json").read_text(encoding="utf-8"))
    assert not ledger["stop_reason"] and len(ledger["records"])==20
    assert all(r["returncode"]==0 for r in ledger["records"])
    assert ledger["charged_worker_seconds"]<=900
    # 4 production paths/class; scaling by T/.1, plus25% reserve for slower
    # population states. This remains a forecast, not an upper cost bound.
    compute_forecast=sum(r["compute_seconds"] for r in timing)*600*4
    adjusted=compute_forecast*1.25
    result=dict(status="DEVELOPMENT CHECK PASS; no scientific production",raw_records=records,
                qualification_tests=len(q["tests"]),marginals=marginals,
                compute_forecast_seconds=compute_forecast,slowdown_adjusted_seconds=adjusted,
                compute_forecast_core_hours=compute_forecast/3600,slowdown_adjusted_core_hours=adjusted/3600,
                eight_worker_compute_wall_hours=adjusted/3600/8,
                proposed_ceiling_core_hours=48,forecast_within_75_percent_ceiling=adjusted<=.75*48*3600,
                development_charged_seconds=ledger["charged_worker_seconds"],
                dependencies="Independent explicit integer reconstruction and merging CDF; shares raw arrays, NumPy/SciPy libraries and registered theoretical laws. Same author, not external replication.")
    (BASE/"independent_check.json").write_text(json.dumps(result,indent=2)+"\n",encoding="utf-8")
    print(json.dumps({k:v for k,v in result.items() if k!="raw_records"},indent=2))


if __name__=="__main__":main()

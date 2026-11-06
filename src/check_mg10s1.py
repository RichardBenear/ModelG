"""MG046 separate integer-ledger/CDF recomputation. No new sampler/evaluator imports."""
from pathlib import Path
import argparse,hashlib,json,math
import numpy as np
from runner import seed_for

ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/"results/MG10S1"
OBS=np.array([.05,.10,.25]);REPS=32768;M=96
NAMES=("total_X","total_B","total_W","R1_count","R3_count","R4_count","node0_X","node0_Y")


def expected_initial(case):
    fixtures=[[[16],[16],[2],[40],[40],[0],[0]],[[2],[12],[3],[8],[1],[0],[0]],
              [[2],[8],[2],[0],[0],[8],[8]],[[10,14,18],[12,20,28],[1,2,4],[20,40,10],[1,8,40],[0,0,0],[0,0,0]]]
    return np.array(fixtures[case],dtype=np.int64)


def integer_audit(data,initial):
    assert np.issubdtype(data.dtype,np.integer) and np.all(data>=0)
    f1,f2,fr,f3,f4,f5,f1r=[data[:,:,i] for i in range(7,14)]
    predicted=np.stack([initial[0]+f1-f2+fr-f1r,initial[1]+f2-fr-f3+f4-f5,
                        initial[2]+f3-f4,initial[3]-f1+f1r,initial[4]-f3,initial[5]+f3,initial[6]+f5],axis=2)
    n=initial.shape[1]
    for s in range(3):
        left,right=data[:,:,14+2*s],data[:,:,15+2*s]
        for node in range(n):predicted[:,:,s,node]+=left[:,:,(node+1)%n]+right[:,:,(node-1)%n]-left[:,:,node]-right[:,:,node]
    assert np.array_equal(predicted,data[:,:,:7])
    assert np.all(np.sum(data[:,:,:7],axis=(2,3))==initial.sum())
    assert np.all(f1r==0)


def measurements(data):
    return np.stack([np.sum(data[:,:,1],axis=-1),np.sum(data[:,:,4],axis=-1),np.sum(data[:,:,6],axis=-1),
                     np.sum(data[:,:,7],axis=-1),np.sum(data[:,:,10],axis=-1),np.sum(data[:,:,11],axis=-1),
                     data[:,:,1,0],data[:,:,2,0]],axis=2)


def distance(a,b):
    points=np.union1d(a,b)
    ca=np.searchsorted(np.sort(a),points,side="right")/len(a)
    cb=np.searchsorted(np.sort(b),points,side="right")/len(b)
    return float(np.max(np.abs(ca-cb)))


def fixtures_check():
    folder=BASE/"qualification"
    for file in folder.glob("development_case*_m*.npz"):
        case=int(file.stem.split("case")[1].split("_")[0])
        with np.load(file) as raw:integer_audit(raw["data"],expected_initial(case))
    saved=json.loads((folder/"synthetic_summary.json").read_text())
    width=math.sqrt(2*math.log(4*96/.01)/REPS)
    with np.load(folder/"synthetic_families.npz") as raw:
        for row in saved:
            differences=np.abs(raw[row["label"]+"_a"]-raw[row["label"]+"_b"])/REPS
            low=np.maximum(0,differences-width);high=np.minimum(1,differences+width)
            coverage=float(np.mean(np.all((low<=row["true_distance"])&(high>=row["true_distance"]),axis=1)))
            qualified=float(np.mean(np.all(high<=.05,axis=1)))
            assert coverage==row["coverage"] and qualified==row["qualification_frequency"]
            assert coverage>=.98
            if row["true_distance"]==0:assert qualified>=.90
            if row["true_distance"]==.10:assert qualified==0
    analytic=json.loads((folder/"A_laws.json").read_text())
    for method in (0,1):
        with np.load(folder/f"A_only_m{method}.npz") as raw:
            data=raw["data"]
            initial=np.array([[0],[0],[0],[5],[0],[0],[0]])
            integer_audit(data,initial)
            for k,t in enumerate(OBS):
                p=math.exp(-2.8*t) if method==0 else (1-2.8*.001)**round(t/.001)
                expected=np.array([sum(math.comb(5,j)*p**j*(1-p)**(5-j) for j in range(c+1)) for c in range(6)])
                actual=np.array([np.mean(data[:,k,3,0]<=c) for c in range(6)])
                error=float(np.max(abs(actual-expected)))
                row=next(r for r in analytic if r["method"]==method and abs(r["time"]-t)<1e-10)
                assert abs(error-row["cdf_error"])<1e-12 and error<=row["bound"]
    with np.load(folder/"one_step_samples.npz") as raw:
        counts,delta=raw["counts"],raw["delta"]
        target=np.array([14,2,1.2,348,396,10.8,0])*.00025
        tolerance=7*np.sqrt((target+target*target)/len(counts))+1e-12
        assert np.all(abs(counts.mean(axis=0)-target)<=tolerance)
        assert np.all((counts[:,3]==0)|(counts[:,3]==1))
        expected=np.stack([counts[:,0]-counts[:,1]+counts[:,2],counts[:,1]-counts[:,2]-counts[:,3]+counts[:,4]-counts[:,5],
                           counts[:,3]-counts[:,4],-counts[:,0],-counts[:,3],counts[:,3],counts[:,5]],axis=1)
        assert np.array_equal(delta,expected)


def check(development=False):
    fixtures_check()
    folder=BASE/("timing" if development else "production")
    observations={};complete=0;expected_reps=512 if development else REPS
    for case in range(4):
        for method in range(4):
            key=f"case{case}_method{method}";meta=folder/(key+".json");file=folder/(key+".npz")
            if not meta.exists() or not file.exists():continue
            record=json.loads(meta.read_text())
            if record["status"]!="COMPLETE":continue
            assert hashlib.sha256(file.read_bytes()).hexdigest()==record["sha256"]
            with np.load(file) as raw:
                data,initial,stats,seeds=raw["data"],raw["initial"],raw["stats"],raw["seeds"]
                assert np.array_equal(initial,expected_initial(case)) and np.array_equal(raw["observations"],OBS)
                assert data.shape==(expected_reps,3,20,initial.shape[1])
                a0=[40,8,8,40][case];b0=[40,1,8,40][case]
                assert np.all(raw["a0"]==a0) and np.all(raw["b0"]==b0)
                integer_audit(data,initial)
                if method:assert np.max(stats[:,2])<=.25+1e-14
                stage=4602 if development else 4603
                expected_seeds=[seed_for(10,case,method,stage,k) for k in range(expected_reps//512)]
                assert np.array_equal(seeds,expected_seeds)
                observations[(case,method)]=measurements(data);complete+=1
    ledger=json.loads((folder/"ledger.json").read_text())
    cost=ledger["prior_charge"]+sum(r["charged_worker_seconds"] for r in ledger["records"])
    assert abs(cost-ledger["charged_worker_seconds"])<1e-8
    if development:
        assert complete==16
        result=dict(status="PASS",groups=16,exact_integer_ledgers=True,analytic_and_synthetic_checks=True,charged_worker_seconds=cost)
    else:
        original=json.loads((folder/"evaluation.json").read_text())
        width=math.sqrt(2*math.log(384/.01)/REPS)
        assert abs(original["simultaneous_band_halfwidth"]-width)<1e-14
        rows=[]
        for case in range(4):
            if (case,0) not in observations or (case,3) not in observations:continue
            for t in range(3):
                for k,name in enumerate(NAMES):
                    d=distance(observations[(case,0)][:,t,k],observations[(case,3)][:,t,k])
                    low,high=max(0,d-width),min(1,d+width)
                    row=next(r for r in original["primary"] if r["case"]==case and r["time"]==OBS[t] and r["observable"]==name)
                    assert abs(d-row["empirical_CDF_distance"])<1e-14
                    assert abs(low-row["lower"])<1e-14 and abs(high-row["upper"])<1e-14
                    rows.append(dict(case=case,time=float(OBS[t]),observable=name,distance=d,lower=low,upper=high))
        if complete!=16 or len(rows)!=96:verdict="UNRESOLVED"
        elif all(row["upper"]<=.05 for row in rows):verdict="QUALIFIED ON TESTED SMALL SYSTEMS"
        elif any(row["lower"]>.05 for row in rows):verdict="NOT QUALIFIED ON TESTED CASES"
        else:verdict="UNRESOLVED"
        assert verdict==original["verdict"]
        result=dict(status="AGREES",verdict=verdict,complete_groups=complete,primary=rows,charged_worker_seconds=cost,
                    dependencies="Original seed helper only; no new sampler/evaluator imports. Integer stoichiometry and sorted empirical-CDF checking independently expressed. Same author, not external replication.")
    (folder/"recomputed.json").write_text(json.dumps(result,indent=2)+"\n",encoding="utf-8")
    print({k:v for k,v in result.items() if k!="primary"})


if __name__=="__main__":
    p=argparse.ArgumentParser();p.add_argument("--development",action="store_true");check(p.parse_args().development)

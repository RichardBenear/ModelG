"""Separate MG042 raw-data recomputation. No model/evaluator imports.

Shares frozen MG2 reference, detector and smoothing only. Ledger equations,
availability, error reductions, classifications and budget summation are separate.
"""
from pathlib import Path
import argparse
import hashlib
import json
import time
import numpy as np
from runner import reference
from measure import detect, smooth_periodic

ROOT=Path(__file__).resolve().parents[1]


def audit_raw(times,z,job):
    c=job["capacity"];a0=c*4200.;b0=c*29*(140/9)*300
    f1,f2,fr,f3,f4,f5,f1r=[z[:,k] for k in range(5,12)]
    hs=z[:,12:15]
    incoming=(np.roll(hs,-1,axis=-1)+np.roll(hs,1,axis=-1))/2
    lap=incoming-hs
    predicted=np.stack([f1-f2+fr+lap[:,0], f2-fr-f3+f4-f5+lap[:,1], f3-f4+lap[:,2],
                        -f1+f1r,-f3,f3,f5],axis=1)
    through=np.stack([f1+f2+fr+f1r,f2+fr+f3+f4+f5,f3+f4,f1+f1r,f3,f3,f5],axis=1)
    through[:,:3]+=hs+incoming
    if job["mode"]=="CHEMOSTAT":
        predicted[:,3]+=f1-f1r
        predicted[:,4]+=f3
        through[:,3]+=np.abs(f1-f1r)
        through[:,4]+=f3
    actual=np.stack([z[:,0]-z[0,0],z[:,1]-z[0,1],z[:,2]-z[0,2],
                     -z[:,3]*a0,-z[:,4]*b0,z[:,8],z[:,10]],axis=1)
    residual=actual-predicted
    per_site=float(np.max(abs(residual)/(1+through)))
    global_error=float(np.max(abs(np.sum(residual,axis=-1))/(1+np.sum(through,axis=-1))))
    if job["mode"]=="CHEMOSTAT": oracle=np.ones(len(times))
    elif job["method"]=="O": oracle=np.exp(-times/(c*300))
    else:
        n={"E0":1500,"E1":3000,"E2":6000}[job["method"]]
        oracle=np.exp(n*times*np.log1p(-1/(n*c*300)))
    availability=float(np.max(abs(1-z[:,3]-oracle[:,None])))
    complete=bool(np.array_equal(times,np.arange(job["T"]+1)) and z.shape==(job["T"]+1,15,job["n"])
                  and np.all(np.isfinite(z)) and np.all(z>=0) and np.all(z[:,3:5]<=1))
    return dict(per_site_scaled_residual=per_site,global_scaled_residual=global_error,
                availability_oracle_error=availability,complete=complete,
                passed=complete and per_site<=1e-8 and global_error<=1e-8 and availability<=1e-8)


def check_development():
    began=time.perf_counter()
    folder=ROOT/"results/MG10/timing"
    ledger=json.loads((folder/"ledger.json").read_text(encoding="utf-8"))
    rows=[]
    for record in ledger["records"]:
        job=record["job"]
        with np.load(folder/(job["id"]+".npz")) as raw:
            audit=audit_raw(raw["times"],raw["state"],job)
        assert audit["passed"],str((job,audit))
        rows.append(dict(id=job["id"],**audit))
    result={"development_recomputation":"PASS","rows":rows,"elapsed_worker_seconds":time.perf_counter()-began,
            "shared_dependencies":"MG2 reference, detector/smoothing only; no model/evaluator imports"}
    (folder/"recomputed.json").write_text(json.dumps(result,indent=2)+"\n",encoding="utf-8")
    print(result)


def check_production():
    folder=ROOT/"results/MG10/production"
    evaluation=json.loads((folder/"evaluation.json").read_text(encoding="utf-8"))
    manifest=json.loads((folder/"manifest.json").read_text(encoding="utf-8"))
    ref=reference(29.,.2)
    baseline=ref["S0"]
    snapshots,counts,audits={},{},[]
    for job in manifest["jobs"]:
        key=job["id"]; file=folder/(key+".npz");meta=folder/(key+".json")
        if not file.exists() or not meta.exists(): continue
        record=json.loads(meta.read_text(encoding="utf-8"))
        assert hashlib.sha256(file.read_bytes()).hexdigest()==record["sha256"]
        with np.load(file) as raw:
            t,z=raw["times"],raw["state"]
        audit=audit_raw(t,z,job);audits.append(dict(id=key,**audit))
        if record["status"]!="COMPLETE" or not audit["passed"]: continue
        snapshots[key]=z[:,:3]
        smoothed=smooth_periodic(z[:,1]-baseline[1],.2)
        counts[key]=[len(detect(v,ref["x"],.2,ref["thr"],ref["lam"])) for v in smoothed]
        assert counts[key]==evaluation["detector_counts"][key]
    def metric(a,b,scale):
        difference=a[50:]-b[50:]
        squared=np.sum(difference*difference,axis=-1)/500
        return float(np.max(np.sqrt(squared)/scale[None,:]))
    recomputed=[]
    for row in evaluation["rows"]:
        cap,start=row["capacity"],row["start"]
        keys=[f"{mode}_{start}_{method}_C{c}" for mode,c in (("FINITE",cap),("CHEMOSTAT",1)) for method in ("E2","O")]
        if not all(k in snapshots for k in keys):
            assert row["verdict"]=="NUMERICALLY UNRESOLVED"
            recomputed.append(dict(capacity=cap,start=start,verdict="NUMERICALLY UNRESOLVED"));continue
        fe,fo,ce,co=[snapshots[k] for k in keys]
        scale=baseline if start=="H" else np.sqrt(np.sum(ref["phi"]**2,axis=1)/500)
        target=0 if start=="H" else 1
        nfe,nfo,nce,nco=[np.array(counts[k][50:]) for k in keys]
        ef,ec=metric(fe,fo,scale),metric(ce,co,scale)
        e,o=metric(fe,ce,scale),metric(fo,co,scale)
        for name,value in (("primary_E2",e),("primary_O",o),("numerical_agreement_finite",ef),("numerical_agreement_control",ec)):
            assert abs(row[name]-value)<1e-10*(1+abs(value)),(name,row[name],value)
        if ef>.001 or ec>.001 or not np.array_equal(nfe,nfo) or not np.array_equal(nce,nco) or not (np.all(nce==target) and np.all(nco==target)):
            verdict="NUMERICALLY UNRESOLVED"
        elif np.all(nfe==target)!=np.all(nfo==target): verdict="UNRESOLVED AT MARGIN"
        elif not np.all(nfe==target): verdict="BASELINE NOT RECOVERED"
        elif e<=.01 and o<=.01: verdict="RECOVERS BASELINE"
        elif e>.01 and o>.01: verdict="BASELINE NOT RECOVERED"
        else: verdict="UNRESOLVED AT MARGIN"
        assert verdict==row["verdict"]
        recomputed.append(dict(capacity=cap,start=start,verdict=verdict,primary_E2=e,primary_O=o))
    if any(not a["passed"] for a in audits): verdict="QUALIFICATION FAILED"
    elif len(snapshots)!=64 or any(r["verdict"] in ("NUMERICALLY UNRESOLVED","UNRESOLVED AT MARGIN") for r in recomputed): verdict="UNRESOLVED"
    elif all(r["verdict"]=="RECOVERS BASELINE" for r in recomputed if r["capacity"] in (100000,1000000)): verdict="RECOVERS BASELINE ON THE TESTED RANGE"
    else: verdict="BASELINE NOT RECOVERED ON THE TESTED RANGE"
    assert verdict==evaluation["verdict"]
    ledger=json.loads((folder/"ledger.json").read_text(encoding="utf-8"))
    cost=sum(x["charged_worker_seconds"] for x in ledger["records"])+ledger["previous_charge"]
    assert abs(cost-ledger["charged_worker_seconds"])<1e-8
    result=dict(recomputation="AGREES",verdict=verdict,rows=recomputed,accounting=audits,
                complete_paths=len(snapshots),charged_worker_seconds=cost,
                source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                shared_dependencies="MG2 reference, detector and smoothing; no model/evaluator imports. Same-author recomputation, not external replication.")
    (folder/"recomputed.json").write_text(json.dumps(result,indent=2)+"\n",encoding="utf-8")
    print(verdict, "paths",len(snapshots),"core-hours",cost/3600)


if __name__=="__main__":
    parser=argparse.ArgumentParser();parser.add_argument("--development",action="store_true")
    if parser.parse_args().development: check_development()
    else: check_production()

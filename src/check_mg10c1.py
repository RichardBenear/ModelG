"""MG044 separate raw-data checking; no new model or evaluator imports.

Shares original MG2 reference/smoothing and MG10 independent ledger helper.
Periodic prominence is recomputed by contour walks instead of scipy.find_peaks.
"""
from pathlib import Path
import argparse,hashlib,json
import numpy as np
from runner import reference
from measure import smooth_periodic
from check_mg10 import audit_raw

ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/"results/MG10C1"
SCALE=np.array([140/9,140/9,29/(140/9),1.,1.])


def contour_count(field,ref):
    phi=field-140/9
    mean=np.mean(phi)
    y=smooth_periodic(phi-mean,.2)+mean
    eta=1e-8*max(1.,140/9,abs(ref["A_s"]))
    n=len(y);candidates=[]
    for i in range(n):
        if not (y[i]<y[(i-1)%n] and y[i]<=y[(i+1)%n] and y[i]<ref["thr"]):continue
        maxima=[]
        for sign in (-1,1):
            level=y[i]
            for step in range(1,n+1):
                v=y[(i+sign*step)%n]
                if v<y[i]:break
                level=max(level,v)
            maxima.append(level)
        if min(maxima)-y[i]>=eta:candidates.append(i)
    kept=[]
    for i in sorted(candidates,key=lambda i:y[i]):
        if all(min(abs(i-j),n-abs(i-j))*.2>=ref["lam"]/2 for j in kept):kept.append(i)
    return len(kept)


def audit_files(folder):
    records=[];data={}
    for meta in Path(folder).glob("*.json"):
        record=json.loads(meta.read_text(encoding="utf-8"))
        if record.get("status")!="COMPLETE" or "job" not in record:continue
        job=record["job"];file=Path(folder)/(job["id"]+".npz")
        assert hashlib.sha256(file.read_bytes()).hexdigest()==record["sha256"]
        with np.load(file) as raw:z={k:raw[k] for k in raw.files}
        assert np.array_equal(z["times"],np.arange(job["T"]+1))
        if job["kind"]=="full":
            audit=audit_raw(z["times"],z["state"],dict(mode="FINITE",capacity=1,method="E2" if job["method"]=="E2" else "O",T=job["T"],n=500))
            assert audit["passed"]
        elif job["kind"].startswith("reduced_"):
            assert np.max(abs(1-z["state"][:,3]-np.exp(-z["times"]/300)))<=1e-10
            assert np.max(z["interpolation_error"])<=1e-10
            # Power coefficients independently evaluated at integer snapshots.
            reconstructed=[]
            for t in z["times"]:
                k=max(0,min(len(z["knots"])-2,int(np.searchsorted(z["knots"],t)-1)))
                u=(t-z["knots"][k])/(z["knots"][k+1]-z["knots"][k])
                reconstructed.append(np.polynomial.polynomial.polyval(u,z["coefficients"][k]))
            assert np.max(abs(np.array(reconstructed)-z["state"])/SCALE)<=1e-10
        else:
            assert np.all(np.isfinite(z["response"]))
            assert np.array_equal(z["modes"],np.arange(1,251))
            # BDF dense evaluation at t=0 may round the supplied unit vector;
            # compare at the already specified linear solver absolute tolerance.
            np.testing.assert_allclose(z["response"][0,1],1,rtol=0,atol=1e-12)
        data[job["id"]]=z;records.append(record)
    ledger=json.loads((Path(folder)/"ledger.json").read_text())
    cost=ledger["prior_charge"]+sum(r["charged_worker_seconds"] for r in ledger["records"])
    assert abs(cost-ledger["charged_worker_seconds"])<1e-8
    return data,records,cost


def check(development=False):
    ref=reference(29.,.2)
    fixtures=json.loads((BASE/"qualification/detector_fixtures.json").read_text())
    with np.load(BASE/"qualification/detector_fixtures.npz") as raw:
        for row in fixtures:
            count=contour_count(raw[row["name"]],ref)
            assert count==row["expected"]==row["actual"],(row,count)
    folder=BASE/("timing" if development else "production")
    raw,records,cost=audit_files(folder)
    if development:
        assert len(raw)==6
        result=dict(status="PASS",raw_files=6,fixtures=len(fixtures),charged_seconds=cost)
        (folder/"recomputed.json").write_text(json.dumps(result,indent=2)+"\n",encoding="utf-8");print(result);return
    original=json.loads((folder/"evaluation.json").read_text())
    if len(raw)!=10:
        assert original["status"]=="INCOMPLETE"
        result=dict(status="AGREES: INCOMPLETE",complete=len(raw),planned=10)
        (folder/"recomputed.json").write_text(json.dumps(result,indent=2)+"\n",encoding="utf-8");print(result);return
    a,b=raw["reduced_primary"]["state"],raw["reduced_reference"]["state"]
    agreement=float(np.max(abs(a-b)/SCALE))
    uniform=raw["full_zero_E2"]["state"][:,:3]
    uniform_error=float(np.max(np.sqrt(np.sum((uniform-a[:,:3,None])**2,axis=2)/500)/SCALE[:3]))
    span=float(np.max(np.ptp(uniform,axis=2)))
    p,q=raw["linear_primary"]["response"],raw["linear_reference"]["response"]
    peaks=np.max(abs(p[:,1]),axis=0);other=np.max(abs(q[:,1]),axis=0)
    selected=int(np.argmax(peaks));mode=selected+1
    log_error=float(np.max(abs(np.log(peaks[peaks>=1])-np.log(other[peaks>=1]))))
    for key,value in (("reduced_normalized_error",agreement),("uniform_E2_normalized_error",uniform_error),("linear_peak_log_gain_error",log_error)):
        assert abs(original[key]-value)<=1e-10*(1+abs(value))
    assert original["selected_mode"]==mode
    window_rows=[]
    wave=np.cos(2*np.pi*mode*np.arange(500)/500)-1j*np.sin(2*np.pi*mode*np.arange(500)/500)
    for label,epsilon in (("eps6",1e-6),("eps8",1e-8)):
        eligible=(abs(p[:,1,selected])>=10)&(epsilon*np.max(abs(p[:,:3,selected]),axis=1)<=1e-3)
        predicted=epsilon*p[:,1,selected]
        for method in ("E2","Oplus"):
            key=f"full_{label}_{method}"
            x=raw[key]["state"][:,1]
            x=(x-np.mean(x,axis=1)[:,None])/(140/9)
            coefficient=np.sum(x*wave,axis=1)/500*(1 if mode==250 else 2)
            relative=abs(coefficient[eligible]-predicted[eligible])/abs(predicted[eligible])
            maximum=float(np.max(relative)) if len(relative) else float("inf")
            expected=next(w for w in original["windows"] if w["id"]==key)
            assert expected["eligible_count"]==int(eligible.sum())
            assert abs(expected["relative_error_max"]-maximum)<=1e-8*(1+abs(maximum))
            window_rows.append(dict(id=key,eligible_count=int(eligible.sum()),relative_error_max=maximum))
    for key,data in raw.items():
        if key.startswith("full_"):
            counts=[contour_count(x,ref) for x in data["state"][:,1]]
            assert counts==original["resolved_minimum_counts"][key],key
    reference_verdict="HOMOGENEOUS REFERENCE RESOLVED" if agreement<=1e-8 and uniform_error<=1e-4 and span==0 else "UNRESOLVED"
    if log_error>.02:amp="UNRESOLVED"
    elif peaks[selected]<100:amp="NO >=100x AMPLIFICATION IN THE TESTED FAMILY/HORIZON"
    elif all(w["eligible_count"]>=3 and w["relative_error_max"]<=.05 for w in window_rows):amp="CONTROLLED AMPLIFICATION SUPPORTED"
    else:amp="UNRESOLVED"
    assert reference_verdict==original["reference_verdict"]
    assert amp==original["amplification_verdict"]
    assert original["detector_verdict"]=="RESOLVED-MINIMUM QUALIFICATION PASS"
    result=dict(status="AGREES",complete=10,reference_verdict=reference_verdict,amplification_verdict=amp,
                detector_verdict="RESOLVED-MINIMUM QUALIFICATION PASS",mode=mode,peak_gain=float(peaks[selected]),
                windows=window_rows,charged_seconds=cost,source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                shared_dependencies="Original reference/smoothing and MG10 independent ledger helper. No new model/evaluator imports; contour-walk detector recomputation. Same author, not external replication.")
    (folder/"recomputed.json").write_text(json.dumps(result,indent=2)+"\n",encoding="utf-8");print(result)


if __name__=="__main__":
    p=argparse.ArgumentParser();p.add_argument("--development",action="store_true")
    check(p.parse_args().development)

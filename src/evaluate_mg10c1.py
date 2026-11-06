"""MG044 prospective follow-up evaluator, distinct from frozen MG10 scoring."""
import json
from pathlib import Path
import numpy as np
from runner import reference
from mg10c1 import BASE,SCALE,resolved_minima


def classify_reference(complete,agreement,uniform_error,uniform_range,ledgers):
    return "HOMOGENEOUS REFERENCE RESOLVED" if complete and agreement<=1e-8 and uniform_error<=1e-4 and uniform_range==0 and ledgers else "UNRESOLVED"


def classify_amplification(complete,peak,linear_ok,windows):
    if not complete or not linear_ok:return "UNRESOLVED"
    if peak<100:return "NO >=100x AMPLIFICATION IN THE TESTED FAMILY/HORIZON"
    return "CONTROLLED AMPLIFICATION SUPPORTED" if len(windows)==4 and all(w["eligible_count"]>=3 and w["relative_error_max"]<=.05 for w in windows) else "UNRESOLVED"


def evaluate(directory=None):
    directory=Path(directory or BASE/"production")
    records={p.stem:json.loads(p.read_text(encoding="utf-8")) for p in directory.glob("*.json") if p.stem not in ("ledger","selection","resolved_manifest","evaluation","recomputed")}
    raw={}
    for key,record in records.items():
        if record.get("status")=="COMPLETE" and (directory/(key+".npz")).exists():
            with np.load(directory/(key+".npz")) as z:raw[key]={k:z[k] for k in z.files}
    required=["reduced_primary","reduced_reference","linear_primary","linear_reference"]+[f"full_{label}_{method}" for label in ("zero","eps6","eps8") for method in ("E2","Oplus")]
    if any(key not in raw for key in required):
        result=dict(status="INCOMPLETE",missing=[k for k in required if k not in raw],planned=10,complete=len(raw),
                    reference_verdict="UNRESOLVED",amplification_verdict="UNRESOLVED")
        (directory/"evaluation.json").write_text(json.dumps(result,indent=2)+"\n",encoding="utf-8");return result
    for key,data in raw.items():
        if not np.array_equal(data["times"],np.arange(301)):raise ValueError(f"Missing endpoint {key}")
    agreement=float(np.max(abs(raw["reduced_primary"]["state"]-raw["reduced_reference"]["state"])/SCALE))
    zero=raw["full_zero_E2"]["state"]
    uniform_error=float(np.max(np.sqrt(np.mean((zero[:,:3]-raw["reduced_primary"]["state"][:,:3,None])**2,axis=2))/SCALE[:3]))
    uniform_range=float(np.max(np.ptp(zero[:,:3],axis=2)))
    primary=raw["linear_primary"]["response"]
    alternate=raw["linear_reference"]["response"]
    peaks=np.max(abs(primary[:,1]),axis=0);refpeaks=np.max(abs(alternate[:,1]),axis=0)
    selected=int(np.argmax(peaks));mode=int(raw["linear_primary"]["modes"][selected])
    log_error=float(np.max(abs(np.log(np.maximum(peaks,1e-300))-np.log(np.maximum(refpeaks,1e-300)))[peaks>=1]))
    windows=[]
    counts={}
    ref=reference(29.,.2)
    for label,epsilon in (("eps6",1e-6),("eps8",1e-8)):
        vector=primary[:,:,selected]
        eligible=(abs(vector[:,1])>=10)&(epsilon*np.max(abs(vector[:,:3]),axis=1)<=1e-3)
        predicted=epsilon*vector[:,1]
        for method in ("E2","Oplus"):
            key=f"full_{label}_{method}"
            x=raw[key]["state"][:,1]
            centered=(x-x.mean(axis=1,keepdims=True))/SCALE[1]
            coefficient=(1 if mode==250 else 2)*np.mean(centered*np.exp(-2j*np.pi*mode*np.arange(500)/500),axis=1)
            error=abs(coefficient[eligible]-predicted[eligible])/abs(predicted[eligible])
            windows.append(dict(id=key,eligible_count=int(eligible.sum()),eligible_times=np.where(eligible)[0].tolist(),
                                relative_error_max=float(np.max(error)) if len(error) else float("inf"),
                                observed_mode_real=coefficient.real.tolist(),observed_mode_imag=coefficient.imag.tolist()))
    for key,data in raw.items():
        if key.startswith("full_"):counts[key]=[len(resolved_minima(x,ref)) for x in data["state"][:,1]]
    ledgers=all(r.get("audit",{}).get("pass",False) for r in records.values() if r.get("job",{}).get("kind")=="full")
    fixtures=json.loads((BASE/"qualification/detector_fixtures.json").read_text())
    detector_pass=all(r["actual"]==r["expected"] for r in fixtures)
    result=dict(status="COMPLETE",planned=10,complete=10,selected_mode=mode,selected_peak_gain=float(peaks[selected]),
                reduced_normalized_error=agreement,uniform_E2_normalized_error=uniform_error,uniform_E2_max_range=uniform_range,
                linear_peak_log_gain_error=log_error,
                reference_verdict=classify_reference(True,agreement,uniform_error,uniform_range,ledgers),
                amplification_verdict=classify_amplification(True,float(peaks[selected]),log_error<=.02,windows),
                detector_verdict="RESOLVED-MINIMUM QUALIFICATION PASS" if detector_pass else "QUALIFICATION FAILED",
                windows=windows,resolved_minimum_counts=counts,
                zero_Oplus_max_active_spatial_range=np.ptp(raw["full_zero_Oplus"]["state"][:,:3],axis=2).max(axis=0).tolist())
    (directory/"evaluation.json").write_text(json.dumps(result,indent=2)+"\n",encoding="utf-8")
    return result


if __name__=="__main__":
    result=evaluate();print({k:v for k,v in result.items() if k not in ("windows","resolved_minimum_counts")})

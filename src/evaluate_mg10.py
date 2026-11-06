"""MG042 evaluator; production criteria fixed by MG040 plus adopted MG041."""
from pathlib import Path
import json
import numpy as np
from mg10_reservoir import CAPACITIES, S0, ROOT, accounting
from runner import reference
from measure import detect, smooth_periodic, fourier_shift


def classify(complete, numerical, control_ok, finite_e, finite_o, err_e, err_o):
    if not complete or not numerical or not control_ok:
        return "NUMERICALLY UNRESOLVED"
    if finite_e != finite_o:
        return "UNRESOLVED AT MARGIN"
    if not finite_e:
        return "BASELINE NOT RECOVERED"
    if err_e <= .01 and err_o <= .01:
        return "RECOVERS BASELINE"
    if err_e > .01 and err_o > .01:
        return "BASELINE NOT RECOVERED"
    return "UNRESOLVED AT MARGIN"


def complete(times, data, status="COMPLETE"):
    return bool(status == "COMPLETE" and np.array_equal(times, np.arange(301))
                and data.shape == (301,15,500) and np.all(np.isfinite(data))
                and np.all(data >= 0) and np.all(data[:,3:5] <= 1))


def run_id(mode, start, method, capacity):
    return f"{mode}_{start}_{method}_C{capacity}"


def detection_counts(data, ref):
    smoothed = smooth_periodic(data[:,1]-S0[1], .2)
    return np.array([len(detect(row, ref["x"], .2, ref["thr"], ref["lam"])) for row in smoothed])


def normalized(a, b, scale):
    return float(np.max(np.sqrt(np.mean((a[50:,:3]-b[50:,:3])**2, axis=2))/scale[None,:]))


def evaluate(directory):
    directory = Path(directory)
    ref = reference(29., .2)
    rows, paths, audits = [], {}, []
    manifest = json.loads((directory / "manifest.json").read_text(encoding="utf-8"))
    for job in manifest["jobs"]:
        key = job["id"]
        file = directory / (key+".npz")
        meta = directory / (key+".json")
        if not file.exists() or not meta.exists():
            continue
        record = json.loads(meta.read_text(encoding="utf-8"))
        with np.load(file) as z:
            times, data = z["times"], z["state"]
        if not complete(times,data,record["status"]):
            continue
        audit = accounting(times, data, job)
        audits.append({"id": key, **audit})
        paths[key] = (data, detection_counts(data,ref))
    for start in ("H", "S"):
        scale = S0 if start == "H" else np.sqrt(np.mean(ref["phi"]**2,axis=1))
        target = 0 if start == "H" else 1
        for cap in CAPACITIES:
            keys = [run_id(mode,start,method,c) for mode,c in (("FINITE",cap),("CHEMOSTAT",1)) for method in ("E2","O")]
            have = all(k in paths for k in keys)
            row = {"capacity":cap,"start":start,"complete":have}
            if not have:
                row["verdict"] = classify(False,False,False,False,False,0,0)
                rows.append(row); continue
            fe,fo,ce,co = [paths[k][0] for k in keys]
            nfe,nfo,nce,nco = [paths[k][1][50:] for k in keys]
            agreement_f,agreement_c = normalized(fe,fo,scale),normalized(ce,co,scale)
            numerical = agreement_f <= .001 and agreement_c <= .001 and np.array_equal(nfe,nfo) and np.array_equal(nce,nco)
            control = bool(np.all(nce==target) and np.all(nco==target))
            det_e,det_o = bool(np.all(nfe==target)),bool(np.all(nfo==target))
            err_e,err_o = normalized(fe,ce,scale),normalized(fo,co,scale)
            row.update(numerical_agreement_finite=agreement_f,numerical_agreement_control=agreement_c,
                       control_valid=control,finite_detector_E2=det_e,finite_detector_O=det_o,
                       primary_E2=err_e,primary_O=err_o,
                       verdict=classify(True,numerical,control,det_e,det_o,err_e,err_o))
            descriptive = {}
            for method in ("E0","E1","E2","O"):
                fk,ck = run_id("FINITE",start,method,cap),run_id("CHEMOSTAT",start,method,1)
                if fk not in paths or ck not in paths:
                    continue
                f,c = paths[fk][0],paths[ck][0]
                diff = f[50:,:3]-c[50:,:3]
                desc = {"primary":normalized(f,c,scale),
                        "max_absolute_rms_by_species":np.sqrt(np.mean(diff**2,axis=2)).max(axis=0).tolist(),
                        "max_absolute_by_species":np.max(np.abs(diff),axis=(0,2)).tolist(),
                        "max_absolute_mean_offset_by_species":np.max(np.abs(diff.mean(axis=2)),axis=0).tolist(),
                        "max_background_corrected_normalized_rms":float(np.max(np.sqrt(np.mean((diff-diff.mean(axis=2,keepdims=True))**2,axis=2))/scale))}
                # Alignment is descriptive only, and only when both paths have one particle.
                if start=="S" and np.all(paths[fk][1][50:]==1) and np.all(paths[ck][1][50:]==1):
                    aligned = []
                    for frow,crow in zip(f[50:],c[50:]):
                        positions = [detect(smooth_periodic(r[1]-S0[1],.2),ref["x"],.2,ref["thr"],ref["lam"])[0] for r in (frow,crow)]
                        d = fourier_shift(frow[:3],.2,positions[0])-fourier_shift(crow[:3],.2,positions[1])
                        aligned.append(np.sqrt(np.mean(d*d,axis=1))/scale)
                    desc["aligned_primary_descriptive"] = float(np.max(aligned))
                descriptive[method] = desc
            row["descriptive"] = descriptive
            rows.append(row)
    valid_audits = len(audits)==64 and all(x["pass"] for x in audits)
    top = [r for r in rows if r["capacity"] in (100000,1000000)]
    if not valid_audits:
        verdict = "QUALIFICATION FAILED" if any(not a["pass"] for a in audits) else "UNRESOLVED"
    elif len(paths)!=64 or any(r["verdict"] in ("NUMERICALLY UNRESOLVED","UNRESOLVED AT MARGIN") for r in rows):
        verdict = "UNRESOLVED"
    elif all(r["verdict"]=="RECOVERS BASELINE" for r in top):
        verdict = "RECOVERS BASELINE ON THE TESTED RANGE"
    else:
        verdict = "BASELINE NOT RECOVERED ON THE TESTED RANGE"
    result = {"verdict":verdict,"complete_trajectories":len(paths),"planned_trajectories":64,
              "rows":rows,"accounting":audits,"detector_counts":{k:v[1].tolist() for k,v in paths.items()}}
    (directory / "evaluation.json").write_text(json.dumps(result,indent=2)+"\n",encoding="utf-8")
    return result


if __name__ == "__main__":
    result = evaluate(ROOT / "results/MG10/production")
    print(result["verdict"])
    for row in result["rows"]:
        print(row["start"],row["capacity"],row["verdict"],row.get("primary_E2"),row.get("primary_O"))

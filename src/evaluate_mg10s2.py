"""MG049 two-family spatial-observable CDF qualification."""
from pathlib import Path
import json
import numpy as np
from mg10s2_measure import NAMES

ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/"results/MG10S2"
R=16384;ALPHA=.005;MARGIN=.08


def distance(a,b):
    points=np.union1d(a,b)
    return float(np.max(abs(np.searchsorted(np.sort(a),points,side="right")/len(a)-np.searchsorted(np.sort(b),points,side="right")/len(b))))


def width(m):return float(2*np.sqrt(np.log(4*m/ALPHA)/(2*R)))


def classify(rows,m,label,invalid=False):
    if invalid:return "QUALIFICATION FAILED"
    if len(rows)!=m:return "UNRESOLVED"
    if all(r["upper"]<=MARGIN for r in rows):return label
    if any(r["lower"]>MARGIN for r in rows):return "NOT QUALIFIED ON TESTED CASES"
    return "UNRESOLVED"


def evaluate():
    folder=BASE/"production"
    registration=json.loads((BASE/"registration.json").read_text())
    grouped={};invalid=False;missing=[]
    for job in registration["jobs"]:
        meta=folder/(job["id"]+".json")
        if not meta.exists():missing.append(job["id"]);continue
        record=json.loads(meta.read_text())
        if record["status"]=="QUALIFICATION FAILED":invalid=True
        if record["status"]!="COMPLETE":missing.append(job["id"]);continue
        values=[]
        for shard in record["shards"]:
            with np.load(folder/shard["file"]) as raw:values.append(raw["measurements"])
        combined=np.concatenate(values)
        if len(combined)!=job["reps"]:missing.append(job["id"]);continue
        grouped.setdefault((job["case"],job["method"]),[]).append(combined)
    groups={key:np.concatenate(value) for key,value in grouped.items() if sum(len(x) for x in value)==R}
    families={"snapshot":[],"short":[]};descriptive=[]
    def compare(a,b,case,method,t,mode,family):
        for k,name in enumerate(NAMES):
            d=distance(a[:,k],b[:,k]);row=dict(case=case,method=method,time=float(t),mode=int(mode),observable=name,distance=d)
            if family:
                half=width(192 if family=="snapshot" else 96)
                row.update(lower=max(0.,d-half),upper=min(1.,d+half));families[family].append(row)
            else:descriptive.append(row)
    for case in range(4):
        for step,t in enumerate((1/6000,1/12000,1/24000)):
            a,b=groups.get((case,2*step)),groups.get((case,2*step+1))
            if a is None or b is None:continue
            for index,mode in enumerate((32,250)):compare(a[:,index],b[:,index],case,step,t,mode,"snapshot")
    for case in (4,5):
        exact=groups.get((case,0))
        if exact is None:continue
        for method in (1,2,3):
            approx=groups.get((case,method))
            if approx is None:continue
            for t_index,t in enumerate((1/8000,1/4000,1/2000)):
                for k,mode in enumerate((1,8)):compare(approx[:,t_index,k],exact[:,t_index,k],case,method,t,mode,"short" if method==3 else None)
    snapshot_label=classify(families["snapshot"],192,"QUALIFIED ON TESTED SNAPSHOT LAWS",invalid)
    short_label=classify(families["short"],96,"QUALIFIED ON TESTED SHORT-TIME SPATIAL LAWS",invalid)
    if invalid:overall="QUALIFICATION FAILED"
    elif snapshot_label=="QUALIFIED ON TESTED SNAPSHOT LAWS" and short_label=="QUALIFIED ON TESTED SHORT-TIME SPATIAL LAWS" and not missing:overall="QUALIFIED WITHIN REGISTERED S2 SCOPE"
    elif "NOT QUALIFIED ON TESTED CASES" in (snapshot_label,short_label):overall="NOT QUALIFIED ON TESTED CASES"
    else:overall="UNRESOLVED"
    result=dict(verdict=overall,snapshot_verdict=snapshot_label,short_verdict=short_label,complete_logical_groups=len(groups),planned_logical_groups=32,
                missing_jobs=missing,alpha_per_family=ALPHA,margin=MARGIN,snapshot_halfwidth=width(192),short_halfwidth=width(96),
                snapshot=families["snapshot"],short=families["short"],coarse_descriptive=descriptive,
                snapshot_max_upper=max((r["upper"] for r in families["snapshot"]),default=None),short_max_upper=max((r["upper"] for r in families["short"]),default=None))
    (folder/"evaluation.json").write_text(json.dumps(result,indent=2)+"\n",encoding="utf-8")
    return result


if __name__=="__main__":
    result=evaluate();print({k:v for k,v in result.items() if k not in ("snapshot","short","coarse_descriptive")})

"""MG046 prespecified simultaneous marginal-CDF equivalence evaluation."""
from pathlib import Path
import json
import numpy as np
import mg10s1 as m

NAMES=("total_X","total_B","total_W","R1_count","R3_count","R4_count","node0_X","node0_Y")
M=96
ALPHA=.01
MARGIN=.05


def empirical_distance(a,b):
    largest=int(max(a.max(),b.max()))
    first=np.bincount(a,minlength=largest+1).cumsum()/len(a)
    second=np.bincount(b,minlength=largest+1).cumsum()/len(b)
    return float(np.max(abs(first-second)))


def half_error():return float(2*np.sqrt(np.log(4*M/ALPHA)/(2*m.R)))


def classify(complete,invalid,bands):
    if invalid:return "QUALIFICATION FAILED"
    if not complete or len(bands)!=M:return "UNRESOLVED"
    if all(row[1]<=MARGIN for row in bands):return "QUALIFIED ON TESTED SMALL SYSTEMS"
    if any(row[0]>MARGIN for row in bands):return "NOT QUALIFIED ON TESTED CASES"
    return "UNRESOLVED"


def evaluate():
    folder=m.BASE/"production";observed={};audits=[];invalid=False
    for case in range(4):
        for method in range(4):
            key=f"case{case}_method{method}";meta=folder/(key+".json");file=folder/(key+".npz")
            if not meta.exists():continue
            record=json.loads(meta.read_text())
            if record["status"]=="QUALIFICATION FAILED":invalid=True
            if record["status"]!="COMPLETE" or not file.exists():continue
            with np.load(file) as raw:data,initial,stats=raw["data"],raw["initial"],raw["stats"]
            if data.shape!=(m.R,3,20,initial.shape[1]):continue
            audit=m.ledger(data,initial)
            audit.update(id=key,max_probability=float(stats[:,2].max()),adapted_steps=float(stats[:,1].sum()))
            if not audit["pass_check"] or (method>0 and stats[:,2].max()>.25+1e-14):invalid=True
            observed[key]=m.observables(data);audits.append(audit)
    rows=[];descriptive=[];bands=[];width=half_error()
    for case in range(4):
        exact=observed.get(f"case{case}_method0")
        for method in (1,2,3):
            approx=observed.get(f"case{case}_method{method}")
            if exact is None or approx is None:continue
            for t in range(3):
                for index,name in enumerate(NAMES):
                    distance=empirical_distance(exact[:,t,index],approx[:,t,index])
                    row=dict(case=case,method=method,time=float(m.OBS[t]),observable=name,empirical_CDF_distance=distance)
                    if method==3:
                        lower,upper=max(0.,distance-width),min(1.,distance+width)
                        row.update(lower=lower,upper=upper)
                        bands.append((lower,upper));rows.append(row)
                    else:descriptive.append(row)
    complete=len(observed)==16 and len(rows)==96
    result=dict(verdict=classify(complete,invalid,bands),complete_groups=len(observed),planned_groups=16,replicates_per_group=m.R,
                alpha=ALPHA,margin=MARGIN,simultaneous_band_halfwidth=width,primary=rows,coarse_descriptive=descriptive,
                max_primary_empirical_distance=max((r["empirical_CDF_distance"] for r in rows),default=None),
                max_primary_upper=max((r["upper"] for r in rows),default=None),audits=audits)
    (folder/"evaluation.json").write_text(json.dumps(result,indent=2)+"\n",encoding="utf-8")
    return result


if __name__=="__main__":
    result=evaluate();print({k:v for k,v in result.items() if k not in ("primary","coarse_descriptive","audits")})

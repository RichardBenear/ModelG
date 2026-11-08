"""MG074 fixed-origin branch measurements and balanced paired-block analysis."""
import os
for key in ("OPENBLAS_NUM_THREADS","MKL_NUM_THREADS","OMP_NUM_THREADS","NUMEXPR_NUM_THREADS"):os.environ[key]="1"
from pathlib import Path
import json,hashlib
import numpy as np
from mg10d1_pool import atomic_json
from mg10e2 import window_average as late_average,balanced_analysis,field_metrics
from mg10p1 import assess
from mg10p1 import continuity
ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/"results/MG10E2"
analyze_blocks=balanced_analysis

def main(base=BASE):
    base=Path(base);r=json.loads((base/"registration.json").read_text());records={};parents=0
    for job in r["parent_jobs"]:
        path=base/"parents"/(job["id"]+".json")
        if path.exists() and json.loads(path.read_text()).get("status")=="COMPLETE":parents+=1
    for job in r["branch_jobs"]:
        path=base/"branches"/(job["id"]+".npz");meta=path.with_suffix(".json")
        if not path.exists():continue
        with np.load(path) as raw:data=raw["data"].copy();times=raw["times"].copy()
        complete=meta.exists() and json.loads(meta.read_text()).get("status")=="COMPLETE" and np.array_equal(times,np.arange(101)/4)
        fields=data[0,:,:3]/1000.;metrics=field_metrics(fields);matches=[[assess(f[1]) for f in fields]];projection=metrics["projection"]
        fields=data[0,:,:3]/1000.;track=None
        if job["start"]:
            first=matches[0][0]
            track=continuity(fields[:,1],times,x0=first["position"] if first["status"]=="REFERENCE_MATCH" else 0.)
        item=dict(job=job,complete=bool(complete),raw_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),times=times.tolist(),
          projection=projection.tolist(),movement_metrics={name:value.tolist() for name,value in metrics.items() if name not in ("projection","field_mean")},field_mean=metrics["field_mean"].tolist(),matches=matches[0],initial_reference_continuity=track,
          source_counts=data[0,:,7].sum(axis=-1).tolist(),external_inputs=data[0,:,20:22].sum(axis=-1).tolist(),
          early_projection=late_average(projection,times,0.,5.).tolist() if len(times)>=21 else None,initial_stocks=data[0,0,3:7].sum(axis=-1).tolist(),final_stocks=data[0,-1,3:7].sum(axis=-1).tolist())
        if complete:
            item.update(late_projection=late_average(projection,times).tolist(),late_physical_profile=late_average(fields,times).tolist())
        records[job["id"]]=item
    effects,summary=analyze_blocks(records);complete=sum(v["complete"] for v in records.values())
    status="COMPLETE RESTORATION ESTIMATION" if parents==32 and complete==128 else "INCOMPLETE RESTORATION ESTIMATION"
    result=dict(status=status,parents_complete=parents,branches_complete=complete,branches_available=len(records),records=records,block_effects=effects,primary_by_step=summary,
      inference="Nominal97.5% t interval per step, approximate95% two-step family, n8 maximum; no nonsignificance/equivalence or physical mortality claim")
    atomic_json(base/"evaluation.json",result)
    full=all(summary[str(k)]["n"]==8 for k in range(2))
    intervals=[summary[str(k)]["ci_nominal"] for k in range(2)]
    direction=full and all(ci is not None and ci[0]>0 for ci in intervals) or full and all(ci is not None and ci[1]<0 for ci in intervals)
    paragraph=("E2 completed its planned comparison of interrupted and restored A supply with seeded and unseeded controls. " if status.startswith("COMPLETE ") else "E2 did not complete its planned restoration comparison; missing branches or parents limit interpretation. ")
    paragraph+=("Both numerical steps show a same-direction seed-history-associated response to restoring local supply under the registered approximate intervals. " if direction else "The registered two-step criterion for direction-consistent pilot evidence is unresolved. ")
    paragraph+="This tests restoration response inside the implemented model; return of a reference signal does not prove full recovery, survival of the same entity, or quantum memory; it does not identify physical particles, equate reference-match loss with death, or derive A/B from a vacuum. Separate raw audit must pass before accepting the result."
    text=f"# MG10-E2 results\n\n{paragraph}\n\nStatus: {status}; PI closure pending.\n\nParents{parents}/32; branches{complete}/128. Primary fixed-origin X late-window restoration-versus-continuous-off seed interaction:\n\n"
    for k in range(2):text+=f"- Step{k}: {summary[str(k)]}\n"
    text+="\nFull profiles, source histories, individual contrasts, P1 diagnostics and incomplete records: evaluation.json. Small-n t approximation, finite horizon, fixed coordinates and background-response limitations apply. No old result rescored.\n"
    (base/"MG10E2_results.md").write_text(text,encoding="utf-8")
if __name__=="__main__":main()

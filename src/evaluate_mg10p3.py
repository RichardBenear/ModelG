"""MG066 descriptive sampled-reference persistence and conditional averaging."""
import os
for key in ("OPENBLAS_NUM_THREADS","MKL_NUM_THREADS","OMP_NUM_THREADS","NUMEXPR_NUM_THREADS"):
    os.environ[key]="1"
from pathlib import Path
import json,hashlib
import numpy as np
from scipy.stats import t as student_t
from mg10d1_pool import atomic_json
from mg10d2 import ledger as integer_ledger
from mg10p1 import assess
from mg10p2 import binomial_interval,reference
from mg10p3 import summarize,OBS

ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/"results/MG10P3"


def jsonable(value):
    if isinstance(value,np.ndarray):return value.tolist()
    if isinstance(value,np.generic):return value.item()
    if isinstance(value,dict):return {k:jsonable(v) for k,v in value.items()}
    if isinstance(value,(list,tuple)):return [jsonable(v) for v in value]
    return value


def describe(values):
    a=np.asarray(values,dtype=float);n=len(a)
    if n==0:return dict(n=0,mean=None,min=None,max=None,ci95_approx=None)
    if not np.all(np.isfinite(a)):raise ValueError("Nonfinite summary; do not drop invalid outputs")
    m=a.mean(axis=0);ci=None
    if n>1:
        half=student_t.ppf(.975,n-1)*a.std(axis=0,ddof=1)/np.sqrt(n);ci=[(m-half).tolist(),(m+half).tolist()]
    return dict(n=n,mean=m.tolist(),min=a.min(axis=0).tolist(),max=a.max(axis=0).tolist(),ci95_approx=ci)


def average_bundle(items,key):
    selected=[p[key] for p in items if p[key] is not None];n=len(selected)
    bundle=dict(n=n,profile_mean=None,stats={},gain=[describe([]) for _ in range(3)])
    if not n:return bundle
    bundle["profile_mean"]=np.mean([v["profile"] for v in selected],axis=0).tolist()
    for metric in ("projection","residual_rms","contrast_rms","mean_individual_residual_rms"):
        bundle["stats"][metric]=describe([v[metric] for v in selected])
    bundle["gain"]=[describe([v["gain"][j] for v in selected if v["gain"][j] is not None]) for j in range(3)]
    return bundle


def evaluate():
    registration=json.loads((BASE/"registration.json").read_text(encoding="utf-8"));ledger=json.loads((BASE/"production/ledger.json").read_text(encoding="utf-8"))
    records={v["id"]:v for v in ledger["records"]};paths={};missing=[];manifest=[]
    for job in registration["jobs"]:
        file=BASE/"production"/(job["id"]+".npz");record=records.get(job["id"],{})
        if not file.exists():missing.append(dict(job=job,reason="No retained checkpoint",last_time=None));continue
        sha=hashlib.sha256(file.read_bytes()).hexdigest()
        with np.load(file) as raw:
            data=raw["data"];times=raw["times"];count=data.shape[1]
            if data.shape!=(1,count,22,500) or not 1<=count<=101 or not np.array_equal(times,OBS[:count]) or int(raw["seed"])!=job["seed"]:raise ValueError("Wrong checkpoint/seed schema")
            if not integer_ledger(data,raw["initial"])["pass_check"]:raise ValueError("Integer accounting")
            successful=record.get("status")=="COMPLETE" and record.get("returncode")==0
            manifest.append(dict(path="production/"+file.name,sha256=sha,checkpoints=count,last_time=float(times[-1]),successful=successful,job=job))
            if not successful or count!=101:
                if successful:raise ValueError("Success metadata disagrees with checkpoint grid")
                missing.append(dict(job=job,reason=record.get("termination_reason",record.get("status","Incomplete")),last_time=float(times[-1])));continue
            if sha!=record["raw_sha256"]:raise ValueError("Raw hash changed")
            result=summarize(data[0],times,job["start"])
            result["observations"]=result["track"]["observations"] if job["start"] else [assess(frame[1]/1000) for frame in data[0,:,:3]]
            result["event20"]=bool(result["track"]["persistent_reference_signal"] if job["start"] else result["any_observed_episode20"])
            paths[job["id"]]=jsonable(dict(job=job,stats=raw["stats"].tolist(),**result))
    groups={}
    for start in range(2):
        for step in range(2):
            selected=[v for v in paths.values() if (v["job"]["start"],v["job"]["step"])==(start,step)];n=len(selected);k=sum(v["event20"] for v in selected)
            reasons={}
            if start:
                for p in selected:
                    reason=p["track"]["censor_reason"] or p["track"]["status"];reasons[reason]=reasons.get(reason,0)+1
            groups[f"s{start}_k{step}"]=dict(start=start,step=step,n=n,expected_n=16,events20=k,fraction=k/n if n else None,
                ci95_descriptive=binomial_interval(k,n) if n else None,
                event_definition="initial reference span>=20" if start else "any distinct observed reference episode>=20",
                matched_through=[p["track"]["matched_through"] for p in selected] if start else None,censor_counts=reasons,
                fixed=average_bundle(selected,"fixed"),aligned=average_bundle(selected,"aligned"),
                alignment_scope="Conditional on initial20-unit X continuity for S; never provided for H")
    complete=len(paths)==64 and len(records)==64 and not ledger["stop_reason"] and all(v.get("status")=="COMPLETE" and v.get("returncode")==0 for v in records.values())
    verdict="COMPLETE TEMPORAL ESTIMATION" if complete else "INCOMPLETE"
    result=dict(verdict=verdict,complete=complete,planned_paths=64,completed_paths=len(paths),paths=paths,groups=groups,missing=missing,
                charged_worker_seconds=ledger["charged_worker_seconds"],stop_reason=ledger["stop_reason"],
                scope="Sampled field continuity, not physical particle lineage/death; aligned averages conditional, no independence-of-frames assumption or rate success threshold")
    atomic_json(BASE/"production/evaluation.json",result);atomic_json(BASE/"production/raw_manifest.json",dict(raw=manifest,retention="All full/partial input archives retained locally; upload verification separate"))
    lines=["# MG10-P3 sampled temporal persistence and coupled-field averaging", "", "**"+verdict+"**. PI closure pending.","",
        f"Completed{len(paths)}/64 registered paths; production{ledger['charged_worker_seconds']/3600:.6f}core-hours under24. Stop reason: {ledger['stop_reason']}.","",
        "Full maintained-supply P17 baseline, H/S, two numerical steps,16 independent paths/cell,T25,101 quarter-unit checkpoints. No favorable persistence-rate threshold; completion is not a physical survival verdict.","",
        "## Descriptive20-unit observed events", "", "|Start|Step|Events/complete|Original n|95% exact interval|Aligned eligible|", "|---|---|---:|---:|---|---:|"]
    for g in groups.values():lines.append(f"|{'S: initial signal' if g['start'] else 'H: any episode'}|{'fine' if g['step'] else 'coarse'}|{g['events20']}/{g['n']}|16|{g['ci95_descriptive']}|{g['aligned']['n']}|")
    lines += ["", "S tracks stop at their first unresolved frame; later episodes cannot restore initial identity. H candidates are distinct maximal episodes, never stitched. A gap/mismatch is unresolved, not physical death. Incomplete trajectories are missing, not counted as failed persistence, even when their prefix includes20 units. Descriptive intervals concern completed paths; incomplete runs may be selected by runtime and do not support the full registered population.","",
        "## Averaging and eligibility", "", "All completed paths have fixed-frame averages over0..20 (81 frames). Only S paths with initial20-unit continuity also have X-phase-aligned averages of the same window. All three fields use X's phase; G/Y are never independently fitted. Signed projections, residuals, mean individual residuals, observed gains and full profiles are in production/evaluation.json. Undefined gain is missing, never zero/infinity. Conditional alignment and coupled noise can bias profiles; no population-wide coupled-coherence decision follows.","",
        "Temporal samples are correlated. No sqrt(81) gain or frame-independent confidence interval is assumed. Approximate non-simultaneous t intervals summarize independent trajectories, with eligible counts explicit. Two finite step sizes are compared descriptively, not declared equivalent. R16 limits rare-event precision; this does not repeat P2's5% match cap. T25 does not answer indefinite survival or depleted/intervention regimes.","",
        "Independent audit, full input retention/upload and PI closure are separate records. Earlier P1/P2/D2 results remain unchanged.",""]
    if paths:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig,ax=plt.subplots(figsize=(8,4),layout="constrained")
        names=[];positions=[]
        for i,g in enumerate(groups.values()):
            if not g["n"]:continue
            lo,hi=g["ci95_descriptive"];v=g["fraction"]
            ax.errorbar(i,v,yerr=np.array([[v-lo],[hi-v]]),fmt="o",capsize=4)
            positions.append(i);names.append(f"{'S initial' if g['start'] else 'H episode'}\n{'fine' if g['step'] else 'coarse'} ({g['n']}/16)")
        ax.set_xticks(positions,names);ax.set_ylim(-.03,1.03);ax.set_ylabel("Observed20-unit event fraction");ax.set_title("Descriptive exact95% intervals; not physical survival");ax.grid(axis="y",alpha=.2)
        fig.savefig(BASE/"MG10P3_fractions.png",dpi=150);plt.close(fig)
        _,template=reference();x=(np.arange(500)-250)*.2
        fig,axes=plt.subplots(2,3,figsize=(12,7),sharex=True,layout="constrained")
        for row,key in enumerate(("fixed","aligned")):
            for j,species in enumerate("GXY"):
                ax=axes[row,j];ax.plot(x,template[j],"k--",label="P17 reference")
                for step in range(2):
                    g=groups[f"s1_k{step}"];bundle=g[key]
                    if bundle["n"]:ax.plot(x,np.array(bundle["profile_mean"])[j],label=f"{'fine' if step else 'coarse'} n={bundle['n']}/16")
                ax.set_xlim(-10,10);ax.set_title(f"{species}: {'all completed S, fixed frame' if not row else 'eligible S, X-aligned'}");ax.set_xlabel("position");ax.set_ylabel("Centered concentration");ax.grid(alpha=.2);ax.legend(fontsize=7)
        fig.savefig(BASE/"MG10P3_profiles.png",dpi=150);plt.close(fig)
        lines += ["![Observed event fractions](MG10P3_fractions.png)","", "![Coupled temporal averages](MG10P3_profiles.png)",""]
    (BASE/"MG10P3_results.md").write_text("\n".join(lines),encoding="utf-8")
    print(f"{verdict}; {len(paths)}/64 paths",flush=True)


if __name__=="__main__":evaluate()

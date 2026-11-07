"""MG062 frozen P2 frequency-family and descriptive observability evaluation."""
import os
for key in ("OPENBLAS_NUM_THREADS","MKL_NUM_THREADS","OMP_NUM_THREADS","NUMEXPR_NUM_THREADS"):
    os.environ[key]="1"
from pathlib import Path
import json,hashlib
import numpy as np
from scipy.stats import t as student_t
from mg10d1_pool import atomic_json
from mg10d2 import ledger as integer_ledger
import mg10p2 as model

ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/"results/MG10P2"


def describe(values):
    a=np.asarray(values,dtype=float);n=len(a)
    if n==0:return dict(n=0,mean=None,min=None,max=None,ci95_approx=None)
    mean=a.mean(axis=0);ci=None
    if n>1:
        half=student_t.ppf(.975,n-1)*a.std(axis=0,ddof=1)/np.sqrt(n)
        ci=[(mean-half).tolist(),(mean+half).tolist()]
    return dict(n=n,mean=mean.tolist(),min=a.min(axis=0).tolist(),max=a.max(axis=0).tolist(),ci95_approx=ci)


def evaluate():
    registration=json.loads((BASE/"registration.json").read_text(encoding="utf-8"))
    ledger=json.loads((BASE/"production/ledger.json").read_text(encoding="utf-8"))
    records={v["id"]:v for v in ledger["records"]};paths={};manifest=[];missing=[]
    for job in registration["jobs"]:
        file=BASE/"production"/(job["id"]+".npz")
        if not file.exists():missing.append(dict(job=job,paths_missing=8));continue
        sha=hashlib.sha256(file.read_bytes()).hexdigest();record=records.get(job["id"],{})
        if record.get("status")=="COMPLETE" and record.get("returncode")==0 and sha!=record["raw_sha256"]:raise ValueError("Raw hash changed")
        with np.load(file) as raw:
            data=raw["data"];n=len(data)
            if not 1<=n<=8 or data.shape!=(n,4,22,500):raise ValueError("Invalid atomic prefix")
            if not np.array_equal(raw["times"],model.OBS) or list(raw["seeds"])!=job["seeds"][:n] or list(raw["replicates"])!=job["replicates"][:n]:raise ValueError("Wrong streams/checkpoints")
            if not integer_ledger(data,raw["initial"])["pass_check"]:raise ValueError("Integer accounting")
            metrics,matches=model.measure(data)
            for index,rep in enumerate(job["replicates"][:n]):
                key=f"c{job['condition']}_s{job['start']}_k{job['step']}_r{rep}"
                if key in paths:raise ValueError("Duplicate path")
                paths[key]=dict(condition=job["condition"],start=job["start"],step=job["step"],replicate=rep,seed=job["seeds"][index],
                    metrics={k:v[index].tolist() for k,v in metrics.items()},matches=matches[index],stats=raw["stats"][index].tolist())
            if n<8:missing.append(dict(job=job,paths_missing=8-n))
        manifest.append(dict(path="production/"+file.name,sha256=sha,paths=n,job=job))
    complete=len(paths)==1024 and len(records)==128 and not ledger["stop_reason"] and all(v.get("status")=="COMPLETE" and v.get("returncode")==0 for v in records.values())
    _,template=model.reference();reference_rms=np.sqrt(np.mean(template**2,axis=1));groups={};primary=[]
    for c in range(2):
        for s in range(2):
            for step in range(2):
                key=f"c{c}_s{s}_k{step}";selected=[v for v in paths.values() if (v["condition"],v["start"],v["step"])==(c,s,step)]
                counts=[sum(v["matches"][t]["status"]=="REFERENCE_MATCH" for v in selected) for t in range(4)]
                n=len(selected);summaries={};proxy=None
                if n:
                    for metric in selected[0]["metrics"]:summaries[metric]=describe([v["metrics"][metric] for v in selected])
                    intervals=[model.binomial_interval(k,n) for k in counts]
                    rates=[k/n for k in counts]
                    if not s:
                        smooth=np.array(summaries["smooth_rms"]["mean"])
                        proxy=[[float(reference_rms[j]/smooth[t,j]) if smooth[t,j]>0 else None for j in range(3)] for t in range(4)]
                else:intervals=None;rates=None
                group=dict(condition=c,start=s,step=step,n=n,expected_n=128,match_counts=counts,match_rates=rates,match_ci95_descriptive=intervals,
                    summaries=summaries,reference_over_mean_H_smooth_rms=proxy)
                if not s:
                    upper=model.upper_bound(counts[-1],n) if n else None
                    item=dict(group=key,matches=counts[-1],n=n,upper=upper,margin=.05,scored=complete,
                              cap_met=bool(upper<=.05) if complete else None)
                    group["primary"]=item;primary.append(item)
                groups[key]=group
    verdict="INCOMPLETE" if not complete else "SHORT-HORIZON UNSEEDED MATCH CAP SUPPORTED ON TESTED CONDITIONS" if all(v["cap_met"] for v in primary) else "MATCH CAP NOT QUALIFIED"
    result=dict(verdict=verdict,complete=complete,planned_paths=1024,completed_paths=len(paths),missing=missing,paths=paths,groups=groups,primary=primary,
        family_alpha=.05,group_alpha=.0125,match_cap=.05,endpoint=.25,observations=model.OBS.tolist(),reference_smooth_rms=reference_rms.tolist(),
        charged_worker_seconds=ledger["charged_worker_seconds"],stop_reason=ledger["stop_reason"],
        interpretation="H is unseeded reference-match frequency, not physical false-particle rate. S retention/GXY observability descriptive; no survival or joint-field qualification.")
    atomic_json(BASE/"production/evaluation.json",result)
    atomic_json(BASE/"production/raw_manifest.json",dict(raw=manifest,retention="Complete raw spatial/event/input checkpoints retained locally; upload verification separate"))
    lines=["# MG10-P2 model-noise observability and short-horizon match calibration", "", "**"+verdict+"**. PI closure pending.","",
        f"{len(paths)}/1024 paths; production charge{ledger['charged_worker_seconds']/3600:.6f}core-hours under4. Stop reason: {ledger['stop_reason']}.","",
        "## Primary unseeded endpoint family", "", "Four one-sided exact Clopper-Pearson upper bounds, groupalpha.0125, simultaneous coverage at least95%; fixed match-frequency cap5%. No primary family is scored for an incomplete batch.","",
        "|Condition|Step|Matches/n at.25|Upper bound|Cap met|", "|---|---|---:|---:|---|"]
    for v in primary:lines.append(f"|{v['group'].split('_s')[0]}|{v['group'].split('_k')[1]}|{v['matches']}/{v['n']}|{v['upper'] if v['upper'] is not None else 'missing'}|{v['cap_met']}|")
    lines += ["", "## Seeded match retention (descriptive)","", "A planted shape need not remain stationary, particularly in the prescribed depleted state. These are not sensitivities to a known surviving particle.","",
        "|Cell|Matches/n at.25|Descriptive95% interval|", "|---|---:|---|"]
    for key,g in groups.items():
        if g["start"]:lines.append(f"|{key}|{g['match_counts'][-1]}/{g['n']}|{g['match_ci95_descriptive'][-1] if g['n'] else 'missing'}|")
    lines += ["", "## Reference contrast / mean unseeded smoothed contrast at.25", "", "This is an observability proxy, not a calibrated physical SNR. TemplateRMS G/X/Y: "+str(reference_rms.tolist())+".","",
        "|Cell|G|X|Y|", "|---|---:|---:|---:|"]
    for key,g in groups.items():
        if not g["start"] and g["n"]:
            values=g["reference_over_mean_H_smooth_rms"][-1];lines.append("|"+key+"|"+"|".join(f"{v:.6g}" if v is not None else "missing" for v in values)+"|")
    lines += ["", "## Scope", "", "FULL is maintained A/B; DEPLETED starts at prescribed82% stocks with original denominators and a rounded effective-rate mass-action background. Seeded initial conditions add the committed G/X/Y reference; they are not assumed stationary depleted solitons. Time.25 is not stationary noise and is below the20-unit persistence criterion. Conditions do not isolate a causal depletion effect.","",
        "P1 is unchanged. G/X/Y signed fixed-phase projections, contrasts and residuals with approximate non-simultaneous t intervals are descriptive. Weak G contrast is not particle absence; no X match is not death. No physical particle lineage, coupled-field qualification, stochastic equivalence or prior-study rescoring is inferred. Full per-path inputs/results and eight-cell statistics are in production/evaluation.json; independent audit and upload verification are separate files.", ""]
    if len(paths):
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig,axes=plt.subplots(1,2,figsize=(11,4),sharey=True,layout="constrained")
        for c,ax in enumerate(axes):
            for s in range(2):
                for k in range(2):
                    g=groups[f"c{c}_s{s}_k{k}"]
                    if g["n"]:ax.plot(model.OBS,g["match_rates"],linestyle="-" if k else "--",marker="o",label=f"{'S' if s else 'H'}, {'fine' if k else 'coarse'}")
            ax.set_title("FULL maintained" if not c else "Prescribed82% finite");ax.set_xlabel("model time");ax.set_ylim(-.03,1.03);ax.grid(alpha=.2);ax.legend(fontsize=8)
        axes[0].set_ylabel("P1 reference-match fraction")
        fig.savefig(BASE/"MG10P2_match_fractions.png",dpi=150);plt.close(fig)
        lines += ["![Descriptive match fractions](MG10P2_match_fractions.png)",""]
    (BASE/"MG10P2_results.md").write_text("\n".join(lines),encoding="utf-8")
    print(f"{verdict}; {len(paths)}/1024 paths",flush=True)


if __name__=="__main__":evaluate()

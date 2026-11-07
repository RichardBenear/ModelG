"""MG055 descriptive D2 evaluation; no particle or equivalence verdict."""
import os
for key in ("OPENBLAS_NUM_THREADS","MKL_NUM_THREADS","OMP_NUM_THREADS","NUMEXPR_NUM_THREADS"):
    os.environ[key]="1"
from pathlib import Path
import json
import numpy as np
from scipy.stats import t as student_t
from bounded_pool import atomic_json
import mg10d2 as model

ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/"results/MG10D2"


def describe(values):
    values=np.asarray(values,dtype=float);n=len(values)
    if not n:return dict(n=0,mean=None,min=None,max=None,ci95_approx=None)
    mean=values.mean(axis=0)
    interval=None
    if n>=2:
        half=student_t.ppf(.975,n-1)*values.std(axis=0,ddof=1)/np.sqrt(n)
        interval=[(mean-half).tolist(),(mean+half).tolist()]
    return dict(n=n,mean=mean.tolist(),min=values.min(axis=0).tolist(),max=values.max(axis=0).tolist(),ci95_approx=interval)


def evaluate():
    registration=json.loads((BASE/"registration.json").read_text(encoding="utf-8"))
    ledger=json.loads((BASE/"production/ledger.json").read_text(encoding="utf-8"))
    records={r["id"]:r for r in ledger["records"]};paths={};manifest=[];excluded=[]
    for job in registration["jobs"]:
        item=records.get(job["id"])
        if not item or item.get("status")!="COMPLETE" or item.get("returncode")!=0:
            excluded.append(dict(job=job,reason="Not completed successfully; kept in denominator"));continue
        with np.load(BASE/"production"/item["raw"]) as raw:
            if raw["data"].shape!=(1,61,22,500):raise ValueError("Incomplete checkpoint dimensions")
            if not model.ledger(raw["data"],raw["initial"])["pass_check"]:raise ValueError("Integer ledger")
            metrics=model.measurements(raw["data"],model.fixture(job["arm"],job["start"]))
            paths[job["id"]]=dict(job=job,metrics={k:v[0].tolist() for k,v in metrics.items()},stats=raw["stats"].tolist())
        manifest.append(dict(path="production/"+item["raw"],sha256=item["raw_sha256"],job=job))
    groups={}
    for arm in range(2):
        for start in range(2):
            for step in range(2):
                selected=[p for p in paths.values() if (p["job"]["arm"],p["job"]["start"],p["job"]["step"])==(arm,start,step)]
                key=f"a{arm}_s{start}_k{step}";summary={}
                if selected:
                    for metric in selected[0]["metrics"]:
                        values=np.array([p["metrics"][metric] for p in selected])
                        summary[metric]=dict(time_series=describe(values),endpoint=describe(values[:,-1]),
                                             time_average=describe(np.trapezoid(values,dx=1,axis=1)/60))
                groups[key]=dict(arm=model.ARMS[arm],start="S" if start else "H",step=step,n=len(selected),expected_n=4,summary=summary)
    contrasts={}
    pairs=[]
    for arm in range(2):
        for start in range(2):pairs.append((f"step_a{arm}_s{start}",f"a{arm}_s{start}_k1",f"a{arm}_s{start}_k0"))
    for arm in range(1,2):
        for start in range(2):
            for step in range(2):pairs.append((f"supply_a{arm}_s{start}_k{step}",f"a{arm}_s{start}_k{step}",f"a0_s{start}_k{step}"))
    for key,left,right in pairs:
        l,r=groups[left],groups[right];differences={}
        if l["n"] and r["n"]:
            for metric in l["summary"]:
                differences[metric]={label:(np.array(l["summary"][metric][label]["mean"])-np.array(r["summary"][metric][label]["mean"])).tolist() for label in ("endpoint","time_average")}
        contrasts[key]=dict(left=left,right=right,n_left=l["n"],n_right=r["n"],mean_difference=differences,interpretation="Descriptive independent-sample difference, no equivalence decision")
    complete=len(paths)==32 and not ledger["stop_reason"]
    verdict="COMPLETE DESCRIPTIVE PILOT" if complete else "BUDGET-SHORT / INCOMPLETE" if "BUDGET" in str(ledger["stop_reason"]) else "INCOMPLETE"
    result=dict(verdict=verdict,planned_paths=32,completed_paths=len(paths),excluded=excluded,paths=paths,groups=groups,contrasts=contrasts,
                charged_worker_seconds=ledger["charged_worker_seconds"],stop_reason=ledger["stop_reason"],
                numerical_equivalence="NOT TESTED: step contrasts descriptive, four independent paths/group",
                particle_survival="NOT MEASURED",ci_scope="Approximate non-simultaneous t intervals; n=4 and nonnormality limit coverage")
    atomic_json(BASE/"production/evaluation.json",result)
    atomic_json(BASE/"production/raw_manifest.json",dict(raw=manifest,retention="All local raw and progress checkpoints retained; upload verification is separate"))
    # Report every completed path's traces and mean/range through JSON; figures
    # are an additional readable view, never the sole decision inputs.
    figures=[]
    if paths:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        colors=plt.get_cmap("tab10").colors
        for start in range(2):
            fig,axes=plt.subplots(4,2,figsize=(12,12),sharex=True,layout="constrained")
            panels=[("field_mean",i,f"{s} spatial mean") for i,s in enumerate("GXY")]
            panels += [("field_rms",i,f"{s} RMS spatial contrast") for i,s in enumerate("GXY")]
            panels += [("availability_mean",i,f"{s} mean availability") for i,s in enumerate("AB")]
            for ax,(metric,component,title) in zip(axes.flat,panels):
                for arm in range(2):
                    for step in range(2):
                        chosen=[p for p in paths.values() if (p["job"]["arm"],p["job"]["start"],p["job"]["step"])==(arm,start,step)]
                        if not chosen:continue
                        values=np.array([p["metrics"][metric] for p in chosen])[:,:,component]
                        style="-" if step else "--"
                        for v in values:ax.plot(range(61),v,color=colors[arm],alpha=.15,linewidth=.6,linestyle=style)
                        ax.plot(range(61),values.mean(axis=0),color=colors[arm],linestyle=style,label=model.ARMS[arm] if step else None)
                        if step:ax.fill_between(range(61),values.min(axis=0),values.max(axis=0),color=colors[arm],alpha=.07)
                ax.set_title(title);ax.grid(alpha=.2);ax.set_xlabel("model time")
            axes[0,0].legend(fontsize=8)
            fig.suptitle(f"MG10-D2 {'seeded' if start else 'homogeneous'}: individual paths + means; solid fine / dashed coarse")
            name=f"MG10D2_{'S' if start else 'H'}_fields.png";fig.savefig(BASE/name,dpi=150);plt.close(fig);figures.append(name)
    report=["# MG10-D2 finite-supply estimation pilot — results", "",f"**{verdict}**. Completed {len(paths)}/32 registered trajectories. PI closure pending.","",
            f"Production charge: {ledger['charged_worker_seconds']/3600:.6f} local core-hours; ceiling24. Stop reason: {ledger['stop_reason']}.","",
            "Two constructed-supply controls, N500/h.2/Omega1000, T60, H/S starts, four independent paths per arm/start/step. Numerical steps1/24000 and1/48000. Exact accounting includes maintained A/B external inputs.","",
            "## Fine-step endpoint estimates", "", "Means across completed independent paths; A/B are fractional availability. Incomplete denominators are explicit. All per-path traces, mean/ranges, approximate t intervals, endpoint and time-average summaries for both steps are in production/evaluation.json.","",
            "|Arm|Start|n/4|Mean X|X contrast RMS|Mean A|Mean B|", "|---|---|---:|---:|---:|---:|---:|"]
    for group in groups.values():
        if group["step"]!=1:continue
        if not group["n"]:report.append(f"|{group['arm']}|{group['start']}|0|—|—|—|—|");continue
        get=lambda metric,index:group["summary"][metric]["endpoint"]["mean"][index]
        report.append(f"|{group['arm']}|{group['start']}|{group['n']}|{get('field_mean',1):.6g}|{get('field_rms',1):.6g}|{get('availability_mean',0):.8g}|{get('availability_mean',1):.8g}|")
    report += ["", "## Limits and checks", "", "This is descriptive estimation. Four paths do not establish 1% recovery, stochastic equivalence, precise survival rates or discovery. A non-significant difference is not equivalence. RMS contrast does not establish localized particles or their identity. T60 does not answer T300/T500 behavior. Reservoirs and external inputs remain constructed resources, not a physical vacuum or quantum entanglement.","",
               "Student-t intervals are approximate and non-simultaneous, with limited reliability for n4/nonnormal populations. Step contrasts use independent samples and do not test equivalence. Numerical sensitivity may remain unresolved. Original MG10/S1/S2 outcomes remain unchanged.","",
               "The same-author independent raw checker writes production/independent_check.json after this report. Inspect that file before treating the accounting/recomputation as verified. Raw upload and PI closure are separate steps.", ""]
    report += [f"![Field histories]({name})\n" for name in figures]
    (BASE/"MG10D2_results.md").write_text("\n".join(report),encoding="utf-8")
    print(f"{verdict}; {len(paths)}/32 paths; {ledger['charged_worker_seconds']/3600:.5f} core-hours",flush=True)


if __name__=="__main__":evaluate()

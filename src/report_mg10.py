"""MG043 post-run descriptive diagnostics/figure. No new trajectories or rescoring."""
from pathlib import Path
import json
import numpy as np
from measure import smooth_periodic

ROOT=Path(__file__).resolve().parents[1]


def run():
    folder=ROOT/"results/MG10/production"
    evaluation=json.loads((folder/"evaluation.json").read_text(encoding="utf-8"))
    diagnostics={"scope":"POST_HOC_EXISTING_RAW_DATA_ONLY; no new run or rescoring", "rows":{}}
    curves={}
    for method in ("E2","O"):
        key=f"FINITE_H_{method}_C1"
        with np.load(folder/(key+".npz")) as raw:
            z,t=raw["state"],raw["times"]
        ranges=np.ptp(z[:,:3],axis=2)
        nonuniform=np.where(ranges.max(axis=1)>1e-6)[0]
        counts=np.array(evaluation["detector_counts"][key])
        exact_uniform=(ranges.max(axis=1)==0)
        artifacts=np.where(exact_uniform & (counts>0))[0]
        samples=[]
        for i in artifacts[:3]:
            smoothed=smooth_periodic(z[i,1]-140/9,.2)
            samples.append(dict(time=int(t[i]),raw_X_spatial_range=float(ranges[i,1]),
                                smoothed_X_spatial_range=float(np.ptp(smoothed)),detector_count=int(counts[i])))
        diagnostics["rows"][key]=dict(first_range_above_1e_minus6=int(t[nonuniform[0]]) if len(nonuniform) else None,
                                      max_spatial_range_by_species=ranges.max(axis=0).tolist(),
                                      scoring_count_range=[int(counts[50:].min()),int(counts[50:].max())],
                                      exactly_uniform_snapshots_with_nonzero_counts=len(artifacts),artifact_examples=samples)
        curves[method]=(t,ranges[:,1])
    (ROOT/"results/MG10/posthoc_diagnostics.json").write_text(json.dumps(diagnostics,indent=2)+"\n",encoding="utf-8")
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(1,2,figsize=(11,4.5),layout="constrained")
    for start,color,label in (("H","#2874a6","Homogeneous"),("S","#a04000","Soliton")):
        rows=[r for r in evaluation["rows"] if r["start"]==start]
        x=[r["capacity"] for r in rows]
        for method,style in (("E2","o-"),("O","x--")):
            axes[0].loglog(x,[r["primary_"+method] for r in rows],style,color=color,label=f"{label}, {method}",markersize=5)
    axes[0].axhline(.01,color="black",linestyle=":",label="1% recovery margin")
    axes[0].set(xlabel="Reservoir capacity C",ylabel="Maximum normalized active-field RMS error",
                title="High-capacity recovery; batch UNRESOLVED")
    axes[0].legend(fontsize=8);axes[0].grid(alpha=.2)
    for method,(t,span) in curves.items():
        axes[1].semilogy(t,np.maximum(span,1e-15),label=method)
    axes[1].set(xlabel="Physical time",ylabel="Spatial range of X (plot floor 1e-15)",
                title="C=1 homogeneous start: numerical disagreement")
    axes[1].legend();axes[1].grid(alpha=.2)
    fig.savefig(ROOT/"results/MG10/MG10_summary.png",dpi=170)
    plt.close(fig)
    print(json.dumps(diagnostics,indent=2))


if __name__=="__main__": run()

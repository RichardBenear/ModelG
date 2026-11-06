"""MG051 descriptive figures/moments from frozen S2 results, no new runs."""
from pathlib import Path
import json
import numpy as np

ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/"results/MG10S2"


def run():
    folder=BASE/"production";evaluation=json.loads((folder/"evaluation.json").read_text())
    labels=("H, Omega100","S, Omega100","H, Omega1000","S, Omega1000","SSA, Omega100","SSA, Omega1000")
    rows=[]
    for case in range(6):
        source=evaluation["snapshot"] if case<4 else evaluation["short"]
        case_rows=[r for r in source if r["case"]==case]
        worst=max(case_rows,key=lambda r:r["upper"])
        rows.append(dict(case=case,label=labels[case],distance=worst["distance"],lower=worst["lower"],upper=worst["upper"]))
    moments=[]
    for case in range(4):
        for method in range(6):
            record=json.loads((folder/f"snapshot_c{case}_m{method}_b0.json").read_text())
            values=[];q=q0=None
            for part in record["shards"]:
                with np.load(folder/part["file"]) as raw:
                    values.append(raw["measurements"])
                    if q is None:q,q0=raw["q_snapshot"],raw["q_instant"]
            empirical=np.concatenate(values).mean(axis=0)
            for k,mode in enumerate((32,250)):
                diagonal=np.diag(q[k]).real
                errors=abs(empirical[k,:3]-diagonal)/diagonal
                pairs=((0,1),(0,2),(1,2))
                cross=[abs(empirical[k,3+j]-q[k,a,b].real)/np.sqrt(diagonal[a]*diagonal[b]) for j,(a,b) in enumerate(pairs)]
                diag0=np.diag(q0[k]).real
                defect=float(np.max(abs(q[k]-q0[k])/np.sqrt(diag0[:,None]*diag0[None,:])))
                moments.append(dict(case=case,method=method,mode=mode,empirical_diagonal=empirical[k,:3].tolist(),
                                    exact_snapshot_diagonal=diagonal.tolist(),max_relative_diagonal_difference=float(errors.max()),
                                    max_normalized_cross_difference=float(max(cross)),snapshot_vs_instantaneous_normalized_difference=defect))
    result=dict(scope="POSTHOC_DESCRIPTION_ONLY; no new criterion or trajectory",case_bounds=rows,moments=moments,
                max_descriptive_diagonal_difference=max(r["max_relative_diagonal_difference"] for r in moments),
                max_descriptive_cross_difference=max(r["max_normalized_cross_difference"] for r in moments))
    (BASE/"descriptive_summary.json").write_text(json.dumps(result,indent=2)+"\n",encoding="utf-8")
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(1,2,figsize=(12,4.6),layout="constrained")
    x=np.arange(6);d=np.array([r["distance"] for r in rows]);low=np.array([r["lower"] for r in rows]);high=np.array([r["upper"] for r in rows])
    axes[0].errorbar(x[:4],d[:4],yerr=np.vstack([d[:4]-low[:4],high[:4]-d[:4]]),fmt="o",capsize=4,label="Full snapshot laws")
    axes[0].errorbar(x[4:],d[4:],yerr=np.vstack([d[4:]-low[4:],high[4:]-d[4:]]),fmt="o",capsize=4,label="Short SSA comparison")
    axes[0].axhline(.08,color="#a93226",linestyle="--",label="Registered margin")
    axes[0].set(xticks=x,xticklabels=["H\n100","S\n100","H\n1000","S\n1000","SSA\n100","SSA\n1000"],
                ylim=(-.002,.09),ylabel="Absolute CDF distance",xlabel="Initial state / Omega",
                title="All registered bounds meet the margin")
    axes[0].legend(fontsize=8);axes[0].grid(alpha=.2)
    axes[0].text(.02,.94,"Dots: empirical distance; whiskers: simultaneous bounds",transform=axes[0].transAxes,fontsize=8)
    for case in (4,5):
        maxima=[]
        for method in (1,2,3):
            source=evaluation["short"] if method==3 else evaluation["coarse_descriptive"]
            maxima.append(max(r["distance"] for r in source if r["case"]==case and r["method"]==method))
        axes[1].plot([1/6000,1/12000,1/24000],maxima,"o-",label=f"Omega{100 if case==4 else 1000}")
    axes[1].set_xscale("log");axes[1].set(xlabel="Maximum requested leap (checkpoint splitting applies)",ylabel="Maximum empirical CDF distance",
                                         title="Short-time refinement (coarser steps descriptive)",ylim=(0,.035))
    axes[1].legend();axes[1].grid(alpha=.2)
    fig.savefig(BASE/"MG10S2_summary.png",dpi=170);plt.close(fig)
    print(dict(case_bounds=rows,max_descriptive_diagonal_difference=result["max_descriptive_diagonal_difference"],max_descriptive_cross_difference=result["max_descriptive_cross_difference"]))


if __name__=="__main__":run()

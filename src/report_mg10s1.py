"""MG048 descriptive plot from registered ensemble decision records only."""
from pathlib import Path
import json
import numpy as np

ROOT=Path(__file__).resolve().parents[1]


def run():
    base=ROOT/"results/MG10S1"
    evaluation=json.loads((base/"production/evaluation.json").read_text())
    names=("Closed","Scarce B","Depleted","Transport")
    rows=[]
    for case in range(4):
        row=dict(case=case,name=names[case])
        for method in (1,2,3):
            source=evaluation["primary"] if method==3 else evaluation["coarse_descriptive"]
            row[f"distance_{method}"]=max(r["empirical_CDF_distance"] for r in source if r["case"]==case and r["method"]==method)
        row["upper"]=max(r["upper"] for r in evaluation["primary"] if r["case"]==case)
        rows.append(row)
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(1,2,figsize=(11,4.4),layout="constrained")
    x=np.arange(4);d=np.array([r["distance_3"] for r in rows]);upper=np.array([r["upper"] for r in rows])
    axes[0].errorbar(x,d,yerr=np.vstack([d,upper-d]),fmt="o",capsize=5,color="#185a7d",label="Finest step: worst marginal per case")
    axes[0].axhline(.05,color="#a93226",linestyle="--",label="Registered equivalence margin")
    axes[0].set(xticks=x,xticklabels=names,ylim=(-.001,.055),ylabel="Absolute CDF distance",
                title="All simultaneous upper bounds meet the margin")
    axes[0].legend(fontsize=8);axes[0].grid(alpha=.2)
    axes[0].text(.02,.02,"Dots: empirical distance\nWhiskers: simultaneous 99% bounds",transform=axes[0].transAxes,fontsize=8)
    for row in rows:
        axes[1].plot([.001,.0005,.00025],[row[f"distance_{i}"] for i in (1,2,3)],"o-",label=row["name"])
    axes[1].set_xscale("log")
    axes[1].set(xlabel="Maximum requested leap (adaptive safety applies)",ylabel="Maximum empirical CDF distance",
                title="Step refinement: descriptive, not an order fit",ylim=(0,.022))
    axes[1].legend(fontsize=8);axes[1].grid(alpha=.2)
    fig.savefig(base/"MG10S1_summary.png",dpi=170);plt.close(fig)
    (base/"descriptive_summary.json").write_text(json.dumps(rows,indent=2)+"\n",encoding="utf-8")
    print(rows)


if __name__=="__main__":run()

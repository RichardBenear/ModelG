"""MG045 descriptive figure from frozen follow-up raw data; no new runs."""
from pathlib import Path
import json
import numpy as np

ROOT=Path(__file__).resolve().parents[1]


def run():
    base=ROOT/"results/MG10C1";folder=base/"production"
    evaluation=json.loads((folder/"evaluation.json").read_text())
    with np.load(folder/"reduced_primary.npz") as raw:t,background=raw["times"],raw["state"]
    with np.load(folder/"full_zero_E2.npz") as raw:euler=raw["state"][:,1]
    with np.load(folder/"full_zero_Oplus.npz") as raw:ode=raw["state"][:,1]
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(1,2,figsize=(11.5,4.6),layout="constrained")
    axes[0].plot(t,background[:,1],color="#174a70",label="Uniform reduced ODE")
    axes[0].plot(t[::12],euler.mean(axis=1)[::12],"o",color="#174a70",markersize=3,label="Uniform Euler")
    axes[0].fill_between(t,ode.min(axis=1),ode.max(axis=1),alpha=.22,color="#b35415",label="Unseeded BDF: spatial range")
    axes[0].set(xlabel="Physical time",ylabel="X concentration",title="Exact uniform evolution versus numerical seeding")
    axes[0].legend(fontsize=8);axes[0].grid(alpha=.2)
    with np.load(folder/"linear_primary.npz") as raw:
        mode=evaluation["selected_mode"];v=raw["response"][:,1,mode-1]
    for label,epsilon,color in (("eps6",1e-6,"#b35415"),("eps8",1e-8,"#174a70")):
        window=next(w for w in evaluation["windows"] if w["id"]==f"full_{label}_E2")
        end=max(window["eligible_times"])
        axes[1].semilogy(t[:end+1],abs(epsilon*v[:end+1]),"--",color=color,label=f"Linear, epsilon={epsilon:g}")
        for method in ("E2","Oplus"):
            with np.load(folder/f"full_{label}_{method}.npz") as raw:x=raw["state"][:,1]
            x=(x-x.mean(axis=1,keepdims=True))/(140/9)
            amplitude=abs(2*np.mean(x*np.exp(-2j*np.pi*mode*np.arange(500)/500),axis=1))
            if method=="E2":axes[1].semilogy(t[:91:4],amplitude[:91:4],"o",color=color,markersize=3)
            else:axes[1].semilogy(t[:91],amplitude[:91],color=color,alpha=.7)
    axes[1].set(xlabel="Physical time",ylabel="Selected X mode amplitude / X0",
                title="Specified disturbances follow the linear prediction",ylim=(1e-10,2))
    axes[1].legend(fontsize=8);axes[1].grid(alpha=.2)
    axes[1].text(.02,.02,"Solid: BDF; dots: Euler\nDashed prediction stops at the linear-window bound",transform=axes[1].transAxes,fontsize=8)
    fig.savefig(base/"MG10C1_summary.png",dpi=170);plt.close(fig)
    diagnostics={"scope":"POSTHOC_DESCRIPTION_OF_REGISTERED_RAW_DATA_ONLY","selected_mode":mode,
                 "selected_wavelength":100/mode,"unseeded_BDF_range_onset_above_1e_minus6":int(np.where(np.ptp(ode,axis=1)>1e-6)[0][0]),
                 "uniform_Euler_exact_spatial_range":float(np.ptp(euler,axis=1).max()),
                 "uniform_Euler_resolved_counts_max":max(evaluation["resolved_minimum_counts"]["full_zero_E2"]),
                 "formal_peak_gain_not_a_physical_amplitude":evaluation["selected_peak_gain"]}
    (base/"descriptive_diagnostics.json").write_text(json.dumps(diagnostics,indent=2)+"\n",encoding="utf-8")
    print(diagnostics)


if __name__=="__main__":run()

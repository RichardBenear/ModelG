"""MG066 separate episodes, direct interpolation and temporal measurements."""
from pathlib import Path
from functools import lru_cache
import numpy as np
from mg10p2_oracle import tools,circular
from check_mg10p3_preparation import independent_episodes

ROOT=Path(__file__).resolve().parents[1]


@lru_cache(maxsize=1)
def averaging_tools():
    _,smooth=tools()
    with np.load(ROOT/"results/MG003/ref_profiles.npz") as raw:profile=raw["h0.2_b29.0"][1:]
    template=np.array([smooth(f) for f in profile]);norm2=np.sum(template**2,axis=1)
    offsets=(np.arange(500)[:,None]-np.arange(500)[None,:])%500
    return smooth,template,norm2,offsets


def average(fields,phases=None):
    smooth,template,norm2,offsets=averaging_tools()
    values=np.array([[smooth(f) for f in frame] for frame in fields])
    if phases is not None:
        shifted=[]
        for frame,p in zip(values,phases):
            d=(np.arange(500)+p/.2+250)%500-250
            weights=np.ones(500);mask=abs(d)>1e-12;weights[mask]=np.sin(np.pi*d[mask])/(500*np.tan(np.pi*d[mask]/500))
            matrix=weights[offsets];shifted.append(np.array([matrix@f for f in frame]))
        values=np.array(shifted)
    projection=np.sum(values*template[None],axis=-1)/norm2[None]
    individual=np.sqrt(np.mean((values-projection[:,:,None]*template[None])**2,axis=-1))
    mean=values.mean(axis=0);coeff=np.sum(mean*template,axis=1)/norm2
    residual=np.sqrt(np.mean((mean-coeff[:,None]*template)**2,axis=1));mean_individual=individual.mean(axis=0)
    return dict(profile=mean,projection=coeff,residual_rms=residual,contrast_rms=np.sqrt(np.mean(mean**2,axis=1)),
        mean_individual_residual_rms=mean_individual,gain=[float(mean_individual[j]/residual[j]) if residual[j]>1e-12 else None for j in range(3)],
        individual_projection=projection,individual_residual_rms=individual,frames=len(fields),alignment="X_phase" if phases is not None else "fixed_frame")


def summarize(data,times,start):
    fields=data[:,:3].astype(float)/1000;assess,_=tools();observations=[assess(f[1]) for f in fields]
    runs=independent_episodes(observations,times);track=None
    if start:
        track=dict(status="INITIAL_UNRESOLVED",matched_through=None,visible_duration=None,persistent_reference_signal=False,censor_reason=None,observations=observations)
        if observations[0]["status"]=="REFERENCE_MATCH" and circular(observations[0]["position"],0)<=.1+1e-10:
            run=runs[0];reason=run["termination"]
            reason=None if reason=="END_OF_RECORD" else "MOTION_AMBIGUITY" if reason=="MOTION_SPLIT" else reason
            track.update(status="CONTINUOUS" if reason is None else "CENSORED",matched_through=run["last"],visible_duration=run["duration"],persistent_reference_signal=run["persistent"],censor_reason=reason)
    fixed=average(fields[:81]);aligned=None
    if start and track["persistent_reference_signal"]:aligned=average(fields[:81],[v["position"] for v in observations[:81]])
    return dict(track=track,observations=observations,episodes=runs,any_observed_episode20=any(v["persistent"] for v in runs),fixed=fixed,aligned=aligned,
                event20=bool(track["persistent_reference_signal"] if start else any(v["persistent"] for v in runs)))

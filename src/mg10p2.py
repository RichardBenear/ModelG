"""MG061 prepared states and observability metrics; P2 production is draft."""
from pathlib import Path
import numpy as np
from scipy.stats import beta
from mg10p1 import centered_smooth,assess
from mg10d2 import fixture as old_fixture

ROOT=Path(__file__).resolve().parents[1]
OBS=np.array([0.,.05,.1,.25])
REPS=128
CHUNK=8


def reference():
    with np.load(ROOT/"results/MG003/ref_profiles.npz") as raw:profile=raw["h0.2_b29.0"][1:].copy()
    return profile,np.array([centered_smooth(f) for f in profile])


def fixture(condition,start):
    if condition not in (0,1) or start not in (0,1):raise ValueError("Invalid condition/start")
    case=old_fixture(0,0)
    if condition:
        case["initial"][3:5]=np.rint(.82*case["initial"][3:5]).astype(np.int64)
        case["clamp_a"]=False;case["clamp_b"]=False
    a=14*case["initial"][3,0]/case["a0"][0]
    b=29*case["initial"][4,0]/case["b0"][0]
    background=np.array([a/.9,a/.9,b/(a/.9)])
    fields=np.repeat(background[:,None],500,axis=1)
    if start:fields+=reference()[0]
    case["initial"][:3]=np.rint(1000*fields).astype(np.int64)
    case.update(condition=condition,start=start,background=background,arm="FULL_CHEMO" if condition==0 else "PRESCRIBED_DEPLETED_FINITE")
    return case


def measure(data):
    _,template=reference();norm2=np.sum(template**2,axis=-1)
    fields=data[:,:,:3]/1000.
    offset=fields-fields[:,:,:,:1];mean=fields[:,:,:,0]+offset.mean(axis=-1)
    centered=offset-offset.mean(axis=-1,keepdims=True)
    smoothed=np.empty_like(centered)
    for r in range(len(fields)):
        for t in range(fields.shape[1]):
            for s in range(3):smoothed[r,t,s]=centered_smooth(fields[r,t,s])
    projection=np.sum(smoothed*template[None,None],axis=-1)/norm2[None,None]
    residual=smoothed-projection[:,:,:,None]*template[None,None]
    matches=[[assess(frame[1]) for frame in path] for path in fields]
    return dict(field_mean=mean,contrast_rms=np.sqrt(np.mean(centered**2,axis=-1)),smooth_rms=np.sqrt(np.mean(smoothed**2,axis=-1)),
                projection=projection,residual_rms=np.sqrt(np.mean(residual**2,axis=-1))),matches


def upper_bound(k,n,alpha=.0125):
    if n<1 or k<0 or k>n:raise ValueError("Invalid binomial counts")
    return 1. if k==n else float(beta.ppf(1-alpha,k+1,n-k))


def binomial_interval(k,n,alpha=.05):
    if n<1 or k<0 or k>n:raise ValueError("Invalid binomial counts")
    return [0. if k==0 else float(beta.ppf(alpha/2,k,n-k+1)),1. if k==n else float(beta.ppf(1-alpha/2,k+1,n-k))]

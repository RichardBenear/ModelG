"""MG062 direct instrument oracle, adapted from qualified P1 dense/direct checker.

Does not import P1/P2 measurement or classification functions.
"""
from pathlib import Path
from functools import lru_cache
import json
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
def circular(a,b):return min((a-b)%100.,(b-a)%100.)

@lru_cache(maxsize=1)
def tools():
    registration=json.loads((ROOT/"results/MG10P1/registration.json").read_text(encoding="utf-8"))
    p=registration["parameters"];n=p["N"];h=p["h"];length=p["L"]
    k=2*np.pi*np.fft.fftfreq(n,d=h)
    kernel=np.fft.ifft(np.exp(-.5*(p["sigma"]*k)**2)).real
    convolution=kernel[(np.arange(n)[:,None]-np.arange(n)[None,:])%n]
    def smooth(field):
        z=field-field[0];z=z-np.sum(z)/n
        y=convolution@z
        return y-np.sum(y)/n
    with np.load(ROOT/"results/MG003/ref_profiles.npz") as raw:template=smooth(raw["h0.2_b29.0"][2])
    basis=np.array([np.roll(template,shift) for shift in range(n)])
    norm2=np.sum(template**2);grid=(np.arange(n)-n//2)*h
    def assess(field):
        y=smooth(field);energy=float(np.sum(y**2))
        floor=1e-8*max(1.,140/9,float(np.max(np.abs(template))))
        if np.sqrt(energy/n)<=floor:return dict(status="FLAT",position=None,amplitude=0.,error=None,local_energy=None)
        scores=basis@y;peak=float(np.max(scores));best=int(np.where(scores>=peak-1e-12*max(1.,abs(peak)))[0][0]);a=float(scores[best]/norm2)
        result=dict(status="NO_REFERENCE_MATCH",position=None,amplitude=a,error=None,local_energy=None)
        if min(abs(a-p["amp_min"]),abs(a-p["amp_max"]))<=p["border"]:result["status"]="BORDERLINE";return result
        if a<=p["amp_min"] or a>=p["amp_max"]:return result
        far=np.array([circular(i*h,best*h)>p["ambiguity_radius"] for i in range(n)])
        if max(scores[far])>=p["ambiguity_ratio"]*scores[best]:result["status"]="AMBIGUOUS_MATCH";return result
        residual=y-a*basis[best];error=float(np.sqrt(np.sum(residual**2))/(a*np.sqrt(norm2)))
        position=(best*h+length/2)%length-length/2
        local=np.array([circular(v,position)<10 for v in grid]);fraction=float(np.sum(y[local]**2)/energy)
        result.update(error=error,local_energy=fraction)
        if abs(error-p["error_max"])<=p["border"] or abs(fraction-p["energy_min"])<=p["border"]:result["status"]="BORDERLINE";return result
        if error>p["error_max"] or fraction<p["energy_min"]:return result
        c0,c1,c2=scores[(best-1)%n],scores[best],scores[(best+1)%n];curvature=c0+c2-2*c1
        delta=.5*(c0-c2)/curvature if curvature<0 else 0.
        delta=min(.5,max(-.5,float(delta)))
        result.update(status="REFERENCE_MATCH",position=float(((best+delta)*h+length/2)%length-length/2))
        return result
    return assess,smooth


def independent_measure(data):
    assess,smooth=tools()
    with np.load(ROOT/"results/MG003/ref_profiles.npz") as raw:profiles=raw["h0.2_b29.0"][1:]
    template=np.array([smooth(f) for f in profiles]);norm2=np.sum(template**2,axis=1)
    fields=data[:,:,:3].astype(float)/1000
    smoothed=np.array([[np.array([smooth(f) for f in frame]) for frame in path] for path in fields])
    means=fields.sum(axis=-1)/500
    centered=fields-fields.mean(axis=-1,keepdims=True)
    projection=np.sum(smoothed*template[None,None],axis=-1)/norm2[None,None]
    return dict(field_mean=means,contrast_rms=np.sqrt(np.mean(centered**2,axis=-1)),smooth_rms=np.sqrt(np.mean(smoothed**2,axis=-1)),
                projection=projection,residual_rms=np.sqrt(np.mean((smoothed-projection[:,:,:,None]*template[None,None])**2,axis=-1))),[[assess(frame[1]) for frame in path] for path in fields]

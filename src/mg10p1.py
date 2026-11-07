"""MG058 single reference-like X motif and conservative sampled continuity.

An instrument qualification, not particle counting, lineage or physical survival.
"""
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
N=500;H=.2;L=100.;SIGMA=.4
AMP_MIN=.5;AMP_MAX=1.5;ERROR_MAX=.5;ENERGY_MIN=.8
AMBIG_RADIUS=2.;AMBIG_RATIO=.9;BORDER=1e-10
MAX_SPEED=1.5;MAX_GAP=1.;PERSIST=20.


def distance(a,b):return abs((a-b+L/2)%L-L/2)


def centered_smooth(field):
    offset=np.asarray(field,dtype=float)-field[0]
    offset-=offset.mean()
    k=2*np.pi*np.fft.fftfreq(N,d=H)
    result=np.fft.ifft(np.fft.fft(offset)*np.exp(-.5*(SIGMA*k)**2)).real
    return result-result.mean()


def reference():
    with np.load(ROOT/"results/MG003/ref_profiles.npz") as raw:profile=raw["h0.2_b29.0"][2].copy()
    return profile,centered_smooth(profile)


def assess(field):
    field=np.asarray(field,dtype=float)
    if field.shape!=(N,) or not np.all(np.isfinite(field)) or np.any(field<0):raise ValueError("Invalid concentration field")
    _,template=reference();smooth=centered_smooth(field)
    norm2=float(template@template);energy=float(smooth@smooth)
    floor=1e-8*max(1.,140/9,float(np.max(np.abs(template))))
    if np.sqrt(energy/N)<=floor:return dict(status="FLAT",position=None,amplitude=0.,error=None,local_energy=None)
    scores=np.fft.ifft(np.fft.fft(smooth)*np.conj(np.fft.fft(template))).real
    peak_score=float(np.max(scores))
    best=int(np.flatnonzero(scores>=peak_score-1e-12*max(1.,abs(peak_score)))[0]);amplitude=float(scores[best]/norm2)
    result=dict(status="NO_REFERENCE_MATCH",position=None,amplitude=amplitude,error=None,local_energy=None)
    if abs(amplitude-AMP_MIN)<=BORDER or abs(amplitude-AMP_MAX)<=BORDER:
        result["status"]="BORDERLINE";return result
    if not AMP_MIN<amplitude<AMP_MAX:return result
    shifts=np.arange(N)*H
    far=np.minimum((shifts-best*H)%L,(best*H-shifts)%L)>AMBIG_RADIUS
    if float(np.max(scores[far]))>=AMBIG_RATIO*float(scores[best]):
        result["status"]="AMBIGUOUS_MATCH";return result
    aligned=np.roll(template,best)
    error=float(np.linalg.norm(smooth-amplitude*aligned)/(amplitude*np.sqrt(norm2)))
    grid=(np.arange(N)-N//2)*H
    position=((best*H+L/2)%L)-L/2
    local=np.minimum((grid-position)%L,(position-grid)%L)<10.
    fraction=float(np.sum(smooth[local]**2)/energy)
    result.update(error=error,local_energy=fraction)
    if abs(error-ERROR_MAX)<=BORDER or abs(fraction-ENERGY_MIN)<=BORDER:
        result["status"]="BORDERLINE";return result
    if error>ERROR_MAX or fraction<ENERGY_MIN:return result
    left,peak,right=scores[(best-1)%N],scores[best],scores[(best+1)%N]
    curvature=left-2*peak+right
    subpixel=float(.5*(left-right)/curvature) if curvature<0 else 0.
    subpixel=max(-.5,min(.5,subpixel))
    result.update(status="REFERENCE_MATCH",position=float(((best+subpixel)*H+L/2)%L-L/2))
    return result


def continuity(fields,times,x0=0.):
    times=np.asarray(times,dtype=float);fields=np.asarray(fields,dtype=float)
    if times.ndim!=1 or not len(times) or not np.all(np.isfinite(times)) or np.any(np.diff(times)<=0) or fields.shape!=(len(times),N) or not np.isfinite(x0):raise ValueError("Invalid sequence/times/initial position")
    observations=[assess(f) for f in fields]
    result=dict(status="INITIAL_UNRESOLVED",matched_through=None,visible_duration=None,persistent_reference_signal=False,
                censor_reason=None,observations=observations,identity_scope="Sampled field continuity only; not token/particle lineage or physical survival")
    if observations[0]["status"]!="REFERENCE_MATCH" or distance(observations[0]["position"],x0)>H/2+BORDER:return result
    last=0;reason=None
    for index in range(1,len(times)):
        if times[index]-times[index-1]>MAX_GAP+BORDER:reason="SAMPLING_GAP";break
        if observations[index]["status"]!="REFERENCE_MATCH":reason=observations[index]["status"];break
        if distance(observations[index]["position"],observations[index-1]["position"])>MAX_SPEED*(times[index]-times[index-1])+BORDER:
            reason="MOTION_AMBIGUITY";break
        last=index
    duration=float(times[last]-times[0])
    result.update(status="CONTINUOUS" if reason is None else "CENSORED",matched_through=float(times[last]),visible_duration=duration,
                  persistent_reference_signal=duration>=PERSIST-BORDER,censor_reason=reason)
    return result

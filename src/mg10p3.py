"""MG065 P3 temporal/averaging measurements; scientific protocol is draft."""
import numpy as np
from mg10p1 import assess,continuity,distance,MAX_GAP,MAX_SPEED,BORDER,PERSIST,centered_smooth
from mg10p2 import reference

OBS=np.arange(101,dtype=float)/4


def episodes(observations,times):
    times=np.asarray(times,dtype=float)
    if len(observations)!=len(times) or not len(times) or not np.all(np.isfinite(times)) or np.any(np.diff(times)<=0):raise ValueError("Invalid episode observations/times")
    out=[];first=None;last=None
    def close(reason):
        if first is not None:
            duration=float(times[last]-times[first])
            out.append(dict(start=float(times[first]),last=float(times[last]),duration=duration,persistent=duration>=PERSIST-BORDER,termination=reason))
    for i,item in enumerate(observations):
        if item["status"]!="REFERENCE_MATCH":
            close(item["status"]);first=last=None;continue
        if first is not None:
            dt=times[i]-times[last]
            reason="SAMPLING_GAP" if dt>MAX_GAP+BORDER else "MOTION_SPLIT" if distance(item["position"],observations[last]["position"])>MAX_SPEED*dt+BORDER else None
            if reason:close(reason);first=last=None
        if first is None:first=i
        last=i
    close("END_OF_RECORD")
    return out


def align(field,position):
    return np.fft.ifft(np.fft.fft(field)*np.exp(2j*np.pi*np.fft.fftfreq(500)*position/.2)).real


def average(fields,phases=None):
    fields=np.asarray(fields,dtype=float)
    if fields.ndim!=3 or fields.shape[1:]!=(3,500) or not len(fields) or not np.all(np.isfinite(fields)):raise ValueError("Invalid averaging fields")
    if phases is not None and (len(phases)!=len(fields) or not np.all(np.isfinite(phases))):raise ValueError("Invalid phases")
    _,template=reference();norm2=np.sum(template**2,axis=1)
    smooth=np.array([[centered_smooth(f) for f in frame] for frame in fields])
    if phases is not None:smooth=np.array([[align(f,p) for f in frame] for frame,p in zip(smooth,phases)])
    projection=np.sum(smooth*template[None],axis=-1)/norm2[None]
    residual=smooth-projection[:,:,None]*template[None]
    individual=np.sqrt(np.mean(residual**2,axis=-1))
    averaged=smooth.mean(axis=0);amplitude=np.sum(averaged*template,axis=1)/norm2
    error=np.sqrt(np.mean((averaged-amplitude[:,None]*template)**2,axis=1))
    mean_individual=individual.mean(axis=0)
    gain=[float(mean_individual[j]/error[j]) if error[j]>1e-12 else None for j in range(3)]
    return dict(profile=averaged,projection=amplitude,residual_rms=error,contrast_rms=np.sqrt(np.mean(averaged**2,axis=1)),
                mean_individual_residual_rms=mean_individual,gain=gain,individual_projection=projection,individual_residual_rms=individual,
                frames=len(fields),alignment="X_phase" if phases is not None else "fixed_frame")


def summarize(data,times,start):
    times=np.asarray(times,dtype=float)
    if not np.array_equal(times,OBS):raise ValueError("Full0..25 grid required; incomplete paths are not complete estimates")
    fields=np.asarray(data)[:,:3]/1000
    window=times<=20+BORDER
    if not np.array_equal(times[window],np.arange(81)/4):raise ValueError("Incomplete fixed0..20 window")
    if start:
        track=continuity(fields[:,1],times,x0=0.)
        observations=track["observations"]
    else:track=None;observations=[assess(f[1]) for f in fields]
    runs=episodes(observations,times)
    fixed=average(fields[window]);aligned=None
    if start and track["persistent_reference_signal"]:
        aligned=average(fields[window],[v["position"] for v in observations[:81]])
    return dict(track=track,episodes=runs,any_observed_episode20=any(r["persistent"] for r in runs),fixed=fixed,aligned=aligned,
                identity_scope="Sampled reference-field continuity, not physical particle/token lineage or death")

"""MG073 E2 prospective schedules and field-response measurements."""
import numpy as np
from mg10e1 import activity as e1_activity,fork_case,ledger,primary_interval
from mg10e1_stream import advance,seed_rng
from mg10p1 import centered_smooth
from mg10p2 import reference

ARMS=("UNCHANGED","LOCAL_CONTINUOUS_OFF","LOCAL_OFF_RESTORE","REMOTE_OFF_RESTORE")
OBS=np.arange(101,dtype=float)/4
SWITCH=5.
def activity(arm,t,switch=SWITCH):
    if arm not in range(4) or not np.isfinite(t) or t<0 or not np.isfinite(switch) or switch<0:raise ValueError("Invalid schedule")
    if arm==0 or (arm in (2,3) and t>=switch):return np.ones(500)
    return e1_activity(2 if arm==3 else 1)
def advance_interval(case,state,flux,hops,external,tau,start,target,arm,switch=SWITCH):
    if not np.isfinite(start) or not np.isfinite(target) or start<0 or target<start or not np.isfinite(tau) or tau<=0:raise ValueError("Invalid interval")
    # Split even an observation interval that straddles the intervention event.
    boundaries=[start]+([switch] if start<switch<target and arm in (2,3) else [])+[target]
    stats=np.zeros(3)
    for left,right in zip(boundaries[:-1],boundaries[1:]):
        alpha=activity(arm,left,switch)
        d=advance(state,flux,hops,external,case["a0"],case["b0"],case["omega"],case["diff"],np.ones(7),tau,left,right,case["clamp_a"],case["clamp_b"],alpha)
        stats[:2]+=d[:2];stats[2]=max(stats[2],d[2])
    return stats
def simulate(case,tau,times,seed,arm,switch=SWITCH):
    times=np.asarray(times,dtype=float)
    if times.ndim!=1 or not len(times) or times[0]!=0 or not np.all(np.isfinite(times)) or np.any(np.diff(times)<=0):raise ValueError("Invalid checkpoints")
    activity(arm,0,switch)
    state=case["initial"].copy();n=state.shape[1]
    if n!=500:raise ValueError("E2 full lattice required")
    flux=np.zeros((7,n),dtype=np.int64);hops=np.zeros((6,n),dtype=np.int64);external=np.zeros((2,n),dtype=np.int64)
    out=np.zeros((1,len(times),22,n),dtype=np.int64);stats=np.zeros((1,3));out[0,0,:7]=state
    advance_interval(case,state.copy(),flux.copy(),hops.copy(),external.copy(),tau,0,0,arm,switch);seed_rng(seed)
    for i in range(1,len(times)):
        d=advance_interval(case,state,flux,hops,external,tau,times[i-1],times[i],arm,switch);stats[0,:2]+=d[:2];stats[0,2]=max(stats[0,2],d[2]);out[0,i]=np.concatenate([state,flux,hops,external])
    return out,stats
def window_average(values,times,left=20.,right=25.):
    values=np.asarray(values,dtype=float);times=np.asarray(times,dtype=float)
    if values.shape[0]!=len(times) or not np.all(np.isfinite(values)) or not np.all(np.isfinite(times)) or np.any(np.diff(times)<=0):raise ValueError("Invalid window values")
    select=(times>=left)&(times<=right)
    if not np.array_equal(times[select],np.arange(round(left*4),round(right*4)+1)/4) or right<=left:raise ValueError("Incomplete window")
    return np.trapezoid(values[select],x=times[select],axis=0)/(right-left)
def paired_effects(late):
    late=np.asarray(late,dtype=float)
    if late.shape!=(2,4,3) or not np.all(np.isfinite(late)):raise ValueError("Invalid paired block")
    excess=late[1]-late[0]
    return dict(primary_restore_continuous=excess[2]-excess[1],restore_sham=excess[2]-excess[0],remote_restore_sham=excess[3]-excess[0],restore_remote=excess[2]-excess[3],seed_excess_by_arm=excess)
def balanced_analysis(records):
    effects={};summary={}
    for k in range(2):
        values=[]
        for b in range(8):
            keys=[f"a{a}_s{s}_k{k}_b{b}" for s in range(2) for a in range(4)]
            if not all(key in records and records[key]["complete"] for key in keys):continue
            late=np.array([[records[f"a{a}_s{s}_k{k}_b{b}"]["late_projection"] for a in range(4)] for s in range(2)])
            result=paired_effects(late);effects[f"k{k}_b{b}"]={name:value.tolist() for name,value in result.items()};values.append(result["primary_restore_continuous"][1])
        summary[str(k)]=primary_interval(values)
    return effects,summary
def field_metrics(fields):
    fields=np.asarray(fields,dtype=float)
    if fields.ndim!=3 or fields.shape[1:]!=(3,500) or not np.all(np.isfinite(fields)) or np.any(fields<0):raise ValueError("Invalid physical fields")
    _,template=reference();norm=np.sum(template**2,axis=-1)
    smooth=np.array([[centered_smooth(f) for f in frame] for frame in fields])
    fixed=np.sum(smooth*template[None],axis=-1)/norm[None]
    positions=[];aligned=[];errors=[]
    for frame in smooth:
        scores=np.fft.ifft(np.fft.fft(frame[1])*np.conj(np.fft.fft(template[1]))).real
        best=int(np.flatnonzero(scores>=scores.max()-1e-12*max(1.,abs(scores.max())))[0])
        shifted=np.roll(template,best,axis=-1);coefficient=np.sum(frame*shifted,axis=-1)/norm
        positions.append(float((best*.2+50)%100-50));aligned.append(coefficient)
        errors.append(np.sqrt(np.mean((frame-coefficient[:,None]*shifted)**2,axis=-1)))
    return dict(projection=fixed,field_mean=fields.mean(axis=-1),smooth_rms=np.sqrt(np.mean(smooth**2,axis=-1)),x_phase_grid=np.array(positions),
                x_aligned_projection=np.array(aligned),x_aligned_residual_rms=np.array(errors))

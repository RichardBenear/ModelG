"""MG055 prospective D2 source correction; original D1 source remains frozen.

Only finite-A sampling delegates to a stable rare-event binomial algorithm.
Low stocks retain the original RNG calls; no Poisson approximation is made.
"""
from pathlib import Path
import numpy as np
import numba as nb
from mg10s1 import _exits4, STOICH
from mg10d2_source import finite_source
from mg10s2 import fixture as s2_fixture

ARMS=("CHEMO","HIGH","BOTH","A_ONLY","B_ONLY")
ROOT=Path(__file__).resolve().parents[1]

def fixture(arm,start):
    if arm not in range(5) or start not in range(2):raise ValueError("Invalid arm/start")
    case=s2_fixture(2+start)
    if arm==1:
        case["initial"][3:5]*=1000000
        case["a0"]*=1000000;case["b0"]*=1000000
    case.update(arm=ARMS[arm],clamp_a=arm in (0,4),clamp_b=arm in (0,3))
    return case

@nb.njit(cache=True)
def _controlled_leap(state,flux,hops,a0,b0,omega,diff,mask,requested,clamp_a,clamp_b,external):
    n=state.shape[1];maxrate=0.
    for i in range(n):
        x,b=state[1,i],state[4,i]
        r4=mask[4]*x*(x-1)/(omega*omega) if x>=2 else 0.
        rate3=mask[3]*29.*x if clamp_b else mask[3]*29.*x*b/b0[i]
        supply=x if clamp_b else min(x,b)
        joint=rate3/supply if supply>0 else 0.
        maxrate=max(maxrate,0. if clamp_a else mask[0]*omega*14./a0[i],mask[1]+2*diff[0],
                    mask[2]*.1+mask[5]*.9+(mask[3]*29. if clamp_b else mask[3]*29.*b/b0[i])+2*diff[1],r4+2*diff[2],joint)
    dt=min(requested,.2/maxrate) if maxrate>0 else requested
    delta=np.zeros_like(state);largest=0.
    for i in range(n):
        left,right=(i-1)%n,(i+1)%n
        g,x,y,a,b=state[0,i],state[1,i],state[2,i],state[3,i],state[4,i]
        pa=0. if clamp_a else mask[0]*omega*14.*dt/a0[i]
        source=np.random.poisson(mask[0]*omega*14.*dt) if clamp_a else finite_source(a,pa)
        gx,gl,gr,_=_exits4(g,mask[1]*dt,diff[0]*dt,diff[0]*dt,0.)
        rate3=mask[3]*29.*x if clamp_b else mask[3]*29.*x*b/b0[i]
        m=x if clamp_b else min(x,b);joint=rate3*dt/m if m>0 else 0.
        conversion=np.random.binomial(m,joint)
        p3x=rate3*dt/x if x>0 else 0.
        denominator=1-p3x
        pg=mask[2]*.1*dt/denominator;pw=mask[5]*.9*dt/denominator
        ph=diff[1]*dt/denominator
        xg,xw,xl,xr=_exits4(x-conversion,pg,pw,ph,ph)
        py=mask[4]*x*(x-1)*dt/(omega*omega) if x>=2 else 0.
        yx,yl,yr,_=_exits4(y,py,diff[2]*dt,diff[2]*dt,0.)
        delta[0,i]+=source-gx-gl-gr+xg
        delta[1,i]+=gx-xg-xw-xl-xr-conversion+yx
        delta[2,i]+=conversion-yx-yl-yr
        if clamp_a:external[0,i]+=source
        else:delta[3,i]-=source
        if clamp_b:external[1,i]+=conversion
        else:delta[4,i]-=conversion
        delta[5,i]+=conversion;delta[6,i]+=xw
        delta[0,left]+=gl;delta[0,right]+=gr
        delta[1,left]+=xl;delta[1,right]+=xr
        delta[2,left]+=yl;delta[2,right]+=yr
        flux[0,i]+=source;flux[1,i]+=gx;flux[2,i]+=xg
        flux[3,i]+=conversion;flux[4,i]+=yx;flux[5,i]+=xw
        hops[0,i]+=gl;hops[1,i]+=gr;hops[2,i]+=xl;hops[3,i]+=xr;hops[4,i]+=yl;hops[5,i]+=yr
        largest=max(largest,pa,(mask[1]+2*diff[0])*dt,py+2*diff[2]*dt,joint,pg+pw+2*ph)
    state+=delta
    if np.any(state<0):raise ValueError("Negative count")
    if largest>.25+1e-14:raise ValueError("Safety bound violated")
    return dt,largest


@nb.njit(cache=True)
def controlled_chunk(initial,a0,b0,omega,diff,mask,tau,observations,reps,seed,clamp_a,clamp_b):
    np.random.seed(seed)
    n=initial.shape[1];out=np.zeros((reps,len(observations),22,n),dtype=np.int64);stats=np.zeros((reps,3))
    for r in range(reps):
        state=initial.copy();flux=np.zeros((7,n),dtype=np.int64);hops=np.zeros((6,n),dtype=np.int64)
        external=np.zeros((2,n),dtype=np.int64)
        t=0.;leaps=0;adapted=0;maximum=0.
        for k,target in enumerate(observations):
            while t<target:
                requested=min(tau,target-t)
                dt,prob=_controlled_leap(state,flux,hops,a0,b0,omega,diff,mask,requested,clamp_a,clamp_b,external)
                if dt<=0:raise ValueError("Nonpositive step")
                adapted+=int(dt<requested*(1-1e-12));leaps+=1;maximum=max(maximum,prob)
                if dt==target-t:t=target
                else:t+=dt
            out[r,k,:7]=state;out[r,k,7:14]=flux;out[r,k,14:20]=hops;out[r,k,20:22]=external
        stats[r,0]=leaps;stats[r,1]=adapted;stats[r,2]=maximum
    return out,stats


def simulate(case,tau,observations,reps,seed,mask=None):
    mask=np.ones(7) if mask is None else np.asarray(mask,dtype=float)
    obs=np.asarray(observations,dtype=float)
    if not np.isfinite(tau) or tau<=0 or reps<=0:raise ValueError("Invalid step/replicates")
    if obs.ndim!=1 or not len(obs) or not np.all(np.isfinite(obs)) or obs[0]<0 or np.any(np.diff(obs)<=0):raise ValueError("Invalid checkpoints")
    initial=case["initial"]
    if initial.dtype!=np.int64 or initial.shape[0]!=7 or np.any(initial<0):raise ValueError("Invalid counts")
    for name in ("a0","b0"):
        if case[name].shape!=(initial.shape[1],) or not np.all(np.isfinite(case[name])) or np.any(case[name]<=0):raise ValueError("Invalid denominators")
    if not np.isfinite(case["omega"]) or case["omega"]<=0 or mask.shape!=(7,) or not np.all(np.isfinite(mask)) or np.any(mask<0):raise ValueError("Invalid rates")
    if case["diff"].shape!=(3,) or not np.all(np.isfinite(case["diff"])) or np.any(case["diff"]<0):raise ValueError("Invalid diffusion")
    return controlled_chunk(initial,case["a0"],case["b0"],case["omega"],case["diff"],mask,tau,obs,reps,seed,case["clamp_a"],case["clamp_b"])


def ledger(data,initial):
    expected=np.einsum("sr,btrn->btsn",STOICH,data[:,:,7:14])
    expected[:,:,3:5]+=data[:,:,20:22]
    for species in range(3):
        left,right=data[:,:,14+2*species],data[:,:,15+2*species]
        expected[:,:,species]+=np.roll(left,-1,axis=-1)+np.roll(right,1,axis=-1)-left-right
    residual=int(np.max(np.abs(data[:,:,:7]-initial[None,None]-expected)))
    total=int(np.max(np.abs(data[:,:,:7].sum(axis=(2,3))-initial.sum()-data[:,:,20:22].sum(axis=(2,3)))))
    return dict(max_integer_residual=residual,max_total_input_residual=total,nonnegative=bool(np.all(data>=0)),pass_check=bool(residual==0 and total==0 and np.all(data>=0)))


def measurements(data,case):
    active=data[:,:,:3]/case["omega"]
    offset=active-active[:,:,:,:1]
    offset_mean=offset.mean(axis=-1)
    means=active[:,:,:,0]+offset_mean
    rms=np.sqrt(np.mean((offset-offset_mean[:,:,:,None])**2,axis=-1))
    ratios=data[:,:,3:5]/np.stack([case["a0"],case["b0"]])[None,None]
    return dict(field_mean=means,field_rms=rms,availability_mean=ratios.mean(axis=-1),availability_min=ratios.min(axis=-1),availability_max=ratios.max(axis=-1),availability_sd=ratios.std(axis=-1),reaction_totals=data[:,:,7:14].sum(axis=-1),external_totals=data[:,:,20:22].sum(axis=-1))

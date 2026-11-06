"""MG046 P17 finite-reactant stochastic qualification kernels; integer state.

Snapshot rows: seven populations, seven reaction counts, six directed-hop counts.
The exact SSA propensity implementation does not call leap hazard/draw helpers.
"""
from pathlib import Path
import numpy as np
import numba as nb

ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/"results/MG10S1"
OBS=np.array([.05,.10,.25])
R=32768
CHUNK=512
TAUS=(0.,.001,.0005,.00025)
STOICH=np.array([[1,-1,1,0,0,0,-1],[0,1,-1,-1,1,-1,0],[0,0,0,1,-1,0,0],
                 [-1,0,0,0,0,0,1],[0,0,0,-1,0,0,0],[0,0,0,1,0,0,0],[0,0,0,0,0,1,0]],dtype=np.int64)


def fixture(index):
    if index==0:
        name="closed";omega=1.;state=np.array([[16],[16],[2],[40],[40],[0],[0]],dtype=np.int64);a0=b0=np.array([40.]);diff=np.zeros(3)
    elif index==1:
        name="scarce_B";omega=1.;state=np.array([[2],[12],[3],[8],[1],[0],[0]],dtype=np.int64);a0=np.array([8.]);b0=np.array([1.]);diff=np.zeros(3)
    elif index==2:
        name="depleted";omega=1.;state=np.array([[2],[8],[2],[0],[0],[8],[8]],dtype=np.int64);a0=b0=np.array([8.]);diff=np.zeros(3)
    elif index==3:
        name="transport";omega=2.;state=np.array([[10,14,18],[12,20,28],[1,2,4],[20,40,10],[1,8,40],[0,0,0],[0,0,0]],dtype=np.int64);a0=b0=np.full(3,40.);diff=np.array([1.,1.,12.])
    else:raise ValueError(index)
    return dict(name=name,state=state,a0=a0.copy(),b0=b0.copy(),omega=omega,diff=diff)


@nb.njit(cache=True)
def _exits4(n,p0,p1,p2,p3):
    if min(p0,p1,p2,p3)<0 or p0+p1+p2+p3>1:raise ValueError("Invalid exit probabilities")
    if n==0:return 0,0,0,0
    k0=np.random.binomial(n,p0)
    if k0==n:return k0,0,0,0
    k1=np.random.binomial(n-k0,p1/(1-p0))
    if k0+k1==n:return k0,k1,0,0
    k2=np.random.binomial(n-k0-k1,p2/(1-p0-p1))
    if k0+k1+k2==n:return k0,k1,k2,0
    k3=np.random.binomial(n-k0-k1-k2,p3/(1-p0-p1-p2))
    return k0,k1,k2,k3


@nb.njit(cache=True)
def _leap(state,flux,hops,a0,b0,omega,diff,mask,requested):
    n=state.shape[1];maxrate=0.
    for i in range(n):
        x,b=state[1,i],state[4,i]
        r4=mask[4]*x*(x-1)/(omega*omega) if x>=2 else 0.
        rate3=mask[3]*29.*x*b/b0[i]
        joint=rate3/min(x,b) if min(x,b)>0 else 0.
        maxrate=max(maxrate,mask[0]*omega*14./a0[i],mask[1]+2*diff[0],
                    mask[2]*.1+mask[5]*.9+mask[3]*29.*b/b0[i]+2*diff[1],r4+2*diff[2],joint)
    dt=min(requested,.2/maxrate) if maxrate>0 else requested
    delta=np.zeros_like(state);largest=0.
    for i in range(n):
        left,right=(i-1)%n,(i+1)%n
        g,x,y,a,b=state[0,i],state[1,i],state[2,i],state[3,i],state[4,i]
        pa=mask[0]*omega*14.*dt/a0[i]
        source=np.random.binomial(a,pa)
        gx,gl,gr,_=_exits4(g,mask[1]*dt,diff[0]*dt,diff[0]*dt,0.)
        rate3=mask[3]*29.*x*b/b0[i]
        m=min(x,b);joint=rate3*dt/m if m>0 else 0.
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
        delta[3,i]-=source;delta[4,i]-=conversion
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
def leap_chunk(initial,a0,b0,omega,diff,mask,tau,observations,reps,seed):
    np.random.seed(seed)
    n=initial.shape[1];out=np.zeros((reps,len(observations),20,n),dtype=np.int64);stats=np.zeros((reps,3))
    for r in range(reps):
        state=initial.copy();flux=np.zeros((7,n),dtype=np.int64);hops=np.zeros((6,n),dtype=np.int64)
        t=0.;leaps=0;adapted=0;maximum=0.
        for k,target in enumerate(observations):
            while t<target:
                requested=min(tau,target-t)
                dt,prob=_leap(state,flux,hops,a0,b0,omega,diff,mask,requested)
                if dt<=0:raise ValueError("Nonpositive step")
                adapted+=int(dt<requested*(1-1e-12));leaps+=1;maximum=max(maximum,prob)
                if dt==target-t:t=target
                else:t+=dt
            out[r,k,:7]=state;out[r,k,7:14]=flux;out[r,k,14:20]=hops
        stats[r,0]=leaps;stats[r,1]=adapted;stats[r,2]=maximum
    return out,stats


@nb.njit(cache=True)
def _ssa_clock(pop,a0,b0,omega,diff,mask):
    # Separate direct propensity implementation; no leap helpers or factorization.
    n=pop.shape[1];rates=np.zeros(13*n)
    for node in range(n):
        j=13*node;G,X,Y,A,B=pop[0,node],pop[1,node],pop[2,node],pop[3,node],pop[4,node]
        rates[j]=mask[0]*14*omega*A/a0[node]
        rates[j+1]=mask[1]*G
        rates[j+2]=mask[2]*X/10
        rates[j+3]=mask[3]*29*X*B/b0[node]
        rates[j+4]=mask[4]*X*(X-1)*Y/(omega*omega) if X>=2 else 0.
        rates[j+5]=mask[5]*9*X/10
        for species in range(3):
            rates[j+7+2*species]=diff[species]*pop[species,node]
            rates[j+8+2*species]=diff[species]*pop[species,node]
    total=rates.sum()
    if total==0:return np.inf,-1
    waiting=np.random.exponential(1/total)
    mark=np.random.random()*total;accumulated=0.
    for event in range(len(rates)):
        accumulated+=rates[event]
        if mark<accumulated:return waiting,event
    raise ValueError("SSA event selection failed")


@nb.njit(cache=True)
def ssa_chunk(initial,a0,b0,omega,diff,mask,observations,reps,seed):
    np.random.seed(seed)
    n=initial.shape[1];out=np.zeros((reps,len(observations),20,n),dtype=np.int64);stats=np.zeros((reps,3))
    for r in range(reps):
        pop=initial.copy();reactions=np.zeros((7,n),dtype=np.int64);hops=np.zeros((6,n),dtype=np.int64)
        future,event=_ssa_clock(pop,a0,b0,omega,diff,mask);events=0
        for k,target in enumerate(observations):
            while future<=target:
                node=event//13;channel=event%13
                if channel==0:pop[3,node]-=1;pop[0,node]+=1
                elif channel==1:pop[0,node]-=1;pop[1,node]+=1
                elif channel==2:pop[1,node]-=1;pop[0,node]+=1
                elif channel==3:pop[4,node]-=1;pop[1,node]-=1;pop[2,node]+=1;pop[5,node]+=1
                elif channel==4:pop[2,node]-=1;pop[1,node]+=1
                elif channel==5:pop[1,node]-=1;pop[6,node]+=1
                elif channel==6:raise ValueError("R1re must be zero at P17")
                else:
                    species=(channel-7)//2;direction=(channel-7)%2
                    other=(node-1)%n if direction==0 else (node+1)%n
                    pop[species,node]-=1;pop[species,other]+=1
                    hops[2*species+direction,node]+=1
                if channel<7:reactions[channel,node]+=1
                if np.any(pop<0):raise ValueError("Negative SSA state")
                waiting,event=_ssa_clock(pop,a0,b0,omega,diff,mask)
                future+=waiting;events+=1
            out[r,k,:7]=pop;out[r,k,7:14]=reactions;out[r,k,14:20]=hops
        stats[r,0]=events
    return out,stats


def simulate(case,method,reps,seed,mask=None):
    mask=np.ones(7) if mask is None else np.asarray(mask,dtype=float)
    if method not in (0,1,2,3) or reps<=0:raise ValueError("Invalid method/replicate count")
    if np.any(case["state"]<0) or np.any(case["a0"]<=0) or np.any(case["b0"]<=0):raise ValueError("Invalid population/reference stock")
    if case["omega"]<=0 or np.any(case["diff"]<0) or np.any(mask<0):raise ValueError("Invalid rate parameters")
    if method==0:
        return ssa_chunk(case["state"],case["a0"],case["b0"],case["omega"],case["diff"],mask,OBS,reps,seed)
    return leap_chunk(case["state"],case["a0"],case["b0"],case["omega"],case["diff"],mask,TAUS[method],OBS,reps,seed)


def ledger(data,initial):
    actual=data[:,:,:7]-initial[None,None]
    expected=np.einsum("sr,btrn->btsn",STOICH,data[:,:,7:14])
    for s in range(3):
        left,right=data[:,:,14+2*s],data[:,:,15+2*s]
        expected[:,:,s]+=np.roll(left,-1,axis=-1)+np.roll(right,1,axis=-1)-left-right
    error=int(np.max(np.abs(actual-expected)))
    total=int(np.max(np.abs(data[:,:,:7].sum(axis=(2,3))-initial.sum())))
    return dict(max_integer_residual=error,max_total_count_error=total,nonnegative=bool(np.all(data>=0)),pass_check=bool(error==0 and total==0 and np.all(data>=0)))


def observables(data):
    return np.stack([data[:,:,1].sum(axis=-1),data[:,:,4].sum(axis=-1),data[:,:,6].sum(axis=-1),
                     data[:,:,7].sum(axis=-1),data[:,:,10].sum(axis=-1),data[:,:,11].sum(axis=-1),
                     data[:,:,1,0],data[:,:,2,0]],axis=-1)

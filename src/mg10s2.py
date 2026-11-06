"""MG049 actual-spacing/density fixtures and independent snapshot-law oracle.

The S1 leap/SSA sources are unchanged. Oracle sampling uses NumPy multinomials
and binomial A rather than the S1 joint/conditional exit helper.
"""
from pathlib import Path
import hashlib
import numpy as np
import mg10s1 as s1

ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/"results/MG10S2"
TAUS=(1/6000,1/12000,1/24000)
OBS=np.array([1/8000,1/4000,1/2000])
REPS=16384
S0=np.array([140/9,140/9,29/(140/9)])
DIFF=np.array([25.,25.,300.])


def fixture(case):
    full=case<4
    if case not in range(6):raise ValueError(case)
    omega=(100.,1000.)[case//2] if full else (100.,1000.)[case-4]
    n=500 if full else 16
    initial=np.zeros((7,n),dtype=np.int64)
    concentrations=np.repeat(S0[:,None],n,axis=1)
    if full and case%2:
        with np.load(ROOT/"results/MG003/ref_profiles.npz") as raw:concentrations+=raw["h0.2_b29.0"][1:]
    initial[:3]=np.rint(omega*concentrations).astype(np.int64)
    initial[3]=round(omega*14*300);initial[4]=round(omega*29*S0[1]*300)
    return dict(initial=initial,a0=initial[3].astype(float),b0=initial[4].astype(float),omega=omega,diff=DIFF.copy(),
                modes=np.array([32,250] if full else [1,8],dtype=int),start="S" if full and case%2 else "H")


def conditional_mean(case,tau):
    n=case["initial"][:3].astype(float);g,x,y=n
    nonlinear=x*(x-1)*y/case["omega"]**2
    lap=np.roll(n,1,axis=1)+np.roll(n,-1,axis=1)-2*n
    drift=np.stack([case["omega"]*14-g+.1*x,g-30*x+nonlinear,29*x-nonlinear])
    return n+tau*(drift+case["diff"][:,None]*lap)


def reconstruct(initial,events):
    event=events.astype(np.int64)
    f1,f2,fr,f3,f4,f5,f1r=[event[:,i] for i in range(7)]
    populations=np.stack([initial[0]+f1-f2+fr-f1r,initial[1]+f2-fr-f3+f4-f5,
                          initial[2]+f3-f4,initial[3]-f1+f1r,initial[4]-f3,initial[5]+f3,initial[6]+f5],axis=1)
    for species in range(3):
        left,right=event[:,7+2*species],event[:,8+2*species]
        populations[:,species]+=np.roll(left,-1,axis=-1)+np.roll(right,1,axis=-1)-left-right
    return populations


def population_hash(populations):
    return hashlib.sha256(np.asarray(populations,dtype="<i8").tobytes(order="C")).hexdigest()


def snapshot(case,step,sampler,reps,seed):
    tau=TAUS[step];initial=case["initial"];n=initial.shape[1]
    if sampler==0:
        data,stats=s1.leap_chunk(initial,case["a0"],case["b0"],case["omega"],case["diff"],np.ones(7),tau,np.array([tau]),reps,seed)
        if not s1.ledger(data,initial)["pass_check"]:raise ValueError("Snapshot ledger")
        if not np.all(stats[:,0]==1) or np.any(stats[:,1]!=0):raise ValueError("Unexpected snapshot safety adaptation")
        events=np.ascontiguousarray(data[:,0,7:20]);actual=np.ascontiguousarray(data[:,0,:7])
    else:
        rng=np.random.default_rng(seed)
        g,x,y=initial[:3]
        if not np.all(initial[4]>=x):raise ValueError("Oracle requires B surplus")
        pg=np.tile([tau,25*tau,25*tau,1-51*tau],(n,1))
        px=np.tile([.1*tau,29*tau,.9*tau,25*tau,25*tau,1-80*tau],(n,1))
        py=np.stack([x*(x-1)/case["omega"]**2*tau,np.full(n,300*tau),np.full(n,300*tau),
                     1-(x*(x-1)/case["omega"]**2+600)*tau],axis=1)
        if min(pg.min(),px.min(),py.min())<0:raise ValueError("Oracle probability")
        G=rng.multinomial(g,pg,size=(reps,n))
        X=rng.multinomial(x,px,size=(reps,n))
        Y=rng.multinomial(y,py,size=(reps,n))
        A=rng.binomial(initial[3],case["omega"]*14*tau/case["a0"],size=(reps,n))
        events=np.zeros((reps,13,n),dtype=np.int64)
        events[:,0]=A;events[:,1]=G[:,:,0];events[:,2]=X[:,:,0];events[:,3]=X[:,:,1]
        events[:,4]=Y[:,:,0];events[:,5]=X[:,:,2]
        events[:,7:9]=G[:,:,1:3].transpose(0,2,1)
        events[:,9:11]=X[:,:,3:5].transpose(0,2,1)
        events[:,11:13]=Y[:,:,1:3].transpose(0,2,1)
        actual=reconstruct(initial,events)
        stats=np.zeros((reps,3));stats[:,0]=1
    if not np.array_equal(reconstruct(initial,events),actual) or np.any(actual<0):raise ValueError("Reconstruction")
    if np.any(actual.sum(axis=(1,2))!=initial.sum()):raise ValueError("Global count conservation")
    if np.min(events)<0 or np.max(events)>np.iinfo(np.uint32).max:raise ValueError("Packing range")
    return events.astype(np.uint32),actual,stats


def short(case,method,reps,seed):
    if method==0:
        data,stats=s1.ssa_chunk(case["initial"],case["a0"],case["b0"],case["omega"],case["diff"],np.ones(7),OBS,reps,seed)
    else:
        data,stats=s1.leap_chunk(case["initial"],case["a0"],case["b0"],case["omega"],case["diff"],np.ones(7),TAUS[method-1],OBS,reps,seed)
    if not s1.ledger(data,case["initial"])["pass_check"]:raise ValueError("Short ledger")
    return data,stats


def covariance(case,tau,modes,instantaneous=False):
    """Exact same-mode snapshot covariance or CTMC instantaneous noise matrix."""
    initial=case["initial"];n=initial.shape[1];omega=case["omega"]
    eye=np.eye(3);result=[]
    for mode in modes:
        theta=2*np.pi*mode/n;left,right=np.exp(1j*theta)-1,np.exp(-1j*theta)-1
        total=np.zeros((3,3),dtype=complex)
        for i in range(n):
            g,x,y,A=initial[:4,i]
            pa=omega*14*tau/case["a0"][i]
            total+=A*pa*(1 if instantaneous else 1-pa)*np.outer(eye[0],eye[0])
            groups=[(g,np.array([tau,25*tau,25*tau]),[-eye[0]+eye[1],eye[0]*left,eye[0]*right]),
                    (x,np.array([.1,29,.9,25,25])*tau,[eye[0]-eye[1],-eye[1]+eye[2],-eye[1],eye[1]*left,eye[1]*right]),
                    (y,np.array([x*(x-1)/omega**2,300,300])*tau,[eye[1]-eye[2],eye[2]*left,eye[2]*right])]
            for count,probabilities,vectors in groups:
                vectors=np.array(vectors);mean=probabilities@vectors
                matrix=sum(p*np.outer(v,v.conj()) for p,v in zip(probabilities,vectors))
                if not instantaneous:matrix-=np.outer(mean,mean.conj())
                total+=count*matrix
        result.append(total/(n*omega*tau))
    return np.array(result)

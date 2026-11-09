"""MG100 conditional full-lattice pilot; frozen G4 law, streamed clock/counters."""
import os
for key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[key]='1'
from pathlib import Path
import time
import numpy as np
import numba as nb
from scipy.integrate import solve_ivp
from scipy.stats import t as student
from mg11g4 import step,arrays,validate,MAX_STEPS
from mg11g3 import factory as parent_factory,donor_path,rhs,pattern
from mg11g2 import initial
from mg11g2_pattern_v2 import metrics
from runner import seed_for

ROOT=Path(__file__).resolve().parents[1];OBS=np.arange(61)/4;TAUS=(1/24000,1/48000);OMEGA=1000

def contexts():return [dict(origin=o,start=s,arm=a) for o in range(2) for s in range(2) for a in range(2)]
def parent_job(j):return dict(id=f'o{j["origin"]}_s{j["start"]}_v1_a{j["arm"]}',origin=j['origin'],start=j['start'],setting=1,arm=j['arm'])
def factory(j):
    case,ga,gb,k=parent_factory(parent_job(j));case=dict(case,initial=np.rint(OMEGA*case['initial']).astype(np.int64),a0=np.rint(OMEGA*case['a0']).astype(np.int64),b0=np.rint(OMEGA*case['b0']).astype(np.int64),omega=float(OMEGA));validate(case,ga,gb,k,np.ones(7));return case,ga,gb,k
def jobs():
    controls=[dict(c,id=f'd_o{c["origin"]}_s{c["start"]}_a{c["arm"]}_f{f}',kind='deterministic',finite=f) for c in contexts() for f in range(2)]
    branches=[]
    for v in range(2):
        for b in range(4):
            for c in contexts():
                cell=4*c['origin']+2*c['start']+c['arm'];branches.append(dict(c,id=f'p_v{v}_b{b}_o{c["origin"]}_s{c["start"]}_a{c["arm"]}',kind='stochastic',setting=v,block=b,seed=seed_for(11,cell,v,10003,b),warmup_seed=seed_for(11,cell,v,10003,b,1)))
    return controls+branches

@nb.njit(cache=True)
def segment(pop,rx,hop,ef,a0,b0,omega,diff,mask,tau,t0,targets,seed,src,dst,w,k,degree,offsets,ids,stats):
    if seed>=0:np.random.seed(seed)
    n=pop.shape[1];out=np.empty((len(targets),20,n),np.int64);env=np.empty((len(targets),2,src.shape[1]),np.int64);t=t0
    for j,target in enumerate(targets):
        while t<target:
            if stats[0]>=MAX_STEPS:raise ValueError('Per-path leap cap')
            requested=min(tau,target-t);dt,bound=step(pop,rx,hop,ef,a0,b0,omega,diff,mask,requested,src,dst,w,k,degree,offsets,ids)
            if dt<=0:raise ValueError('Nonpositive numerical step')
            t=target if dt==target-t else t+dt;stats[0]+=1;stats[1]+=int(dt<requested*(1-1e-12));stats[2]=max(stats[2],bound)
        out[j,:7]=pop;out[j,7:14]=rx;out[j,14:]=hop;env[j]=ef
    return out,env

def streamed(case,ga,gb,k,times,seed,tau,deadline,checkpoint=None,chunk=4):
    times=np.asarray(times,dtype=float)
    if times.ndim!=1 or len(times)<2 or times[0]!=0 or np.any(~np.isfinite(times)) or np.any(np.diff(times)<=0) or not isinstance(seed,(int,np.integer)) or not 0<=seed<2**32 or not np.isfinite(tau) or tau<=0 or not isinstance(chunk,int) or chunk<1:raise ValueError('Invalid stream bounds')
    validate(case,ga,gb,k,np.ones(7));src,dst,w,degree,offsets,ids=arrays(ga,gb,k);n=ga['n'];pop=case['initial'].copy();rx=np.zeros((7,n),np.int64);hop=np.zeros((6,n),np.int64);ef=np.zeros((2,len(ga['src'])),np.int64);stats=np.zeros(3)
    data=np.zeros((1,len(times),20,n),np.int64);env=np.zeros((1,len(times),2,len(ga['src'])),np.int64);data[0,0,:7]=pop;compute=0.;io=0.
    def save(right):
        nonlocal io
        if checkpoint:
            clock=time.perf_counter();checkpoint(times[:right],(data[:,:right],env[:,:right],stats[None].copy()));io+=time.perf_counter()-clock
    save(1)
    for left in range(1,len(times),chunk):
        if time.perf_counter()>deadline:raise TimeoutError('G5 stream deadline')
        right=min(len(times),left+chunk);clock=time.perf_counter();z,e=segment(pop,rx,hop,ef,case['a0'],case['b0'],case['omega'],case['diff'],np.ones(7),tau,times[left-1],times[left:right],seed if left==1 else -1,src,dst,w,k,degree,offsets,ids,stats);compute+=time.perf_counter()-clock;data[0,left:right]=z;env[0,left:right]=e;save(right)
    return (data,env,stats[None]),dict(compute_seconds=compute,io_seconds=io)

def concentration_case(case):return dict(case,initial=case['initial']/OMEGA,a0=case['a0']/OMEGA,b0=case['b0']/OMEGA,omega=1.)
def drift(t,y,case,ga,gb,k,finite):
    out=rhs(t,y,case,ga,gb,k)
    if finite:
        n=ga['n'];x=y[n:2*n];yy=y[2*n:3*n];delta=x*(np.maximum(x-1/OMEGA,0)-x)*yy;out[n:2*n]+=delta;out[2*n:3*n]-=delta;out[11*n:12*n]+=delta
    return out
def deterministic(j,times,deadline):
    integer,ga,gb,k=factory(j);case=concentration_case(integer)
    def fun(t,y):
        if time.perf_counter()>deadline:raise TimeoutError('G5 deterministic deadline')
        return drift(t,y,case,ga,gb,k,j['finite'])
    sol=solve_ivp(fun,(0,float(times[-1])),initial(case,ga),t_eval=times,method='BDF',jac_sparsity=pattern(ga,gb),rtol=1e-11,atol=1e-13,max_step=.05)
    if not sol.success or sol.y.shape[1]!=len(times):raise RuntimeError('G5 deterministic incomplete')
    return sol.y.T,case,ga,gb,k

def projection(pop,times):return np.asarray(metrics(pop,times,full=False)['projection'])
def late_value(pop,times):
    times=np.asarray(times);keep=times>=10
    if not np.array_equal(times,OBS):raise ValueError('Incomplete primary observations')
    return float(np.trapezoid(projection(pop,times)[keep,1],times[keep])/5)
def interval(values,confidence=.975):
    x=np.asarray(values,float)
    if x.shape!=(4,) or np.any(~np.isfinite(x)):raise ValueError('Four balanced blocks required')
    mean=float(x.mean());sd=float(x.std(ddof=1));half=float(student.ppf((1+confidence)/2,3)*sd/2);return dict(n=4,mean=mean,sd=sd,lower=mean-half,upper=mean+half)
def summarize(records):
    samples=[];out=[]
    for v in range(2):
        vals=[]
        for b in range(4):
            keys=[(v,b,c['origin'],c['start'],c['arm']) for c in contexts()]
            if all(key in records for key in keys):vals.append(sum((1 if s==1 else -1)*(1 if a==1 else -1)*records[(v,b,o,s,a)] for _,_,o,s,a in keys)/2)
        samples.append(vals);out.append(interval(vals) if len(vals)==4 else dict(n=len(vals),values=vals,inference=None))
    agreement=None
    if all(len(x)==4 for x in samples):
        x,y=map(np.asarray,samples);a=x.var(ddof=1)/4;b=y.var(ddof=1)/4;se=np.sqrt(a+b);df=(a+b)**2/(a*a/3+b*b/3) if se else 6.;d=float(y.mean()-x.mean());half=float(student.ppf(.95,df)*se);agreement=dict(mean=d,lower=d-half,upper=d+half,df=float(df),within_margin=bool(d-half>=-.02 and d+half<=.02))
    return dict(settings=out,refinement=agreement,direction_consistent=bool(agreement and agreement['within_margin'] and all(v['lower']>0 for v in out)))

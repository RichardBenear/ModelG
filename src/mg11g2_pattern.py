"""MG088 deterministic pattern-response design and sparse full-ledger adapter.

Frozen N1 continuum law and Step1 flux layout retained. No stochastic sampler.
"""
import os
for key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[key]='1'
from pathlib import Path
import time
import numpy as np
from scipy.integrate import BDF
from scipy.sparse import coo_matrix
from mg11n1 import graph
from mg11g2 import rhs,initial,unpack,sizes
from mg10p1 import centered_smooth,continuity
ROOT=Path(__file__).resolve().parents[1]
N=500;S0=np.array([140/9,140/9,29/(140/9)])
OBS=np.arange(101,dtype=float)/4
CONFIGS=('MAINTAINED','ISOLATED_C1','LOCAL_C1','CROSS_RING_C1','ISOLATED_C100','LOCAL_C100','CROSS_RING_C100')
SETTINGS=(dict(rtol=1e-8,atol=1e-10,max_step=.1),dict(rtol=1e-10,atol=1e-12,max_step=.05))

def reference():
    with np.load(ROOT/'results/MG003/ref_profiles.npz') as raw:phi=raw['h0.2_b29.0'][1:].copy()
    return phi,np.array([centered_smooth(f) for f in phi])

def environment(kind,n=N):
    if kind=='NONE':return graph(n,[])
    if n<8 or n%4:raise ValueError('Environmental test ring requires N divisible by4 and >=8')
    edges=[(i,(i+1)%n,1.) for i in range(n)]
    if kind=='LOCAL':edges +=[(i,(i+3)%n,1.) for i in range(0,n,2)]
    elif kind=='CROSS_RING':edges +=[(i,i+n//2,1.) for i in range(n//2)]
    else:raise ValueError('Unknown environment')
    return graph(n,edges)

def factory(config,start,n=N):
    if config not in CONFIGS or start not in (0,1):raise ValueError('Unknown condition/start')
    if start and n!=N:raise ValueError('Seed requires unchanged500node reference')
    capacity=100 if config.endswith('C100') else 1
    kind='LOCAL' if config.startswith('LOCAL') else 'CROSS_RING' if config.startswith('CROSS_RING') else 'NONE'
    g=environment(kind,n);a0=np.full(n,4200.*capacity);b0=np.full(n,29*(140/9)*300*capacity);pop=np.zeros((7,n));pop[:3]=S0[:,None]
    if start:pop[:3]+=reference()[0]
    pop[3]=a0;pop[4]=b0
    case=dict(initial=pop,a0=a0,b0=b0,omega=1.,diff=np.array([25.,25.,300.]))
    return case,g,np.full(2,1/3) if kind!='NONE' else np.zeros(2),'CHEMOSTAT' if config=='MAINTAINED' else 'FINITE'

def validate(case,g,k,mode):
    n=g['n'];pop=case['initial']
    if pop.shape!=(7,n) or not np.all(np.isfinite(pop)) or np.any(pop<0):raise ValueError('Invalid deterministic populations')
    if case['omega']!=1. or mode not in ('FINITE','CHEMOSTAT'):raise ValueError('Concentration adapter only')
    for key in ('a0','b0'):
        if case[key].shape!=(n,) or np.any(~np.isfinite(case[key])) or np.any(case[key]<=0):raise ValueError('Invalid references')
    if case['diff'].shape!=(3,) or np.any(~np.isfinite(case['diff'])) or np.any(case['diff']<0) or np.asarray(k).shape!=(2,) or np.any(~np.isfinite(k)) or np.any(np.asarray(k)<0):raise ValueError('Invalid transport')
    if mode=='CHEMOSTAT' and (not np.array_equal(pop[3],case['a0']) or not np.array_equal(pop[4],case['b0']) or np.ptp(case['a0']) or np.ptp(case['b0'])):raise ValueError('Uniform maintained stocks required')

def sparsity(g):
    """Conservative structural pattern; columns only physical populations.

    Extra local dependencies are allowed; missing dependencies are tested by
    independent finite differences before use. Ledger rows have no feedback.
    """
    n=g['n'];e=len(g['src']);rows=[];cols=[]
    def add(row,col):rows.append(row);cols.append(col)
    for i in range(n):
        for s in range(7):
            for c in range(5):add(s*n+i,c*n+i);add(7*n+s*n+i,c*n+i)
        for s in range(3):
            for j in ((i-1)%n,(i+1)%n):add(s*n+i,s*n+j)
            for d in range(2):add(14*n+(2*s+d)*n+i,s*n+i)
        for s in range(2):
            for c in range(5):add(20*n+2*e+s*n+i,c*n+i)
    for s in range(2):
        for edge,(i,j) in enumerate(zip(g['src'],g['dst'])):
            add((3+s)*n+i,(3+s)*n+i);add((3+s)*n+j,(3+s)*n+i);add(20*n+s*e+edge,(3+s)*n+i)
            # Maintained stocks would need to cancel incoming/outgoing flow.
            add(20*n+2*e+s*n+i,(3+s)*n+i);add(20*n+2*e+s*n+j,(3+s)*n+i)
    return coo_matrix((np.ones(len(rows),dtype=bool),(rows,cols)),shape=(sum(sizes(g)),sum(sizes(g)))).tocsr()

def solve_pattern(case,g,k,mode,times,setting,deadline,checkpoint=None):
    validate(case,g,k,mode);times=np.asarray(times,dtype=float)
    if times.ndim!=1 or len(times)<2 or times[0]!=0 or np.any(~np.isfinite(times)) or np.any(np.diff(times)<=0):raise ValueError('Invalid observations')
    controls=SETTINGS[setting];start=time.perf_counter()
    def bounded(t,y):
        if time.perf_counter()>deadline:raise TimeoutError('Per-job development/production deadline')
        return rhs(t,y,case,g,k,mode)
    homogeneous=bool(np.all(case['initial']==case['initial'][:,:1]))
    if homogeneous:
        # Exact invariant subspace, not artificial stabilization of perturbed data.
        # Two additional integrals generate each directed environmental flux.
        small=dict(case,initial=case['initial'][:,:1],a0=case['a0'][:1],b0=case['b0'][:1]);empty=graph(1,[])
        def fun(t,y):
            if time.perf_counter()>deadline:raise TimeoutError('Per-job deadline')
            return np.r_[rhs(t,y[:22],small,empty,np.zeros(2),mode),y[3:5]]
        y0=np.r_[initial(small,empty),[0.,0.]];jac=None
        def expand(y):
            pop,rx,hops,env,imports=unpack(y[:22],empty);return np.r_[np.repeat(pop,g['n'],axis=1).ravel(),np.repeat(rx,g['n'],axis=1).ravel(),np.repeat(hops,g['n'],axis=1).ravel(),(k[:,None]*g['weight'][None,:]*y[22:24,None]).ravel(),np.repeat(imports,g['n'],axis=1).ravel()]
    else:fun=bounded;y0=initial(case,g);jac=sparsity(g);expand=lambda y:y.copy()
    solver=BDF(fun,0.,y0,float(times[-1]),jac_sparsity=jac,**controls);saved=[expand(y0)];next_index=1
    while solver.status=='running':
        solver.step()
        if solver.status=='failed':raise RuntimeError('BDF failed')
        if next_index<len(times) and times[next_index]<=solver.t:
            dense=solver.dense_output()
            while next_index<len(times) and times[next_index]<=solver.t:
                saved.append(expand(dense(times[next_index])));next_index+=1
                if checkpoint:checkpoint(times[:next_index],np.asarray(saved))
    if len(saved)!=len(times) or not np.all(np.isfinite(saved)):raise RuntimeError('Incomplete/nonfinite record')
    return dict(times=times,state=np.asarray(saved),elapsed_seconds=time.perf_counter()-start,homogeneous_reduction=homogeneous,nfev=solver.nfev,njev=solver.njev)

def metrics(pop,times,full=True):
    fields=pop[:,:3];template=reference()[1];smoothed=np.array([[centered_smooth(v) for v in frame] for frame in fields]);projection=np.sum(smoothed*template[None],axis=-1)/np.sum(template*template,axis=-1)[None]
    residual=smoothed-projection[:,:,None]*template[None];out=dict(projection=projection.tolist(),field_mean=fields.mean(axis=-1).tolist(),residual_rms=np.sqrt(np.mean(residual*residual,axis=-1)).tolist(),contrast_rms=np.sqrt(np.mean(smoothed*smoothed,axis=-1)).tolist())
    if full:
        if not np.array_equal(times,OBS):raise ValueError('Full observation grid required')
        out['late_projection']=(np.trapezoid(projection[80:],x=times[80:],axis=0)/5).tolist();out['reference_continuity']=continuity(fields[:,1],times,x0=0.)
    return out

def jobs():return [dict(id=f'c{c}_s{s}_v{v}',config=config,condition=c,start=s,setting=v) for c,config in enumerate(CONFIGS) for s in (0,1) for v in (0,1)]

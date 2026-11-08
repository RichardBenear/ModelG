"""MG087 Step1 deterministic coupling/ledgers, reusing frozen N1 continuum law.

No new stochastic kernel. Concentration mean field is NOT finite-count SSA drift.
Ledger state: populations7N, reactions7N, directed active hops6N,
directed environmental hops2E, external imports2N. Equal-volume symmetric graph.
"""
import time
import numpy as np
from scipy.integrate import solve_ivp
from mg11n1 import continuum_rhs,validate

def sizes(g):
    n=g['n'];e=len(g['src']);return (7*n,7*n,6*n,2*e,2*n)

def unpack(y,g):
    n=g['n'];e=len(g['src']);parts=np.split(np.asarray(y),np.cumsum(sizes(g))[:-1])
    return tuple(a.reshape(shape) for a,shape in zip(parts,((7,n),(7,n),(6,n),(2,e),(2,n))))

def initial(case,g):
    y=np.zeros(sum(sizes(g)));y[:7*g['n']]=case['initial'].ravel()/case['omega'];return y

def configuration(case,g,kappa,mode):
    if mode not in ('FINITE','CHEMOSTAT'):raise ValueError('Unknown coupling mode')
    validate(case,g,np.asarray(kappa,dtype=float),np.ones(7))
    # Uniform reference stocks avoid signed local withdrawals needed to maintain
    # heterogeneous clamped baths against transport; no such boundary is assumed.
    if mode=='CHEMOSTAT':
        for s,key in ((3,'a0'),(4,'b0')):
            if not np.all(case[key]==case[key][0]) or not np.array_equal(case['initial'][s],case[key]):raise ValueError('Chemostat requires uniform fully available stocks')

def rhs(t,y,case,g,kappa,mode='FINITE'):
    pop,_,_,_,_=unpack(y,g);n=g['n'];G,X,Y,A,B=pop[:5]
    flux=np.array([14*A/(case['a0']/case['omega']),G,.1*X,29*X*B/(case['b0']/case['omega']),X*X*Y,.9*X,np.zeros(n)])
    active=np.repeat(case['diff'],2)[:,None]*np.repeat(pop[:3],2,axis=0)
    environment=np.asarray(kappa)[:,None]*g['weight'][None,:]*pop[3:5,g['src']]
    drift=continuum_rhs(t,pop.ravel(),case,g,kappa).reshape(7,n);external=np.zeros((2,n))
    if mode=='CHEMOSTAT':external=-drift[3:5].copy();drift[3:5]=0
    return np.concatenate((drift.ravel(),flux.ravel(),active.ravel(),environment.ravel(),external.ravel()))

def reconstruct(case,g,y):
    pop,rx,active,env,external=unpack(y,g);out=case['initial']/case['omega']
    # Explicit reaction increments rather than the primary continuum function.
    f0,f1,f2,f3,f4,f5,f6=rx
    out=out+np.array([f0-f1+f2-f6,f1-f2-f3+f4-f5,f3-f4,-f0+f6,-f3,f3,f5])
    for species in range(3):
        left,right=active[2*species:2*species+2];out[species]+=np.roll(left,-1)+np.roll(right,1)-left-right
    for species in range(2):
        np.add.at(out[3+species],g['src'],-env[species]);np.add.at(out[3+species],g['dst'],env[species]);out[3+species]+=external[species]
    return out

def solve(case,g,kappa,times,mode='FINITE',rtol=1e-9,atol=1e-11,deadline=None):
    kappa=np.asarray(kappa,dtype=float);configuration(case,g,kappa,mode);times=np.asarray(times,dtype=float)
    if times.ndim!=1 or len(times)<2 or times[0]!=0 or np.any(~np.isfinite(times)) or np.any(np.diff(times)<=0):raise ValueError('Invalid observation grid')
    def bounded(t,y):
        if deadline is not None and time.perf_counter()>deadline:raise TimeoutError('Qualification elapsed cap')
        return rhs(t,y,case,g,kappa,mode)
    sol=solve_ivp(bounded,(0,float(times[-1])),initial(case,g),t_eval=times,method='BDF',rtol=rtol,atol=atol,max_step=.05)
    if not sol.success or sol.y.shape[1]!=len(times) or not np.all(np.isfinite(sol.y)):raise RuntimeError('Incomplete/nonfinite coupling integration')
    # Never clip numerical negatives or normalize away a ledger discrepancy.
    return dict(times=sol.t,state=sol.y.T,mode=mode,nfev=sol.nfev)

def well_mixed_rhs(t,y,case):
    """Shared finite bath after fast mixing, not a chemostat.

    Local active species stay on the same ring. A/B per-node concentration equals
    its spatial mean. Reaction demand is aggregated over all contacts.
    """
    n=case['initial'].shape[1];z=y.reshape(7,n);G,X,Y,A,B=z[:5]
    a=A.mean();b=B.mean();f0=14*a/(case['a0']/case['omega']);f3=29*X*b/(case['b0']/case['omega']);f4=X*X*Y
    out=np.array([f0-G+.1*X,G-X-f3+f4,f3-f4,np.full(n,-f0.mean()),np.full(n,-f3.mean()),f3,.9*X])
    for species in range(3):out[species]+=case['diff'][species]*(np.roll(z[species],1)+np.roll(z[species],-1)-2*z[species])
    return out.ravel()

def mixed_initial(case):
    z=case['initial']/case['omega'];z[3:5]=z[3:5].mean(axis=1,keepdims=True);return z

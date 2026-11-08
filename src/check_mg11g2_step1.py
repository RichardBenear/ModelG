"""MG087 separate population ODE, ledger incidence and qualification audit."""
import os
for key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[key]='1'
from pathlib import Path
import time,json,hashlib
import numpy as np
from scipy.integrate import solve_ivp
from mg11g2_cases import fixture,TIMES
from mg10d1_pool import atomic_json
ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/'results/MG11G2/step1'

def reference(t,flat,case,g,kappa,mode,mixed=False):
    n=g['n'];z=flat.reshape(7,n);out=np.zeros((7,n));a0=case['a0']/case['omega'];b0=case['b0']/case['omega']
    # Individual node equations and explicit directed incidence, no primary RHS.
    for i in range(n):
        gg,x,y,a,b,_,_=z[:,i];a=z[3].mean() if mixed else a;b=z[4].mean() if mixed else b
        source=14*a/a0[i];convert=29*b*x/b0[i];feedback=x*x*y
        out[:,i]=[source-gg+.1*x,gg-.1*x-convert+feedback-.9*x,convert-feedback,-source,-convert,convert,.9*x]
        for s in range(3):
            for j in ((i-1)%n,(i+1)%n):out[s,i]+=case['diff'][s]*(z[s,j]-z[s,i])
    if mixed:
        out[3]=out[3].mean();out[4]=out[4].mean()
    else:
        for s in range(2):
            for i,j,w in zip(g['src'],g['dst'],g['weight']):flow=kappa[s]*w*z[3+s,i];out[3+s,i]-=flow;out[3+s,j]+=flow
    if mode=='CHEMOSTAT':out[3:5]=0
    return out.ravel()

def integrate(case,g,kappa,mode,deadline,mixed=False):
    z=case['initial']/case['omega']
    if mixed:z[3:5]=np.broadcast_to(z[3:5].mean(axis=1,keepdims=True),z[3:5].shape)
    def bounded(t,y):
        if time.perf_counter()>deadline:raise TimeoutError('Qualification elapsed cap')
        return reference(t,y,case,g,kappa,mode,mixed)
    sol=solve_ivp(bounded,(0,1),z.ravel(),method='Radau',t_eval=TIMES,rtol=1e-11,atol=1e-12,max_step=.025)
    assert sol.success and sol.y.shape==(7*g['n'],4) and np.all(np.isfinite(sol.y));return sol.y.T.reshape(4,7,g['n'])

def incidence(case,g,y):
    n=g['n'];e=len(g['src']);pop=y[:7*n].reshape(7,n);rx=y[7*n:14*n].reshape(7,n);active=y[14*n:20*n].reshape(6,n);env=y[20*n:20*n+2*e].reshape(2,e);imports=y[-2*n:].reshape(2,n)
    matrix=np.array([[1,-1,1,0,0,0,-1],[0,1,-1,-1,1,-1,0],[0,0,0,1,-1,0,0],[-1,0,0,0,0,0,1],[0,0,0,-1,0,0,0],[0,0,0,1,0,0,0],[0,0,0,0,0,1,0]],dtype=float)
    predicted=case['initial']/case['omega']+matrix@rx
    for s in range(3):
        for i in range(n):
            for d,j in enumerate(((i-1)%n,(i+1)%n)):predicted[s,i]-=active[2*s+d,i];predicted[s,j]+=active[2*s+d,i]
    for s in range(2):
        for edge,(i,j) in enumerate(zip(g['src'],g['dst'])):predicted[3+s,i]-=env[s,edge];predicted[3+s,j]+=env[s,edge]
    predicted[3:5]+=imports;return pop,predicted,imports

def main(base=BASE,deadline=None):
    started=time.perf_counter();deadline=deadline or started+300;base=Path(base);plan=json.loads((base/'plan.json').read_text(encoding='utf-8'));summary=json.loads((base/'qualification.json').read_text(encoding='utf-8'));checks=[]
    for name,sha in plan['source_hashes'].items():assert hashlib.sha256((ROOT/name).read_text(encoding='utf-8').replace('\r\n','\n').encode()).hexdigest()==sha
    for job in plan['jobs']:
        path=base/(job['id']+'.npz');case,g=fixture(job['fixture'],job['capacity']);k=np.array(job['kappa']);scale=np.maximum(1,(case['initial']/case['omega']).max(axis=1))
        with np.load(path) as raw:
            np.testing.assert_array_equal(raw['times'],TIMES)
            for name,value in dict(initial=case['initial'],a0=case['a0'],b0=case['b0'],omega=case['omega'],diff=case['diff'],src=g['src'],dst=g['dst'],weight=g['weight'],kappa=k).items():np.testing.assert_array_equal(raw[name],value)
            states=raw['state'];assert states.shape==(4,22*g['n']+2*len(g['src']))
        reference_states=integrate(case,g,k,job['mode'],deadline);populations=states[:,:7*g['n']].reshape(4,7,g['n']);error=float(np.max(abs(populations-reference_states)/scale[None,:,None]));assert error<=2e-6
        ledger_error=0.;balance=0.;minimum=float(states.min());assert minimum>=-1e-10
        for y in states:
            pop,predicted,imports=incidence(case,g,y);ledger_error=max(ledger_error,float(np.max(abs(pop-predicted)/scale[:,None])));balance=max(balance,abs(float(pop.sum()-case['initial'].sum()/case['omega']-imports.sum()))/max(1,case['initial'].sum()/case['omega']))
            if job['mode']=='FINITE':assert not np.any(imports)
        assert ledger_error<=1e-8 and balance<=1e-8 and np.min(np.diff(states[:,7*g['n']:],axis=0))>=-1e-10
        np.savez_compressed(base/(job['id']+'_reference.npz'),times=TIMES,population=reference_states)
        checks.append(dict(job=job['id'],sha256=hashlib.sha256(path.read_bytes()).hexdigest(),reference_error=error,ledger_error=ledger_error,balance_error=balance,minimum=minimum))
    # Independently repeat both limiting comparisons using separate references.
    c,g=fixture('mixing');mixed=integrate(c,g,[0,0],'FINITE',deadline,mixed=True);sc=np.array([140/9,140/9,29/(140/9),40,80,20,20]);mixerrors=[]
    for k in (1,10,100,1000):
        path=base/f'mixing_k{k}.npz'
        with np.load(path) as raw:pop=raw['state'][:,:7*g['n']].reshape(4,7,g['n'])
        mixerrors.append(float(np.max(abs(pop[1:]-mixed[1:])/sc[None,:,None])))
    np.testing.assert_allclose(mixerrors,summary['mixing_errors'],rtol=1e-6,atol=2e-8);assert np.all(np.diff(mixerrors)<0) and mixerrors[-1]<=.005
    assert np.max(abs(mixed.sum(axis=(1,2))-mixed[0].sum()))<1e-8
    c,g=fixture('chemostat');control=integrate(c,g,[.7,.7],'CHEMOSTAT',deadline);caperrors=[]
    for capacity in (1,10,100):
        with np.load(base/f'capacity_c{capacity}.npz') as raw:active=raw['state'][:,:3*g['n']].reshape(4,3,g['n'])
        caperrors.append(float(np.max(abs(active-control[:,:3])/sc[None,:3,None])))
    np.testing.assert_allclose(caperrors,summary['capacity_errors'],rtol=1e-6,atol=2e-8);assert np.all(np.diff(caperrors)<0) and caperrors[-1]<=.001
    atomic_json(base/'independent_check.json',dict(passed=True,raw_checks=checks,mixing_errors=mixerrors,capacity_errors=caperrors,elapsed_seconds=time.perf_counter()-started,dependencies='Same author/shared fixture constants/NumPy/SciPy. Separately written population equations, directed incidence ledger, Radau integration and limiting comparisons; primary frozen N1 continuum law not called. Not external replication.'))
    return checks
if __name__=='__main__':main()

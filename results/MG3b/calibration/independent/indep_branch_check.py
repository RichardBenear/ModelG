# Independent check (no calib_mg3b imports): S=0.073, 96 runs/arm, r=0.71, branch 'tau' Wald slope coverage
import numpy as np
from scipy import stats, optimize
rng=np.random.default_rng(777)
L=100; Om=np.array([130,120,110,100.]); S=0.073; k=4.6e-5*np.exp(-S*(Om-120)); ps=[0.15,0.35,0.33,0.31]
def run(n,lam,p):
    st=rng.random(n)<p; t=rng.exponential(1/lam,n); ev=(~st)&(t<=430); return ev.sum(), np.where(st,0,np.minimum(t,430)).sum()
def cell_stop(lam,p):
    N=E=0
    for b in range(4):
        n,e=run(16,lam,p); N+=n;E+=e
        if N>=20: break
    return N,E
def cp(n1,E1,n2,E2,lev):
    ci=stats.binomtest(n2,n1+n2).proportion_ci(lev,'exact'); f=lambda p: np.inf if p>=1 else p/(1-p)*E1/E2
    return f(ci.low),f(ci.high)
hits=[];R=6000
for _ in range(R):
    n1,E1=run(96,k[1]*L,ps[1]); n2,E2=run(96,0.71*k[1]*L,ps[1])
    lo,hi=cp(n1,E1,n2,E2,0.90)
    if not (lo>=0.5 and hi<=2): continue          # keep EQUIVALENT -> branch tau
    n=[0]*4;E=[0]*4
    for j in (0,2,3): n[j],E[j]=cell_stop(k[j]*L,ps[j])
    n[1],E[1]=n1,E1
    n=np.array(n,float);E=np.array(E); x=Om-Om.mean()
    nll=lambda th: -(n*(th[0]+th[1]*x+np.log(L*E))-L*E*np.exp(th[0]+th[1]*x)).sum()
    res=optimize.minimize(nll,[np.log(n.sum()/(L*E.sum())),0],method='BFGS')
    mu=L*E*np.exp(res.x[0]+res.x[1]*x); I=np.array([[mu.sum(),(mu*x).sum()],[(mu*x).sum(),(mu*x*x).sum()]])
    se=np.sqrt(np.linalg.inv(I)[1,1]); s=res.x[1]
    hits.append((s-1.96*se<=-S<=s+1.96*se, -S>s+1.96*se, -S<s-1.96*se))
h=np.array(hits); print('n branch',len(h),'cov %.4f miss_hi %.4f miss_lo %.4f'%tuple(h.mean(0)), 'MC se %.4f'%np.sqrt(.95*.05/len(h)))

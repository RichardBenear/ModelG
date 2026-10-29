"""Independent synthetic check of interval coverage after selecting/reusing a timestep arm.
No lattice simulation. Reimplements the exponential simulator and confidence intervals.
"""
from pathlib import Path
import json
import numpy as np
from scipy.stats import beta,chi2
ROOT=Path(__file__).parent
def cell(rng,lam,ps,N):
    startup=rng.random((N,64))<ps
    wait=rng.exponential(1/lam,(N,64))
    n=((~startup)&(wait<=430)).sum(1)
    E=np.where(startup,0,np.minimum(wait,430)).sum(1)
    return n,E
def ratio_ci(a,b,level):
    n1,E1=a;n2,E2=b;alpha=1-level;n=n1+n2
    lo=np.where(n2>0,beta.ppf(alpha/2,np.maximum(n2,1),n1+1),0)
    hi=np.where(n1>0,beta.ppf(1-alpha/2,n2+1,np.maximum(n1,1)),1)
    with np.errstate(divide='ignore'):
        return lo/(1-lo)*E1/E2,hi/(1-hi)*E1/E2
def decision(a,b):
    lo,hi=ratio_ci(a,b,.90);eq=(lo>=.5)&(hi<=2)
    lo,hi=ratio_ci(a,b,.95)
    return np.where(eq,0,np.where((lo>1)|(hi<1),1,2))
def check(cfg,R=100000):
    rng=np.random.default_rng(np.random.SeedSequence(20261030,spawn_key=(cfg['index'],)))
    total={b:dict(n=0,covered=0,miss_below=0,miss_above=0) for b in ['tau','tau2_inc','tau2_1b','stop']}
    ntrigger=0;falseeq=0
    lam=.0046
    for i in range(0,R,10000):
        N=min(10000,R-i)
        a=cell(rng,lam,cfg['ps1'],N);b=cell(rng,lam*cfg['r'],cfg['ps2'],N)
        c=cell(rng,lam*cfg['r']*cfg['r2'],cfg['ps2'],N)
        d=decision(a,b);e=decision(b,c)
        trigger=d==1;ntrigger+=int(trigger.sum());falseeq+=int(np.sum(trigger&(e==0)))
        branch=np.where(d==0,'tau',np.where(d==2,'tau2_inc',np.where(e==0,'tau2_1b','stop')))
        n=np.where(d==0,a[0],b[0]);E=np.where(d==0,a[1],b[1]);true=np.where(d==0,lam,lam*cfg['r'])
        lo=np.where(n>0,chi2.ppf(.025,2*np.maximum(n,1))/2,0)/E
        hi=chi2.ppf(.975,2*n+2)/2/E
        for name,v in total.items():
            s=branch==name;v['n']+=int(s.sum());v['covered']+=int(np.sum(s&(lo<=true)&(hi>=true)))
            v['miss_below']+=int(np.sum(s&(true<lo)));v['miss_above']+=int(np.sum(s&(true>hi)))
    for v in total.values():
        v['coverage']=v['covered']/v['n'] if v['n'] else None
        v['miss_lo']=v['miss_below']/v['n'] if v['n'] else None
        v['miss_hi']=v['miss_above']/v['n'] if v['n'] else None
    return dict(cfg=cfg,R=R,branches=total,n_part1b=ntrigger,n_part1b_equivalent=falseeq,part1b_equiv_given_trigger=falseeq/ntrigger if ntrigger else None,part1b_trigger_and_equiv=falseeq/R)
if __name__=='__main__':
    cfgs=[dict(index=i,r=r,r2=r2,ps1=.35,ps2=ps2) for i,(r,r2,ps2) in enumerate([(1,1,.35),(.71,1,.35),(.5,.5,.35),(.71,.5,.35),(.5,2,.35),(.5,.5,.05)])]
    out=[]
    for cfg in cfgs:
        row=check(cfg);out.append(row);print(json.dumps(row),flush=True)
    (ROOT/'selection_check.json').write_text(json.dumps(out,indent=2))

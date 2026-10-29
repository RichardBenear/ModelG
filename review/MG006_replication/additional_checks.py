"""Independent mean-map control and quantitative design diagnostics."""
import json,math
from pathlib import Path
import numpy as np
from scipy import stats
from numba import njit
from rescore import rows,prim,load
ROOT=Path(__file__).parent; INPUT=ROOT/'inputs'
@njit
def evolve(G,X,Y,steps=750000,dt=6.666666666666667e-4):
    g=G.copy();x=X.copy();y=Y.copy();N=len(g)
    for _ in range(steps):
        for j in range(N):
            k=(j-1)%N;l=(j+1)%N;v=X[j]*X[j]*Y[j]
            g[j]=G[j]+dt*((G[k]-2*G[j]+G[l])/0.04-G[j]+.1*X[j]+14)
            x[j]=X[j]+dt*((X[k]-2*X[j]+X[l])/0.04+G[j]-30*X[j]+v)
            y[j]=Y[j]+dt*(12*(Y[k]-2*Y[j]+Y[l])/0.04+29*X[j]-v)
        G,g=g,G;X,x=x,X;Y,y=y,Y
    return max(np.max(np.abs(G-140/9)),np.max(np.abs(X-140/9)),np.max(np.abs(Y-29/(140/9))))
if __name__=='__main__':
    a=np.full(500,140/9);c=np.full(500,29/(140/9));rng=np.random.default_rng(0)
    exact=evolve(a.copy(),a.copy(),c.copy())
    pert=evolve(a+1e-3*rng.standard_normal(500),a+1e-3*rng.standard_normal(500),c.copy())
    # exp(-lambda*50)-exp(-lambda*150) maximum at log(3)/100
    h=math.log(3)/100;p=math.exp(-50*h)-math.exp(-150*h)
    # Check that shared outcomes / attribution are benign in the stored Stage B records.
    rr=[r for o in [150,120,100] for r in rows(f'B_{o}')]
    hist={str(o):[dict(rep=r['rep'],first=r['first']['t_start'] if r['first'] else None,births=r['births']) for r in rows(f'B_{o}')] for o in [150,120,100]}
    # Validity counterexample for the claimed exact Poisson count law: one risk subject, no finite censoring.
    # Event count is identically one, exposure is exponentially distributed, hence conditional N|E is not Poisson.
    # Directly measure the resulting Garwood coverage, which is conservative in this example, not exact 95%.
    e=rng.exponential(size=1000000);lo=stats.chi2.ppf(.025,2)/2;hi=stats.chi2.ppf(.975,4)/2
    cov=np.mean((lo/e<=1)&(hi/e>=1))
    # List the three b=29.3 post-burn-in event waits; their entire exposure precedes t=200.
    waits=[r['first']['t_start']-50 for r in rows('C4_b29.3_120') if prim(r)['event']]
    out=dict(N2=dict(exact=exact,perturbed=pert),stageA=dict(max_probability=p,hazard_at_max=h,four_of_four=p**4,max_chance_one_four_of_four_12_independent_cells=1-(1-p**4)**12),shared_stageB=sum(b['shared'] for r in rr for b in r['births']),persistent_near_false=sum(b['cls']=='persistent' and b['near'] is False for r in rr for b in r['births']),C4_high_b_waits=waits,garwood_counterexample_coverage=cov,garwood_counterexample_analytic_coverage=math.exp(-lo)-math.exp(-hi),birth_history=hist)
    (ROOT/'additional_checks.json').write_text(json.dumps(out,indent=2))
    print(json.dumps({k:v for k,v in out.items() if k!='birth_history'},indent=2))

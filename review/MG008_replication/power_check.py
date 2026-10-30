"""Independent Monte Carlo check of the reported stationarity-power labels.
Uses inverse integrated hazards and an independently computed exposure median.
"""
from pathlib import Path
import json
import numpy as np
from scipy import stats
ROOT=Path(__file__).parent
def power(lam,ps,s,fac,fixed,R=10000,seed=1):
    rng=np.random.default_rng(seed);start=rng.random((R,64))<ps;e=rng.exponential(size=(R,64))
    u=np.where(e<lam*s,e/lam,s+(e-lam*s)/(fac*lam));event=(u<=430)&~start
    exposure=np.where(start,0,np.minimum(u,430))
    if fixed:used=np.ones((R,64),bool)
    else:
        hit=event.reshape(R,4,16).sum(2).cumsum(1)>=20
        batches=np.where(hit.any(1),hit.argmax(1)+1,4)
        used=np.arange(64)[None,:]<16*batches[:,None]
    exposure*=used;event&=used
    lo=np.zeros(R);hi=np.full(R,430.)
    for _ in range(45):
        mid=(lo+hi)/2;below=np.minimum(exposure,mid[:,None]).sum(1)<exposure.sum(1)/2
        lo=np.where(below,mid,lo);hi=np.where(below,hi,mid)
    split=(lo+hi)/2;n=event.sum(1);ne=(event&(u<split[:,None])).sum(1)
    p=np.minimum(1,2*stats.binom.cdf(np.minimum(ne,n-ne),n,.5))
    return float(np.mean(p<.01))
if __name__=='__main__':
    r=json.loads((ROOT/'independent_rescore.json').read_text());out={}
    for i,(om,c) in enumerate(r['part2'].items()):
        out[om]={}
        for j,fac in enumerate([2,.5]):
            p=power(c['events']/c['exposure'],c['excluded']/c['runs'],c['stationarity']['split']-50,fac,om=='120',seed=20261030+i*10+j)
            out[om][str(fac)]=dict(power=p,R=10000,mc_se=float(np.sqrt(p*(1-p)/10000)))
    (ROOT/'independent_power.json').write_text(json.dumps(out,indent=2));print(json.dumps(out,indent=2))

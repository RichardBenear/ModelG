"""MG128 count-space accepted-event strict-tick prototype; not qualified."""
import numpy as np
from mg12h1_reference import transitions

def channels(available,config):
    """Return rule, output increment, reserved snapshot counts, total rate."""
    a=np.asarray(available)
    if a.dtype.kind not in 'iu' or a.ndim!=2 or a.shape[0]!=7 or np.any(a<0):raise ValueError('Invalid counts')
    for rule,dest,rate in transitions(a,*config):
        change=np.asarray(dest,dtype=np.int64).reshape(a.shape)-a
        reserved=np.maximum(-change,0)
        if rule=='R4':
            site=int(np.flatnonzero(change[2])[0]);reserved[1,site]=2
        if not np.isfinite(rate) or rate<=0 or np.any(reserved>a):raise ValueError('Invalid channel')
        yield rule,change,reserved,rate

def apply(available,output,channel):
    rule,change,reserved,rate=channel
    if np.any(reserved>available):raise ValueError('Stale channel')
    remaining=available-reserved;result=output+change
    if np.any(result<0) or result.sum()!=output.sum():raise ValueError('Invalid update')
    return remaining,result


def tick(snapshot,config,delta,rng,max_events=100000):
    """Accepted-event race; injected RNG, separate available/output populations."""
    import math
    if not math.isfinite(delta) or delta<=0:raise ValueError('Invalid tick')
    available=np.array(snapshot,copy=True);output=np.array(snapshot,copy=True);clock=0.;trace=[]
    while True:
        options=list(channels(available,config));total=sum(c[3] for c in options)
        if not total:break
        wait=float(rng.exponential(1/total))
        if not math.isfinite(wait) or wait<=0:raise ValueError('Invalid exponential draw')
        clock+=wait
        if clock>=delta:break
        if len(trace)>=max_events:raise ValueError('Accepted event cap')
        u=float(rng.random())
        if not 0<=u<1:raise ValueError('Invalid uniform draw')
        threshold=u*total;acc=0.;chosen=None
        for c in options:
            acc+=c[3]
            if threshold<acc:chosen=c;break
        if chosen is None:raise ValueError('Channel selection failed')
        available,output=apply(available,output,chosen)
        trace.append(dict(time=clock,rule=chosen[0],increment=chosen[1].tolist(),reserved=chosen[2].tolist()))
    return dict(counts=output,available=available,trace=trace)

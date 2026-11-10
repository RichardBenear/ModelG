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

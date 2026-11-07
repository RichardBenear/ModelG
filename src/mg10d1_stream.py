"""MG053 resumable-in-process segments of the qualified D1 event sequence."""
import numpy as np
import numba as nb
from mg10d1 import _controlled_leap


@nb.njit(cache=True)
def seed_rng(seed):
    np.random.seed(seed)


@nb.njit(cache=True)
def advance(state,flux,hops,external,a0,b0,omega,diff,mask,tau,start,target,clamp_a,clamp_b):
    t=start;leaps=0;adapted=0;maximum=0.
    while t<target:
        requested=min(tau,target-t)
        dt,prob=_controlled_leap(state,flux,hops,a0,b0,omega,diff,mask,requested,clamp_a,clamp_b,external)
        if dt<=0:raise ValueError("Nonpositive step")
        adapted+=int(dt<requested*(1-1e-12));leaps+=1;maximum=max(maximum,prob)
        if dt==target-t:t=target
        else:t+=dt
    return np.array([float(leaps),float(adapted),maximum])

"""MG055 stable finite-binomial source at huge n/tiny p, not Poisson sampling."""
import numpy as np
import numba as nb


@nb.njit(cache=True)
def geometric_binomial(n,p):
    # The spacings between successes in independent Bernoulli trials are
    # geometric. Count successes whose cumulative spacings fit in n trials.
    # log1p avoids forming rounded 1-p. Compare before converting to int so a
    # spacing beyond the remaining trial count cannot overflow an integer.
    if n<0 or not 0<=p<=1:raise ValueError("Invalid binomial parameters")
    if n==0 or p==0:return 0
    if p==1:return n
    remaining=n;successes=0;denominator=np.log1p(-p)
    while remaining>0:
        failures=np.floor(np.log1p(-np.random.random())/denominator)
        if failures>=remaining:break
        remaining-=int(failures)+1
        successes+=1
    return successes


@nb.njit(cache=True)
def finite_source(n,p):
    if n<0 or not 0<=p<=1:raise ValueError("Invalid finite source parameters")
    if n<=1000000000:return np.random.binomial(n,p)
    if p==0:return 0
    if p==1:return n
    if n>9007199254740992 or p>=1e-6 or n*p>30:
        raise ValueError("Outside qualified huge-stock rare-source domain")
    return geometric_binomial(n,p)

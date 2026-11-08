"""MG083 degree/distance-matched environmental geometry-supply proposal.
Uses frozen N1 exact SSA; no new reaction or transport kernel.
"""
import numpy as np
from scipy.linalg import expm
from scipy.stats import beta
from mg11n1 import graph
SOURCE=0
CONTACT=3
TIMES=np.array([0.,.5,1.,2.,4.])
TOPOLOGIES=('TRIANGULAR_PRISM','COMPLETE_BIPARTITE_3_3')
def topology(index):
    if index==0:edges=[(0,1),(1,2),(0,2),(3,4),(4,5),(3,5),(0,3),(1,4),(2,5)]
    elif index==1:edges=[(i,j) for i in range(3) for j in range(3,6)]
    else:raise ValueError('Invalid topology')
    return graph(6,[(i,j,1.) for i,j in edges])
def cell(resource,network):
    if resource not in (0,1):raise ValueError('Invalid resource')
    g=topology(network);pop=np.zeros((7,6),dtype=np.int64);mask=np.zeros(7)
    if resource==0:pop[3,SOURCE]=1;mask[0]=1
    else:pop[4,SOURCE]=1;pop[1,CONTACT]=1;mask[3]=1
    case=dict(initial=pop,a0=np.full(6,28.),b0=np.full(6,36.25),omega=1.,diff=np.zeros(3))
    return case,g,np.full(2,1/3),TIMES.copy(),mask
def generator(case,g,k,resource):
    n=g['n'];size=2*n if resource==0 else n+1;Q=np.zeros((size,size))
    for i,j,w in zip(g['src'],g['dst'],g['weight']):Q[i,j]+=k[resource]*w
    if resource==0:
        rates=14*case['omega']/case['a0']
        for i in range(n):Q[i,n+i]=rates[i]
    else:
        for i in range(n):Q[i,n]=29*case['initial'][1,i]/case['b0'][i]
    Q[np.arange(size),np.arange(size)]=-Q.sum(axis=1);return Q
def exact_probabilities(case,g,k,resource,t):
    p=np.zeros(2*g['n'] if resource==0 else g['n']+1);p[int(case['initial'][3+resource].argmax())]=1
    return p@expm(generator(case,g,k,resource)*t)
def delivery(data,resource):
    # Count actual contact events, not a detector/particle classification.
    return data[:,:,7 if resource==0 else 10,CONTACT]
def source_seed_case(resource,network):return resource*2+network
def binomial_interval(k,n,alpha=.025):
    if not isinstance(k,(int,np.integer)) or not isinstance(n,(int,np.integer)) or n<0 or not 0<=k<=n:raise ValueError('Invalid count')
    if n==0:return None
    return [0. if k==0 else float(beta.ppf(alpha/2,k,n-k+1)),1. if k==n else float(beta.ppf(1-alpha/2,k+1,n-k))]
def difference_interval(k_prism,n_prism,k_bipartite,n_bipartite):
    a=binomial_interval(k_prism,n_prism);b=binomial_interval(k_bipartite,n_bipartite)
    if a is None or b is None:return dict(mean=None,ci=None)
    return dict(mean=k_bipartite/n_bipartite-k_prism/n_prism,ci=[b[0]-a[1],b[1]-a[0]],scope='At least95% joint two-group coverage from two97.5% two-sided exact binomial intervals; B primaryT4 only, no full secondary family coverage')

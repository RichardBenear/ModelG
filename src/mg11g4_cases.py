"""MG096 fixed small integer fixtures and separately labeled performance input."""
import numpy as np
from mg11n1 import graph
TAUS=(0.,.001,.0005,.00025,0.)
R=16384;CHUNK=2048

def fixture(c):
    if c==0:
        pop=np.array([[16],[16],[2],[40],[40],[0],[0]],np.int64);a0=b0=np.array([40.]);om=1.;diff=np.zeros(3);ga=gb=graph(1,[]);times=np.array([0.,.05,.1]);k=np.zeros(2)
    elif c==1:
        pop=np.array([[2,3,1],[4,5,2],[1,2,1],[8,1,4],[1,2,0],[0,0,0],[0,0,0]],np.int64);a0=np.full(3,8.);b0=np.full(3,4.);om=2.;diff=np.array([.2,.3,.4]);ga=graph(3,[(0,1,1.),(1,2,.5)]);gb=graph(3,[(0,2,.75),(2,1,1.25)]);times=np.array([0.,.05,.1]);k=np.array([.4,.7])
    elif c==2:
        pop=np.array([[10,14,18],[12,20,28],[1,2,4],[20,40,10],[1,8,40],[0,0,0],[0,0,0]],np.int64);a0=b0=np.full(3,40.);om=2.;diff=np.array([25.,25.,300.]);ga=graph(3,[(0,1,1.),(1,2,.5)]);gb=graph(3,[(0,2,.75),(2,1,1.25)]);times=np.array([0.,.01,.02]);k=np.array([.4,.7])
    else:raise ValueError('Unknown fixture')
    return dict(initial=pop,a0=a0.copy(),b0=b0.copy(),omega=om,diff=diff),ga,gb,k,times,np.ones(7)

def one_step(large=False):
    pop=np.zeros((7,2),np.int64);pop[3,0]=100000000 if large else 1
    if not large:pop[4,0]=1;pop[1,1]=1
    case=dict(initial=pop,a0=np.full(2,28.),b0=np.full(2,36.25),omega=1.,diff=np.zeros(3));g=graph(2,[(0,1,1.)]);mask=np.zeros(7)
    if not large:mask[[0,3]]=1
    return case,g,g,np.ones(2),np.array([0.,.00001 if large else .1]),mask

def benchmark(cross=False):
    from mg11g2_pattern_v2 import factory as deterministic
    case,ga,k,_=deterministic('LOCAL_C1',1);_,gb,_,_=deterministic('CROSS_RING_C1' if cross else 'LOCAL_C1',1)
    initial=np.rint(case['initial']*1000).astype(np.int64);return dict(initial=initial,a0=np.rint(case['a0']*1000),b0=np.rint(case['b0']*1000),omega=1000.,diff=case['diff'].copy()),ga,gb,k,np.array([0.,.001,.002]),np.ones(7)

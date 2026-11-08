"""MG087 fixed small nonlinear coupling fixtures; no pattern/parameter search."""
import numpy as np
from mg11n1 import graph

TIMES=np.array([0.,.25,.5,1.])
def fixture(name,capacity=1):
    if name in ('legacy','capacity','chemostat'):
        n=6;omega=9.;pop=np.zeros((7,n),dtype=np.int64);pop[:3]=np.array([[150,160,140,170,150,160],[140,160,150,170,140,160],[20,15,21,18,22,16]])
        a0=np.full(n,37800.*capacity);b0=np.full(n,1218000.*capacity);pop[3]=a0.astype(np.int64);pop[4]=b0.astype(np.int64);diff=np.array([25.,25.,300.]);g=graph(n,[(i,(i+1)%n,1.) for i in range(n)])
    elif name=='mixing':
        n=3;omega=10.;pop=np.array([[140,160,150],[140,150,160],[20,18,19],[200,400,600],[600,800,1000],[0,0,0],[0,0,0]],dtype=np.int64);a0=np.full(n,400.);b0=np.full(n,800.);diff=np.array([.2,.3,.4]);g=graph(n,[(0,1,1.),(1,2,.7)])
    else:raise ValueError('Unknown fixture')
    return dict(initial=pop,a0=a0,b0=b0,omega=omega,diff=diff),g

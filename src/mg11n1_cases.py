"""MG079 prospective fixed network qualification cells, not selected by outcomes."""
import numpy as np
from mg11n1 import graph
def cell(index):
    mask=np.zeros(7);diff=np.zeros(3)
    if index in (0,1,2):
        n=3;pop=np.zeros((7,n),dtype=np.int64);g=graph(n,[(0,1,1.),(1,2,.5)]);a0=np.array([20.,28.,35.]);b0=np.full(n,36.25);omega=1.;times=np.array([0.,.25,.5]);kappa=np.array([.4,.7])
        if index==0:name='transport';pop[3,0]=1;pop[4,2]=1
        elif index==1:name='A_source';pop[3,0]=1;mask[0]=1
        else:name='B_contact';pop[4,2]=1;pop[1,0]=1;mask[3]=1;times=np.array([0.,.5,1.])
    elif index==3:
        name='full_tiny';n=2;pop=np.array([[2,2],[2,1],[1,0],[4,2],[1,3],[0,0],[0,0]],dtype=np.int64);g=graph(n,[(0,1,.7)]);a0=np.full(n,10.);b0=np.full(n,5.);omega=2.;diff=np.array([.2,.2,1.]);mask[:]=1;times=np.array([0.,.05,.1]);kappa=np.array([.4,.6])
    else:raise ValueError(index)
    return dict(name=name,initial=pop,a0=a0,b0=b0,omega=omega,diff=diff),g,kappa,times,mask

"""MG106 independent count generator and finite-state master equation.

No primary matcher/commit/executor calls. Counts collapse exchangeable identities.
"""
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import expm_multiply
def transitions(pop,omega,refs,links,diffusion,coupling,mask):
    n=pop.shape[1];out=[]
    def add(rule,i,changes,rate):
        if rate<=0:return
        other=pop.copy()
        for s,j,amount in changes:other[s,j]+=amount
        if np.any(other<0):raise ValueError('Reference invalid transition')
        out.append((rule,tuple(other.ravel()),float(rate)))
    for i in range(n):
        g,x,y,a,b,z,w=map(int,pop[:,i]);a0,b0=refs[i]
        rates=[14*omega*a/a0,g,.1*x,29*b*x/b0,x*(x-1)*y/omega**2,.9*x]
        increments=[[(3,i,-1),(0,i,1)],[(0,i,-1),(1,i,1)],[(1,i,-1),(0,i,1)],[(4,i,-1),(1,i,-1),(2,i,1),(5,i,1)],[(2,i,-1),(1,i,1)],[(1,i,-1),(6,i,1)]]
        for rule,changes,rate,on in zip(('R1','R2','R2r','R3','R4','R5'),increments,rates,mask):
            if on:add(rule,i,changes,rate)
    for kind,i,j,weight in links:
        for s in range(5):
            if kind=='ActiveLink' and s<3:rate=diffusion[s]*weight*pop[s,i]
            elif kind=='ALink' and s==3:rate=coupling[0]*weight*pop[s,i]
            elif kind=='BLink' and s==4:rate=coupling[1]*weight*pop[s,i]
            else:continue
            add('H_'+('G','X','Y','A','B')[s],i,[(s,i,-1),(s,j,1)],rate)
    return out
def master(initial,omega,refs,links,diffusion,coupling,mask,horizon,max_states=2048):
    states=[tuple(initial.ravel())];indices={states[0]:0};edges=[];left=0
    while left<len(states):
        pop=np.asarray(states[left],np.int64).reshape(initial.shape);total=0.
        for rule,dest,rate in transitions(pop,omega,refs,links,diffusion,coupling,mask):
            if dest not in indices:
                if len(states)>=max_states:raise ValueError('Reference state-space cap')
                indices[dest]=len(states);states.append(dest)
            edges.append((left,indices[dest],rate));total+=rate
        edges.append((left,left,-total));left+=1
    row,col,data=zip(*edges) if edges else ([],[],[]);q=coo_matrix((data,(row,col)),shape=(len(states),len(states))).tocsr();p0=np.zeros(len(states));p0[0]=1;prob=expm_multiply(q.T*horizon,p0);assert prob.min()>=-1e-12 and abs(prob.sum()-1)<1e-10
    return states,prob,q

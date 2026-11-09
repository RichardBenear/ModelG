"""MG092 separately assembled two-resource graph incidence and Radau reference."""
import time
import numpy as np
from scipy.sparse import coo_matrix
from scipy.integrate import solve_ivp
from mg11g2_pattern_v2_oracle import population_rhs,population_pattern,laplacians,dense_projection,late

def independent_solve(case,ga,gb,k,times,deadline):
    n=ga['n'];mats=laplacians(ga);LB,_=laplacians(gb);base=population_pattern(ga);rows=np.r_[4*n+gb['src'],4*n+gb['dst']];cols=np.tile(4*n+gb['src'],2);jac=(base+coo_matrix((np.ones(len(rows),dtype=bool),(rows,cols)),shape=base.shape)).tocsr().astype(bool)
    def fun(t,y):
        if time.perf_counter()>deadline:raise TimeoutError('G3 Radau deadline')
        out=population_rhs(t,y,case,ga,[k[0],0.],'FINITE',mats).reshape(7,n);out[4]+=k[1]*(LB@y.reshape(7,n)[4]);return out.ravel()
    result=solve_ivp(fun,(0,float(times[-1])),case['initial'].ravel(),t_eval=times,method='Radau',jac_sparsity=jac,rtol=1e-11,atol=1e-12,max_step=.025)
    if not result.success or result.y.shape!=(7*n,len(times)):raise RuntimeError('G3 reference incomplete')
    return result.y.T.reshape(len(times),7,n)

def accounting(state,case,ga,gb):
    n=ga['n'];e=len(ga['src']);t=len(state);pop=state[:,:7*n].reshape(t,7,n);rx=state[:,7*n:14*n].reshape(t,7,n);hop=state[:,14*n:20*n].reshape(t,6,n);env=state[:,20*n:20*n+2*e].reshape(t,2,e);imports=state[:,-2*n:].reshape(t,2,n)
    F=rx;pred=case['initial'][None]+np.stack((F[:,0]-F[:,1]+F[:,2]-F[:,6],F[:,1]-F[:,2]-F[:,3]+F[:,4]-F[:,5],F[:,3]-F[:,4],-F[:,0]+F[:,6],-F[:,3],F[:,3],F[:,5]),axis=1)
    nodes=np.arange(n);src=np.tile(nodes,2);dest=np.r_[(nodes-1)%n,(nodes+1)%n];idx=np.arange(2*n);D=coo_matrix((np.r_[np.ones(2*n),-np.ones(2*n)],(np.r_[dest,src],np.r_[idx,idx])),shape=(n,2*n)).tocsr()
    for s in range(3):pred[:,s]+=(D@hop[:,2*s:2*s+2].reshape(t,2*n).T).T
    for s,g in enumerate((ga,gb)):
        I=coo_matrix((np.r_[np.ones(e),-np.ones(e)],(np.r_[g['dst'],g['src']],np.tile(np.arange(e),2))),shape=(n,e)).tocsr();pred[:,3+s]+=(I@env[:,s].T).T+imports[:,s]
    scales=np.maximum(1,case['initial'].max(axis=1));local=float(np.max(abs(pop-pred)/scales[None,:,None]));total=float(np.max(abs(pop.sum(axis=(1,2))-case['initial'].sum()))/case['initial'].sum());minimum=float(state.min());increment=float(np.diff(state[:,7*n:],axis=0).min()) if t>1 else 0.
    return dict(local_residual=local,total_residual=total,minimum=minimum,minimum_ledger_increment=increment,passed=bool(local<=1e-7 and total<=1e-10 and minimum>=-1e-9 and increment>=-1e-9 and not np.any(imports)))

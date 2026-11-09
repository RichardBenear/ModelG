"""MG088 separate sparse population equations, incidence and real-space metrics."""
import numpy as np
from scipy.sparse import coo_matrix,diags
from scipy.integrate import solve_ivp
from mg11g2_pattern import S0,N,OBS

def laplacians(g):
    n=g['n'];src=g['src'];dst=g['dst'];w=g['weight'];L=coo_matrix((np.r_[w,-w],(np.r_[dst,src],np.r_[src,src])),shape=(n,n)).tocsr()
    nodes=np.arange(n);active=coo_matrix((np.tile([1.,1.,-2.],n),(np.repeat(nodes,3),np.stack(((nodes-1)%n,(nodes+1)%n,nodes),axis=1).ravel())),shape=(n,n)).tocsr();return L,active

def population_rhs(t,y,case,g,k,mode,matrices=None):
    L,D=laplacians(g) if matrices is None else matrices;n=g['n'];z=y.reshape(7,n);gg,x,yy,a,b=z[:5];out=np.zeros_like(z)
    source=14*a/case['a0'];conversion=29*x*b/case['b0'];feedback=x*x*yy
    out[0]=source-gg+.1*x;out[1]=gg-.1*x-conversion+feedback-.9*x;out[2]=conversion-feedback;out[3]=-source;out[4]=-conversion;out[5]=conversion;out[6]=.9*x
    for s in range(3):out[s]+=case['diff'][s]*(D@z[s])
    for s in range(2):out[3+s]+=k[s]*(L@z[3+s])
    if mode=='CHEMOSTAT':out[3:5]=0
    return out.ravel()

def population_pattern(g):
    n=g['n'];L,D=laplacians(g);rows=[];cols=[]
    for i in range(n):
        for s in range(7):
            for c in range(5):rows.append(s*n+i);cols.append(c*n+i)
    for s,matrix in [(s,D) for s in range(3)]+[(s,L) for s in (3,4)]:
        q=matrix.tocoo();rows.extend((s*n+q.row).tolist());cols.extend((s*n+q.col).tolist())
    return coo_matrix((np.ones(len(rows),dtype=bool),(rows,cols)),shape=(7*n,7*n)).tocsr()

def independent_solve(case,g,k,mode,times,deadline):
    import time
    matrices=laplacians(g)
    def fun(t,y):
        if time.perf_counter()>deadline:raise TimeoutError('Independent solve deadline')
        return population_rhs(t,y,case,g,k,mode,matrices)
    result=solve_ivp(fun,(0,float(times[-1])),case['initial'].ravel(),t_eval=times,method='Radau',jac_sparsity=population_pattern(g),rtol=1e-11,atol=1e-12,max_step=.025)
    if not result.success or result.y.shape!=(7*g['n'],len(times)) or not np.all(np.isfinite(result.y)):raise RuntimeError('Independent Radau failed')
    return result.y.T.reshape(len(times),7,g['n'])

def incidence_check(state,case,g):
    n=g['n'];e=len(g['src']);t=len(state);pop=state[:,:7*n].reshape(t,7,n);rx=state[:,7*n:14*n].reshape(t,7,n);hops=state[:,14*n:20*n].reshape(t,6,n);env=state[:,20*n:20*n+2*e].reshape(t,2,e);imports=state[:,20*n+2*e:].reshape(t,2,n)
    matrix=np.array([[1,-1,1,0,0,0,-1],[0,1,-1,-1,1,-1,0],[0,0,0,1,-1,0,0],[-1,0,0,0,0,0,1],[0,0,0,-1,0,0,0],[0,0,0,1,0,0,0],[0,0,0,0,0,1,0]],dtype=float)
    predicted=case['initial'][None]+np.einsum('sr,trn->tsn',matrix,rx)
    nodes=np.arange(n);src=np.tile(nodes,2);dst=np.r_[(nodes-1)%n,(nodes+1)%n];index=np.arange(2*n);A=coo_matrix((np.r_[np.ones(2*n),-np.ones(2*n)],(np.r_[dst,src],np.r_[index,index])),shape=(n,2*n)).tocsr()
    B=coo_matrix((np.r_[np.ones(e),-np.ones(e)],(np.r_[g['dst'],g['src']],np.tile(np.arange(e),2))),shape=(n,e)).tocsr()
    for s in range(3):predicted[:,s]+=(A@hops[:,2*s:2*s+2].reshape(t,2*n).T).T
    for s in range(2):predicted[:,3+s]+=(B@env[:,s].T).T+imports[:,s]
    scales=np.maximum(1,case['initial'].max(axis=1));local=float(np.max(abs(pop-predicted)/scales[None,:,None]));total=float(np.max(abs(pop.sum(axis=(1,2))-case['initial'].sum()-imports.sum(axis=(1,2))))/case['initial'].sum());minimum=float(state.min());monotone=float(np.diff(state[:,7*n:],axis=0).min()) if t>1 else 0.
    return dict(local_residual=local,total_residual=total,minimum=minimum,minimum_ledger_increment=monotone,passed=bool(local<=1e-7 and total<=1e-10 and minimum>=-1e-9 and monotone>=-1e-9))

def dense_projection(pop,phi):
    n=pop.shape[-1];frequency=2*np.pi*np.fft.fftfreq(n,d=.2);kernel=np.fft.ifft(np.exp(-.5*(.4*frequency)**2)).real;matrix=kernel[(np.arange(n)[:,None]-np.arange(n)[None,:])%n]
    template=phi-phi.mean(axis=-1,keepdims=True);template=template@matrix.T;fields=pop[:,:3];centered=fields-fields.mean(axis=-1,keepdims=True);smooth=centered@matrix.T;scores=np.einsum('tsn,sn->ts',smooth,template)/np.sum(template*template,axis=-1)[None];return scores

def late(values,times):
    index=np.flatnonzero((times>=20)&(times<=25))
    if len(index)!=21 or not np.array_equal(times[index],np.arange(80,101)/4):raise ValueError('Incomplete late window')
    values=values[index];return np.sum((values[:-1]+values[1:])*.125,axis=0)/5

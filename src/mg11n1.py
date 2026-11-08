"""MG079 fixed equal-volume environmental network; exact small-system CTMC oracle.
No replenishment, recycling or growing geometry. Historical kernels unchanged.
"""
from pathlib import Path
import numpy as np
import numba as nb
from mg10s1 import STOICH,_ssa_clock as legacy_clock
ROOT=Path(__file__).resolve().parents[1]

def graph(n,edges):
    if not isinstance(n,(int,np.integer)) or n<1:raise ValueError('Invalid graph size')
    seen=set();directed=[]
    for i,j,w in edges:
        if not isinstance(i,(int,np.integer)) or not isinstance(j,(int,np.integer)) or not 0<=i<n or not 0<=j<n or i==j or not np.isfinite(w) or w<=0:raise ValueError('Invalid edge')
        key=tuple(sorted((i,j)))
        if key in seen:raise ValueError('Duplicate undirected edge')
        seen.add(key);directed.extend(((i,j,float(w)),(j,i,float(w))))
    directed.sort()
    return dict(n=n,src=np.array([i for i,j,w in directed],dtype=np.int64),dst=np.array([j for i,j,w in directed],dtype=np.int64),weight=np.array([w for i,j,w in directed],dtype=float))

def validate(case,g,kappa,mask):
    pop=case['initial'];n=g['n']
    if pop.shape!=(7,n) or pop.dtype!=np.int64 or np.any(pop<0) or sum(map(int,pop.ravel()))>2**53:raise ValueError('Invalid or excessive counts')
    for key in ('a0','b0'):
        if case[key].shape!=(n,) or not np.all(np.isfinite(case[key])) or np.any(case[key]<=0):raise ValueError('Invalid fixed denominators')
    if not np.isfinite(case['omega']) or case['omega']<=0 or case['diff'].shape!=(3,) or np.any(case['diff']<0) or not np.all(np.isfinite(case['diff'])):raise ValueError('Invalid active rates')
    if kappa.shape!=(2,) or np.any(kappa<0) or not np.all(np.isfinite(kappa)) or mask.shape!=(7,) or np.any(mask<0) or not np.all(np.isfinite(mask)):raise ValueError('Invalid network/reaction rates')
    e=len(g['src'])
    if g['src'].dtype!=np.int64 or g['dst'].dtype!=np.int64 or g['src'].shape!=(e,) or g['dst'].shape!=(e,) or g['weight'].shape!=(e,):raise ValueError('Invalid edge arrays')
    edges={}
    for i,j,w in zip(g['src'],g['dst'],g['weight']):
        if not 0<=i<n or not 0<=j<n or i==j or not np.isfinite(w) or w<=0 or (int(i),int(j)) in edges:raise ValueError('Invalid directed graph')
        edges[int(i),int(j)]=w
    if any(edges.get((j,i))!=w for (i,j),w in edges.items()):raise ValueError('Only symmetric equal-volume networks qualified')

def simulate(case,g,kappa,times,reps,seed,mask=None,max_events=1000000):
    kappa=np.asarray(kappa,dtype=float);mask=np.ones(7) if mask is None else np.asarray(mask,dtype=float);times=np.asarray(times,dtype=float)
    validate(case,g,kappa,mask)
    if not isinstance(seed,(int,np.integer)) or not 0<=seed<2**32:raise ValueError('Invalid seed')
    if times.ndim!=1 or not len(times) or times[0]<0 or not np.all(np.isfinite(times)) or np.any(np.diff(times)<=0) or not isinstance(reps,(int,np.integer)) or reps<1 or not isinstance(max_events,(int,np.integer)) or max_events<1:raise ValueError('Invalid run bounds')
    return network_chunk(case['initial'],case['a0'],case['b0'],case['omega'],case['diff'],mask,times,reps,seed,g['src'],g['dst'],g['weight'],kappa,max_events)

def continuum_rhs(t,flat,case,g,kappa):
    """Concentration mean-field equations, distinct from finite-count R4 drift."""
    z=np.asarray(flat).reshape(7,g['n']);G,X,Y,A,B=z[:5];a0=case['a0']/case['omega'];b0=case['b0']/case['omega']
    flux=np.array([14*A/a0,G,.1*X,29*X*B/b0,X*X*Y,.9*X,np.zeros(g['n'])]);out=STOICH@flux
    for s in range(3):out[s]+=case['diff'][s]*(np.roll(z[s],1)+np.roll(z[s],-1)-2*z[s])
    for s in range(2):
        flow=kappa[s]*g['weight']*z[3+s,g['src']]
        np.add.at(out[3+s],g['src'],-flow);np.add.at(out[3+s],g['dst'],flow)
    return out.ravel()

@nb.njit(cache=True)
def _network_clock(pop,a0,b0,omega,diff,mask,src,dst,weight,kappa):
    if len(src)==0 or (kappa[0]==0 and kappa[1]==0):return legacy_clock(pop,a0,b0,omega,diff,mask)
    # Separate direct propensity implementation; no leap helpers or factorization.
    n=pop.shape[1];rates=np.zeros(13*n+2*len(src))
    for node in range(n):
        j=13*node;G,X,Y,A,B=pop[0,node],pop[1,node],pop[2,node],pop[3,node],pop[4,node]
        rates[j]=mask[0]*14*omega*A/a0[node]
        rates[j+1]=mask[1]*G
        rates[j+2]=mask[2]*X/10
        rates[j+3]=mask[3]*29*X*B/b0[node]
        rates[j+4]=mask[4]*X*(X-1)*Y/(omega*omega) if X>=2 else 0.
        rates[j+5]=mask[5]*9*X/10
        for species in range(3):
            rates[j+7+2*species]=diff[species]*pop[species,node]
            rates[j+8+2*species]=diff[species]*pop[species,node]
    for s in range(2):
        for edge in range(len(src)):rates[13*n+s*len(src)+edge]=kappa[s]*weight[edge]*pop[3+s,src[edge]]
    total=rates.sum()
    if not np.isfinite(total):raise ValueError("Nonfinite total event rate")
    if total==0:return np.inf,-1
    waiting=np.random.exponential(1/total)
    mark=np.random.random()*total;accumulated=0.
    for event in range(len(rates)):
        accumulated+=rates[event]
        if mark<accumulated:return waiting,event
    raise ValueError("SSA event selection failed")

@nb.njit(cache=True)
def network_chunk(initial,a0,b0,omega,diff,mask,observations,reps,seed,src,dst,weight,kappa,max_events):
    np.random.seed(seed)
    n=initial.shape[1];out=np.zeros((reps,len(observations),20,n),dtype=np.int64);stats=np.zeros((reps,3));environment=np.zeros((reps,len(observations),2,len(src)),dtype=np.int64)
    for r in range(reps):
        pop=initial.copy();reactions=np.zeros((7,n),dtype=np.int64);hops=np.zeros((6,n),dtype=np.int64)
        ef=np.zeros((2,len(src)),dtype=np.int64)
        future,event=_network_clock(pop,a0,b0,omega,diff,mask,src,dst,weight,kappa);events=0
        for k,target in enumerate(observations):
            while future<=target:
                if events>=max_events:raise ValueError("Per-trajectory exact-event budget reached")
                if event>=13*n:
                    q=event-13*n;species=q//len(src);edge=q%len(src)
                    pop[3+species,src[edge]]-=1;pop[3+species,dst[edge]]+=1;ef[species,edge]+=1
                else:
                    node=event//13;channel=event%13
                    if channel==0:pop[3,node]-=1;pop[0,node]+=1
                    elif channel==1:pop[0,node]-=1;pop[1,node]+=1
                    elif channel==2:pop[1,node]-=1;pop[0,node]+=1
                    elif channel==3:pop[4,node]-=1;pop[1,node]-=1;pop[2,node]+=1;pop[5,node]+=1
                    elif channel==4:pop[2,node]-=1;pop[1,node]+=1
                    elif channel==5:pop[1,node]-=1;pop[6,node]+=1
                    elif channel==6:raise ValueError("R1re must be zero at P17")
                    else:
                        species=(channel-7)//2;direction=(channel-7)%2
                        other=(node-1)%n if direction==0 else (node+1)%n
                        pop[species,node]-=1;pop[species,other]+=1
                        hops[2*species+direction,node]+=1
                    if channel<7:reactions[channel,node]+=1
                if np.any(pop<0):raise ValueError("Negative SSA state")
                waiting,event=_network_clock(pop,a0,b0,omega,diff,mask,src,dst,weight,kappa)
                future+=waiting;events+=1
            out[r,k,:7]=pop;out[r,k,7:14]=reactions;out[r,k,14:20]=hops;environment[r,k]=ef
        stats[r,0]=events
    return out,environment,stats

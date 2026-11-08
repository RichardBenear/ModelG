"""MG080 checkpointable ensembles; same held-stream law as MG079 samplers.
Numba seed set once; Python oracle retains one Generator across blocks.
Historical samplers unchanged; new stream functions copy their path loops.
"""
import numpy as np
import numba as nb
from mg11n1 import _network_clock

@nb.njit(cache=True)
def stream_chunk(initial,a0,b0,omega,diff,mask,observations,reps,seed,src,dst,weight,kappa,max_events):
    if seed>=0:np.random.seed(seed)
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


def oracle_stream(case,g,kappa,times,reps,rng,mask,max_events=1000000):
    """Python/Generator independent of production clock, STOICH and update code."""
    n=g['n'];e=len(g['src']);out=np.zeros((reps,len(times),20,n),dtype=np.int64);env=np.zeros((reps,len(times),2,e),dtype=np.int64);stats=np.zeros((reps,3))
    for r in range(reps):
        p=case['initial'].copy();rx=np.zeros((7,n),dtype=np.int64);hop=np.zeros((6,n),dtype=np.int64);ef=np.zeros((2,e),dtype=np.int64);clock=0.;events=0
        def next_event():
            rates=[];labels=[]
            for i in range(n):
                G,X,Y,A,B=map(int,p[:5,i]);values=(14*case['omega']*A/case['a0'][i],float(G),X/10,29*X*B/case['b0'][i],X*(X-1)*Y/case['omega']**2 if X>=2 else 0.,.9*X)
                for c,rate in enumerate(values):rates.append(rate*mask[c]);labels.append(('r',i,c))
                for s in range(3):
                    for direction in range(2):rates.append(case['diff'][s]*p[s,i]);labels.append(('h',i,2*s+direction))
            for s in range(2):
                for edge,(i,j,w) in enumerate(zip(g['src'],g['dst'],g['weight'])):rates.append(kappa[s]*w*p[3+s,i]);labels.append(('e',edge,s))
            total=sum(rates)
            if total==0:return np.inf,None
            wait=rng.exponential(1/total);pick=rng.random()*total;cumulative=0.
            for rate,label in zip(rates,labels):
                cumulative+=rate
                if pick<cumulative:return wait,label
            raise ValueError('Oracle mark selection')
        wait,label=next_event();future=clock+wait
        for k,target in enumerate(times):
            while future<=target:
                if events>=max_events:raise ValueError('Oracle event cap')
                kind,i,c=label
                if kind=='r':
                    if c==0:p[3,i]-=1;p[0,i]+=1
                    elif c==1:p[0,i]-=1;p[1,i]+=1
                    elif c==2:p[1,i]-=1;p[0,i]+=1
                    elif c==3:p[4,i]-=1;p[1,i]-=1;p[2,i]+=1;p[5,i]+=1
                    elif c==4:p[2,i]-=1;p[1,i]+=1
                    elif c==5:p[1,i]-=1;p[6,i]+=1
                    rx[c,i]+=1
                elif kind=='h':
                    s=c//2;other=(i-1)%n if c%2==0 else (i+1)%n;p[s,i]-=1;p[s,other]+=1;hop[c,i]+=1
                else:s=c;p[3+s,g['src'][i]]-=1;p[3+s,g['dst'][i]]+=1;ef[s,i]+=1
                if np.any(p<0):raise ValueError('Oracle negative state')
                events+=1;clock=future;wait,label=next_event();future=clock+wait
            out[r,k,:7]=p;out[r,k,7:14]=rx;out[r,k,14:20]=hop;env[r,k]=ef
        stats[r,0]=events
    return out,env,stats

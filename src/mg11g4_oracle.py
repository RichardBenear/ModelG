"""MG096 independent Python/Generator two-graph exact SSA, literal explicit updates."""
import numpy as np

def direct_ssa(case,ga,gb,kappa,times,reps,seed,mask,max_events=1000000,_rng=None):
    """Python/Generator independent of production clock, STOICH and update code."""
    rng=np.random.default_rng(seed) if _rng is None else _rng;n=ga['n'];e=len(ga['src']);out=np.zeros((reps,len(times),20,n),dtype=np.int64);env=np.zeros((reps,len(times),2,e),dtype=np.int64);stats=np.zeros((reps,3))
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
                for edge,(i,j,w) in enumerate(zip((ga,gb)[s]['src'],(ga,gb)[s]['dst'],(ga,gb)[s]['weight'])):rates.append(kappa[s]*w*p[3+s,i]);labels.append(('e',edge,s))
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
                else:s=c;p[3+s,(ga,gb)[s]['src'][i]]-=1;p[3+s,(ga,gb)[s]['dst'][i]]+=1;ef[s,i]+=1
                if np.any(p<0):raise ValueError('Oracle negative state')
                events+=1;clock=future;wait,label=next_event();future=clock+wait
            out[r,k,:7]=p;out[r,k,7:14]=rx;out[r,k,14:20]=hop;env[r,k]=ef
        stats[r,0]=events
    return out,env,stats

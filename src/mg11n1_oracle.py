"""MG079 independent tiny-network direct SSA and linear matrix reference."""
import numpy as np
from scipy.linalg import expm
def linear_law(case,g,kappa,kind,t):
    n=g['n'];w=np.zeros((n,n))
    for i,j,value in zip(g['src'],g['dst'],g['weight']):w[i,j]=value
    L=w-np.diag(w.sum(axis=1))
    if kind=='transport':
        a=np.zeros(n);b=np.zeros(n);a[np.flatnonzero(case['initial'][3])[0]]=1;b[np.flatnonzero(case['initial'][4])[0]]=1
        return np.outer(a@expm(kappa[0]*L*t),b@expm(kappa[1]*L*t)).ravel()
    if kind=='A_source':
        lam=14*case['omega']/case['a0'];Q=np.zeros((2*n,2*n));Q[:n,:n]=kappa[0]*L-np.diag(lam);Q[np.arange(n),n+np.arange(n)]=lam
        p=np.zeros(2*n);p[np.flatnonzero(case['initial'][3])[0]]=1;return p@expm(Q*t)
    if kind=='B_contact':
        lam=29*case['initial'][1]/case['b0'];Q=np.zeros((n+1,n+1));Q[:n,:n]=kappa[1]*L-np.diag(lam);Q[:n,n]=lam
        p=np.zeros(n+1);p[np.flatnonzero(case['initial'][4])[0]]=1;return p@expm(Q*t)
    raise ValueError(kind)
def categories(data,kind):
    pop=data[:,:,:7];n=pop.shape[-1]
    if kind=='transport':
        assert np.all(pop[:,:,3].sum(axis=-1)==1) and np.all(pop[:,:,4].sum(axis=-1)==1)
        return pop[:,:,3].argmax(axis=-1)*n+pop[:,:,4].argmax(axis=-1)
    if kind=='A_source':
        assert np.all(pop[:,:,3].sum(axis=-1)+pop[:,:,0].sum(axis=-1)==1)
        return np.where(pop[:,:,3].sum(axis=-1)>0,pop[:,:,3].argmax(axis=-1),n+pop[:,:,0].argmax(axis=-1))
    if kind=='B_contact':
        assert np.all(pop[:,:,4].sum(axis=-1)+data[:,:,10].sum(axis=-1)==1)
        return np.where(pop[:,:,4].sum(axis=-1)>0,pop[:,:,4].argmax(axis=-1),n)
    raise ValueError(kind)
def direct_ssa(case,g,kappa,times,reps,seed,mask,max_events=1000000):
    """Python/Generator independent of production clock, STOICH and update code."""
    rng=np.random.default_rng(seed);n=g['n'];e=len(g['src']);out=np.zeros((reps,len(times),20,n),dtype=np.int64);env=np.zeros((reps,len(times),2,e),dtype=np.int64);stats=np.zeros((reps,3))
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

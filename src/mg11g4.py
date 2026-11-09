"""MG096 bounded synchronous network leap; finite P17 resources and two graphs.

All draws read one snapshot; reserve joint B/X conversion before remaining
conditional exits. Products/arrivals cannot act in the same numerical step.
Continuous-time approximation, not a fundamental-tick or token-lineage law.
"""
import numpy as np
import numba as nb
from mg10s1 import _exits4,_leap as legacy_leap,_ssa_clock as legacy_clock
from mg11n1 import validate as legacy_validate
from mg10d2_source import finite_source
MAX_COUNT=1000000000
MAX_TOTAL=1000000000000
MAX_STEPS=2000000

def validate(case,ga,gb,k,mask):
    if ga['n']!=gb['n'] or len(ga['src'])!=len(gb['src']):raise ValueError('Equal edge-array sizes required')
    k=np.asarray(k,dtype=float);mask=np.asarray(mask,dtype=float)
    legacy_validate(case,ga,k,mask);legacy_validate(case,gb,k,mask)
    if np.max(case['initial'])>MAX_COUNT or sum(map(int,case['initial'].ravel()))>MAX_TOTAL:raise ValueError('Outside qualified count domain')

@nb.njit(cache=True)
def _hop(available,probability,species,node,src,dst,weight,kappa,degree,offsets,edge_ids,delta,ef):
    if probability==0 or available==0:return
    left=np.random.binomial(available,probability);remaining_weight=degree[species,node]
    for slot in range(offsets[species,node],offsets[species,node+1]):
        edge=edge_ids[species,slot]
        w=kappa[species]*weight[species,edge]
        if w==0:continue
        count=left if slot==offsets[species,node+1]-1 else np.random.binomial(left,min(1.,w/remaining_weight))
        left-=count;remaining_weight-=w;delta[3+species,node]-=count;delta[3+species,dst[species,edge]]+=count;ef[species,edge]+=count
    if left!=0:raise ValueError('Environmental multinomial allocation failed')

@nb.njit(cache=True)
def step(pop,rx,hops,ef,a0,b0,omega,diff,mask,requested,src,dst,weight,kappa,degree,offsets,edge_ids):
    if np.max(pop)>MAX_COUNT:raise ValueError('Path left qualified per-species count domain')
    if len(src[0])==0 or (kappa[0]==0 and kappa[1]==0):return legacy_leap(pop,rx,hops,a0,b0,omega,diff,mask,requested)
    n=pop.shape[1];maximum=0.
    for i in range(n):
        x,b=pop[1,i],pop[4,i];rate3=mask[3]*29.*x*b/b0[i];joint=rate3/min(x,b) if min(x,b)>0 else 0.;r4=mask[4]*x*(x-1)/(omega*omega) if x>=2 else 0.
        maximum=max(maximum,mask[0]*omega*14/a0[i]+degree[0,i],mask[1]+2*diff[0],mask[2]*.1+mask[5]*.9+mask[3]*29*b/b0[i]+2*diff[1],r4+2*diff[2],joint,(rate3/b if b>0 else 0.)+degree[1,i])
    if not np.isfinite(maximum):raise ValueError('Nonfinite leap hazard')
    dt=min(requested,.2/maximum) if maximum>0 else requested;delta=np.zeros_like(pop);largest=0.
    for i in range(n):
        left,right=(i-1)%n,(i+1)%n;g,x,y,a,b=pop[0,i],pop[1,i],pop[2,i],pop[3,i],pop[4,i]
        pa=mask[0]*omega*14*dt/a0[i];source=finite_source(a,pa)
        gx,gl,gr,_=_exits4(g,mask[1]*dt,diff[0]*dt,diff[0]*dt,0.)
        rate3=mask[3]*29*x*b/b0[i];m=min(x,b);joint=rate3*dt/m if m else 0.;conversion=np.random.binomial(m,joint)
        p3x=rate3*dt/x if x else 0.;p3b=rate3*dt/b if b else 0.;denominator=1-p3x
        pg=mask[2]*.1*dt/denominator;pw=mask[5]*.9*dt/denominator;ph=diff[1]*dt/denominator
        xg,xw,xl,xr=_exits4(x-conversion,pg,pw,ph,ph)
        py=mask[4]*x*(x-1)*dt/(omega*omega) if x>=2 else 0.;yx,yl,yr,_=_exits4(y,py,diff[2]*dt,diff[2]*dt,0.)
        delta[0,i]+=source-gx-gl-gr+xg;delta[1,i]+=gx-xg-xw-xl-xr-conversion+yx;delta[2,i]+=conversion-yx-yl-yr;delta[3,i]-=source;delta[4,i]-=conversion;delta[5,i]+=conversion;delta[6,i]+=xw
        delta[0,left]+=gl;delta[0,right]+=gr;delta[1,left]+=xl;delta[1,right]+=xr;delta[2,left]+=yl;delta[2,right]+=yr
        rx[0,i]+=source;rx[1,i]+=gx;rx[2,i]+=xg;rx[3,i]+=conversion;rx[4,i]+=yx;rx[5,i]+=xw
        hops[0,i]+=gl;hops[1,i]+=gr;hops[2,i]+=xl;hops[3,i]+=xr;hops[4,i]+=yl;hops[5,i]+=yr
        ah=degree[0,i]*dt/(1-pa);bh=degree[1,i]*dt/(1-p3b)
        _hop(a-source,ah,0,i,src,dst,weight,kappa,degree,offsets,edge_ids,delta,ef);_hop(b-conversion,bh,1,i,src,dst,weight,kappa,degree,offsets,edge_ids,delta,ef)
        largest=max(largest,pa+degree[0,i]*dt,(mask[1]+2*diff[0])*dt,py+2*diff[2]*dt,joint,pg+pw+2*ph,ah,bh)
    pop+=delta
    if np.any(pop<0) or largest>.25+1e-14:raise ValueError('Negative state or leap probability bound')
    return dt,largest

@nb.njit(cache=True)
def leap_chunk(initial,a0,b0,omega,diff,mask,tau,observations,reps,seed,src,dst,weight,kappa,degree,offsets,edge_ids,max_steps):
    if seed>=0:np.random.seed(seed)
    n=initial.shape[1];e=src.shape[1];out=np.zeros((reps,len(observations),20,n),np.int64);env=np.zeros((reps,len(observations),2,e),np.int64);stats=np.zeros((reps,3))
    for r in range(reps):
        p=initial.copy();rx=np.zeros((7,n),np.int64);hop=np.zeros((6,n),np.int64);ef=np.zeros((2,e),np.int64);t=0.;leaps=0;adapted=0;largest=0.
        for j,target in enumerate(observations):
            while t<target:
                if leaps>=max_steps:raise ValueError('Per-path leap cap')
                requested=min(tau,target-t);dt,bound=step(p,rx,hop,ef,a0,b0,omega,diff,mask,requested,src,dst,weight,kappa,degree,offsets,edge_ids)
                if dt<=0:raise ValueError('Nonpositive numerical step')
                t=target if dt==target-t else t+dt;leaps+=1;adapted+=int(dt<requested*(1-1e-12));largest=max(largest,bound)
            out[r,j,:7]=p;out[r,j,7:14]=rx;out[r,j,14:20]=hop;env[r,j]=ef
        stats[r,0]=leaps;stats[r,1]=adapted;stats[r,2]=largest
    return out,env,stats

def arrays(ga,gb,k):
    src=np.stack([ga['src'],gb['src']]);dst=np.stack([ga['dst'],gb['dst']]);w=np.stack([ga['weight'],gb['weight']]);degree=np.zeros((2,ga['n']))
    for s in range(2):np.add.at(degree[s],src[s],k[s]*w[s])
    offsets=np.zeros((2,ga['n']+1),np.int64);ids=np.empty_like(src)
    for s in range(2):
        ids[s]=np.argsort(src[s],kind='stable');offsets[s,1:]=np.cumsum(np.bincount(src[s],minlength=ga['n']))
    return src,dst,w,degree,offsets,ids

def simulate(case,ga,gb,k,times,reps,seed,tau=0.,mask=None,max_events=1000000,max_steps=MAX_STEPS):
    k=np.asarray(k,dtype=float);mask=np.ones(7) if mask is None else np.asarray(mask,dtype=float);times=np.asarray(times,dtype=float);validate(case,ga,gb,k,mask)
    if times.ndim!=1 or not len(times) or times[0]<0 or np.any(~np.isfinite(times)) or np.any(np.diff(times)<=0) or not isinstance(reps,(int,np.integer)) or reps<1 or not isinstance(seed,(int,np.integer)) or not 0<=seed<2**32 or not np.isfinite(tau) or tau<0 or not isinstance(max_steps,(int,np.integer)) or not 1<=max_steps<=MAX_STEPS or not isinstance(max_events,(int,np.integer)) or max_events<1:raise ValueError('Invalid execution bounds')
    src,dst,w,degree,offsets,edge_ids=arrays(ga,gb,k)
    if np.any(~np.isfinite(degree)):raise ValueError('Nonfinite graph exit hazard')
    if tau==0:return exact_chunk(case['initial'],case['a0'],case['b0'],case['omega'],case['diff'],mask,times,reps,seed,src,dst,w,k,max_events)
    return leap_chunk(case['initial'],case['a0'],case['b0'],case['omega'],case['diff'],mask,tau,times,reps,seed,src,dst,w,k,degree,offsets,edge_ids,max_steps)

@nb.njit(cache=True)
def single_steps(initial,a0,b0,omega,diff,mask,dt,reps,seed,src,dst,w,k,degree,offsets,edge_ids):
    np.random.seed(seed);n=initial.shape[1];e=src.shape[1];data=np.zeros((reps,2,20,n),np.int64);env=np.zeros((reps,2,2,e),np.int64);stats=np.zeros((reps,3))
    for r in range(reps):
        p=initial.copy();rx=np.zeros((7,n),np.int64);hop=np.zeros((6,n),np.int64);ef=np.zeros((2,e),np.int64);data[r,0,:7]=p;actual,bound=step(p,rx,hop,ef,a0,b0,omega,diff,mask,dt,src,dst,w,k,degree,offsets,edge_ids)
        if actual!=dt:raise ValueError('Single-step fixture unexpectedly adaptive')
        data[r,1,:7]=p;data[r,1,7:14]=rx;data[r,1,14:20]=hop;env[r,1]=ef;stats[r,0]=1;stats[r,2]=bound
    return data,env,stats


@nb.njit(cache=True)
def _exact_clock(pop,a0,b0,omega,diff,mask,src,dst,weight,kappa):
    if src.shape[1]==0 or (kappa[0]==0 and kappa[1]==0):return legacy_clock(pop,a0,b0,omega,diff,mask)
    # Separate direct propensity implementation; no leap helpers or factorization.
    n=pop.shape[1];rates=np.zeros(13*n+2*src.shape[1])
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
        for edge in range(src.shape[1]):rates[13*n+s*src.shape[1]+edge]=kappa[s]*weight[s,edge]*pop[3+s,src[s,edge]]
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
def exact_chunk(initial,a0,b0,omega,diff,mask,observations,reps,seed,src,dst,weight,kappa,max_events):
    if seed>=0:np.random.seed(seed)
    n=initial.shape[1];out=np.zeros((reps,len(observations),20,n),dtype=np.int64);stats=np.zeros((reps,3));environment=np.zeros((reps,len(observations),2,src.shape[1]),dtype=np.int64)
    for r in range(reps):
        pop=initial.copy();reactions=np.zeros((7,n),dtype=np.int64);hops=np.zeros((6,n),dtype=np.int64)
        ef=np.zeros((2,src.shape[1]),dtype=np.int64)
        future,event=_exact_clock(pop,a0,b0,omega,diff,mask,src,dst,weight,kappa);events=0
        for k,target in enumerate(observations):
            while future<=target:
                if events>=max_events:raise ValueError("Per-trajectory exact-event budget reached")
                if event>=13*n:
                    q=event-13*n;species=q//src.shape[1];edge=q%src.shape[1]
                    pop[3+species,src[species,edge]]-=1;pop[3+species,dst[species,edge]]+=1;ef[species,edge]+=1
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
                waiting,event=_exact_clock(pop,a0,b0,omega,diff,mask,src,dst,weight,kappa)
                future+=waiting;events+=1
            out[r,k,:7]=pop;out[r,k,7:14]=reactions;out[r,k,14:20]=hops;environment[r,k]=ef
        stats[r,0]=events
    return out,environment,stats

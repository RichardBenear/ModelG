"""MG135 compact strict-tick prototype. Frozen qualified engines unchanged."""
from dataclasses import dataclass
import math,time
import numpy as np

@dataclass(frozen=True)
class Channel:
    rule:str
    reserve:tuple # (flattened coordinate, count)
    change:tuple
    coefficient:float
    factors:tuple # (coordinate, falling-factorial order)
    def rate(self,flat):
        value=self.coefficient
        for coordinate,order in self.factors:
            n=int(flat[coordinate])
            for k in range(order):value*=max(0,n-k)
        if not math.isfinite(value) or value<0:raise ValueError('Invalid hazard')
        return value

class Model:
    def __init__(self,n,config):
        omega,refs,links,diff,coupling,mask=config
        if type(n) is not int or n<1 or len(refs)!=n or len(mask)!=6 or any(type(v) not in (bool,np.bool_) for v in mask) or not math.isfinite(omega) or omega<=0:raise ValueError('Invalid configuration')
        if len(diff)!=3 or len(coupling)!=2 or any(not math.isfinite(v) or v<0 for v in tuple(diff)+tuple(coupling)) or any(len(v)!=2 or any(not math.isfinite(x) or x<=0 for x in v) for v in refs):raise ValueError('Invalid coefficients')
        channels=[]
        def add(rule,reserve,change,coefficient,factors):
            if coefficient>0:channels.append(Channel(rule,tuple(reserve),tuple(change),float(coefficient),tuple(factors)))
        for i in range(n):
            g,x,y,a,b,z,w=[s*n+i for s in range(7)];a0,b0=refs[i]
            definitions=[('R1',[(a,1)],[(a,-1),(g,1)],14*omega/a0,[(a,1)]),('R2',[(g,1)],[(g,-1),(x,1)],1.,[(g,1)]),('R2r',[(x,1)],[(x,-1),(g,1)],.1,[(x,1)]),('R3',[(b,1),(x,1)],[(b,-1),(x,-1),(y,1),(z,1)],29/b0,[(b,1),(x,1)]),('R4',[(x,2),(y,1)],[(y,-1),(x,1)],1/omega**2,[(x,2),(y,1)]),('R5',[(x,1)],[(x,-1),(w,1)],.9,[(x,1)])]
            for on,d in zip(mask,definitions):
                if on:add(*d)
        for kind,i,j,weight in links:
            if kind not in ('ActiveLink','ALink','BLink') or not 0<=i<n or not 0<=j<n or i==j or not math.isfinite(weight) or weight<=0:raise ValueError('Invalid link')
            for s in (range(3) if kind=='ActiveLink' else [3 if kind=='ALink' else 4]):
                rate=diff[s] if s<3 else coupling[s-3];source=s*n+i;dest=s*n+j;add('H_'+('G','X','Y','A','B')[s],[(source,1)],[(source,-1),(dest,1)],rate*weight,[(source,1)])
        self.n=n;self.channels=tuple(channels);deps=[[] for _ in range(7*n)]
        for index,c in enumerate(channels):
            for coordinate,_ in c.factors:deps[coordinate].append(index)
        self.dependencies=tuple(tuple(v) for v in deps)
    def validate(self,pop):
        a=np.asarray(pop)
        if a.dtype.kind!='i' or a.shape!=(7,self.n) or np.any(a<0):raise ValueError('Invalid signed integer counts')
        if sum(map(int,a.ravel()))>np.iinfo(np.int64).max:raise ValueError('Population overflow')
        return a.astype(np.int64,copy=True).ravel()
    def workspace(self,pop):return Workspace(self,self.validate(pop))

class Tree:
    def __init__(self,rates):
        self.values=np.asarray(rates,float).copy();self.tree=np.zeros(len(rates)+1)
        for i,value in enumerate(self.values):
            j=i+1
            while j<len(self.tree):self.tree[j]+=value;j+=j&-j
    def total(self):
        value=0.;i=len(self.values)
        while i:value+=self.tree[i];i-=i&-i
        return value
    def update(self,index,value):
        delta=value-self.values[index];self.values[index]=value;i=index+1
        while i<len(self.tree):self.tree[i]+=delta;i+=i&-i
    def select(self,u):
        total=self.total()
        if not 0<=u<1 or total<=0:raise ValueError('Invalid selection')
        target=u*total;index=0;bit=1<<(len(self.values).bit_length()-1)
        while bit:
            candidate=index+bit
            if candidate<len(self.tree) and self.tree[candidate]<=target:target-=self.tree[candidate];index=candidate
            bit>>=1
        if index>=len(self.values) or self.values[index]<=0:raise ValueError('Numerical tree selection')
        return index

class Workspace:
    def __init__(self,model,flat):
        self.model=model;self.available=flat.copy();self.output=flat.copy();self.tree=Tree([c.rate(flat) for c in model.channels]);self.positive=int(np.count_nonzero(self.tree.values))
    def accept(self,index):
        channel=self.model.channels[index]
        if self.tree.values[index]<=0 or any(self.available[i]<count for i,count in channel.reserve):raise ValueError('Unavailable event')
        touched=set()
        for i,count in channel.reserve:self.available[i]-=count;touched.update(self.model.dependencies[i])
        for i,change in channel.change:self.output[i]+=change
        for i in touched:
            value=self.model.channels[i].rate(self.available);self.positive+=int(value>0)-int(self.tree.values[i]>0);self.tree.update(i,value)
        if self.positive==0:self.tree=Tree(np.zeros(len(self.tree.values)))
        return channel

def tick(model,pop,delta,rng,deadline=math.inf,max_events=1000000,trace=True):
    if not math.isfinite(delta) or delta<=0:raise ValueError('Invalid tick')
    work=model.workspace(pop);clock=0.;events=[];count=0
    while work.positive:
        if time.perf_counter()>deadline:raise TimeoutError('Sparse tick deadline')
        total=work.tree.total()
        if not math.isfinite(total) or total<=0:raise ValueError('Numerical total')
        wait=float(rng.exponential(1/total))
        if not math.isfinite(wait) or wait<=0:raise ValueError('Invalid mark')
        clock+=wait
        if clock>=delta:break
        if count>=max_events:raise ValueError('Event cap')
        index=work.tree.select(float(rng.random()));c=work.accept(index);count+=1
        if trace:events.append(dict(time=clock,channel=index,rule=c.rule,reserve=c.reserve,change=c.change))
    return dict(counts=work.output.reshape(7,model.n),available=work.available.reshape(7,model.n),events=count,trace=events)

"""MG106 explicit typed-incidence reference engine, finite fixed substrate.

Identity-rich matches and event traces; no aggregate reaction kernel calls.
Continuous-time SSA and frozen-proposal strict ticks are distinct executors.
"""
from dataclasses import dataclass,replace
from itertools import combinations
import math
import numpy as np
SPECIES=('G','X','Y','A','B','Z','W');RULES=('R1','R2','R2r','R3','R4','R5')

@dataclass(frozen=True)
class Relation:
    kind:str
    roles:tuple
    weight:float=1.

@dataclass(frozen=True)
class State:
    sites:tuple
    relations:tuple
    omega:float
    references:tuple # (site,original_a_count,original_b_count)
    diffusion:tuple=(1.,1.,12.)
    coupling:tuple=(1/3,1/3)
    enabled:tuple=(True,)*6
    def validate(self):
        if not self.sites or len(set(self.sites))!=len(self.sites):raise ValueError('Duplicate/empty sites')
        if not math.isfinite(self.omega) or self.omega<=0 or len(self.diffusion)!=3 or len(self.coupling)!=2 or any(not math.isfinite(x) or x<0 for x in self.diffusion+self.coupling):raise ValueError('Invalid rates')
        if len(self.enabled)!=6 or any(type(x) is not bool for x in self.enabled):raise ValueError('Invalid mask')
        refs={r[0]:r[1:] for r in self.references}
        if len(refs)!=len(self.references) or set(refs)!=set(self.sites) or any(len(r)!=3 or any(not math.isfinite(x) or x<=0 for x in r[1:]) for r in self.references):raise ValueError('Invalid fixed references')
        occupied=set();links=set()
        for r in self.relations:
            if not math.isfinite(r.weight) or r.weight<=0:raise ValueError('Invalid relation weight')
            if r.kind=='Occupies':
                if len(r.roles)!=3:raise ValueError('Occupies arity')
                token,site,species=r.roles
                if token in occupied or token in self.sites or site not in self.sites or species not in SPECIES or r.weight!=1:raise ValueError('Invalid occupation')
                occupied.add(token)
            elif r.kind in ('ActiveLink','ALink','BLink'):
                if len(r.roles)!=2 or r.roles[0] not in self.sites or r.roles[1] not in self.sites or r.roles[0]==r.roles[1] or (r.kind,r.roles) in links:raise ValueError('Invalid/duplicate link')
                links.add((r.kind,r.roles))
            else:raise ValueError('Unknown relation sort')
        return self
    def occupations(self):return {r.roles[0]:(r.roles[1],r.roles[2]) for r in self.relations if r.kind=='Occupies'}
    def counts(self):
        pop=np.zeros((7,len(self.sites)),dtype=np.int64);index={s:i for i,s in enumerate(self.sites)}
        for site,species in self.occupations().values():pop[SPECIES.index(species),index[site]]+=1
        return pop

@dataclass(frozen=True)
class Event:
    rule:str
    participants:frozenset
    before:tuple # (token,site,species) for ALL participants including catalysts
    after:tuple # replacement Occupies relations for writable participants
    rate:float
    def key(self):return self.rule,self.participants,frozenset(self.before),frozenset(self.after)

def matches(state):
    state.validate();occ=state.occupations();out=[];refs={r[0]:r[1:] for r in state.references}
    def add(rule,tokens,changes,rate):
        if rate>0:out.append(Event(rule,frozenset(tokens),tuple((t,*occ[t]) for t in tokens),tuple((t,*changes[t]) for t in changes),float(rate)))
    for site in state.sites:
        by={s:[t for t,v in occ.items() if v==(site,s)] for s in SPECIES};a0,b0=refs[site]
        for s,rule,target,rate,enabled in [('A','R1','G',14*state.omega/a0,state.enabled[0]),('G','R2','X',1.,state.enabled[1]),('X','R2r','G',.1,state.enabled[2]),('X','R5','W',.9,state.enabled[5])]:
            if enabled:
                for t in by[s]:add(rule,[t],{t:(site,target)},rate)
        if state.enabled[3]:
            for b in by['B']:
                for x in by['X']:add('R3',[b,x],{b:(site,'Z'),x:(site,'Y')},29/b0)
        if state.enabled[4]:
            for x1,x2 in combinations(by['X'],2):
                for y in by['Y']:add('R4',[x1,x2,y],{y:(site,'X')},2/state.omega**2)
    for link in state.relations:
        if link.kind=='Occupies':continue
        source,dest=link.roles
        for token,(site,species) in occ.items():
            if site!=source:continue
            if link.kind=='ActiveLink' and species in SPECIES[:3]:rate=state.diffusion[SPECIES.index(species)]*link.weight
            elif link.kind=='ALink' and species=='A':rate=state.coupling[0]*link.weight
            elif link.kind=='BLink' and species=='B':rate=state.coupling[1]*link.weight
            else:continue
            add('H_'+species,[token],{token:(dest,species)},rate)
    if any(not math.isfinite(e.rate) for e in out):raise ValueError('Nonfinite propensity')
    return out

def commit(state,events):
    occ=state.occupations();reserved=set();new=dict(occ);enabled={e.key():e for e in matches(state)}
    for event in events:
        if event.key() not in enabled or event!=enabled[event.key()] or event.rate<=0 or event.participants & reserved:raise ValueError('Stale/invalid/conflicting event')
        if set(t for t,_,_ in event.before)!=set(event.participants) or not set(t for t,_,_ in event.after)<=event.participants:raise ValueError('Invalid event roles')
        for token,site,species in event.before:
            if occ.get(token)!=(site,species):raise ValueError('Stale event')
        reserved.update(event.participants)
        for token,site,species in event.after:new[token]=(site,species)
    relations=tuple(Relation('Occupies',(t,*v)) for t,v in new.items())+tuple(r for r in state.relations if r.kind!='Occupies');result=replace(state,relations=relations).validate()
    if set(result.occupations())!=set(occ):raise ValueError('Token identity accounting')
    return result

def select(state,proposals,delta):
    if not math.isfinite(delta) or delta<=0:raise ValueError('Invalid tick')
    enabled={e.key():e for e in matches(state)};seen=set();marked=[]
    for event,time in proposals:
        if event.key() not in enabled or event!=enabled[event.key()] or event.key() in seen or not math.isfinite(time) or time<0:raise ValueError('Invalid/duplicate proposal')
        seen.add(event.key())
        if time<delta:marked.append((time,event))
    times=[t for t,_ in marked]
    if len(set(times))!=len(times):raise ValueError('Exact proposal-time tie')
    accepted=[];trace=[];owners={}
    for time,event in sorted(marked,key=lambda pair:pair[0]):
        blockers=tuple(sorted({owners[t] for t in event.participants if t in owners}))
        ok=not blockers;index=len(trace);trace.append(dict(time=float(time),rule=event.rule,participants=tuple(event.participants),before=event.before,after=event.after,accepted=ok,blocked_by=blockers))
        if ok:
            accepted.append(event)
            for token in event.participants:owners[token]=index
    return commit(state,accepted),trace

def tick(state,delta,rng):
    events=matches(state);proposals=[(event,float(rng.exponential(1/event.rate))) for event in events];return select(state,proposals,delta)

def simulate(initial,times,seed,delta=None,max_steps=100000,deadline=math.inf):
    import time
    times=np.asarray(times,dtype=float)
    if times.ndim!=1 or len(times)<2 or times[0]!=0 or np.any(~np.isfinite(times)) or np.any(np.diff(times)<=0) or type(seed) is not int or not 0<=seed<2**32 or type(max_steps) is not int or max_steps<1 or (delta is not None and (not math.isfinite(delta) or delta<=0)):raise ValueError('Invalid simulation bounds')
    state=initial.validate();rng=np.random.default_rng(seed);clock=0.;steps=0;trace=[];saved=[state.counts()];frames=[state];parents={};boundary=1
    while boundary<len(times):
        if time.perf_counter()>deadline:raise TimeoutError('Development deadline')
        if steps>=max_steps:raise ValueError('Event/tick cap')
        if delta is None:
            events=matches(state);total=sum(e.rate for e in events)
            if not total:
                while boundary<len(times):saved.append(state.counts());frames.append(state);boundary+=1
                break
            wait=float(rng.exponential(1/total));future=clock+wait
            while boundary<len(times) and times[boundary]<future:saved.append(state.counts());frames.append(state);boundary+=1
            if boundary==len(times):break
            pick=rng.random()*total;acc=0.;chosen=None
            for event in events:
                acc+=event.rate
                if pick<acc:chosen=event;break
            if chosen is None:raise ValueError('Exact event selection failed')
            clock=future;state=commit(state,[chosen]);events_trace=[dict(time=clock,rule=chosen.rule,participants=tuple(chosen.participants),before=chosen.before,after=chosen.after,accepted=True,blocked_by=())]
        else:
            target=min(clock+delta,float(times[boundary]));state,events_trace=tick(state,target-clock,rng)
            for row in events_trace:row['time']+=clock
            clock=target
        step_start=len(trace)
        for row in events_trace:
            row.update(step=steps,parents=tuple(sorted({parents[t] for t in row['participants'] if t in parents})))
            row['blocked_by']=tuple(step_start+x for x in row['blocked_by'])
            if row['accepted']:
                # Tick parents refer to pre-tick history; same-tick events share no tokens.
                for token in row['participants']:parents[token]=len(trace)
            trace.append(row)
        steps+=1
        while boundary<len(times) and times[boundary]<=clock:saved.append(state.counts());frames.append(state);boundary+=1
    return dict(counts=np.asarray(saved),frames=frames,trace=trace,steps=steps)

def relabel(state,mapping):
    if set(mapping)!=set(state.sites)|set(state.occupations()) or len(set(mapping.values()))!=len(mapping):raise ValueError('Invalid relabeling')
    relations=[]
    for r in state.relations:
        roles=(mapping[r.roles[0]],mapping[r.roles[1]],r.roles[2]) if r.kind=='Occupies' else tuple(mapping[x] for x in r.roles);relations.append(Relation(r.kind,roles,r.weight))
    return replace(state,sites=tuple(mapping[s] for s in state.sites),relations=tuple(relations),references=tuple((mapping[s],a,b) for s,a,b in state.references)).validate()

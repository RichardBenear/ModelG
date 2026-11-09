"""MG110 fixed multi-incident cells; explicit channel provenance, H1 frozen."""
from dataclasses import dataclass,replace
import math,time
import numpy as np
from mg12h1 import State,Relation,Event,matches as old_matches,SPECIES
CELLS=('ActiveCell','ACell','BCell')
@dataclass(frozen=True)
class HyperState(State):
    def validate(self):
        plain=replace(self,relations=tuple(r for r in self.relations if r.kind not in CELLS));State.validate(plain);seen=set()
        for r in self.relations:
            if r.kind not in CELLS:continue
            key=r.kind,frozenset(r.roles)
            if len(r.roles)<2 or len(set(r.roles))!=len(r.roles) or not set(r.roles)<=set(self.sites) or key in seen or not math.isfinite(r.weight) or r.weight<=0:raise ValueError('Invalid/duplicate cell')
            seen.add(key)
        return self
@dataclass(frozen=True)
class CellEvent(Event):
    cell:tuple=()
    def key(self):return super().key()+(self.cell[0],frozenset(self.cell[1]))
def matches(state):
    state.validate();out=old_matches(replace(state,relations=tuple(r for r in state.relations if r.kind not in CELLS)));occ=state.occupations()
    for cell in state.relations:
        if cell.kind not in CELLS:continue
        for token,(site,sp) in occ.items():
            if site not in cell.roles:continue
            if cell.kind=='ActiveCell' and sp in SPECIES[:3]:rate=state.diffusion[SPECIES.index(sp)]
            elif cell.kind=='ACell' and sp=='A':rate=state.coupling[0]
            elif cell.kind=='BCell' and sp=='B':rate=state.coupling[1]
            else:continue
            if rate==0:continue
            for dest in cell.roles:
                if dest!=site:out.append(CellEvent('H_'+sp,frozenset([token]),((token,site,sp),),((token,dest,sp),),rate*cell.weight/(len(cell.roles)-1),(cell.kind,cell.roles)))
    if len({e.key() for e in out})!=len(out):raise ValueError('Duplicate event channel key')
    return out
def commit(state,events):
    available={e.key():e for e in matches(state)};occ=state.occupations();reserved=set()
    for e in events:
        if e.key() not in available or available[e.key()]!=e or e.participants & reserved:raise ValueError('Invalid/stale/conflicting match')
        reserved.update(e.participants)
        for token,site,sp in e.after:occ[token]=(site,sp)
    out=replace(state,relations=tuple(Relation('Occupies',(t,*v)) for t,v in occ.items())+tuple(r for r in state.relations if r.kind!='Occupies')).validate();assert set(out.occupations())==set(state.occupations());return out
def select(state,proposals,delta):
    if not math.isfinite(delta) or delta<=0:raise ValueError('Invalid tick')
    available={e.key():e for e in matches(state)};seen=set();marked=[]
    for e,t in proposals:
        if e.key() not in available or available[e.key()]!=e or e.key() in seen or not math.isfinite(t) or t<0:raise ValueError('Invalid proposal')
        seen.add(e.key())
        if t<delta:marked.append((t,e))
    if len({t for t,_ in marked})!=len(marked):raise ValueError('Mark tie')
    used=set();events=[];trace=[]
    for t,e in sorted(marked,key=lambda pair:pair[0]):
        ok=not bool(used & e.participants);trace.append(dict(time=float(t),rule=e.rule,before=e.before,after=e.after,participants=tuple(e.participants),cell=e.cell if isinstance(e,CellEvent) else None,accepted=ok))
        if ok:used.update(e.participants);events.append(e)
    return commit(state,events),trace
def simulate(initial,times,seed,delta=None,deadline=math.inf,max_steps=100000):
    times=np.asarray(times,float)
    if type(seed) is not int or not 0<=seed<2**32 or times.ndim!=1 or len(times)<2 or times[0]!=0 or np.any(~np.isfinite(times)) or np.any(np.diff(times)<=0) or (delta is not None and (not math.isfinite(delta) or delta<=0)) or type(max_steps) is not int or max_steps<1:raise ValueError('Invalid execution bounds')
    state=initial.validate();rng=np.random.default_rng(seed);clock=0.;frame=1;steps=0;saved=[state.counts()];trace=[]
    while frame<len(times):
        if time.perf_counter()>deadline:raise TimeoutError('H2 development deadline')
        if steps>=max_steps:raise ValueError('Event/tick cap')
        events=matches(state)
        if delta is None:
            total=sum(e.rate for e in events)
            if not total:
                while frame<len(times):saved.append(state.counts());frame+=1
                break
            future=clock+rng.exponential(1/total)
            while frame<len(times) and times[frame]<future:saved.append(state.counts());frame+=1
            if frame==len(times):break
            pick=rng.random()*total;acc=0.;chosen=None
            for e in events:
                acc+=e.rate
                if pick<acc:chosen=e;break
            if chosen is None:raise ValueError('Mark selection')
            clock=float(future);state=commit(state,[chosen]);rows=[dict(time=clock,rule=chosen.rule,before=chosen.before,after=chosen.after,participants=tuple(chosen.participants),cell=chosen.cell if isinstance(chosen,CellEvent) else None,accepted=True)]
        else:
            target=min(clock+delta,float(times[frame]));state,rows=select(state,[(e,float(rng.exponential(1/e.rate))) for e in events],target-clock)
            for row in rows:row['time']+=clock
            clock=target
        for row in rows:row['step']=steps
        trace.extend(rows);steps+=1
        while frame<len(times) and times[frame]<=clock:saved.append(state.counts());frame+=1
    return dict(counts=np.asarray(saved),trace=trace,final=state)

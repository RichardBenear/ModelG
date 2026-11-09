"""MG105 static typed-match contract checks; no RNG or time integration."""
from itertools import combinations,permutations
from pathlib import Path
from fractions import Fraction as F
import json,hashlib
ROOT=Path(__file__).resolve().parents[1]
SPECIES=('G','X','Y','A','B','Z','W')
def matches(state):
    by={s:[i for i,v in state.items() if v==s] for s in SPECIES};out=[]
    def add(rule,participants,changes,rate):out.append((rule,frozenset(participants),changes,rate))
    for s,rule,target,rate in [('A','R1','G',F(14,5)),('G','R2','X',F(1)),('X','R2r','G',F(1,10)),('X','R5','W',F(9,10))]:
        for i in by[s]:add(rule,[i],{i:target},rate)
    for b in by['B']:
        for x in by['X']:add('R3',[b,x],{b:'Z',x:'Y'},F(29,4))
    for x1,x2 in combinations(by['X'],2):
        for y in by['Y']:add('R4',[x1,x2,y],{y:'X'},F(1,2))
    return out
def apply(state,event):
    assert set(event[2])<=event[1];new=dict(state);new.update(event[2]);return new
def normalize(event):return event[0],event[1],frozenset(event[2].items()),event[3]
def main():
    # omega=2,a0=10,b0=4; at most seven explicit tokens, one local site.
    states=[dict(enumerate(values)) for values in [('A','G','X','X','Y','B','W'),('X','Y'),('X','X','X','Y','Y'),('Z','W')]];count=0;relabels=0
    increments={'R1':(1,0,0,-1,0,0,0),'R2':(-1,1,0,0,0,0,0),'R2r':(1,-1,0,0,0,0,0),'R3':(0,-1,1,0,-1,1,0),'R4':(0,1,-1,0,0,0,0),'R5':(0,-1,0,0,0,0,1)}
    for state in states:
        events=matches(state);n={s:list(state.values()).count(s) for s in SPECIES};expected={'R1':F(14*n['A'],5),'R2':F(n['G']),'R2r':F(n['X'],10),'R3':F(29*n['B']*n['X'],4),'R4':F(n['X']*(n['X']-1)*n['Y'],4),'R5':F(9*n['X'],10)}
        for rule,total in expected.items():assert sum((e[3] for e in events if e[0]==rule),F(0))==total;count+=1
        for event in events:
            new=apply(state,event);assert set(new)==set(state);assert tuple(list(new.values()).count(s)-n[s] for s in SPECIES)==increments[event[0]];count+=1
            if event[0]=='R4':assert len(event[1])==3 and len(event[2])==1;count+=1
        for a,b in combinations(events,2):
            if a[1].isdisjoint(b[1]):assert apply(apply(state,a),b)==apply(apply(state,b),a);count+=1
        # Exhaust all 7! identifier permutations on the seven-token fixture.
        for perm in permutations(state):
            mapping=dict(zip(state,perm));inverse={v:k for k,v in mapping.items()};other={mapping[i]:s for i,s in state.items()};back={(e[0],frozenset(inverse[i] for i in e[1]),frozenset((inverse[i],s) for i,s in e[2].items()),e[3]) for e in matches(other)};assert back=={normalize(e) for e in events};relabels+=1
    state=states[0];events=matches(state);r4=next(e for e in events if e[0]=='R4');cat=next(i for i in r4[1] if state[i]=='X');exit_event=next(e for e in events if e[0]=='R5' and cat in e[1]);assert not r4[1].isdisjoint(exit_event[1]);assert cat not in r4[2] # unchanged catalyst still reserves token
    # Original match set excludes an X newly produced by G->X until a new snapshot.
    g=next(i for i,s in state.items() if s=='G');assert all(g not in e[1] for e in events if e[0]=='R4');after=apply(state,next(e for e in events if e[0]=='R2'));assert any(g in e[1] for e in matches(after) if e[0]=='R4')
    payload=dict(passed=True,scope='static local reaction schemas only; no substrate engine, RNG, transport or temporal qualification',exact_rational_checks=count,bijective_relabelings=relabels,catalyst_overlap_rejected=True,new_product_excluded_from_old_snapshot=True,shared_dependencies='schema code; independently written count formulas and integer increments; no independent dynamics engine',source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),contract_sha256=hashlib.sha256((ROOT/'spec/MG12H0_contract.md').read_bytes()).hexdigest());p=ROOT/'results/MG12H0/contract_checks.json';p.parent.mkdir(parents=True,exist_ok=True);p.write_text(json.dumps(payload,indent=2)+'\n',encoding='utf-8');print(json.dumps(payload))
if __name__=='__main__':main()

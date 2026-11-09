"""MG106 prospective bounded explicit-engine development and separate trace audit."""
import os
for key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[key]='1'
from pathlib import Path
import json,time,hashlib,argparse,traceback,subprocess,math
from dataclasses import replace
from collections import defaultdict
import numpy as np
from mg12h1 import State,Relation,SPECIES,matches,commit,select,simulate,relabel
from mg12h1_reference import transitions,master
from runner import seed_for
from mg10d1_pool import atomic_json
ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/'results/MG12H1/development';TIMES=np.array([0.,.05,.1]);R=256;STEPS=(None,.02,.01,.005)
FILES=('src/mg12h1.py','src/mg12h1_reference.py','src/prepare_mg12h1.py','spec/MG12H0_contract.md','spec/MG12H1_development.md','src/runner.py')
def hashes():return {name:hashlib.sha256((ROOT/name).read_text(encoding='utf-8').replace('\r\n','\n').encode()).hexdigest() for name in FILES}
def read(p):return json.loads(Path(p).read_text(encoding='utf-8'))
def fixture(c):
    sites=tuple(f's{i}' for i in range(2 if c==0 else 3 if c==2 else 1));links=[]
    if c==0:tokens=[(0,'G')];mask=(False,True,True,False,False,True);diff=(.7,.7,.7);pairs=[('ActiveLink',0,1,1.)]
    elif c==1:tokens=[(0,s) for s in ('A','G','X','X','Y','B')];mask=(True,)*6;diff=(0.,0.,0.);pairs=[]
    elif c==2:tokens=[(0,'A'),(1,'A'),(0,'B'),(2,'B')];mask=(False,)*6;diff=(0.,0.,0.);pairs=[('ALink',0,1,1.),('ALink',1,2,.5),('BLink',0,2,.75),('BLink',2,1,1.25)]
    elif c==3:tokens=[(0,s) for s in ('X','X','Y','Y')];mask=(False,False,False,False,True,False);diff=(0.,0.,0.);pairs=[]
    else:raise ValueError('Unknown fixture')
    for kind,a,b,w in pairs:links.extend([Relation(kind,(sites[a],sites[b]),w),Relation(kind,(sites[b],sites[a]),w)])
    state=State(sites,tuple(Relation('Occupies',(f't{i}',sites[s],species)) for i,(s,species) in enumerate(tokens))+tuple(links),2.,tuple((s,10.,4.) for s in sites),diff,(.4,.7),mask).validate();return state
def config(state):
    idx={s:i for i,s in enumerate(state.sites)};refs=[r[1:] for r in state.references];links=[(r.kind,idx[r.roles[0]],idx[r.roles[1]],r.weight) for r in state.relations if r.kind!='Occupies'];return state.omega,refs,links,state.diffusion,state.coupling,state.enabled
def from_counts(template,pop):
    rel=[r for r in template.relations if r.kind!='Occupies'];i=0
    for s,species in enumerate(SPECIES):
        for j,site in enumerate(template.sites):
            for _ in range(int(pop[s,j])):rel.append(Relation('Occupies',(f'n{i}',site,species)));i+=1
    return replace(template,relations=tuple(rel)).validate()
def jobs():return [dict(id=f'c{c}_m{m}',cell=c,method=m,seeds=[seed_for(12,c,m,10600,r) for r in range(R)]) for c in range(4) for m in range(4)]
def plan():
    if (BASE/'plan.json').exists():raise ValueError('Already planned')
    from prepare_mg11g4 import collect_seeds
    js=jobs();values=[s for j in js for s in j['seeds']];assert len(set(values))==len(values);prior=set();scan=subprocess.run(['rg','-i','-l','-g','*.json','"[^"\\n]*seed[^"\\n]*"','results'],cwd=ROOT,capture_output=True,text=True)
    if scan.returncode not in (0,1):raise ValueError(scan.stderr)
    for name in scan.stdout.splitlines():
        if not name.replace('\\','/').startswith('results/MG12H1/'):collect_seeds(read(ROOT/name),prior)
    assert not set(values)&prior
    atomic_json(BASE/'plan.json',dict(stage='MG106',authorization='go ahead with the engine development',scope='bounded development only; not formal production',jobs=js,R=R,times=TIMES.tolist(),steps=STEPS,source_hashes=hashes(),seed_collisions=0,prior_seed_values=len(prior),ceiling_elapsed_seconds=300,workers=1,threads_per_library=1,exact_DKW_alpha_family=.01,exact_DKW_comparisons=4))
def trace_audit(initial,result,delta):
    # No primary commit or match calls: literal role/state updates and identity ledger.
    state=dict(initial.occupations());original=set(state);last_time=0.;step=-1;reserved=set();base={};owners={};parents={};accepted=0
    links={(r.kind,*r.roles) for r in initial.relations if r.kind!='Occupies'}
    for i,row in enumerate(result['trace']):
        assert row['time']>=last_time and row['time']<=TIMES[-1];last_time=row['time']
        if row['step']!=step:assert row['step']>step;step=row['step'];reserved=set();base=dict(state);owners={}
        tokens=set(row['participants']);before={t:(s,sp) for t,s,sp in row['before']};after={t:(s,sp) for t,s,sp in row['after']};assert set(before)==tokens and all(base.get(t)==v for t,v in before.items())
        assert tuple(sorted({parents[t] for t in tokens if t in parents}))==tuple(row['parents'])
        expected_blockers=tuple(sorted({owners[t] for t in tokens if t in owners}));assert tuple(row['blocked_by'])==expected_blockers
        if not row['accepted']:assert delta is not None and tokens & reserved;continue
        assert not tokens & reserved;rule=row['rule'];expected={}
        if rule.startswith('H_'):
            assert len(tokens)==1;token=next(iter(tokens));site,species=before[token];newsite,newsp=after[token];kind='ActiveLink' if species in SPECIES[:3] else 'ALink' if species=='A' else 'BLink';assert species==newsp and (kind,site,newsite) in links;expected=after
        elif rule=='R3':
            assert sorted(sp for _,sp in before.values())==['B','X'] and len({s for s,_ in before.values()})==1
            expected={t:(s,'Z' if sp=='B' else 'Y') for t,(s,sp) in before.items()}
        elif rule=='R4':
            assert sorted(sp for _,sp in before.values())==['X','X','Y'] and len({s for s,_ in before.values()})==1;expected={t:(s,'X') for t,(s,sp) in before.items() if sp=='Y'}
        else:
            source,dest={'R1':('A','G'),'R2':('G','X'),'R2r':('X','G'),'R5':('X','W')}[rule];assert len(tokens)==1;token=next(iter(tokens));site,species=before[token];assert species==source;expected={token:(site,dest)}
        assert after==expected;state.update(after);reserved.update(tokens)
        for token in tokens:owners[token]=i;parents[token]=i
        assert set(state)==original;accepted+=1
    assert state==result['frames'][-1].occupations();assert np.all(result['counts'].sum(axis=(1,2))==len(original));return accepted
def algebra():
    c=fixture(1);es=matches(c);r4=next(e for e in es if e.rule=='R4');x=next(t for t,_,s in r4.before if s=='X');loss=next(e for e in es if e.rule=='R5' and x in e.participants);out,trace=select(c,[(r4,.001),(loss,.002)],.01);assert [r['accepted'] for r in trace]==[True,False]
    disjoint=next((a,b) for a in es for b in es if not a.participants & b.participants);out,trace=select(c,[(disjoint[0],.001),(disjoint[1],.002)],.01);assert all(r['accepted'] for r in trace)
    feedback=fixture(3);events=matches(feedback);out,trace=select(feedback,[(e,.001+i*.001) for i,e in enumerate(events)],.02);assert sum(r['accepted'] for r in trace)==1
    # Transport into a new site cannot immediately create a second hop in old matches.
    moving=fixture(0);hop=next(e for e in matches(moving) if e.rule=='H_G');assert all(e.before[0][1]=='s0' for e in matches(moving));assert commit(moving,[hop]).occupations()['t0'][0]=='s1'
    bad=0
    calls=[lambda:select(c,[(r4,.001),(loss,.001)],.01),lambda:select(c,[(r4,.001),(r4,.002)],.01),lambda:commit(c,[r4,loss]),lambda:commit(c,[replace(r4,rate=999)]),lambda:replace(c,omega=-1).validate(),lambda:replace(c,relations=c.relations+(c.relations[0],)).validate(),lambda:simulate(c,[0,0],1),lambda:simulate(c,TIMES,-1),lambda:simulate(c,TIMES,1,delta=0)]
    for call in calls:
        try:call()
        except ValueError:bad+=1
        else:raise AssertionError('Invalid input accepted')
    mapping={key:f'new_{i*7+3}' for i,key in enumerate(c.sites+tuple(c.occupations()))};inverse={v:k for k,v in mapping.items()};other=relabel(c,mapping)
    back={(e.rule,frozenset(inverse[t] for t in e.participants),frozenset((inverse[t],inverse[s],sp) for t,s,sp in e.after),e.rate) for e in matches(other)};original={(e.rule,e.participants,frozenset(e.after),e.rate) for e in es};assert back==original
    # Explicitly corresponding marks, rather than relying on matching enumeration/RNG order.
    a=next(e for e in matches(other) if e.rule=='R4');out1,_=select(c,[(r4,.001)],.01);out2,_=select(other,[(a,.001)],.01);assert {inverse[t]:(inverse[s],sp) for t,(s,sp) in out2.occupations().items()}==out1.occupations()
    return dict(invalid_rejected=bad,catalyst_conflict=True,shared_site_disjoint=True,feedback_capacity=True,arrival_snapshot=True,relabeling_and_corresponding_marks=True)
def run():
    p=read(BASE/'plan.json')
    if hashes()!=p['source_hashes'] or (BASE/'attempt.json').exists():raise ValueError('Changed source or already attempted')
    subprocess.run(['git','ls-files','--error-unmatch','results/MG12H1/development/plan.json'],cwd=ROOT,check=True,stdout=subprocess.DEVNULL);start=time.perf_counter();deadline=start+300;atomic_json(BASE/'attempt.json',dict(status='RUNNING'));checks=algebra();refs={};generator_states=0
    for c in range(4):
        state=fixture(c);args=config(state);states,prob,q=master(state.counts(),*args,.1);refs[c]=(states,prob);n=len(state.sites)
        for key in states:
            if time.perf_counter()>deadline:raise TimeoutError('Development deadline')
            pop=np.asarray(key,np.int64).reshape(7,n);actual=defaultdict(float)
            explicit=from_counts(state,pop)
            for event in matches(explicit):actual[event.rule,tuple(commit(explicit,[event]).counts().ravel())]+=event.rate
            expected=defaultdict(float)
            for rule,dest,rate in transitions(pop,*args):expected[rule,dest]+=rate
            assert set(actual)==set(expected)
            for k,v in actual.items():np.testing.assert_allclose(v,expected[k],rtol=1e-12,atol=1e-12)
            generator_states+=1
        np.savez_compressed(BASE/f'c{c}_reference.npz',states=np.array(states),probabilities=prob);atomic_json(BASE/f'c{c}_reference.json',dict(states=len(states),row_sum_error=float(abs(q.sum(axis=1)).max())))
    results=[];epsilon=math.sqrt(math.log(2*4/.01)/(2*R))
    for j in p['jobs']:
        state=fixture(j['cell']);tau=STEPS[j['method']];states,prob=refs[j['cell']];lookup={key:i for i,key in enumerate(states)};samples=[];traces=[];accepted=0;clock=time.perf_counter()
        for seed in j['seeds']:
            result=simulate(state,TIMES,seed,tau,deadline=deadline);accepted+=trace_audit(state,result,tau);samples.append(result['counts']);traces.append(result['trace'])
        data=np.asarray(samples);idx=np.array([lookup[tuple(z[-1].ravel())] for z in data]);distance=float(np.max(abs(np.cumsum(np.bincount(idx,minlength=len(states)))/R-np.cumsum(prob))))
        np.savez_compressed(BASE/(j['id']+'.npz'),counts=data,state_index=idx,seeds=np.asarray(j['seeds'],np.uint32),times=TIMES);atomic_json(BASE/(j['id']+'_traces.json'),dict(job=j,traces=traces));meta=dict(job=j,elapsed_seconds=time.perf_counter()-clock,joint_state_index_CDF_distance=distance,exact_DKW_epsilon=epsilon,exact_DKW_supported=bool(distance<=epsilon) if tau is None else None,accepted_events_audited=accepted);atomic_json(BASE/(j['id']+'.json'),meta);results.append(meta)
        # Held replay uses an already registered seed and must reproduce identities/trace.
        replay=simulate(state,TIMES,j['seeds'][0],tau,deadline=deadline);np.testing.assert_array_equal(replay['counts'],data[0]);assert replay['trace']==traces[0]
        if tau is None:assert distance<=epsilon,'Exact-law smoke check failed'
    elapsed=time.perf_counter()-start;assert elapsed<300 and hashes()==p['source_hashes'];atomic_json(BASE/'results.json',dict(passed=True,scope='development only; strict tick distributions descriptive; formal qualification pending',source_hashes=hashes(),checks=checks,generator_states=generator_states,jobs=results,elapsed_seconds=elapsed,R=R,development_paths=4096,held_replays=16,exact_DKW_alpha_family=.01,exact_DKW_epsilon=epsilon));atomic_json(BASE/'attempt.json',dict(status='COMPLETE',elapsed_seconds=elapsed))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--plan',action='store_true');p.add_argument('--run',action='store_true');a=p.parse_args()
    if a.plan:plan()
    elif a.run:
        try:run()
        except BaseException:atomic_json(BASE/'failure.json',dict(error=traceback.format_exc(),source_hashes=hashes()));raise
    else:p.error('Choose action')

"""MG114 mixed incidence kinetics and conditional ring-drift development."""
import os
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[key]='1'
from pathlib import Path
from dataclasses import replace,asdict
from collections import defaultdict
import json,time,math,hashlib,subprocess,argparse,traceback
import numpy as np
from mg12h2 import HyperState,Relation,matches,commit,simulate
from mg12h1 import SPECIES
from mg12h1_reference import transitions,master
from runner import seed_for
from mg10d1_pool import atomic_json
ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/'results/MG12H3/bridge';TIMES=np.array([0.,.05,.1]);R=128;TAUS=(None,.005,.001)
FILES=('src/prepare_mg12h3_bridge.py','src/check_mg12h3_bridge.py','spec/MG12H3_bridge_development.md','src/mg12h1.py','src/mg12h2.py','src/mg12h1_reference.py','src/mg11g2.py','src/mg11n1.py','src/mg11g2_pattern_v2.py','src/runner.py')
def read(p):return json.loads(Path(p).read_text(encoding='utf-8'))
def hashes():return {name:hashlib.sha256((ROOT/name).read_text(encoding='utf-8').replace('\r\n','\n').encode()).hexdigest() for name in FILES}
def fixture(c):
    n=2 if c==3 else 3;species=[('A','X'),('B','X'),('X','X','Y'),('A','B','X','Y')][c];relations=[Relation('Occupies',(f't{i}',0,sp)) for i,sp in enumerate(species)];roles=tuple(range(n));relations.append(Relation('ActiveCell',roles,1.))
    if c==2:relations.append(Relation('ActiveCell',(0,1),.5))
    if c in (0,3):relations.append(Relation('ACell',roles,1.))
    if c in (1,3):relations.append(Relation('BCell',roles,1.))
    return HyperState(tuple(range(n)),tuple(relations),2.,tuple((i,10.,4.) for i in range(n)),(.2,.3,.4),(.4,.7),(True,)*6).validate()
def compile_links(state):
    # Incidence-block algebra, independent of engine event enumeration.
    kinds={'ActiveCell':'ActiveLink','ACell':'ALink','BCell':'BLink'};out=[];n=len(state.sites)
    for kind,target in kinds.items():
        matrix=np.zeros((n,n))
        for cell in state.relations:
            if cell.kind!=kind:continue
            ids=list(cell.roles);matrix[np.ix_(ids,ids)]+=cell.weight*(np.ones((len(ids),len(ids)))-np.eye(len(ids)))/(len(ids)-1)
        out.extend((target,int(i),int(j),float(matrix[i,j])) for i,j in zip(*np.where(matrix>0)))
    return out
def config(state):return state.omega,[r[1:] for r in state.references],compile_links(state),state.diffusion,state.coupling,state.enabled
def instantiate(state,pop):
    relations=[r for r in state.relations if r.kind!='Occupies'];index=0
    for sp,species in enumerate(SPECIES):
        for site in state.sites:
            for _ in range(int(pop[sp,site])):relations.append(Relation('Occupies',(f'p{index}',site,species)));index+=1
    return replace(state,relations=tuple(relations)).validate()
def jobs():return [dict(id=f'c{c}_m{m}',cell=c,method=m,seeds=[seed_for(12,c,m,11400,r) for r in range(R)]) for c in range(4) for m in range(3)]
def plan():
    if (BASE/'plan.json').exists():raise ValueError('Already planned')
    from prepare_mg11g4 import collect_seeds
    prior=set();scan=subprocess.run(['rg','-i','-l','-g','*.json','"[^"\\n]*seed[^"\\n]*"','results'],cwd=ROOT,capture_output=True,text=True)
    if scan.returncode not in (0,1):raise ValueError(scan.stderr)
    for name in scan.stdout.splitlines():
        if not name.replace('\\','/').startswith('results/MG12H3/'):collect_seeds(read(ROOT/name),prior)
    used=set(prior);collisions=[];js=jobs()
    for j in js:
        for r,seed in enumerate(j['seeds']):
            nonce=0
            while seed in used:collisions.append(dict(job=j['id'],replicate=r,nonce=nonce,rejected=seed));nonce+=1;seed=seed_for(12,j['cell'],j['method'],11400,r,nonce)
            j['seeds'][r]=seed;used.add(seed)
    atomic_json(BASE/'plan.json',dict(stage='MG114',authorization='I sign off. Then proceed to doing H3.',scope='mixed/ring-limit bounded development, no particle run',source_hashes=hashes(),fixtures=[asdict(fixture(c)) for c in range(4)],jobs=js,R=R,taus=TAUS,times=TIMES.tolist(),ceiling_elapsed_seconds=300,workers=1,threads_per_library=1,collisions_resolved=collisions,prior_seed_values=len(prior)))
def ring_drift():
    from mg11g2_pattern_v2 import environment
    from mg11g2 import initial,rhs
    n=8;g=environment('LOCAL',n);concentration=np.zeros((7,n));concentration[0]=1;concentration[1]=1+np.arange(n)%2;concentration[2]=1;concentration[3]=2;concentration[4]=3;diff=np.array([.2,.3,.4]);k=np.array([.4,.7]);links=[]
    for i in range(n):
        for j in ((i-1)%n,(i+1)%n):links.append(('ActiveLink',i,j,1.))
    for kind in ('ALink','BLink'):links.extend((kind,int(i),int(j),float(w)) for i,j,w in zip(g['src'],g['dst'],g['weight']))
    rows=[]
    for omega in (1,10,100,1000):
        pop=(concentration*omega).astype(np.int64);refs=[(10*omega,4*omega)]*n;drift=np.zeros_like(concentration)
        for rule,dest,rate in transitions(pop,omega,refs,links,diff,k,(True,)*6):drift+=(np.asarray(dest).reshape(7,n)-pop)*rate/omega
        case=dict(initial=pop,omega=float(omega),a0=np.full(n,10*omega),b0=np.full(n,4*omega),diff=diff);continuum=rhs(0,initial(case,g),case,g,k,'FINITE')[:7*n].reshape(7,n);correction=np.zeros_like(drift);correction[1]=-concentration[1]*concentration[2]/omega;correction[2]=-correction[1];np.testing.assert_allclose(drift,continuum+correction,atol=1e-11,rtol=1e-11);assert abs(drift.sum())<1e-10
        rows.append(dict(omega=omega,max_continuum_difference=float(abs(drift-continuum).max()),max_corrected_residual=float(abs(drift-continuum-correction).max()),global_drift=float(drift.sum())))
    return rows
def run():
    p=read(BASE/'plan.json')
    if p['source_hashes']!=hashes() or (BASE/'attempt.json').exists():raise ValueError('Changed source or existing attempt')
    subprocess.run(['git','ls-files','--error-unmatch','results/MG12H3/bridge/plan.json'],cwd=ROOT,check=True,stdout=subprocess.DEVNULL);start=time.perf_counter();deadline=start+300;atomic_json(BASE/'attempt.json',dict(status='RUNNING'));ring=ring_drift();references={};generator=[]
    for c in range(4):
        state=fixture(c);args=config(state);states,prob,q=master(state.counts(),*args,.1,max_states=8192);n=len(state.sites)
        for key in states:
            if time.perf_counter()>deadline:raise TimeoutError('Combined generator deadline')
            pop=np.asarray(key,np.int64).reshape(7,n);explicit=instantiate(state,pop);actual=defaultdict(float)
            for event in matches(explicit):actual[event.rule,tuple(commit(explicit,[event]).counts().ravel())]+=event.rate
            expected=defaultdict(float)
            for rule,dest,rate in transitions(pop,*args):expected[rule,dest]+=rate
            assert set(actual)==set(expected)
            for key,rate in expected.items():np.testing.assert_allclose(actual[key],rate,rtol=1e-12,atol=1e-12)
        references[c]=(states,prob);generator.append(dict(cell=c,states=len(states),max_row_sum=float(abs(q.sum(axis=1)).max())));np.savez_compressed(BASE/f'c{c}_reference.npz',states=np.array(states),probabilities=prob)
    from check_mg12h3_bridge import audit_trace
    rows=[];accepted=0;epsilon=math.sqrt(math.log(800)/(2*R))
    for job in p['jobs']:
        state=fixture(job['cell']);tau=TAUS[job['method']];states,prob=references[job['cell']];lookup={key:i for i,key in enumerate(states)};data=[];traces=[];clock=time.perf_counter()
        for seed in job['seeds']:
            result=simulate(state,TIMES,seed,tau,deadline);accepted+=audit_trace(state,result['counts'],result['trace'],job['method']);data.append(result['counts']);traces.append(result['trace'])
        data=np.asarray(data);indices=np.array([lookup[tuple(x[-1].ravel())] for x in data]);distance=float(max(abs(np.cumsum(np.bincount(indices,minlength=len(states)))/R-np.cumsum(prob))));np.savez_compressed(BASE/(job['id']+'.npz'),counts=data,seeds=np.asarray(job['seeds'],np.uint32),times=TIMES);atomic_json(BASE/(job['id']+'_traces.json'),dict(job=job,traces=traces));rows.append(dict(id=job['id'],CDF_distance=distance,elapsed_seconds=time.perf_counter()-clock))
        if job['method']==0:assert distance<=epsilon
    elapsed=time.perf_counter()-start;assert elapsed<300 and hashes()==p['source_hashes'];atomic_json(BASE/'results.json',dict(passed=True,source_hashes=hashes(),generator=generator,ring_drift=ring,jobs=rows,paths=1536,accepted_events_audited=accepted,exact_DKW_epsilon=epsilon,elapsed_seconds=elapsed,scope='development: combined generator/ring conditional drift only, no particle trajectory or formal strict-tick accuracy'));atomic_json(BASE/'attempt.json',dict(status='COMPLETE'))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--plan',action='store_true');p.add_argument('--run',action='store_true');a=p.parse_args()
    if a.plan:plan()
    elif a.run:
        try:run()
        except BaseException:atomic_json(BASE/'failure.json',dict(error=traceback.format_exc()));raise
    else:p.error('Choose action')

"""MG110 fixed-cell development, independent matrix/reference/dimension checks."""
import os
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[key]='1'
from pathlib import Path
from dataclasses import replace,asdict
from itertools import product
import json,math,time,hashlib,subprocess,argparse,traceback
import numpy as np
from scipy.linalg import expm
from mg12h2 import HyperState,Relation,matches,commit,select,simulate
from mg12h1 import SPECIES,matches as old_matches,State,relabel
from runner import seed_for
from mg10d1_pool import atomic_json
ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/'results/MG12H2/development';TIMES=np.array([0.,.05,.1]);R=256
FILES=('src/mg12h2.py','src/prepare_mg12h2.py','spec/MG12H2_development.md','src/mg12h1.py','src/runner.py')
def read(p):return json.loads(Path(p).read_text(encoding='utf-8'))
def hashes():return {n:hashlib.sha256((ROOT/n).read_text(encoding='utf-8').replace('\r\n','\n').encode()).hexdigest() for n in FILES}
def fixture(c):
    n=4;cells=[];coords=None;d=0
    if c==0:cells=[('ActiveCell',(i,(i+1)%4),1.) for i in range(4)]
    elif c==1:n=3;cells=[('ActiveCell',(0,1,2),1.)]
    elif c==2:cells=[('ActiveCell',(0,1,2),1.),('ActiveCell',(0,1,2,3),.6)]
    elif c in (3,4):
        d=c-1;coords=list(product(range(3),repeat=d));n=len(coords);idx={x:i for i,x in enumerate(coords)};seen=set()
        for i,x in enumerate(coords):
            for axis in range(d):
                y=list(x);y[axis]=(y[axis]+1)%3;j=idx[tuple(y)];key=frozenset((i,j))
                if key not in seen:cells.append(('ActiveCell',(i,j),1.));seen.add(key)
    elif c==5:cells=[('ACell',(0,1,2),1.),('ACell',(1,2,3),.5),('BCell',(0,2,3),.75),('BCell',(0,1,2,3),1.25)]
    else:raise ValueError('Unknown fixture')
    tokens=[Relation('Occupies',('a',0,'A')),Relation('Occupies',('b',3,'B'))] if c==5 else [Relation('Occupies',('g',0,'G'))];state=HyperState(tuple(range(n)),tuple(tokens)+tuple(Relation(kind,roles,w) for kind,roles,w in cells),2.,tuple((i,10.,4.) for i in range(n)),(1.,1.,1.),(.4,.7),(False,)*6).validate();return state,coords,d
def matrices(state):
    n=len(state.sites);out={};index={s:i for i,s in enumerate(state.sites)}
    for sp in ('G','A','B'):
        q=np.zeros((n,n));kind='ActiveCell' if sp=='G' else 'ACell' if sp=='A' else 'BCell';rate=state.diffusion[0] if sp=='G' else state.coupling[0 if sp=='A' else 1]
        for rel in state.relations:
            if rel.kind!=kind:continue
            ids=np.array([index[s] for s in rel.roles]);block=np.ones((len(ids),len(ids)))-np.eye(len(ids));q[np.ix_(ids,ids)]+=rate*rel.weight*block/(len(ids)-1)
        q-=np.diag(q.sum(axis=1));out[sp]=q
    return out
def reference(state):
    mats=matrices(state);tokens=tuple(state.occupations());species=[state.occupations()[t][1] for t in tokens];locations=list(product(state.sites,repeat=len(tokens)));states=[]
    for loc in locations:
        p=np.zeros((7,len(state.sites)),np.int64)
        for sp,site in zip(species,loc):p[SPECIES.index(sp),site]+=1
        states.append(p)
    if len(tokens)==1:q=mats[species[0]]
    else:q=np.kron(mats[species[0]],np.eye(len(state.sites)))+np.kron(np.eye(len(state.sites)),mats[species[1]])
    initial=locations.index(tuple(state.occupations()[t][0] for t in tokens));prob=np.array([expm(q*t)[initial] for t in TIMES[1:]]);assert prob.min()>=-1e-12 and max(abs(prob.sum(axis=1)-1))<1e-11;return np.array(states),prob,q
def jobs():return [dict(id=f'c{c}_m{m}',cell=c,method=m,seeds=[seed_for(12,c,m,11000,r) for r in range(R)]) for c in range(6) for m in range(2)]
def plan():
    if (BASE/'plan.json').exists():raise ValueError('Already planned')
    from prepare_mg11g4 import collect_seeds
    prior=set();scan=subprocess.run(['rg','-i','-l','-g','*.json','"[^"\\n]*seed[^"\\n]*"','results'],cwd=ROOT,capture_output=True,text=True)
    if scan.returncode not in (0,1):raise ValueError(scan.stderr)
    for name in scan.stdout.splitlines():
        if not name.replace('\\','/').startswith('results/MG12H2/'):collect_seeds(read(ROOT/name),prior)
    js=jobs();used=set(prior);collisions=[]
    for j in js:
        for r,seed in enumerate(j['seeds']):
            nonce=0
            while seed in used:collisions.append(dict(job=j['id'],replicate=r,nonce=nonce,rejected=seed));nonce+=1;seed=seed_for(12,j['cell'],j['method'],11000,r,nonce)
            j['seeds'][r]=seed;used.add(seed)
    atomic_json(BASE/'plan.json',dict(stage='MG110',scope='H2bounded development only',source_hashes=hashes(),jobs=js,fixtures=[asdict(fixture(c)[0]) for c in range(6)],R=R,prior_seed_values=len(prior),collision_resolutions=collisions,ceiling_elapsed_seconds=300,workers=1,threads_per_library=1))
def audit(initial,result,method):
    occ=dict(initial.occupations());frame=0;last=-1.;step=-1;reserved=set();base={};accepted=0
    def check(index):
        p=np.zeros_like(result['counts'][index])
        for site,sp in occ.values():p[SPECIES.index(sp),site]+=1
        np.testing.assert_array_equal(p,result['counts'][index])
    for row in result['trace']:
        while frame<3 and TIMES[frame]<row['time']:check(frame);frame+=1
        assert last<=row['time']<=.1;last=row['time']
        if row['step']!=step:step=row['step'];reserved=set();base=dict(occ)
        tokens=set(row['participants']);before={t:(s,sp) for t,s,sp in row['before']};after={t:(s,sp) for t,s,sp in row['after']};assert len(tokens)==1 and set(before)==set(after)==tokens and all(base[t]==v for t,v in before.items())
        t=next(iter(tokens));source,sp=before[t];dest,newsp=after[t];assert sp==newsp and source!=dest;kind,roles=row['cell'];assert any(r.kind==kind and set(r.roles)==set(roles) for r in initial.relations);assert source in roles and dest in roles
        if row['accepted']:assert not tokens & reserved;occ.update(after);reserved.update(tokens);accepted+=1
        else:assert method==1 and tokens & reserved
    while frame<3:check(frame);frame+=1
    assert np.all(result['counts'].sum(axis=(1,2))==len(initial.occupations()));return accepted
def run():
    p=read(BASE/'plan.json')
    if p['source_hashes']!=hashes() or (BASE/'attempt.json').exists():raise ValueError('Changed source or existing attempt')
    subprocess.run(['git','ls-files','--error-unmatch','results/MG12H2/development/plan.json'],cwd=ROOT,check=True,stdout=subprocess.DEVNULL);start=time.perf_counter();deadline=start+300;atomic_json(BASE/'attempt.json',dict(status='RUNNING'));checks=[];refs={};states_checked=0
    for c in range(6):
        state,coords,d=fixture(c);states,prob,q=reference(state);refs[c]=(states,prob);mats=matrices(state)
        for mat in mats.values():assert abs(mat.sum(axis=1)).max()<1e-12 and abs(mat.sum(axis=0)).max()<1e-12 and abs(mat-mat.T).max()<1e-12
        for idx,pop in enumerate(states):
            if time.perf_counter()>deadline:raise TimeoutError('H2development cap')
            loc=[]
            for token,(old,sp) in state.occupations().items():loc.append((token,int(np.flatnonzero(pop[SPECIES.index(sp)])[0]),sp))
            explicit=replace(state,relations=tuple(Relation('Occupies',x) for x in loc)+tuple(r for r in state.relations if r.kind!='Occupies'));rates=np.zeros(len(states));lookup={tuple(v.ravel()):i for i,v in enumerate(states)}
            for event in matches(explicit):rates[lookup[tuple(commit(explicit,[event]).counts().ravel())]]+=event.rate
            rates[idx]=-rates.sum();np.testing.assert_allclose(rates,q[idx],rtol=1e-12,atol=1e-12);states_checked+=1
        if coords is not None:
            mode=np.cos(2*np.pi*np.array(coords)[:,0]/3);np.testing.assert_allclose(expm(mats['G']*.1)@mode,np.exp(-3*.1)*mode,atol=1e-11,rtol=0)
        # Independent weighted pairwise compilation; generator equivalence, not byte replay.
        links=[]
        for sp,kind in [('G','ActiveLink'),('A','ALink'),('B','BLink')]:
            scale=1 if sp=='G' else state.coupling[0 if sp=='A' else 1]
            for i,j in zip(*np.where(mats[sp]>0)):links.append(Relation(kind,(int(i),int(j)),float(mats[sp][i,j]/scale)))
        old=State(state.sites,tuple(r for r in state.relations if r.kind=='Occupies')+tuple(links),state.omega,state.references,state.diffusion,state.coupling,state.enabled).validate();a=sum(e.rate for e in old_matches(old));b=sum(e.rate for e in matches(state));assert abs(a-b)<1e-12
        checks.append(dict(cell=c,states=len(states),pairwise_initial_exit=a,grid_dimension=d if coords else None));np.savez_compressed(BASE/f'c{c}_reference.npz',states=states,probabilities=prob,generator=q)
    dimensions=[]
    for d in (1,2,3):
        balls=[]
        for radius in (4,8):
            direct=sum(sum(abs(x) for x in coords)<=radius for coords in product(range(-radius,radius+1),repeat=d));formula=sum(2**k*math.comb(d,k)*math.comb(radius,k) for k in range(d+1));assert direct==formula;balls.append(direct)
        dimensions.append(dict(known_dimension=d,ball_r4=balls[0],ball_r8=balls[1],finite_scale_slope=math.log(balls[1]/balls[0])/math.log(2)))
    # Overlapping cells must keep distinct channels sharing the same move.
    mixed=fixture(2)[0];es=matches(mixed);assert len(es)==5 and abs(sum(e.rate for e in es)-1.6)<1e-12;offers=[(e,.001+i*.001) for i,e in enumerate(es)];out,tr=select(mixed,offers,.01);assert sum(x['accepted'] for x in tr)==1
    mapping={key:f'label_{i}' for i,key in enumerate(mixed.sites+tuple(mixed.occupations()))};inverse={v:k for k,v in mapping.items()};changed=relabel(mixed,mapping)
    original={(e.rule,e.participants,frozenset(e.after),e.rate,e.cell[0],frozenset(e.cell[1])) for e in matches(mixed)}
    restored={(e.rule,frozenset(inverse[t] for t in e.participants),frozenset((inverse[t],inverse[site],sp) for t,site,sp in e.after),e.rate,e.cell[0],frozenset(inverse[site] for site in e.cell[1])) for e in matches(changed)};assert original==restored
    first=es[0];corresponding=next(e for e in matches(changed) if frozenset((inverse[t],inverse[site],sp) for t,site,sp in e.after)==frozenset(first.after) and frozenset(inverse[site] for site in e.cell[1])==frozenset(first.cell[1]));left,_=select(mixed,[(first,.001)],.01);right,_=select(changed,[(corresponding,.001)],.01);assert left.occupations()=={inverse[t]:(inverse[site],sp) for t,(site,sp) in right.occupations().items()}
    shared=replace(mixed,relations=mixed.relations+(Relation('Occupies',('other',0,'G')),));events=matches(shared);one=next(e for e in events if 'g' in e.participants);two=next(e for e in events if 'other' in e.participants and e.cell==one.cell);_,rows=select(shared,[(one,.001),(two,.002)],.01);assert all(row['accepted'] for row in rows)
    invalid=0
    for relation in (Relation('ActiveCell',(0,0,1)),Relation('ActiveCell',(0,)),Relation('ActiveCell',(0,9)),Relation('ActiveCell',(0,1),-1)):
        try:replace(mixed,relations=mixed.relations+(relation,)).validate()
        except ValueError:invalid+=1
        else:raise AssertionError('Invalid cell accepted')
    results=[];epsilon=math.sqrt(math.log(12/.01)/(2*R));accepted=0
    for job in p['jobs']:
        state=fixture(job['cell'])[0];states,prob=refs[job['cell']];lookup={tuple(v.ravel()):i for i,v in enumerate(states)};data=[];traces=[];clock=time.perf_counter();tau=None if job['method']==0 else .005
        for seed in job['seeds']:
            result=simulate(state,TIMES,seed,tau,deadline);accepted+=audit(state,result,job['method']);data.append(result['counts']);traces.append(result['trace'])
        data=np.array(data);idx=np.array([lookup[tuple(v[-1].ravel())] for v in data]);distance=float(max(abs(np.cumsum(np.bincount(idx,minlength=len(states)))/R-np.cumsum(prob[1]))));np.savez_compressed(BASE/(job['id']+'.npz'),counts=data,seeds=np.asarray(job['seeds'],np.uint32),times=TIMES);atomic_json(BASE/(job['id']+'_traces.json'),dict(job=job,traces=traces));results.append(dict(id=job['id'],elapsed_seconds=time.perf_counter()-clock,CDF_distance=distance))
        if job['method']==0:assert distance<=epsilon
    elapsed=time.perf_counter()-start;assert elapsed<300 and hashes()==p['source_hashes'];atomic_json(BASE/'results.json',dict(passed=True,source_hashes=hashes(),generator_states=states_checked,checks=checks,dimensions=dimensions,overlap_channels_preserved=True,invalid_rejected=invalid,paths=3072,accepted_events_audited=accepted,jobs=results,exact_DKW_epsilon=epsilon,elapsed_seconds=elapsed,scope='development only; pairwise-equivalent transport, known-dimension calibration, no particle or emergent-space result'));atomic_json(BASE/'attempt.json',dict(status='COMPLETE'))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--plan',action='store_true');p.add_argument('--run',action='store_true');a=p.parse_args()
    if a.plan:plan()
    elif a.run:
        try:run()
        except BaseException:atomic_json(BASE/'failure.json',dict(error=traceback.format_exc()));raise
    else:p.error('Choose action')

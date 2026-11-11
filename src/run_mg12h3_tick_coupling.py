"""MG156 bounded development; freeze plan/source before every stochastic attempt."""
import os
for key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[key]='1'
from pathlib import Path
from collections import defaultdict
from fractions import Fraction
import json,math,time,hashlib,subprocess,traceback,sys,platform,argparse
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import expm_multiply
from mg12h3_sparse import Model,tick
from mg12h3_tick_coupling import pair,replay
from prepare_mg12h3_bridge import fixture,config,instantiate,ROOT,seed_for
from mg12h2 import matches
from mg12h1 import SPECIES
from mg10e2c1_pool import atomic_json
BASE=ROOT/'results/MG12H3/tick_coupling_development'
FILES=['src/mg12h3_tick_coupling.py','src/run_mg12h3_tick_coupling.py','spec/MG12H3_tick_coupling_development.md','src/mg12h3_sparse.py','src/prepare_mg12h3_bridge.py','src/mg12h1.py','src/mg12h2.py','src/runner.py','src/mg12h3_finer_precision.py','src/mg12h3_contrast_estimator.py','src/run_mg12h3_sparse_combined.py']
INPUTS=['results/MG12H3/ring_feasibility/o0_a0_Omega10.npz']+[f'results/MG12H3/bridge/c{c}_reference.npz' for c in range(4)]
def read(p):return json.loads(Path(p).read_text())
def hashes():
    out={p:hashlib.sha256((ROOT/p).read_text().replace('\r\n','\n').encode()).hexdigest() for p in FILES}
    out.update({p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in INPUTS});return out
def plan():
    assert not (BASE/'plan.json').exists();prior=set()
    def collect(x):
        if isinstance(x,dict):
            for v in x.values():collect(v)
        elif isinstance(x,list):
            for v in x:collect(v)
        elif type(x) is int and 0<=x<2**32:prior.add(x)
    for p in (ROOT/'results').rglob('*.json'):
        if 'plan' in p.name or 'registration' in p.name:collect(read(p))
    used=set(prior);collisions=[]
    def seeds(stage,c,n):
        out=[]
        for r in range(n):
            nonce=0;value=seed_for(12,c,0,stage,r,nonce)
            while value in used:
                collisions.append([stage,c,r,nonce,value]);nonce+=1;value=seed_for(12,c,0,stage,r,nonce)
            used.add(value);out.append(value)
        return out
    s=dict(regression=[seeds(15600,c,32) for c in range(4)],marginal=[seeds(15601,c,2048) for c in range(4)],transport=seeds(15602,0,4096),actual_paired=[seeds(15603,c,4) for c in range(2)],actual_single=[seeds(15604,c,8) for c in range(2)])
    atomic_json(BASE/'plan.json',dict(stage='MG156',authorization='Looks good to me. Go ahead and prepare and then run the tests.',scope='bounded development, not production qualification',source_hashes=hashes(),seeds=s,prior_seed_values=len(prior),collisions_resolved=collisions,stage_A_seconds=60,stage_B_seconds=60,workers=1,threads=1,environment=dict(python=sys.version,numpy=np.__version__,platform=platform.platform()),storage='all final samples/histograms/exact probabilities/seeds; selected full joint traces and boundary records'))
def exact(initial,delta,deadline):
    # Identity-level frozen matches, not count channels or the coupling sampler.
    occ=initial.occupations();tokens=list(occ);lookup={t:i for i,t in enumerate(tokens)};n=len(initial.sites)
    start=((1<<len(tokens))-1,tuple(occ[t][0]*7+SPECIES.index(occ[t][1]) for t in tokens))
    marks=[]
    for e in matches(initial):
        mask=sum(1<<lookup[t] for t in e.participants);after=tuple((lookup[t],site*7+SPECIES.index(sp)) for t,site,sp in e.after);marks.append((mask,after,e.rate))
    states=[start];index={start:0};rows=[];cols=[];rates=[];cursor=0
    while cursor<len(states):
        if time.perf_counter()>deadline:raise TimeoutError('Exact generator deadline')
        available,output=states[cursor];total=0.
        for mask,after,rate in marks:
            if available&mask!=mask:continue
            dest=list(output)
            for t,value in after:dest[t]=value
            key=(available^mask,tuple(dest))
            if key not in index:
                if len(states)>=16384:raise ValueError('Exact augmented state cap')
                index[key]=len(states);states.append(key)
            rows.append(cursor);cols.append(index[key]);rates.append(rate);total+=rate
        rows.append(cursor);cols.append(cursor);rates.append(-total);cursor+=1
    q=coo_matrix((rates,(rows,cols)),shape=(len(states),len(states))).tocsr();p0=np.zeros(len(states));p0[0]=1.;p=expm_multiply(q.T*delta,p0)
    assert p.min()>-1e-12 and abs(p.sum()-1)<1e-10
    projected=defaultdict(float)
    for (_,out),mass in zip(states,p):
        count=np.zeros((7,n),np.int64)
        for code in out:site,sp=divmod(code,7);count[sp,site]+=1
        projected[tuple(count.ravel())]+=float(mass)
    support=np.asarray(list(projected),np.int64);prob=np.asarray(list(projected.values()))
    return dict(support=support,probabilities=prob,remaining=np.asarray([s[0] for s in states],np.uint64),outputs=np.asarray([s[1] for s in states]),q_rows=np.asarray(rows),q_cols=np.asarray(cols),q_rates=np.asarray(rates),augmented_probabilities=p)
def save_pair(name,result):
    np.savez_compressed(BASE/(name+'.npz'),**{k:v for k,v in result.items() if k not in ('trace','boundaries')})
    if result['boundaries']:atomic_json(BASE/(name+'_trace.json'),dict(trace=result['trace'],boundaries=result['boundaries']))
class Script:
    def __init__(self,wait=.001):self.wait=wait
    def exponential(self,scale):return self.wait
    def random(self):return .5
def simple(n,mask,diff=(0.,0.,0.),links=()):return Model(n,(1.,[(10.,4.)]*n,list(links),diff,(0.,0.),mask))
def scripts(deadline):
    pop=np.zeros((7,1),np.int64);pop[1,0]=2;pop[2,0]=2;m=simple(1,(False,False,False,False,True,False))
    r=pair(m,pop,pop,(.01,.005),[0,.01],0,deadline,rng=Script(),sentinel=True);assert [x['mask'] for x in r['trace']]==[3,2];assert tuple(r['counts'][:,-1,1,0])==(3,4);replay(m,pop,pop,(.01,.005),[0,.01],r);save_pair('script_R4',r)
    left=np.zeros((7,1),np.int64);left[3,0]=1;right=np.zeros_like(left);right[0,0]=1;m=simple(1,(True,False,False,False,False,False));r=pair(m,left,right,(.01,.005),[0,.01],0,deadline,rng=Script(),sentinel=True);assert [x['mask'] for x in r['trace']]==[1];replay(m,left,right,(.01,.005),[0,.01],r);save_pair('script_one_sided',r)
    r=pair(m,left,left,(.01,.005),[0,.01],0,deadline,rng=Script(.005),sentinel=True);assert r['events']==0;replay(m,left,left,(.01,.005),[0,.01],r);save_pair('script_boundary',r)
    empty=np.zeros_like(left);r=pair(m,empty,empty,(.01,.005),[0,.01],0,deadline,sentinel=True);assert r['events']==0;replay(m,empty,empty,(.01,.005),[0,.01],r);save_pair('script_empty',r)
    for kw,exc in ((dict(deadline=time.perf_counter()-1),TimeoutError),(dict(max_events=0),ValueError)):
        try:pair(m,left,left,(.01,.005),[0,.01],0,rng=Script(),**kw)
        except exc:pass
        else:raise AssertionError('Missing bound guard')
    return 6
def stage_A(p,deadline):
    algebra=0
    for c in range(4):
        m=Model(len(fixture(c).sites),config(fixture(c)))
        with np.load(ROOT/INPUTS[c+1]) as z:states=z['states']
        for i,a in enumerate(states):
            ar=np.asarray([ch.rate(a) for ch in m.channels]);br=np.asarray([ch.rate(states[(i+1)%len(states)]) for ch in m.channels]);common=np.minimum(ar,br)
            np.testing.assert_allclose(common+(ar-common),ar);np.testing.assert_allclose(common+(br-common),br);np.testing.assert_allclose(common+(ar-common)+(br-common),np.maximum(ar,br));algebra+=1
    assert algebra==896;scripted=scripts(deadline);regression=[]
    for c,seeds in enumerate(p['seeds']['regression']):
        state=fixture(c);pop=state.counts();m=Model(len(state.sites),config(state));finals=[]
        for r,seed in enumerate(seeds):
            result=pair(m,pop,pop,(.01,.01),[0,.1],seed,deadline,sentinel=True);rng=np.random.default_rng(seed);single=pop.copy();marks=[]
            for step in range(10):
                out=tick(m,single,.01,rng,deadline);single=out['counts'];marks.extend((step*.01+e['time'],e['channel']) for e in out['trace'])
            np.testing.assert_array_equal(result['counts'][0,-1],single);np.testing.assert_array_equal(result['counts'][0],result['counts'][1]);assert result['rng_state']==json.dumps(rng.bit_generator.state);assert [x['channel'] for x in result['trace']]==[x[1] for x in marks];np.testing.assert_allclose([x['time'] for x in result['trace']],[x[0] for x in marks],rtol=0,atol=1e-12);replay(m,pop,pop,(.01,.01),[0,.1],result);finals.append(single)
            if r==0:save_pair(f'regression_c{c}_r0',result)
        np.savez_compressed(BASE/f'regression_c{c}.npz',counts=finals,seeds=seeds);regression.append(len(seeds))
    marginals=[]
    for c,seeds in enumerate(p['seeds']['marginal']):
        state=fixture(c);left=state.counts();right=left.copy();right[1,0]-=1;right[1,1]+=1;m=Model(len(state.sites),config(state));refs=[]
        for s,pop in enumerate((left,right)):
            ref=exact(instantiate(state,pop),.1,deadline);np.savez_compressed(BASE/f'exact_c{c}_s{s}.npz',**ref);refs.append(ref)
        data=[];hist=[]
        for r,seed in enumerate(seeds):
            result=pair(m,left,right,(.1,.1),[0,.1],seed,deadline,sentinel=r==0);data.append(result['counts'][:,-1]);hist.append(result['histogram'])
            if r==0:replay(m,left,right,(.1,.1),[0,.1],result);save_pair(f'marginal_c{c}_r0',result)
        data=np.asarray(data);np.savez_compressed(BASE/f'marginal_c{c}.npz',counts=data,histograms=hist,seeds=seeds)
        for s,ref in enumerate(refs):
            lookup={tuple(key):i for i,key in enumerate(ref['support'])};ids=[lookup[tuple(x.ravel())] for x in data[:,s]];emp=np.bincount(ids,minlength=len(lookup))/len(seeds);tv=float(abs(emp-ref['probabilities']).sum()/2);bound=math.sqrt((len(lookup)*math.log(2)+math.log(8/.001))/(2*len(seeds)));marginals.append(dict(cell=c,side=s,support=len(lookup),TV=tv,bound=bound));assert tv<=bound
    m=simple(2,(False,)*6,(80.,0.,0.),(('ActiveLink',0,1,1.),('ActiveLink',1,0,1.)));pop=np.zeros((7,2),np.int64);pop[0,0]=1;seeds=p['seeds']['transport'];data=[];hist=[]
    for r,seed in enumerate(seeds):
        result=pair(m,pop,pop,(.005,.002),[0,.01,.02],seed,deadline,sentinel=r==0);data.append(result['counts']);hist.append(result['histogram'])
        if r==0:replay(m,pop,pop,(.005,.002),[0,.01,.02],result);save_pair('transport_r0',result)
    data=np.asarray(data);np.savez_compressed(BASE/'transport.npz',counts=data,histograms=hist,seeds=seeds);transport=[];bound=math.sqrt(math.log(8/.001)/(2*len(seeds)))
    for s,delta in enumerate((.005,.002)):
        for frame,T in enumerate((.01,.02),1):
            expected=(1-(2*math.exp(-80*delta)-1)**round(T/delta))/2;observed=float(data[:,s,frame,0,1].mean());transport.append(dict(side=s,time=T,expected=expected,observed=observed,bound=bound));assert abs(observed-expected)<=bound
    return dict(passed=True,algebra_states=algebra,scripted=scripted,identical_tick_paths=sum(regression),exact_marginals=marginals,transport=transport,exact_pairs=8192,transport_pairs=4096)
def stage_B(p,deadline):
    from mg12h3_finer_precision import setup
    initial,contrast=setup();model=contrast.model;rows=[]
    for c,deltas in enumerate(((.0001,.00002),(.00005,.00002))):
        paired=[];independent=[];runtimes=[]
        for r,seed in enumerate(p['seeds']['actual_paired'][c]):
            clock=time.perf_counter();result=pair(model,initial,initial,deltas,[0,.0001],seed,deadline,sentinel=r==0);elapsed=time.perf_counter()-clock;save_pair(f'actual_c{c}_r{r}',result)
            if r==0:replay(model,initial,initial,deltas,[0,.0001],result)
            final=result['counts'][:,-1];paired.append([contrast.direct(x) for x in final]);runtimes.append(dict(kind='paired',replicate=r,seconds=elapsed,events=result['events']))
            sides=[]
            for s,delta in enumerate(deltas):
                seed=p['seeds']['actual_single'][c][2*r+s];rng=np.random.default_rng(seed);pop=initial.copy();events=0;clock=time.perf_counter()
                for step in range(round(.0001/delta)):
                    out=tick(model,pop,delta,rng,deadline,trace=False);pop=out['counts'];events+=out['events']
                elapsed=time.perf_counter()-clock;sides.append(contrast.direct(pop));np.savez_compressed(BASE/f'actual_single_c{c}_r{r}_s{s}.npz',counts=pop,seed=seed,events=events,seconds=elapsed);runtimes.append(dict(kind='single',replicate=r,side=s,seconds=elapsed,events=events))
            independent.append(sides)
        paired=np.asarray(paired);independent=np.asarray(independent);np.savez_compressed(BASE/f'actual_contrasts_c{c}.npz',paired=paired,independent=independent,weights=contrast.weights)
        pv=np.var(paired[:,0]-paired[:,1],axis=0,ddof=1);iv=np.var(independent[:,0]-independent[:,1],axis=0,ddof=1);rows.append(dict(comparison=c,deltas=deltas,paired_difference_variance=pv.tolist(),independent_difference_variance=iv.tolist(),ratio=(pv/iv).tolist(),paired_mean_difference=np.mean(paired[:,0]-paired[:,1],axis=0).tolist(),runtimes=runtimes))
    return dict(passed=True,comparisons=rows,scope='R4 descriptive microtiming only; no variance benefit or tick equivalence verdict')
def run(stage):
    p=read(BASE/'plan.json');assert p['source_hashes']==hashes();path=BASE/f'stage_{stage}_attempt.json';assert not path.exists()
    subprocess.run(['git','ls-files','--error-unmatch',str((BASE/'plan.json').relative_to(ROOT))],cwd=ROOT,check=True,stdout=subprocess.DEVNULL)
    if stage=='B':assert read(BASE/'stage_A_results.json')['passed']
    start=time.perf_counter();deadline=start+60;atomic_json(path,dict(status='RUNNING',stage=stage))
    try:
        result=(stage_A if stage=='A' else stage_B)(p,deadline);elapsed=time.perf_counter()-start;assert elapsed<=60;assert hashes()==p['source_hashes'];result.update(elapsed_seconds=elapsed,source_hashes=hashes());atomic_json(BASE/f'stage_{stage}_results.json',result);atomic_json(path,dict(status='COMPLETE',elapsed_seconds=elapsed));print(json.dumps(result))
    except BaseException:
        atomic_json(BASE/f'stage_{stage}_failure.json',dict(error=traceback.format_exc(),elapsed_seconds=time.perf_counter()-start));atomic_json(path,dict(status='FAILED'));raise
if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--plan',action='store_true');parser.add_argument('--stage',choices=['A','B']);args=parser.parse_args()
    if args.plan:plan()
    elif args.stage:run(args.stage)
    else:parser.error('Choose plan or stage')

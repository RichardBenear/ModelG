"""MG083 prospectively bounded geometry-delivery development, not production."""
import os
for key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[key]='1'
from pathlib import Path
import json,hashlib,argparse,sys,platform,time,traceback
from mg10d1_pool import atomic_json,run
ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/'results/MG11G1/preparation'
FILES=('src/mg11g1.py','src/mg11g1_sample.py','src/prepare_mg11g1.py','src/check_mg11g1_preparation.py','spec/MG11G1_proposal.md',
       'src/mg11n1.py','src/mg11n1_stream.py','src/check_mg11n1_preparation.py','src/mg10s1.py','src/runner.py','src/mg10d1_pool.py')
def read(p):return json.loads(Path(p).read_text(encoding='utf-8'))
def digest(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def hashes():return {name:hashlib.sha256((ROOT/name).read_text(encoding='utf-8').replace('\r\n','\n').encode()).hexdigest() for name in FILES}
def candidate_jobs():
    from runner import seed_for
    return [dict(id=f'r{r}_g{g}_c{c}',resource=r,network=g,chunk=c,reps=4096,seed=seed_for(11,2*r+g,0,8303,c),warmup_seed=seed_for(11,2*r+g,0,8303,c,1)) for r in range(2) for g in range(2) for c in range(16)]
def plan():
    import numpy,scipy,numba
    from runner import seed_for
    from mg11g1 import cell,exact_probabilities,CONTACT,difference_interval
    if (BASE/'preparation_plan.json').exists():raise ValueError('Already planned; preserve/version changes')
    jobs=candidate_jobs();values=[j[key] for j in jobs for key in ('seed','warmup_seed')]
    dev=[seed_for(11,r,0,8300,0) for r in range(2)]+[seed_for(11,c,1,8300,0,seg) for c in range(4) for seg in (0,1)]+[seed_for(11,c,0,8301,0,seg) for c in range(4) for seg in (0,1)]
    assert len(values)==len(set(values))==128 and len(dev)==len(set(dev))==18 and not set(values)&set(dev)
    prior=set()
    def collect(v):
        if isinstance(v,dict):
            for key,x in v.items():
                if 'seed' in key.lower() and isinstance(x,int):prior.add(x)
                else:collect(x)
        elif isinstance(v,list):
            for x in v:collect(x)
    for directory in (ROOT/'results').glob('MG*'):
        if directory.name=='MG11G1':continue
        for path in directory.rglob('*.json'):collect(read(path))
    assert not (set(values)|set(dev))&prior
    predictions={}
    for r in range(2):
        p=[]
        for g in range(2):
            case,graph,k,t,m=cell(r,g);law=exact_probabilities(case,graph,k,r,4.);p.append(float(law[6+CONTACT] if r==0 else law[-1]))
        predictions[str(r)]=dict(probabilities=p,bipartite_minus_prism=p[1]-p[0])
    p=predictions['1']['probabilities'];design=difference_interval(round(65536*p[0]),65536,round(65536*p[1]),65536)
    atomic_json(BASE/'design_calculation.json',dict(scope='Analytical pre-Monte-Carlo design, not observed outcomes. Graphs/T4/rates fixed before calculation; R65536 selected afterward for primary precision.',predictions=predictions,primary_design_rounded_counts=[round(65536*x) for x in p],primary_design_interval=design))
    atomic_json(BASE/'preparation_plan.json',dict(status='G1 DRAFT; DEVELOPMENT ONLY',source_hashes=hashes(),candidate_jobs=jobs,development_seeds=dev,prior_recorded_integer_seed_values=len(prior),collisions=0,
       development_ceiling_seconds=900,candidate_production_ceiling_seconds=3600,primary_time=4,primary_resource='B',paths_per_cell=65536,timing_reps=1024,checkpoint_block=128,
       environment=dict(python=sys.version,platform=platform.platform(),numpy=numpy.__version__,scipy=scipy.__version__,numba=numba.__version__)))
def qualify(directory):
    import numpy as np
    from mg11g1 import cell,topology,generator,exact_probabilities,delivery,CONTACT,binomial_interval,difference_interval
    from mg11n1 import graph,simulate
    from mg11g1_sample import generate,raw_save
    from runner import seed_for
    from check_mg11n1_preparation import check_raw
    tests=[];raw=[];refs={};geometry=[]
    for gidx in range(2):
        g=topology(gidx);W=np.zeros((6,6));W[g['src'],g['dst']]=g['weight'];assert len(g['src'])==18 and np.all(W.sum(axis=1)==3) and W[0,3]==1
        geometry.append(dict(network=gidx,adjacency=W.tolist(),degree=W.sum(axis=1).tolist(),triangles=int(round(np.trace(W@W@W)/6)),source_target_distance=1))
    assert [g['triangles'] for g in geometry]==[2,0]
    for r in range(2):
        a=cell(r,0);b=cell(r,1)
        for key in a[0]:np.testing.assert_array_equal(a[0][key],b[0][key])
        for i in (2,3,4):np.testing.assert_array_equal(a[i],b[i])
    tests.append(dict(name='degree_edge_capacity_distance_stock_contact_matching',passed=True))
    for r in range(2):
        seed=seed_for(11,r,0,8300,0);previous=None
        for gidx in range(2):
            case,g,k,t,m=cell(r,gidx);result=simulate(case,g,[0.,0.],t,8,seed,m)
            if previous is not None:
                for a,b in zip(result,previous):np.testing.assert_array_equal(a,b)
            previous=result;assert not np.any(delivery(result[0],r)) and not np.any(result[1]);p=directory/f'null_r{r}_g{gidx}.npz';raw_save(p,case,g,np.zeros(2),t,m,seed,*result,8);raw.append(check_raw(p))
        tests.append(dict(name=f'zero_transport_target_null_and_graph_identity_r{r}',passed=True))
    perm=np.array([2,5,1,4,0,3]);reference_results=[]
    for r in range(2):
        for gidx in range(2):
            case,g,k,t,m=cell(r,gidx);W=np.zeros((6,6));W[g['src'],g['dst']]=g['weight'];pc={key:(value.copy() if hasattr(value,'copy') else value) for key,value in case.items()};pc['initial'][:,perm]=case['initial'];pc['a0'][perm]=case['a0'];pc['b0'][perm]=case['b0']
            pg=graph(6,[(int(perm[i]),int(perm[j]),float(w)) for i,j,w in zip(g['src'],g['dst'],g['weight']) if i<j]);mapping=np.r_[perm,6+perm] if r==0 else np.r_[perm,6]
            Q=generator(case,g,k,r);Qp=generator(pc,pg,k,r);np.testing.assert_allclose(Qp[np.ix_(mapping,mapping)],Q,atol=1e-15)
            p=exact_probabilities(case,g,k,r,4.);pp=exact_probabilities(pc,pg,k,r,4.);np.testing.assert_allclose(pp[mapping],p,atol=1e-12)
            key=f'r{r}_g{gidx}';refs[key+'_Q']=Q;refs[key+'_Q_permuted']=Qp;refs[key+'_p']=p;refs[key+'_p_permuted']=pp
            fast=exact_probabilities(case,g,np.full(2,1000.),r,4.);expected=(1-np.exp(-2))/6 if r==0 else 1-np.exp(-.8*4/6);value=fast[6+CONTACT] if r==0 else fast[-1];assert abs(value-expected)<.001
            if r==0:
                for target in t:assert abs(exact_probabilities(case,g,k,r,float(target))[6:].sum()-(1-np.exp(-.5*target)))<1e-12
            reference_results.append(dict(resource=r,network=gidx,fast_probability=fast.tolist(),fast_target=float(value),fast_target_expected=float(expected)))
    np.savez_compressed(directory/'reference_inputs.npz',permutation=perm,**refs);atomic_json(directory/'formal_references.json',dict(geometry=geometry,cases=reference_results))
    tests.append(dict(name='simultaneous_role_graph_permutation_invariance',passed=True));tests.append(dict(name='A_total_null_and_linear_fast_mixing_limits',passed=True))
    intervals=[]
    for n,k in ((1,0),(1,1),(10,0),(10,10),(10,5),(65536,16900),(65536,17700)):
        intervals.append(dict(n=n,k=k,interval=binomial_interval(k,n)))
    assert binomial_interval(0,0) is None and difference_interval(0,0,0,1)['ci'] is None
    for k,n in ((-1,1),(2,1),(0,-1)):
        try:binomial_interval(k,n);raise AssertionError('Bad count accepted')
        except ValueError:pass
    atomic_json(directory/'interval_cases.json',intervals);tests.append(dict(name='exact_binomial_boundaries_empty_invalid_inputs',passed=True))
    for r in range(2):
        for gidx in range(2):
            c=2*r+gidx;job=dict(id=f'equivalence_r{r}_g{gidx}',resource=r,network=gidx,reps=32,seed=seed_for(11,c,1,8300,0),warmup_seed=seed_for(11,c,1,8300,0,1));generate(job,directory,block=16)
            case,g,k,t,m=cell(r,gidx);expected=simulate(case,g,k,t,32,job['seed'],m)
            with np.load(directory/(job['id']+'.npz')) as z:
                for name,value in zip(('data','environment','stats'),expected):np.testing.assert_array_equal(z[name],value)
            raw.append(check_raw(directory/(job['id']+'.npz')));raw.append(check_raw(directory/(job['id']+'_warm.npz')))
    tests.append(dict(name='four_actual_writer_monolithic_seed_regressions',passed=True))
    atomic_json(directory/'qualification_results.json',dict(passed=True,tests=tests,raw_checks=raw,scope='Development/formal checks only; not scientific geometry-delivery estimation'))
def worker(job,directory):
    directory=Path(directory)
    if hashes()!=read(BASE/'preparation_plan.json')['source_hashes']:raise ValueError('Frozen source changed')
    if job['kind']=='qualification':qualify(directory);atomic_json(directory/(job['id']+'.json'),dict(status='COMPLETE',job=job));return
    from runner import seed_for
    from mg11g1_sample import generate
    r=job['resource'];g=job['network'];c=2*r+g;actual=dict(job,reps=1024,seed=seed_for(11,c,0,8301,0),warmup_seed=seed_for(11,c,0,8301,0,1));generate(actual,directory)
def execute():
    if (BASE/'qualification/ledger.json').exists():raise ValueError('Already attempted; no automatic retry')
    q=run([dict(id='qualification',kind='qualification')],BASE/'qualification',Path(__file__).resolve(),ROOT,900,workers=1,timeout=180,poll_seconds=.25)
    if q['stop_reason']:raise ValueError('Qualification failed; preserve')
    jobs=[dict(id=f'timing_r{r}_g{g}',kind='timing',resource=r,network=g) for r in range(2) for g in range(2)]
    t=run(jobs,BASE/'timing',Path(__file__).resolve(),ROOT,900,workers=4,prior=q['charged_worker_seconds'],timeout=180,poll_seconds=.25)
    atomic_json(BASE/'execution.json',dict(qualification=q,timing=t))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--plan',action='store_true');p.add_argument('--execute',action='store_true');p.add_argument('--worker',nargs=2);a=p.parse_args()
    if a.plan:plan()
    elif a.execute:execute()
    elif a.worker:
        j=json.loads(a.worker[0]);d=Path(a.worker[1])
        try:worker(j,d)
        except BaseException:atomic_json(d/(j['id']+'.json'),dict(status='FAILED',job=j,error=traceback.format_exc()));raise
    else:p.error('Development only: plan/execute/worker')

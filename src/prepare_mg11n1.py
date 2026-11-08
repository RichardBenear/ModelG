"""MG079 prospective network qualification development; no production CLI."""
import os
for key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[key]='1'
from pathlib import Path
import argparse,json,hashlib,platform,sys,time,traceback
from mg10d1_pool import atomic_json,run
ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/'results/MG11N1/preparation'
FILES=('src/mg11n1.py','src/mg11n1_cases.py','src/mg11n1_oracle.py','src/prepare_mg11n1.py','src/check_mg11n1_preparation.py','spec/MG11N1_proposal.md',
       'src/mg10s1.py','src/mg10_reservoir.py','src/mg10p2.py','src/mg10d2.py','src/mg10s2.py','src/mg10p1.py','src/runner.py','src/mg10d1_pool.py')
def digest(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def hashes():
    out={name:hashlib.sha256((ROOT/name).read_text(encoding='utf-8').replace('\r\n','\n').encode()).hexdigest() for name in FILES}
    out['results/MG003/ref_profiles.npz']=digest(ROOT/'results/MG003/ref_profiles.npz');return out
def candidate_jobs():
    from runner import seed_for
    jobs=[]
    for c in range(4):
        for engine in range(2 if c==3 else 1):
            for chunk in range(4 if c==3 else 8):jobs.append(dict(id=f'c{c}_m{engine}_r{chunk}',cell=c,engine=engine,chunk=chunk,reps=4096,seed=seed_for(11,c,engine,7903,chunk),warmup_seed=seed_for(11,c,engine,7903,chunk,1)))
    return jobs
def plan():
    import numpy,scipy,numba
    from runner import seed_for
    if (BASE/'preparation_plan.json').exists():raise ValueError('Already planned')
    jobs=candidate_jobs();seeds=[j[key] for j in jobs for key in ('seed','warmup_seed')];dev=[seed_for(11,c,0,7900,0) for c in range(5)]+[seed_for(11,c,m,7901,0,seg) for c in range(4) for m in range(2 if c==3 else 1) for seg in (0,1)]
    assert len(jobs)==32 and len(set(seeds))==64 and len(set(dev))==15 and not set(seeds)&set(dev)
    prior=set()
    def collect(v):
        if isinstance(v,dict):
            for key,value in v.items():
                if 'seed' in key.lower() and isinstance(value,int):prior.add(value)
                else:collect(value)
        elif isinstance(v,list):
            for x in v:collect(x)
    for directory in (ROOT/'results').glob('MG*'):
        if directory.name=='MG11N1':continue
        for path in directory.rglob('*.json'):collect(json.loads(path.read_text(encoding='utf-8')))
    assert not (set(seeds)|set(dev))&prior
    atomic_json(BASE/'preparation_plan.json',dict(status='DRAFT NETWORK QUALIFICATION; DEVELOPMENT ONLY',source_hashes=hashes(),candidate_jobs=jobs,development_seeds=dev,prior_recorded_seed_count=len(prior),collisions=0,
      development_ceiling_seconds=900,candidate_production_cap_seconds=3600,environment=dict(python=sys.version,platform=platform.platform(),numpy=numpy.__version__,scipy=scipy.__version__,numba=numba.__version__)))
def save_raw(path,case,g,kappa,times,mask,seed,result):
    import numpy as np
    data,environment,stats=result
    np.savez_compressed(path,data=data,environment=environment,stats=stats,initial=case['initial'],a0=case['a0'],b0=case['b0'],omega=case['omega'],diff=case['diff'],src=g['src'],dst=g['dst'],weight=g['weight'],kappa=kappa,times=times,mask=mask,seed=seed)
def qualify(directory):
    import numpy as np
    from scipy.linalg import expm
    from mg11n1 import graph,simulate,continuum_rhs
    from mg11n1_cases import cell
    from mg11n1_oracle import linear_law
    from mg10s1 import fixture,ssa_chunk
    from runner import seed_for
    from check_mg11n1_preparation import check_raw
    tests=[];raw=[]
    for i in range(4):
        c=fixture(i);case=dict(initial=c['state'],a0=c['a0'],b0=c['b0'],omega=c['omega'],diff=c['diff']);g=graph(case['initial'].shape[1],[]);times=np.array([0.,.01,.05]);seed=seed_for(11,i,0,7900,0)
        new=simulate(case,g,[0.,0.],times,8,seed);old,stats=ssa_chunk(case['initial'],case['a0'],case['b0'],case['omega'],case['diff'],np.ones(7),times,8,seed)
        np.testing.assert_array_equal(new[0],old);np.testing.assert_array_equal(new[2],stats);assert not np.any(new[1])
        p=directory/f'disconnected_{i}.npz';save_raw(p,case,g,np.zeros(2),times,np.ones(7),seed,new);raw.append(check_raw(p));tests.append(dict(name=f'legacy_finite_identity_{i}',passed=True))
    # Invalid inputs reject before any sampling.
    for edges in ([(0,0,1)],[(0,1,-1)],[(0,1,1),(1,0,1)],[(0,3,1)]):
        try:graph(3,edges);raise AssertionError('Invalid graph accepted')
        except ValueError:pass
    case,g,k,t,m=cell(0)
    for bad in ([-1.,0.],[float('nan'),0.]):
        try:simulate(case,g,bad,t,1,1,m);raise AssertionError('Bad rate accepted')
        except ValueError:pass
    asymmetric=dict(g);asymmetric['dst']=g['dst'].copy();asymmetric['dst'][0]=2
    try:simulate(case,asymmetric,k,t,1,1,m);raise AssertionError('Asymmetric graph accepted')
    except ValueError:pass
    tests.append(dict(name='graph_domain_duplicate_symmetry_rate_rejection',passed=True))
    # Pointwise finite continuum reference, not finite-count moment closure.
    import mg10_reservoir as ref
    old=ref.initial('H',3);old[0]=[15,16,17];old[1]=[14,15,16];old[2]=[1,1.2,1.5];old[3]=[.1,.2,.3];old[4]=[.05,.1,.15];a0,b0=ref.stocks(1)
    z=np.vstack([old[:3],a0*(1-old[3]),b0*(1-old[4]),np.zeros((2,3))]);case=dict(initial=np.zeros((7,3),np.int64),a0=np.full(3,a0),b0=np.full(3,b0),omega=1.,diff=np.array([25.,25.,300.]))
    actual=continuum_rhs(0,z.ravel(),case,g,[0.,0.]).reshape(7,3);r=ref.ode_rhs(0,old.ravel(),1,True).reshape(15,3);expected=np.vstack([r[:3],-a0*r[3],-b0*r[4],r[8],r[10]])
    np.testing.assert_allclose(actual,expected,rtol=1e-12,atol=1e-12)
    np.savez_compressed(directory/'continuum_inputs.npz',old=old,state=z,a0=a0,b0=b0,actual=actual,expected=expected)
    tests.append(dict(name='zero_hop_continuum_finite_reference',passed=True))
    case,g,k,t,m=cell(1);lam=14/case['a0'];predicted=np.r_[np.full(3,np.exp(-lam.mean())/3),lam/lam.sum()*(1-np.exp(-lam.mean()))];mix=[]
    for rate in (1.,10.,100.,1000.):
        law=linear_law(case,g,[rate,0.],'A_source',1.);mix.append(dict(kappa=rate,law=law.tolist(),error=float(abs(law-predicted).max())))
    assert all(mix[i+1]['error']<mix[i]['error'] for i in range(3)) and mix[-1]['error']<.001
    capacity=[]
    for c in (1.,100.,10000.):capacity.append(dict(capacity=c,output=float(14*c*(-np.expm1(-1/c))),relative_error=float(abs(c*(-np.expm1(-1/c))-1))))
    assert all(capacity[i+1]['relative_error']<capacity[i]['relative_error'] for i in range(2)) and capacity[-1]['relative_error']<.0001
    W=np.array([[-1.,1,0,0],[1,-1,0,0],[0,0,-1,1],[0,0,1,-1]])
    law=np.array([1.,0,0,0])@expm(1000*W);np.testing.assert_allclose(law,[.5,.5,0,0],atol=1e-12)
    atomic_json(directory/'linear_limits.json',dict(mixing=mix,well_mixed=predicted.tolist(),capacity=capacity,disconnected_generator=W.tolist(),disconnected_law=law.tolist()))
    tests.append(dict(name='linear_component_fast_mixing_capacity_limits',passed=True))
    # Size/API smoke only: physically tiny horizon at500nodes.
    from mg10p2 import fixture as full
    original_full=full(0,0);case={key:original_full[key] for key in ('initial','a0','b0','omega','diff')};g=graph(500,[(i,(i+1)%500,1.) for i in range(500)]);times=np.array([0.,1e-6]);seed=seed_for(11,4,0,7900,0)
    result=simulate(case,g,[1e-9,1e-9],times,1,seed);p=directory/'full500_smoke.npz';save_raw(p,case,g,np.full(2,1e-9),times,np.ones(7),seed,result);raw.append(check_raw(p));tests.append(dict(name='full500_dimensions_and_ledger_only',passed=True))
    atomic_json(directory/'qualification_results.json',dict(passed=True,tests=tests,raw_checks=raw,scope='Development formal/limiting-case checks, not production distribution verdict'))
def worker(job,directory):
    directory=Path(directory)
    if hashes()!=json.loads((BASE/'preparation_plan.json').read_text(encoding='utf-8'))['source_hashes']:raise ValueError('Frozen source changed')
    if job['kind']=='qualification':qualify(directory);atomic_json(directory/(job['id']+'.json'),dict(status='COMPLETE',job=job));return
    from mg11n1 import simulate
    from mg11n1_oracle import direct_ssa
    from mg11n1_cases import cell
    from runner import seed_for
    from check_mg11n1_preparation import check_raw
    case,g,k,t,m=cell(job['cell']);engine=job['engine'];fn=simulate if engine==0 else direct_ssa;reps=256 if engine==0 else 128
    warm=seed_for(11,job['cell'],engine,7901,0,1);seed=seed_for(11,job['cell'],engine,7901,0)
    result=fn(case,g,k,t,1,warm,m);save_raw(directory/(job['id']+'_warm.npz'),case,g,k,t,m,warm,result)
    clock=time.perf_counter();result=fn(case,g,k,t,reps,seed,m);elapsed=time.perf_counter()-clock
    p=directory/(job['id']+'.npz');save_raw(p,case,g,k,t,m,seed,result);checked=check_raw(p)
    atomic_json(p.with_suffix('.json'),dict(status='COMPLETE',job=job,reps=reps,compute_seconds=elapsed,seed=seed,raw=p.name,sha256=digest(p),raw_check=checked))
def execute():
    if (BASE/'qualification/ledger.json').exists():raise ValueError('Already attempted; no automatic retry')
    q=run([dict(id='qualification',kind='qualification')],BASE/'qualification',Path(__file__).resolve(),ROOT,900,workers=1,timeout=180,poll_seconds=.25)
    if q['stop_reason']:raise ValueError('Qualification failure; preserve')
    jobs=[dict(id=f'timing_c{c}_m{m}',kind='timing',cell=c,engine=m) for c in range(4) for m in range(2 if c==3 else 1)]
    t=run(jobs,BASE/'timing',Path(__file__).resolve(),ROOT,900,workers=5,prior=q['charged_worker_seconds'],timeout=180,poll_seconds=.25)
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

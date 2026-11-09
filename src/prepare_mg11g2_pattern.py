"""MG088 prospective deterministic qualification, actual-writer timing, pipeline fixtures."""
import os
for key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[key]='1'
from pathlib import Path
import argparse,json,hashlib,time,traceback,copy,shutil,subprocess
import numpy as np
from mg11g2_pattern import factory,environment,sparsity,solve_pattern,reference,metrics,OBS,S0,CONFIGS,jobs
from mg11g2 import initial,rhs,reconstruct
from mg11g2_pattern_oracle import population_rhs,independent_solve,incidence_check,dense_projection,late
from run_mg11g2_pattern import generate,raw_save,digest,read
from evaluate_mg11g2_pattern import main as evaluate,decision
from check_mg11g2_pattern import main as audit
from mg10d1_pool import atomic_json
ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/'results/MG11G2/pattern_response/preparation'
FILES=('src/mg11g2_pattern.py','src/mg11g2_pattern_oracle.py','src/prepare_mg11g2_pattern.py','src/evaluate_mg11g2_pattern.py','src/check_mg11g2_pattern.py','src/run_mg11g2_pattern.py','spec/MG11G2_pattern_response.md','src/mg11g2.py','src/mg11n1.py','src/mg10s1.py','src/mg10_reservoir.py','src/mg10p1.py','src/mg10d1_pool.py','src/mg10e2c1_pool.py')
SHORT=np.array([0.,.125,.25])
def hashes():return {name:hashlib.sha256((ROOT/name).read_text(encoding='utf-8').replace('\r\n','\n').encode()).hexdigest() for name in FILES}
def plan():
    if (BASE/'plan.json').exists():raise ValueError('Already planned')
    atomic_json(BASE/'plan.json',dict(stage='MG088',status='PROSPECTIVE DEVELOPMENT ONLY',pi_request='Go ahead with the Pattern-Response study.',jobs=jobs(),short_times=SHORT.tolist(),source_hashes=hashes(),reference_sha256=digest(ROOT/'results/MG003/ref_profiles.npz'),ceiling_elapsed_seconds=900,job_timeout_seconds=120,workers=1,threads_per_library=1,random_draws=0,production_candidate_worker_seconds=7200))

def geometry(g):
    adjacency=[[] for _ in range(g['n'])]
    for i,j in zip(g['src'],g['dst']):adjacency[i].append(j)
    colors={0:0};queue=[0];bipartite=True
    while queue:
        i=queue.pop(0)
        for j in adjacency[i]:
            if j not in colors:colors[j]=1-colors[i];queue.append(j)
            elif colors[j]==colors[i]:bipartite=False
    return dict(degree=[len(v) for v in adjacency],connected=len(colors)==g['n'],bipartite=bipartite)

def formal():
    reports={}
    for kind in ('LOCAL','CROSS_RING'):
        g=environment(kind);v=geometry(g);assert v['connected'] and set(v['degree'])=={3} and len(g['src'])==1500 and np.all(g['weight']==1);assert v['bipartite']==(kind=='LOCAL');reports[kind]=dict(nodes=500,directed_edges=1500,degree=3,connected=True,bipartite=v['bipartite'],outgoing_hazard=1.)
    # Ledger, independent RHS and Jacobian dependencies on a nonuniform N8 field.
    case,g,k,mode=factory('LOCAL_C1',0,n=8);y=initial(case,g);pop=y[:56].reshape(7,8);phase=np.arange(8)
    for s in range(5):pop[s]*=.9+.03*np.sin(phase+s)
    dy=rhs(0,y,case,g,k,mode);scale=np.maximum(1,pop.max(axis=1));ledger=reconstruct(dict(case,initial=np.zeros_like(case['initial'])),g,dy);dpop=dy[:56].reshape(7,8);assert np.max(abs(dpop-ledger)/scale[:,None])<=1e-11;assert abs(dpop.sum())<=1e-10
    other=population_rhs(0,y[:56],case,g,k,mode).reshape(7,8);err=float(np.max(abs(other-dpop)/scale[:,None]));assert err<=1e-11
    pattern=sparsity(g).tocsc();checked=0
    for col in range(56):
        perturbed=y.copy();epsilon=1e-6*max(1,abs(y[col]));perturbed[col]+=epsilon;derivative=(rhs(0,perturbed,case,g,k,mode)-dy)/epsilon;detected=np.flatnonzero(abs(derivative)>1e-8);allowed=pattern.indices[pattern.indptr[col]:pattern.indptr[col+1]];assert set(detected)<=set(allowed);checked+=1
    # Active-ring rotation and graph+role relabeling must commute with the drift.
    permutation=np.roll(np.arange(8),1);inverse=np.argsort(permutation);rotg=dict(g,src=inverse[g['src']],dst=inverse[g['dst']]);rotcase=dict(case,initial=case['initial'][:,permutation],a0=case['a0'][permutation],b0=case['b0'][permutation]);rot=population_rhs(0,pop[:,permutation].ravel(),rotcase,rotg,k,mode).reshape(7,8);np.testing.assert_allclose(rot,other[:,permutation],rtol=1e-10,atol=1e-10)
    from mg10_reservoir import ode_rhs
    c,g0,k0,m0=factory('ISOLATED_C1',0,n=8);old=np.zeros((15,8));old[:3]=pop[:3];old[3]=1-pop[3]/c['a0'];old[4]=1-pop[4]/c['b0'];legacy=ode_rhs(0,old.ravel(),1,True).reshape(15,8);zero=population_rhs(0,pop.ravel(),c,g0,k0,m0).reshape(7,8);np.testing.assert_allclose(zero[:3],legacy[:3],rtol=1e-11,atol=1e-11)
    # Raw reference is a deviation, never an absolute concentration field.
    c,g,k,mode=factory('MAINTAINED',1);np.testing.assert_array_equal(c['initial'][:3],S0[:,None]+reference()[0]);assert c['initial'].min()>=0
    synthetic=np.repeat(c['initial'][None],101,axis=0);alpha=.2+.01*OBS;synthetic[:,:3]=S0[None,:,None]+alpha[:,None,None]*reference()[0][None]+3
    instrument=metrics(synthetic,OBS);np.testing.assert_allclose(instrument['projection'],np.repeat(alpha[:,None],3,axis=1),rtol=1e-10,atol=1e-10);np.testing.assert_allclose(instrument['late_projection'],[.425]*3,rtol=1e-10,atol=1e-10);np.testing.assert_allclose(late(np.c_[OBS,OBS*0+1,OBS*2],OBS),[22.5,1,45],atol=1e-12)
    np.testing.assert_allclose(dense_projection(synthetic,reference()[0]),instrument['projection'],rtol=1e-9,atol=1e-9)
    return dict(graphs=reports,derivative_reference_error=err,jacobian_columns_checked=checked,ledger_rotation_legacy_reference_projection_trapezoid=True)

def archive(base,label):
    for name in ('evaluation.json','independent_check.json','MG11G2_pattern_results.md'):(base/name).rename(base/(label+'_'+name))

def pipeline(deadline):
    base=BASE/'SYNTHETIC_pipeline';directory=base/'production';directory.mkdir(parents=True);atomic_json(base/'registration.json',dict(jobs=jobs(),scope='ALGEBRAIC STATIC LEDGER FIXTURES, NOT SCIENTIFIC SOLUTIONS'))
    for job in jobs():
        if time.perf_counter()>deadline:raise TimeoutError('Development cap')
        case,g,k,mode=factory(job['config'],job['start']);state=np.zeros((101,22*g['n']+2*len(g['src'])));state[:,:7*g['n']]=case['initial'].ravel();path=directory/(job['id']+'.npz');raw_save(path,job,case,g,k,mode,OBS,state,OBS);atomic_json(path.with_suffix('.json'),dict(status='COMPLETE',job=job,sha256=digest(path),scope='ALGEBRAIC ONLY'))
    e=evaluate(base);audit(base,synthetic=True);assert e['complete'] and e['numerical_qualified'] and not e['primary_resolved'];archive(base,'complete')
    # Pure estimand and forced-background invariance, no scientific outcomes.
    records=copy.deepcopy(e['records']);refine=e['refinement']
    for version in (0,1):records[f'c3_s1_v{version}']['metrics']['late_projection'][1]+=.03
    decision1=decision(records,refine);assert decision1['primary_resolved']
    for row in records.values():row['metrics']['late_projection']=(np.array(row['metrics']['late_projection'])+.4).tolist()
    decision2=decision(records,refine);np.testing.assert_allclose(decision1['primary_GXY_by_setting'],decision2['primary_GXY_by_setting'],atol=1e-12)
    # Wrong reported primary must be rejected despite valid raw ledgers.
    wrong=copy.deepcopy(e);wrong['primary_GXY_by_setting'][1][1]+=.1;atomic_json(base/'evaluation.json',wrong)
    try:audit(base,synthetic=True)
    except AssertionError:rejected=True
    else:raise AssertionError('Wrong primary accepted')
    (base/'evaluation.json').rename(base/'wrong_primary_evaluation.json');atomic_json(base/'wrong_primary_rejected.json',dict(rejected=True,raw_accounting_unchanged=True))
    path=directory/'c2_s1_v1.npz';backup=path.with_suffix('.npz.withheld');path.rename(backup);evaluate(base);audit(base,synthetic=True);v=read(base/'evaluation.json');assert not v['complete'] and v['primary_GXY_by_setting'][1] is None and not v['primary_resolved'];archive(base,'missing')
    with np.load(backup) as raw:prefix={key:raw[key] for key in raw.files}
    prefix['times']=prefix['times'][:40];prefix['state']=prefix['state'][:40]
    with path.open('wb') as f:np.savez_compressed(f,**prefix)
    evaluate(base);audit(base,synthetic=True);v=read(base/'evaluation.json');assert v['jobs_complete']==27 and not v['primary_resolved'];archive(base,'prefix');path.rename(path.with_suffix('.npz.prefix_fixture'));backup.rename(path)
    empty=BASE/'SYNTHETIC_empty';empty.mkdir();atomic_json(empty/'registration.json',dict(jobs=jobs()));evaluate(empty);audit(empty,synthetic=True);assert read(empty/'evaluation.json')['jobs_complete']==0
    atomic_json(BASE/'pipeline_qualification.json',dict(passed=True,source_hashes=hashes(),complete_missing_prefix_empty_wrong_primary=True,synthetic_scope='Constant-population, zero-ledger algebraic fixtures; not ODE samples',estimand_background_invariance=True))

def execute():
    started=time.perf_counter();deadline=started+900;plan=read(BASE/'plan.json')
    if hashes()!=plan['source_hashes'] or jobs()!=plan['jobs'] or digest(ROOT/'results/MG003/ref_profiles.npz')!=plan['reference_sha256']:raise ValueError('Source/reference/job change')
    subprocess.run(['git','ls-files','--error-unmatch','results/MG11G2/pattern_response/preparation/plan.json'],cwd=ROOT,check=True,stdout=subprocess.DEVNULL)
    if (BASE/'attempt.json').exists():raise ValueError('Already attempted; no automatic retry')
    atomic_json(BASE/'attempt.json',dict(status='RUNNING DEVELOPMENT ONLY'));formal_results=formal();atomic_json(BASE/'formal.json',formal_results);timings={};rawchecks=[];comparisons=[]
    for job in jobs():
        if time.perf_counter()>deadline:raise TimeoutError('Development cap')
        meta=generate(job,BASE/'timing',times=SHORT,deadline=min(deadline,time.perf_counter()+120),scope='DISCLOSED T.25 DEVELOPMENT, NOT SCIENCE');timings[job['id']]=meta;print('Timed '+job['id'],flush=True)
    for c,config in enumerate(CONFIGS):
        for start in (0,1):
            case,g,k,mode=factory(config,start);reference_pop=independent_solve(case,g,k,mode,SHORT,min(deadline,time.perf_counter()+120));np.savez_compressed(BASE/f'reference_c{c}_s{start}.npz',times=SHORT,population=reference_pop)
            for setting in (0,1):
                id=f'c{c}_s{start}_v{setting}';path=BASE/'timing'/(id+'.npz')
                with np.load(path) as raw:
                    np.testing.assert_array_equal(raw['requested_times'],SHORT);np.testing.assert_array_equal(raw['times'],SHORT);state=raw['state'];pop=state[:,:3500].reshape(3,7,500)
                    for name,value in dict(initial=case['initial'],a0=case['a0'],b0=case['b0'],src=g['src'],dst=g['dst'],weight=g['weight'],kappa=k).items():np.testing.assert_array_equal(raw[name],value)
                error=float(np.max(abs(pop-reference_pop)/np.maximum(1,case['initial'].max(axis=1))[None,:,None]));assert error<=2e-6;v=incidence_check(state,case,g);assert v['passed'];rawchecks.append(dict(file=path.relative_to(ROOT).as_posix(),sha256=digest(path),**v));comparisons.append(dict(job=id,error=error))
                if start==0:np.testing.assert_array_equal(pop,np.broadcast_to(pop[:,:,:1],pop.shape))
                target=case['a0'][None]*np.exp(-14*SHORT[:,None]/case['a0'][None]) if mode=='FINITE' else np.repeat(case['a0'][None],3,axis=0);assert np.max(abs(pop[:,3]-target)/case['a0'][None])<=1e-9
    pipeline(deadline);forecasts={id:meta['non_io_seconds']*100*(3 if meta['job']['start'] else 2)+meta['io_seconds']*450*(3 if meta['job']['start'] else 2)+(15 if meta['job']['start'] else 5) for id,meta in timings.items()};forecast=float(sum(forecasts.values()));elapsed=time.perf_counter()-started;assert elapsed<900 and shutil.disk_usage(ROOT).free>=2*1024**3
    reserve=forecast<=5400;result=dict(passed=True,production_within_reserve=reserve,status='DEVELOPMENT/PREPARATION PASS; SCIENTIFIC BUDGET/CRITERIA APPROVAL PENDING',source_hashes=hashes(),reference_sha256=plan['reference_sha256'],elapsed_seconds=elapsed,formal=formal_results,raw_checks=rawchecks,independent_comparisons=comparisons,forecast_seconds_per_job=forecasts,forecast_worker_seconds=forecast,forecast_core_hours=forecast/3600,eight_worker_wall_minutes=forecast/480,candidate_cap_core_hours=2,raw_array_upper_bytes=316736000,scope='28disclosed T.25 writer paths+14separate short population references; algebraic pipeline fixtures only, no T25science')
    atomic_json(BASE/'preparation_results.json',result);atomic_json(BASE/'attempt.json',dict(status='COMPLETE DEVELOPMENT; no scientific production',elapsed_seconds=elapsed));print(json.dumps({k:v for k,v in result.items() if k in ('passed','production_within_reserve','elapsed_seconds','forecast_worker_seconds','forecast_core_hours','eight_worker_wall_minutes')},indent=2))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--plan',action='store_true');p.add_argument('--execute',action='store_true');a=p.parse_args()
    if a.plan:plan()
    elif a.execute:
        try:execute()
        except BaseException:atomic_json(BASE/'failure.json',dict(status='FAILED; preserve/version, no automatic retry',error=traceback.format_exc()));raise
    else:p.error('Choose plan/execute')

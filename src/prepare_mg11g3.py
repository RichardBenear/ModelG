"""MG092 prospective donor-fork/switch qualification and actual-writer forecast."""
import os
for key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[key]='1'
from pathlib import Path
import argparse,time,json,hashlib,traceback,copy,shutil,subprocess
import numpy as np
from mg11g3 import factory,donor_path,donor_manifest,jobs,rhs,reference_rhs,pattern,OBS,PARENT_INDEX,SWITCH_TIME
from mg11g3_oracle import accounting,independent_solve
from mg11g2 import initial,rhs as old_rhs,reconstruct
from mg11g2_pattern_v2 import factory as old_factory,environment,S0
from run_mg11g3 import generate,raw_save,digest,read
from evaluate_mg11g3 import main as evaluate,decide
from check_mg11g3 import main as audit
from mg10d1_pool import atomic_json
ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/'results/MG11G3/preparation';SHORT=np.array([0.,.125,.25])
FILES=('src/mg11g3.py','src/mg11g3_oracle.py','src/run_mg11g3.py','src/evaluate_mg11g3.py','src/check_mg11g3.py','src/prepare_mg11g3.py','spec/MG11G3_B_transport_switch.md','src/mg11g2.py','src/mg11n1.py','src/mg10s1.py','src/mg11g2_pattern_v2.py','src/mg11g2_pattern_v2_oracle.py','src/mg10p1.py','src/mg10d1_pool.py','src/mg10e2c1_pool.py')
def hashes():return {name:hashlib.sha256((ROOT/name).read_text(encoding='utf-8').replace('\r\n','\n').encode()).hexdigest() for name in FILES}
def plan():
    if (BASE/'plan.json').exists():raise ValueError('Already planned')
    atomic_json(BASE/'plan.json',dict(stage='MG092',pi_request='I sign off. Go ahead with the next stage.',scope='DEVELOPMENT ONLY; scientific budget/criteria approval pending',jobs=jobs(),source_hashes=hashes(),donor_hashes=donor_manifest(),short_times=SHORT.tolist(),ceiling_elapsed_seconds=900,workers=1,threads_per_library=1,random_draws=0,candidate_production_cap_worker_seconds=3600))

def formal():
    # Forks copy ALL physical stores, preserve fixed denominators; counters reset.
    manifest=donor_manifest()
    for job in jobs():
        case,ga,gb,k=factory(job)
        with np.load(donor_path(job)) as raw:np.testing.assert_array_equal(case['initial'].ravel(),raw['state'][40,:3500]);np.testing.assert_array_equal(case['a0'],raw['a0']);np.testing.assert_array_equal(case['b0'],raw['b0'])
        assert np.all(case['initial'][5:]>0) and not np.any(initial(case,ga)[3500:])
        other,*_=factory(job);case['initial'][0,0]+=1;assert case['initial'][0,0]!=other['initial'][0,0]
    assert manifest==donor_manifest()
    case,ga,k,mode=old_factory('LOCAL_C1',0,n=8);gb=environment('CROSS_RING',8);y=initial(case,ga);pop=y[:56].reshape(7,8)
    for s in range(5):pop[s]*=.85+.05*np.sin(np.arange(8)+s)
    actual=rhs(0,y,case,ga,ga,k);expected=old_rhs(0,y,case,ga,k,'FINITE');np.testing.assert_allclose(actual,expected,rtol=1e-12,atol=1e-10)
    switched=rhs(0,y,case,ga,gb,k);np.testing.assert_array_equal(switched[:32],actual[:32]);np.testing.assert_array_equal(switched[40:56],actual[40:56]);np.testing.assert_array_equal(switched[160:184],actual[160:184])
    off=rhs(0,y,case,ga,gb,[k[0],0]);assert not np.any(off[184:208])
    shape=pattern(ga,gb).tocsc();checked=0
    for column in range(56):
        perturbed=y.copy();h=1e-6*max(1,abs(y[column]));perturbed[column]+=h;d=(rhs(0,perturbed,case,ga,gb,k)-switched)/h;nonzero=np.flatnonzero(abs(d)>1e-8);allowed=shape.indices[shape.indptr[column]:shape.indptr[column+1]];assert set(nonzero)<=set(allowed);checked+=1
    # Ledger identity through separate graph incidence at a short Euler algebraic
    # increment: no integration and no clipped negative populations.
    nonuniform=dict(case,initial=pop.copy());start_state=initial(nonuniform,ga);algebraic=np.stack((start_state,start_state+1e-8*rhs(0,start_state,nonuniform,ga,gb,k)));assert accounting(algebraic,nonuniform,ga,gb)['passed']
    separate=reference_rhs(0,y[:56],case,ga,gb,k).reshape(7,8);np.testing.assert_allclose(separate,switched[:56].reshape(7,8),rtol=1e-10,atol=1e-10)
    empty=y.copy();empty[24:40]=0;d=rhs(0,empty,case,ga,gb,k);assert np.all(d[24:40]==0) and np.all(d[56:64]==0) and np.all(d[80:88]==0)
    return dict(donor_files=len(manifest),forks_verified=24,same_graph_identity=True,A_law_and_only_B_intervention=True,B_off_and_empty_stock=True,jacobian_columns_checked=checked,two_graph_ledger_identity=True)

def archive(base,label):
    for name in ('evaluation.json','independent_check.json','MG11G3_results.md'):(base/name).rename(base/(label+'_'+name))

def pipeline(deadline):
    base=BASE/'SYNTHETIC_pipeline';directory=base/'production';directory.mkdir(parents=True);atomic_json(base/'registration.json',dict(jobs=jobs(),scope='ALGEBRAIC STATIC POST-FORK FIXTURES, NOT ODE SAMPLES'))
    for job in jobs():
        if time.perf_counter()>deadline:raise TimeoutError('Development cap')
        case,ga,gb,k=factory(job);state=np.zeros((61,14000));state[:,:3500]=case['initial'].ravel();p=directory/(job['id']+'.npz');raw_save(p,job,case,ga,gb,k,OBS,state,OBS,digest(donor_path(job)));atomic_json(p.with_suffix('.json'),dict(status='COMPLETE',job=job,sha256=digest(p),scope='ALGEBRAIC ONLY'))
    e=evaluate(base);audit(base,synthetic=True);assert e['complete'] and e['numerical_qualified'] and not e['sham_qualified'] and not e['primary_resolved'];archive(base,'complete_failed_shams')
    # Estimand algebra: equal deterministic weighting and background cancellation.
    records=copy.deepcopy(e['records'])
    for job in jobs():records[job['id']]['metrics']['late_projection']=[float(job['start'])]*3;records[job['id']]['metrics']['late_projection'][1]+=.03*job['start']*(job['arm']==1)
    passing_shams=[dict(x,active_error=0.) for x in e['shams']];one=decide(records,e['refinement'],passing_shams);assert one['primary_resolved'];np.testing.assert_allclose(one['primary_GXY_by_setting'],[[0,.03,0]]*2,atol=1e-12)
    for row in records.values():row['metrics']['late_projection']=(np.array(row['metrics']['late_projection'])+.4).tolist()
    two=decide(records,e['refinement'],passing_shams);np.testing.assert_allclose(one['primary_GXY_by_setting'],two['primary_GXY_by_setting'],atol=1e-12)
    wrong=copy.deepcopy(e);wrong['primary_GXY_by_setting'][1][1]+=.1;atomic_json(base/'evaluation.json',wrong)
    try:audit(base,synthetic=True)
    except AssertionError:pass
    else:raise AssertionError('Wrong primary accepted')
    (base/'evaluation.json').rename(base/'wrong_primary_evaluation.json');atomic_json(base/'wrong_primary_rejected.json',dict(rejected=True,raw_accounting_unchanged=True))
    p=directory/'o0_s1_v1_a1.npz';backup=p.with_suffix('.npz.withheld');p.rename(backup);evaluate(base);audit(base,synthetic=True);v=read(base/'evaluation.json');assert not v['complete'] and v['primary_GXY_by_setting'][1] is None and not v['primary_resolved'];archive(base,'missing')
    with np.load(backup) as raw:partial={key:raw[key] for key in raw.files}
    partial['times']=partial['times'][:20];partial['state']=partial['state'][:20]
    with p.open('wb') as f:np.savez_compressed(f,**partial)
    evaluate(base);audit(base,synthetic=True);v=read(base/'evaluation.json');assert v['jobs_complete']==23 and not v['primary_resolved'];archive(base,'prefix');p.rename(p.with_suffix('.npz.prefix_fixture'));backup.rename(p)
    empty=BASE/'SYNTHETIC_empty';empty.mkdir();atomic_json(empty/'registration.json',dict(jobs=jobs()));evaluate(empty);audit(empty,synthetic=True);assert read(empty/'evaluation.json')['jobs_complete']==0
    atomic_json(BASE/'pipeline_qualification.json',dict(passed=True,source_hashes=hashes(),complete_failed_shams_missing_prefix_empty_wrong_primary=True,estimand_equal_weighting_background_invariance=True,synthetic_scope='Static physical stores and zero post ledgers; NOT ODE samples. Sham failure must withhold interpretation.'))

def execute():
    start=time.perf_counter();deadline=start+900;plan=read(BASE/'plan.json')
    if hashes()!=plan['source_hashes'] or donor_manifest()!=plan['donor_hashes']:raise ValueError('Source/donor change')
    subprocess.run(['git','ls-files','--error-unmatch','results/MG11G3/preparation/plan.json'],cwd=ROOT,check=True,stdout=subprocess.DEVNULL)
    if (BASE/'attempt.json').exists():raise ValueError('Already attempted; no retry')
    atomic_json(BASE/'attempt.json',dict(status='RUNNING DEVELOPMENT'));formal_result=formal();atomic_json(BASE/'formal.json',formal_result);timings={};checks=[];comparisons=[];shams=[]
    for job in jobs():
        meta=generate(job,BASE/'timing',times=SHORT,deadline=min(deadline,time.perf_counter()+120),scope='DISCLOSED POSTT.25 DEVELOPMENT, NOT SCIENCE');timings[job['id']]=meta;print('Timed '+job['id'],flush=True)
    for job in jobs():
        case,ga,gb,k=factory(job);ref=independent_solve(case,ga,gb,k,SHORT,min(deadline,time.perf_counter()+120));np.savez_compressed(BASE/(job['id']+'_reference.npz'),times=SHORT,population=ref)
        p=BASE/'timing'/(job['id']+'.npz')
        with np.load(p) as raw:
            state=raw['state'];np.testing.assert_array_equal(raw['initial'],case['initial']);np.testing.assert_array_equal(raw['times'],SHORT);np.testing.assert_array_equal(state[0,:3500],case['initial'].ravel());assert not np.any(state[0,3500:]);pop=state[:,:3500].reshape(3,7,500)
        error=float(np.max(abs(pop-ref)/np.maximum(1,case['initial'].max(axis=1))[None,:,None]));assert error<=2e-6;v=accounting(state,case,ga,gb);assert v['passed'];checks.append(dict(file=p.relative_to(ROOT).as_posix(),sha256=digest(p),**v));comparisons.append(dict(job=job['id'],error=error))
        if job['arm']==2:assert not np.any(state[:,11500:13000])
        if not job['start']:np.testing.assert_array_equal(pop,np.broadcast_to(pop[:,:,:1],pop.shape))
        if job['arm']==job['origin']:
            with np.load(donor_path(job)) as raw:target=raw['state'][41,:1500].reshape(3,500)
            err=float(np.max(abs(pop[-1,:3]-target)/S0[:,None]));assert err<=2e-6;shams.append(dict(job=job['id'],active_error=err))
    pipeline(deadline);forecasts={id:m['non_io_seconds']*60*(3 if m['job']['start'] else 2)+m['io_seconds']*200*(3 if m['job']['start'] else 2)+(15 if m['job']['start'] else 5) for id,m in timings.items()};forecast=sum(forecasts.values());elapsed=time.perf_counter()-start;assert elapsed<900 and shutil.disk_usage(ROOT).free>=2*1024**3
    assert donor_manifest()==plan['donor_hashes'];q=dict(passed=True,status='DEVELOPMENT PASS; concrete scientific criteria/budget/retention approval pending',source_hashes=hashes(),donor_hashes=donor_manifest(),formal=formal_result,raw_checks=checks,independent_comparisons=comparisons,short_sham_comparisons=shams,elapsed_seconds=elapsed,forecast_seconds_per_job=forecasts,forecast_worker_seconds=forecast,forecast_core_hours=forecast/3600,eight_worker_wall_minutes=forecast/480,production_within_reserve=forecast<=2700,candidate_cap_core_hours=1,scope='24short post.25writer branches and24separate Radau references; no full post15scientific branch')
    atomic_json(BASE/'preparation_results.json',q);atomic_json(BASE/'attempt.json',dict(status='COMPLETE DEVELOPMENT',elapsed_seconds=elapsed));print(json.dumps({key:q[key] for key in ('passed','elapsed_seconds','forecast_core_hours','eight_worker_wall_minutes','production_within_reserve')},indent=2))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--plan',action='store_true');p.add_argument('--execute',action='store_true');a=p.parse_args()
    if a.plan:plan()
    elif a.execute:
        try:execute()
        except BaseException:atomic_json(BASE/'failure.json',dict(status='FAILED; preserve/version, no automatic retry',error=traceback.format_exc()));raise
    else:p.error('Choose plan/execute')

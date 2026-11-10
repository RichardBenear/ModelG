"""MG118 separate fixed tolerance diagnostic; original experiment is immutable."""
from prepare_mg12h3_population import *
from scipy.integrate import BDF,Radau
from mg11g2_pattern_v2_oracle import population_rhs,population_pattern,laplacians
from datetime import datetime,timezone
BASE=ROOT/'results/MG12H3/tolerance_check'
OLD=ROOT/'results/MG12H3/population';JOB=next(j for j in jobs() if j['id']=='ISOLATED_C1_s1_v0')
METHODS=(dict(id='bdf_tight',solver='BDF',rtol=1e-11,atol=1e-13,max_step=.05),dict(id='bdf_tighter',solver='BDF',rtol=1e-12,atol=1e-14,max_step=.025),dict(id='radau_reference',solver='Radau',rtol=1e-12,atol=1e-14,max_step=.025))
FILES=('src/run_mg12h3_tolerance.py','src/check_mg12h3_tolerance.py','spec/MG12H3_tolerance_check.md','src/mg12h3_population.py','src/prepare_mg12h3_population.py','src/mg11g2_pattern_v2_oracle.py','src/mg11g2_pattern_v2.py','src/mg11g3.py','src/mg11g2.py','src/mg11n1.py')
def diagnostic_hashes():
    names=set(FILES);names.add(factory(JOB)[5].relative_to(ROOT).as_posix());names.update(p.relative_to(ROOT).as_posix() for p in (OLD/'production').glob('*.npz'));names.update('results/MG12H3/population/'+name for name in ('registration.json','failure.json','ledger.json','stopped_evaluation.json'))
    return {name:hashlib.sha256((ROOT/name).read_bytes() if name.endswith('.npz') else (ROOT/name).read_text(encoding='utf-8').replace('\r\n','\n').encode()).hexdigest() for name in sorted(names)}
def register():
    if (BASE/'registration.json').exists():raise ValueError('Already registered')
    old=read(OLD/'registration.json')
    from prepare_mg12h3_population import hashes as original_hashes
    assert original_hashes()==old['source_hashes']
    atomic_json(BASE/'registration.json',dict(stage='MG118',authorization='Go ahead with the tighter tolerance check',registered=datetime.now(timezone.utc).isoformat(),job=JOB,methods=METHODS,times=(np.arange(101)/4).tolist(),source_original_hashes=diagnostic_hashes(),elapsed_cap=300,workers=1,threads_per_library=1,margin=1e-5,balance_limit=1e-10,public_retention=True,original_verdict='failed/incomplete; unchanged'))
def save_diagnostic(path,method,times,state):
    tmp=path.with_suffix('.npz.tmp')
    with tmp.open('wb') as f:np.savez_compressed(f,state=state,times=times,job_json=json.dumps(JOB),method_json=json.dumps(method))
    os.replace(tmp,path)
def run():
    r=read(BASE/'registration.json')
    if r['source_original_hashes']!=diagnostic_hashes() or (BASE/'run_state.json').exists():raise ValueError('Changed source/original or existing attempt')
    subprocess.run(['git','ls-files','--error-unmatch','results/MG12H3/tolerance_check/registration.json'],cwd=ROOT,check=True,stdout=subprocess.DEVNULL);start=time.perf_counter();deadline=start+300;started=datetime.now(timezone.utc).isoformat();atomic_json(BASE/'run_state.json',dict(status='RUNNING',started=started));case,ga,gb,k,mode,reference,parent=factory(JOB);n=ga['n'];mats=matrices(case,ga,gb);oracle_mats=laplacians(ga);times=np.array(r['times']);records=[]
    with np.load(OLD/'production'/(JOB['id']+'.npz')) as z:
        for frame in (0,50,100):
            y=z['state'][frame];new=rhs(y,case,mats,k,mode)[:7*n];other=population_rhs(0,y[:7*n],case,ga,k,mode,oracle_mats);np.testing.assert_allclose(new,other,atol=1e-9,rtol=1e-12)
    for method in METHODS:
        clock=time.perf_counter();reference_method=method['solver']=='Radau';initial=case['initial'].ravel() if reference_method else np.r_[case['initial'].ravel(),np.zeros(2*n)]
        def fun(t,y):
            if time.perf_counter()>deadline:raise TimeoutError('Tolerance diagnostic cap')
            return population_rhs(t,y,case,ga,k,mode,oracle_mats) if reference_method else rhs(y,case,mats,k,mode)
        solver=(Radau if reference_method else BDF)(fun,0.,initial,25.,jac_sparsity=population_pattern(ga) if reference_method else sparsity(n,mats),rtol=method['rtol'],atol=method['atol'],max_step=method['max_step']);saved=[initial.copy()];index=1;path=BASE/(method['id']+'.npz')
        def output():
            values=np.array(saved)
            if reference_method:values=np.c_[values,np.zeros((len(values),2*n))]
            save_diagnostic(path,method,times[:len(values)],values);return values
        output()
        while solver.status=='running':
            solver.step()
            if solver.status=='failed':raise RuntimeError('Tolerance solver failed')
            while index<len(times) and times[index]<=solver.t:saved.append(solver.dense_output()(times[index]));index+=1;output()
        if index!=101:raise RuntimeError('Incomplete tolerance record')
        state=output();checked=check(JOB,state,times);record=dict(method=method,status='COMPLETE',elapsed_seconds=time.perf_counter()-clock,check=checked,sha256=hashlib.sha256(path.read_bytes()).hexdigest());records.append(record);atomic_json(BASE/(method['id']+'.json'),record);atomic_json(BASE/'ledger.json',dict(records=records,elapsed_seconds=time.perf_counter()-start,elapsed_cap=300))
        assert checked['balance']<=1e-10 and checked['minimum']>=-1e-9,'Diagnostic domain/accounting failure'
        # Accuracy results are scored after all three predeclared diagnostic integrations.
        # A solver/domain failure still stops immediately without retries.
    from check_mg12h3_tolerance import evaluate
    result=evaluate();elapsed=time.perf_counter()-start;assert elapsed<300 and diagnostic_hashes()==r['source_original_hashes'];atomic_json(BASE/'run_state.json',dict(status='DIAGNOSTIC RESULTS/AUDIT READY; original study unchanged',started=started,finished=datetime.now(timezone.utc).isoformat(),elapsed_seconds=elapsed,elapsed_cap=300,supported=result['supported']))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--register',action='store_true');p.add_argument('--run',action='store_true');a=p.parse_args()
    if a.register:register()
    elif a.run:
        try:run()
        except BaseException:atomic_json(BASE/'failure.json',dict(error=traceback.format_exc()));atomic_json(BASE/'run_state.json',dict(status='DIAGNOSTIC STOPPED; data preserved; no retry'));raise
    else:p.error('Choose action')

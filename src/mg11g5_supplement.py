"""MG100b separate drift reference, adversarial checks and secondary report.

Uses already qualified short trajectories; new Radau integrations only post.01.
No new random streams and no change to primary scientific criteria.
"""
from mg11g5 import *
import json,hashlib,argparse,subprocess
from scipy.integrate import solve_ivp
from mg11g3 import pattern as ledger_pattern
from mg11g2_pattern_v2_oracle import population_rhs,population_pattern,laplacians
from mg11g4_io import check_raw,digest,save
from mg10p1 import continuity
from mg10d1_pool import atomic_json
from prepare_mg11g5_v2 import hashes as primary_hashes,read
from evaluate_mg11g5 import evaluate
from check_mg11g5 import audit
BASE=ROOT/'results/MG11G5/preparation_v2';EXTRA=ROOT/'results/MG11G5/supplement'
def hashes():
    return dict(primary_hashes(),**{name:hashlib.sha256((ROOT/name).read_text(encoding='utf-8').replace('\r\n','\n').encode()).hexdigest() for name in ('src/mg11g5_supplement.py','src/run_mg11g5_ready.py')})
def plan():
    if (EXTRA/'plan.json').exists():raise ValueError('Already planned')
    q=read(BASE/'results.json');assert q['passed'] and primary_hashes()==q['source_hashes']
    atomic_json(EXTRA/'plan.json',dict(source_hashes=hashes(),scope='16separate Radau short drift references, invalid input/ledger/primary/metadata tests and secondary reporting; no new RNG',primary_preparation_sha256=digest(BASE/'results.json'),ceiling_elapsed_seconds=300,primary_population_error_bound=2e-6,shared='same initial rounded donors and structural Jacobian; separate population equations and solver'))
def qualify():
    p=read(EXTRA/'plan.json');assert p['source_hashes']==hashes() and not (EXTRA/'attempt.json').exists();subprocess.run(['git','ls-files','--error-unmatch','results/MG11G5/supplement/plan.json'],cwd=ROOT,check=True,stdout=subprocess.DEVNULL);start=time.perf_counter();deadline=start+300;atomic_json(EXTRA/'attempt.json',dict(status='RUNNING'));references=[]
    for j in [j for j in jobs() if j['kind']=='deterministic']:
        integer,ga,gb,k=factory(j);case=concentration_case(integer);LA,active=laplacians(ga);LB,_=laplacians(gb)
        def fun(t,y):
            if time.perf_counter()>deadline:raise TimeoutError('Supplement deadline')
            z=y.reshape(7,500);out=population_rhs(t,y,case,ga,[k[0],0.],'FINITE',(LA,active)).reshape(7,500);out[4]+=k[1]*(LB@z[4])
            if j['finite']:
                # Separate explicit finite-count formula (no primary drift helper).
                old=z[1]**2*z[2];new=z[1]*np.maximum(z[1]-1/1000,0)*z[2];out[1]+=new-old;out[2]-=new-old
            return out.ravel()
        with np.load(BASE/(j['id']+'.npz')) as old:times=old['times'];target=old['state'][:,:3500].reshape(-1,7,500)
        sparsity=ledger_pattern(ga,gb)[:3500,:3500];sol=solve_ivp(fun,(0,.01),case['initial'].ravel(),t_eval=times,method='Radau',jac_sparsity=sparsity,rtol=1e-11,atol=1e-12,max_step=.0025);assert sol.success;actual=sol.y.T.reshape(-1,7,500);scale=np.maximum(1,case['initial'].max(axis=1));error=float(np.max(abs(actual-target)/scale[None,:,None]));assert error<=2e-6
        with (EXTRA/(j['id']+'_reference.npz')).open('wb') as f:np.savez_compressed(f,population=actual,times=times,job_json=json.dumps(j))
        references.append(dict(job=j,max_scaled_population_error=error))
    # Same physical parent across all stochastic settings/blocks/destinations.
    for origin in range(2):
        for s in range(2):
            group=[j for j in jobs() if j['kind']=='stochastic' and j['origin']==origin and j['start']==s];case,*_=factory(group[0])
            for j in group[1:]:
                other,*_=factory(j)
                for key in ('initial','a0','b0'):np.testing.assert_array_equal(case[key],other[key])
    bad_rejected=0;case,ga,gb,k=factory(contexts()[0])
    calls=[lambda:streamed(case,ga,gb,k,[0,.01],-1,TAUS[0],deadline),lambda:streamed(case,ga,gb,k,[0,0],1,TAUS[0],deadline),lambda:streamed(case,ga,gb,k,[0,np.nan],1,TAUS[0],deadline),lambda:streamed(case,ga,gb,k,[0,.01],1,-1,deadline)]
    for call in calls:
        try:call()
        except ValueError:bad_rejected+=1
        else:raise AssertionError('Invalid stream accepted')
    # Deliberately corrupted raw ledger and altered primary JSON must be rejected.
    j=next(j for j in jobs() if j['kind']=='stochastic');source=BASE/'synthetic'/(j['id']+'_complete_reference.npz');bad=EXTRA/'bad_ledger.npz'
    with np.load(source) as raw:
        payload={key:raw[key] for key in raw.files};payload['data']=payload['data'].copy();payload['data'][0,-1,0,0]+=1
    with bad.open('wb') as f:np.savez_compressed(f,**payload)
    try:check_raw(bad)
    except AssertionError:bad_rejected+=1
    else:raise AssertionError('Corrupted ledger accepted')
    altered=read(BASE/'synthetic_prefix.json');altered['settings'][1]['mean']+=1;atomic_json(EXTRA/'bad_primary.json',altered)
    try:audit(BASE/'synthetic',EXTRA/'bad_primary.json',EXTRA/'must_not_pass.json')
    except AssertionError:bad_rejected+=1
    else:raise AssertionError('Wrong primary accepted')
    assert hashes()==p['source_hashes'];elapsed=time.perf_counter()-start;assert elapsed<300;atomic_json(EXTRA/'results.json',dict(passed=True,source_hashes=hashes(),references=references,invalid_cases_rejected=bad_rejected,same_parents=True,elapsed_seconds=elapsed));atomic_json(EXTRA/'attempt.json',dict(status='COMPLETE'))
def report():
    directory=ROOT/'results/MG11G5/production';rows=[];controls={}
    for j in jobs():
        p=directory/(j['id']+'.npz')
        if not p.exists():continue
        with np.load(p) as z:
            times=z['times'];stochastic=j['kind']=='stochastic';pop=z['data'][0,:,:7]/OMEGA if stochastic else z['state'][:,:3500].reshape(-1,7,500)
            rx=z['data'][0,:,7:14]/OMEGA if stochastic else z['state'][:,3500:7000].reshape(-1,7,500);cont=continuity(pop[:,1],10+times);cont.pop('persistent_reference_signal');rows.append(dict(job=j,post_reference_continuity=cont,stocks_AB=pop[:,3:5].sum(axis=-1).tolist(),uptake_AB=rx[:,[0,3]].sum(axis=-1).tolist(),scope='post-switch sampled field matching; no20-unit persistence/lifetime claim'))
            if not stochastic and np.array_equal(times,OBS):controls[j['finite'],j['origin'],j['start'],j['arm']]=late_value(pop,times)
    contrasts={}
    for finite in range(2):
        if all((finite,c['origin'],c['start'],c['arm']) in controls for c in contexts()):contrasts[str(finite)]=sum((1 if c['start'] else -1)*(1 if c['arm'] else -1)*controls[finite,c['origin'],c['start'],c['arm']] for c in contexts())/2
    atomic_json(ROOT/'results/MG11G5/secondary.json',dict(conditional_control_interaction=contrasts,rows=rows,interpretation='control closures are not exact stochastic means; origin also fixes A topology'))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--plan',action='store_true');p.add_argument('--qualify',action='store_true');p.add_argument('--report',action='store_true');a=p.parse_args()
    if a.plan:plan()
    elif a.qualify:qualify()
    elif a.report:report()
    else:p.error('Choose action')

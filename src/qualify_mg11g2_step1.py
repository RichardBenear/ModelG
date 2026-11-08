"""MG087 prospective bounded deterministic Step1 qualification, no sampling."""
import os
for key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[key]='1'
from pathlib import Path
import argparse,time,json,hashlib,traceback,subprocess
import numpy as np
from scipy.integrate import solve_ivp
from mg11g2 import initial,unpack,rhs,reconstruct,solve,configuration,well_mixed_rhs,mixed_initial
from mg11g2_cases import fixture,TIMES
from mg10_reservoir import ode_rhs as legacy_rhs
from mg10d1_pool import atomic_json
ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/'results/MG11G2/step1'
FILES=('src/mg11g2.py','src/mg11g2_cases.py','src/qualify_mg11g2_step1.py','src/check_mg11g2_step1.py','spec/MG11G2_step1.md','src/mg11n1.py','src/mg10_reservoir.py','src/mg10s1.py','src/mg10d1_pool.py')

def hashes():return {name:hashlib.sha256((ROOT/name).read_text(encoding='utf-8').replace('\r\n','\n').encode()).hexdigest() for name in FILES}
def jobs():
    out=[dict(id='zero_network',fixture='legacy',capacity=1,kappa=[0.,0.],mode='FINITE'),dict(id='maintained',fixture='chemostat',capacity=1,kappa=[.7,.7],mode='CHEMOSTAT')]
    out += [dict(id=f'capacity_c{c}',fixture='capacity',capacity=c,kappa=[.7,.7],mode='FINITE') for c in (1,10,100)]
    out += [dict(id=f'mixing_k{k}',fixture='mixing',capacity=1,kappa=[float(k)]*2,mode='FINITE') for k in (1,10,100,1000)]
    return out

def plan():
    if (BASE/'plan.json').exists():raise ValueError('Already planned; preserve/version')
    atomic_json(BASE/'plan.json',dict(stage='MG087',scope='Authorized Step1 deterministic coupling qualification only; no particles/stochastic production',pi_authorization='yes, I like that approach. Go ahead and build the Step 1.',source_hashes=hashes(),jobs=jobs(),times=TIMES.tolist(),elapsed_cap_seconds=300,workers=1,threads_per_library=1,random_draws=0))

def formal():
    metrics={};case,g=fixture('legacy');y=initial(case,g);z=unpack(y,g)[0];z[3]*=.7;z[4]*=.6;n=g['n'];scale=np.maximum(1,z.max(axis=1))
    for name,k,mode in (('closed',[.7,.9],'FINITE'),('zero_network',[0,0],'FINITE')):
        dy=rhs(0,y,case,g,k,mode);dp,rx,hop,env,imports=unpack(dy,g);pred=reconstruct(dict(case,initial=np.zeros_like(case['initial'])),g,dy);err=float(np.max(abs(dp-pred)/scale[:,None]));assert err<=1e-11 and abs(dp.sum()-imports.sum())<=1e-10;metrics[name+'_derivative_ledger_error']=err
    # Legacy state uses depletion fractions, not resource concentrations.
    old=np.zeros((15,n));old[:3]=z[:3];old[3]=1-z[3]/(case['a0']/case['omega']);old[4]=1-z[4]/(case['b0']/case['omega']);d=legacy_rhs(0,old.ravel(),1,True).reshape(15,n);actual=unpack(rhs(0,y,case,g,[0,0]),g)[0];expected=np.zeros((7,n));expected[:3]=d[:3];expected[3]=-d[3]*case['a0']/case['omega'];expected[4]=-d[4]*case['b0']/case['omega'];expected[5]=d[8];expected[6]=d[10];err=float(np.max(abs(actual-expected)/scale[:,None]));assert err<=1e-11;metrics['legacy_nonlinear_drift_error']=err
    case,g=fixture('chemostat');y=initial(case,g);dy=rhs(0,y,case,g,[.7,.7],'CHEMOSTAT');dp,rx,hop,env,imports=unpack(dy,g);old=np.zeros((15,n));old[:3]=unpack(y,g)[0][:3];expected=legacy_rhs(0,old.ravel(),1,False).reshape(15,n)[:3];assert np.max(abs(dp[:3]-expected))<=1e-11 and np.all(dp[3:5]==0) and np.all(imports>=0);assert abs(dp.sum()-imports.sum())<=1e-10
    pred=reconstruct(dict(case,initial=np.zeros_like(case['initial'])),g,dy);assert np.max(abs(dp-pred))<=1e-11;metrics['chemostat_drift_and_imports']=True
    # Quasipositivity of all species at each declared empty boundary.
    boundary_min=0.
    for s in range(7):
        empty=initial(case,g);unpack(empty,g)[0][s]=0;drift=unpack(rhs(0,empty,case,g,[.7,.7]),g)[0];boundary_min=min(boundary_min,float(drift[s].min()));assert drift[s].min()>=-1e-12
    empty=initial(case,g);unpack(empty,g)[0][3:5]=0;dp,rx,*_=unpack(rhs(0,empty,case,g,[0,0]),g);assert not np.any(rx[[0,3]]) and np.all(dp[3:5]==0);metrics['boundary_minimum']=boundary_min;metrics['zero_resource_consumption']=True
    rejected=0
    def reject(call):
        nonlocal rejected
        try:call()
        except ValueError:rejected+=1;return
        raise AssertionError('Invalid input accepted')
    reject(lambda:configuration(case,g,[1,1],'UNKNOWN'));bad=dict(case,initial=case['initial'].copy());bad['initial'][3,0]-=1;reject(lambda:configuration(bad,g,[1,1],'CHEMOSTAT'));reject(lambda:configuration(case,g,[-1,1],'FINITE'));reject(lambda:solve(case,g,[0,0],[0,0]));reject(lambda:solve(case,g,[0,0],[0,float('nan')]))
    metrics['invalid_inputs_rejected']=rejected;return metrics

def execute():
    started=time.perf_counter();deadline=started+300;registration=json.loads((BASE/'plan.json').read_text(encoding='utf-8'))
    if registration['source_hashes']!=hashes() or registration['jobs']!=jobs():raise ValueError('Source/plan changed')
    subprocess.run(['git','ls-files','--error-unmatch','results/MG11G2/step1/plan.json'],cwd=ROOT,check=True,stdout=subprocess.DEVNULL)
    if (BASE/'attempt.json').exists():raise ValueError('Already attempted; no automatic retry')
    atomic_json(BASE/'attempt.json',dict(status='RUNNING',scope=registration['scope']));tests=formal();outputs={}
    for job in registration['jobs']:
        if time.perf_counter()>deadline:raise TimeoutError('Qualification cap')
        case,g=fixture(job['fixture'],job['capacity']);result=solve(case,g,job['kappa'],TIMES,job['mode'],deadline=deadline);outputs[job['id']]=result['state'];path=BASE/(job['id']+'.npz')
        with path.with_suffix('.npz.tmp').open('wb') as f:np.savez_compressed(f,times=result['times'],state=result['state'],initial=case['initial'],a0=case['a0'],b0=case['b0'],omega=case['omega'],diff=case['diff'],src=g['src'],dst=g['dst'],weight=g['weight'],kappa=job['kappa'])
        os.replace(path.with_suffix('.npz.tmp'),path);print('Saved '+job['id'],flush=True)
    case,g=fixture('mixing');scale=np.array([140/9,140/9,29/(140/9),40,80,20,20])
    def bounded_mixed(t,y):
        if time.perf_counter()>deadline:raise TimeoutError('Qualification cap')
        return well_mixed_rhs(t,y,case)
    sol=solve_ivp(bounded_mixed,(0,1),mixed_initial(case).ravel(),t_eval=TIMES,method='BDF',rtol=1e-10,atol=1e-12,max_step=.025);assert sol.success;mixed=sol.y.T.reshape(4,7,g['n']);np.savez_compressed(BASE/'mixed_reference_primary.npz',times=TIMES,population=mixed)
    mixing=[float(np.max(abs(outputs[f'mixing_k{k}'][1:,:7*g['n']].reshape(3,7,g['n'])-mixed[1:])/scale[None,:,None])) for k in (1,10,100,1000)]
    case,g=fixture('capacity');control=outputs['maintained'][:,:3*g['n']].reshape(4,3,g['n']);capacity=[float(np.max(abs(outputs[f'capacity_c{c}'][:,:3*g['n']].reshape(4,3,g['n'])-control)/scale[None,:3,None])) for c in (1,10,100)]
    # Save numerical evidence even when a limiting criterion fails.
    atomic_json(BASE/'qualification.json',dict(passed=False,status='SEPARATE AUDIT PENDING',formal=tests,mixing_errors=mixing,capacity_errors=capacity))
    assert np.all(np.diff(mixing)<0) and mixing[-1]<=.005;assert np.all(np.diff(capacity)<0) and capacity[-1]<=.001
    from check_mg11g2_step1 import main as audit
    checks=audit(BASE,deadline=deadline);elapsed=time.perf_counter()-started;assert elapsed<300
    atomic_json(BASE/'qualification.json',dict(passed=True,status='QUALIFIED WITHIN G2 STEP1 SMALL DETERMINISTIC SCOPE',formal=tests,mixing_errors=mixing,capacity_errors=capacity,elapsed_seconds=elapsed,source_hashes=hashes(),raw_jobs=len(checks)))
    atomic_json(BASE/'attempt.json',dict(status='COMPLETE; no particle/stochastic study run',elapsed_seconds=elapsed))
    text=f'''# MG11-G2 Step1 coupling qualification — MG087

Step1 qualified the coupling between the explicit environmental network and the full nonlinear Model G continuum reactions on the declared small systems. Resource use and transport can now be reconstructed from explicit ledgers, maintained supplies recover the old active equations, and the tested finite-capacity and rapid-mixing limits approach their appropriate references. This prepares a controlled study of how environmental structure affects Model G patterns; it does not yet establish particle persistence, long-horizon or stochastic accuracy, a physical vacuum, or entanglement.

Status: **QUALIFIED WITHIN G2 STEP1 SMALL DETERMINISTIC SCOPE**. Build/qualification authorized; no pattern-response or stochastic production run. All frozen source/spec/plan published before integration; no scientific random draws.

Nine nonlinear coupling trajectories, all population/transfer/reaction/import ledgers retained. Separate population-only Radau references plus explicit incidence reconstruction pass for all nine. Maximum normalized reference error {max(v['reference_error'] for v in checks):.9g}; ledger error {max(v['ledger_error'] for v in checks):.9g}; balance error {max(v['balance_error'] for v in checks):.9g}; minimum recorded component {min(v['minimum'] for v in checks):.9g}. No clipping. Formal zero-stock, boundary, zero-network legacy, maintained drift/import, and five invalid-input tests pass.

Rapid-mixing errors at kappa1/10/100/1000: {mixing}; strictly decreasing, final<=.005. Initial mixing layer excluded from resource-limit comparison (positive checkpoints .25,.5,1 only). Reference is a SHARED FINITE bath that still depletes, not a maintained reservoir. Nonlinear result is restricted to this fixture/horizon, not a general convergence theorem.

Capacity-recovery active errors at C1/10/100: {capacity}; strictly decreasing, final<=.001. Capacity changes reference stocks, not Omega. This is a short-time nonuniform-field test, not recovery of a persistent pattern at arbitrary horizon.

Qualification elapsed {elapsed:.6f}seconds<300, one local process/library threads1. Same-author separate equations/ledger/Radau path shares constants/NumPy/SciPy; no external replication. All arrays, configs, plan, hashes and audit retained here. Historical kernels and studies unchanged. Next scientific design/budget and PI approval remain separate.
'''
    (BASE/'MG11G2_step1_results.md').write_text(text,encoding='utf-8');print(json.dumps(dict(passed=True,elapsed_seconds=elapsed,mixing_errors=mixing,capacity_errors=capacity)),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--plan',action='store_true');p.add_argument('--execute',action='store_true');a=p.parse_args()
    if a.plan:plan()
    elif a.execute:
        try:execute()
        except BaseException:atomic_json(BASE/'failure.json',dict(status='FAILED; preserve/version before correction, no automatic retry',error=traceback.format_exc()));raise
    else:p.error('Choose plan or execute')

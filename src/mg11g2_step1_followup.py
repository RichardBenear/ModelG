"""MG087a disclosed audit of failed mixing range and separate fixed-tolerance extension.

Original source/data/criteria/result remain frozen. No re-run of its trajectories.
"""
import os
for key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[key]='1'
from pathlib import Path
import argparse,time,json,hashlib,traceback
import numpy as np
from mg11g2 import solve
from mg11g2_cases import fixture,TIMES
from check_mg11g2_step1 import integrate,incidence
from qualify_mg11g2_step1 import FILES,hashes
from mg10d1_pool import atomic_json
ROOT=Path(__file__).resolve().parents[1];OLD=ROOT/'results/MG11G2/step1';BASE=ROOT/'results/MG11G2/step1_mixing_extension'
SCALE=np.array([140/9,140/9,29/(140/9),40,80,20,20])
def read(p):return json.loads(Path(p).read_text(encoding='utf-8'))
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def source_hashes():
    out=hashes()
    for name in ('src/mg11g2_step1_followup.py','spec/MG11G2_step1_mixing_extension.md'):out[name]=hashlib.sha256((ROOT/name).read_text(encoding='utf-8').replace('\r\n','\n').encode()).hexdigest()
    return out
def plan():
    if (BASE/'plan.json').exists():raise ValueError('Already planned')
    original=read(OLD/'plan.json');assert hashes()==original['source_hashes'];assert not read(OLD/'qualification.json')['passed']
    atomic_json(BASE/'plan.json',dict(stage='MG087a',scope='Step1 build qualification follow-up only; original .005 limit and failed range remain unchanged',source_hashes=source_hashes(),original_hashes={p.relative_to(ROOT).as_posix():sha(p) for p in OLD.iterdir() if p.is_file()},original_result='FAILED: kappa1000 error.00521369787 exceeds.005',additional_kappa=[10000.,10000.],elapsed_cap_seconds=300,workers=1,threads_per_library=1,random_draws=0))
def verify():
    p=read(BASE/'plan.json');assert source_hashes()==p['source_hashes']
    for name,value in p['original_hashes'].items():assert sha(ROOT/name)==value
    return p
def audit_one(case,g,states,mode,reference):
    n=g['n'];scale=np.maximum(1,(case['initial']/case['omega']).max(axis=1));pop=states[:,:7*n].reshape(4,7,n);err=float(np.max(abs(pop-reference)/scale[None,:,None]));ledger=0.;balance=0.
    for y in states:
        z,pred,imports=incidence(case,g,y);ledger=max(ledger,float(np.max(abs(z-pred)/scale[:,None])));balance=max(balance,abs(float(z.sum()-case['initial'].sum()/case['omega']-imports.sum()))/max(1,case['initial'].sum()/case['omega']))
        if mode=='FINITE':assert not np.any(imports)
    minimum=float(states.min());assert err<=2e-6 and ledger<=1e-8 and balance<=1e-8 and minimum>=-1e-10 and np.diff(states[:,7*n:],axis=0).min()>=-1e-10
    return dict(reference_error=err,ledger_error=ledger,balance_error=balance,minimum=minimum)
def audit_original():
    verify()
    if (BASE/'original_audit.json').exists():raise ValueError('Already audited')
    start=time.perf_counter();deadline=start+300;checks=[]
    for job in read(OLD/'plan.json')['jobs']:
        case,g=fixture(job['fixture'],job['capacity']);path=OLD/(job['id']+'.npz')
        with np.load(path) as raw:
            for name,value in dict(initial=case['initial'],a0=case['a0'],b0=case['b0'],omega=case['omega'],diff=case['diff'],src=g['src'],dst=g['dst'],weight=g['weight'],kappa=job['kappa'],times=TIMES).items():np.testing.assert_array_equal(raw[name],value)
            states=raw['state'];assert states.shape==(4,22*g['n']+2*len(g['src']))
        reference=integrate(case,g,job['kappa'],job['mode'],deadline);v=audit_one(case,g,states,job['mode'],reference);checks.append(dict(job=job['id'],sha256=sha(path),**v));np.savez_compressed(BASE/(job['id']+'_original_reference.npz'),times=TIMES,population=reference)
    case,g=fixture('mixing');mixed=integrate(case,g,[0,0],'FINITE',deadline,mixed=True);mix=[]
    np.savez_compressed(BASE/'mixed_independent_reference.npz',times=TIMES,population=mixed)
    for k in (1,10,100,1000):
        with np.load(OLD/f'mixing_k{k}.npz') as raw:pop=raw['state'][:,:7*g['n']].reshape(4,7,g['n'])
        mix.append(float(np.max(abs(pop[1:]-mixed[1:])/SCALE[None,:,None])))
    case,g=fixture('chemostat');maintained=integrate(case,g,[.7,.7],'CHEMOSTAT',deadline);capacity=[]
    for c in (1,10,100):
        with np.load(OLD/f'capacity_c{c}.npz') as raw:active=raw['state'][:,:3*g['n']].reshape(4,3,g['n'])
        capacity.append(float(np.max(abs(active-maintained[:,:3])/SCALE[None,:3,None])))
    q=read(OLD/'qualification.json');np.testing.assert_allclose(mix,q['mixing_errors'],rtol=1e-6,atol=2e-8);np.testing.assert_allclose(capacity,q['capacity_errors'],rtol=1e-6,atol=2e-8);assert np.all(np.diff(mix)<0) and mix[-1]>.005 and np.all(np.diff(capacity)<0) and capacity[-1]<=.001
    assert np.max(abs(mixed.sum(axis=(1,2))-mixed[0].sum()))<1e-8
    elapsed=time.perf_counter()-start;assert elapsed<300
    atomic_json(BASE/'original_audit.json',dict(audit_passed=True,original_qualification_passed=False,original_failed_criterion='Highest registered mixing rate does not reach fixed .005 tolerance',raw_checks=checks,mixing_errors=mix,capacity_errors=capacity,elapsed_seconds=elapsed,dependencies='Same author/shared fixture constants/NumPy/SciPy; separate population equations, directed incidence and Radau references. No external replication; original result remains failed.'))
    print(json.dumps(dict(original_audit_passed=True,original_mixing_criterion_passed=False,elapsed_seconds=elapsed)))
def extension():
    p=verify();assert read(BASE/'original_audit.json')['audit_passed']
    if (BASE/'extension_attempt.json').exists():raise ValueError('Already attempted; no retry')
    start=time.perf_counter();deadline=start+300;atomic_json(BASE/'extension_attempt.json',dict(status='RUNNING'));case,g=fixture('mixing');k=p['additional_kappa'];sol=solve(case,g,k,TIMES,deadline=deadline);states=sol['state']
    np.savez_compressed(BASE/'mixing_k10000.npz',times=TIMES,state=states,initial=case['initial'],a0=case['a0'],b0=case['b0'],omega=case['omega'],diff=case['diff'],src=g['src'],dst=g['dst'],weight=g['weight'],kappa=k)
    reference=integrate(case,g,k,'FINITE',deadline);v=audit_one(case,g,states,'FINITE',reference);np.savez_compressed(BASE/'mixing_k10000_reference.npz',times=TIMES,population=reference)
    with np.load(OLD/'mixed_reference_primary.npz') as z:mixed=z['population']
    with np.load(BASE/'mixed_independent_reference.npz') as z:independent_mixed=z['population']
    error=float(np.max(abs(states[1:,:7*g['n']].reshape(3,7,g['n'])-mixed[1:])/SCALE[None,:,None]));independent_error=float(np.max(abs(reference[1:]-independent_mixed[1:])/SCALE[None,:,None]));np.testing.assert_allclose(error,independent_error,rtol=1e-6,atol=2e-8)
    elapsed=time.perf_counter()-start;result=dict(passed=bool(error<=.005 and error<read(OLD/'qualification.json')['mixing_errors'][-1]),original_qualification_passed=False,kappa=k,mixing_error=error,independent_mixing_error=independent_error,unchanged_margin=.005,raw_check=v,sha256=sha(BASE/'mixing_k10000.npz'),elapsed_seconds=elapsed,scope='One additional prospectively frozen rate; original registered range FAILED and not rescored. Small deterministic fixture only, no pattern/stochastic production.')
    atomic_json(BASE/'extension_result.json',result);assert elapsed<300 and result['passed'];atomic_json(BASE/'extension_attempt.json',dict(status='COMPLETE',elapsed_seconds=elapsed));print(json.dumps(result,indent=2))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--plan',action='store_true');p.add_argument('--audit-original',action='store_true');p.add_argument('--extension',action='store_true');a=p.parse_args()
    if a.plan:plan()
    else:
        try:
            if a.audit_original:audit_original()
            elif a.extension:extension()
            else:p.error('Choose action')
        except BaseException:atomic_json(BASE/'followup_failure.json',dict(status='FAILED; retain/version, no automatic retry',error=traceback.format_exc()));raise

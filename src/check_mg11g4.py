"""MG096 separate raw/config/seed/CDF-family/ledger audit."""
import os
for key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[key]='1'
from pathlib import Path
import json,math
import numpy as np
from mg11g4_cases import fixture,R,TAUS
from mg11g4_io import check_raw,jobs
from mg10d1_pool import atomic_json
ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/'results/MG11G4'
def read(p):return json.loads(Path(p).read_text(encoding='utf-8'))
def equal(a,b):np.testing.assert_allclose(a,b,rtol=1e-10,atol=1e-11)
def main(base=BASE,synthetic=False):
    base=Path(base);r=read(base/'registration.json');e=read(base/'evaluation.json');assert r['jobs']==jobs();checks=[];groups={};full_jobs={}
    if not synthetic:
        from prepare_mg11g4 import hashes
        assert hashes()==r['source_hashes']
    for job in jobs():
        p=base/'production'/(job['id']+'.npz');case,ga,gb,k,t,m=fixture(job['cell'])
        if not p.exists():full_jobs[job['id']]=False;continue
        v=check_raw(p);checks.append(v);assert v['sha256']==e['records'][job['id']]['sha256']
        with np.load(p) as raw:
            for name,value in dict(initial=case['initial'],a0=case['a0'],b0=case['b0'],omega=case['omega'],diff=case['diff'],src=ga['src'],dst=ga['dst'],weight=ga['weight'],srcB=gb['src'],dstB=gb['dst'],weightB=gb['weight'],kappa=k,times=t,mask=m,engine=job['engine'],tau=TAUS[job['engine']],seed=job['seed'],requested_reps=job['reps']).items():np.testing.assert_array_equal(raw[name],value)
            data=raw['data'][:,:,:7].copy();groups.setdefault((job['cell'],job['engine']),[]).append(data)
        meta=read(p.with_suffix('.json')) if p.with_suffix('.json').exists() else {};ok=len(data)==job['reps'] and meta.get('status')=='COMPLETE';full_jobs[job['id']]=ok;assert e['records'][job['id']]['n']==len(data) and e['records'][job['id']]['complete']==ok
        if ok:assert meta['job']==job and meta['sha256']==v['sha256']
        warm=p.with_name(p.stem+'_warm.npz');checks.append(check_raw(warm))
        with np.load(warm) as raw:
            assert len(raw['data'])==1 and raw['requested_reps']==1 and raw['seed']==job['warmup_seed'] and raw['engine']==job['engine']
            for name,value in dict(initial=case['initial'],a0=case['a0'],b0=case['b0'],omega=case['omega'],diff=case['diff'],src=ga['src'],dst=ga['dst'],weight=ga['weight'],srcB=gb['src'],dstB=gb['dst'],weightB=gb['weight'],kappa=k,times=t,mask=m,tau=TAUS[job['engine']]).items():np.testing.assert_array_equal(raw[name],value)
    groups={key:np.concatenate(value) for key,value in groups.items()};full=all(full_jobs.values())
    for c in range(3):
        for m in range(5):
            key=(c,m);n=len(groups[key]) if key in groups else 0;ok=n==R and all(full_jobs[j['id']] for j in jobs() if (j['cell'],j['engine'])==key);assert e['groups'][str(key)]==dict(n=n,complete=ok)
    assert full==e['complete']==e['family_complete'] and e['jobs_complete']==sum(full_jobs.values()) and e['paths_available']==sum(len(x) for x in groups.values());eps=math.sqrt(math.log(1568000)/(2*R));seen=set();passed=[]
    assert len(e['formal'])==196 and len(e['coarse_medium_descriptive'])==196
    for formal,rows in ((True,e['formal']),(False,e['coarse_medium_descriptive'])):
        for row in rows:
            c,m,ref=row['cell'],row['engine'],row['reference_engine'];case,ga,gb,k,t,mask=fixture(c);index=int(np.flatnonzero(t==row['time'])[0]);assert index>0 and 0<=row['species']<7 and 0<=row['node']<ga['n'];assert (m,ref) in (((3,0),(0,4)) if formal else ((1,0),(2,0)));identity=(formal,c,m,ref,row['time'],row['species'],row['node']);assert identity not in seen;seen.add(identity);x=groups.get((c,m));y=groups.get((c,ref))
            assert row['n']==(len(x) if x is not None else 0) and row['n_reference']==(len(y) if y is not None else 0)
            if x is None or y is None or not len(x) or not len(y):assert row['distance'] is None;continue
            a=x[:,index,row['species'],row['node']];b=y[:,index,row['species'],row['node']];grid=np.union1d(a,b);d=float(np.max(abs(np.searchsorted(np.sort(a),grid,side='right')/len(a)-np.searchsorted(np.sort(b),grid,side='right')/len(b))));equal(d,row['distance'])
            if formal:
                if full:equal(eps,row['epsilon_each']);equal(d+2*eps,row['upper']);assert row['passed']==bool(d+2*eps<=.08);passed.append(d+2*eps<=.08)
                else:assert row['upper'] is None and row['passed'] is None and row['epsilon_each'] is None
    # Empty/missing comparisons must also suppress their entire formal family.
    if not full:assert all(x['passed'] is None and x['upper'] is None for x in e['formal'])
    qualified=full and r['formal_passed'] and len(passed)==196 and all(passed);expected='QUALIFIED WITHIN REGISTERED SMALL-SYSTEM STOCHASTIC SCOPE' if qualified else 'NOT QUALIFIED UNDER REGISTERED G4 CRITERIA' if full else 'INCOMPLETE G4 STOCHASTIC QUALIFICATION';assert e['status']==expected and e['scientific_pass']==qualified;charge=None
    if not synthetic:
        l=read(base/'production/ledger.json');charge=l['charged_worker_seconds'];assert l['ceiling']==7200 and l['prior_charge']==0 and charge<=7200;equal(charge,math.fsum(x['charged_worker_seconds'] for x in l['records']))
        if full:assert not l['stop_reason'] and len(l['records'])==120 and all(x['status']=='COMPLETE' and x['returncode']==0 for x in l['records'])
    atomic_json(base/'independent_check.json',dict(passed=True,status=expected,raw_checks=checks,formal_comparisons=196,charged_worker_seconds=charge,dependencies='Same author/shared law/constants/libraries and explicit two-graph integer reconstruction (no kernel update helper). Separate sorted/searchsorted CDFs and family/budget arithmetic. Python exact SSA independent of leap/Numba exact clocks. Not external replication or full temporal/large-lattice qualification.'));return checks
if __name__=='__main__':main()

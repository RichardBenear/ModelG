"""MG080 separate raw reconstruction, joint-generator action and sorted CDF audit."""
import os
for key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[key]='1'
from pathlib import Path
import json,hashlib,math
import numpy as np
from check_mg11n1_preparation import check_raw,categorical,probability
from mg11n1_cases import cell
from mg10d1_pool import atomic_json
ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/'results/MG11N1'
def cdf(sample,grid):return np.searchsorted(np.sort(sample),grid,side='right')/len(sample)
def equal(a,b):np.testing.assert_allclose(a,b,rtol=1e-10,atol=1e-11)
def main(base=BASE,synthetic=False):
    base=Path(base);r=json.loads((base/'registration.json').read_text());e=json.loads((base/'evaluation.json').read_text());checks=[];groups={};complete={}
    from prepare_mg11n1 import candidate_jobs
    assert r['jobs']==candidate_jobs() and len(e['analytic'])==6 and len(e['two_sample'])==28
    if not synthetic:
        for name,sha in r['source_hashes'].items():
            p=ROOT/name;content=p.read_bytes() if name.endswith('.npz') else p.read_text(encoding='utf-8').replace('\r\n','\n').encode();assert hashlib.sha256(content).hexdigest()==sha
        for name,sha in r['formal_evidence'].items():assert hashlib.sha256((BASE/'preparation'/name).read_bytes()).hexdigest()==sha
    for job in r['jobs']:
        path=base/'production'/(job['id']+'.npz');key=(job['cell'],job['engine']);case,g,k,t,m=cell(job['cell'])
        if not path.exists():complete[job['id']]=False;continue
        v=check_raw(path);checks.append(v);assert v['sha256']==e['records'][job['id']]['sha256']
        with np.load(path) as raw:
            for name,expected in dict(initial=case['initial'],a0=case['a0'],b0=case['b0'],omega=case['omega'],diff=case['diff'],src=g['src'],dst=g['dst'],weight=g['weight'],kappa=k,times=t,mask=m).items():np.testing.assert_array_equal(raw[name],expected)
            assert raw['seed']==job['seed'] and raw['requested_reps']==job['reps'] and raw['reps_completed']==len(raw['data']) and len(raw['data'])<=job['reps']
            data=raw['data'].copy();groups.setdefault(key,[]).append(data)
        meta=path.with_suffix('.json');ok=len(data)==job['reps'] and meta.exists() and json.loads(meta.read_text()).get('status')=='COMPLETE';complete[job['id']]=ok;assert ok==e['records'][job['id']]['complete']
        if ok:
            actual_meta=json.loads(meta.read_text());assert actual_meta['sha256']==v['sha256']
            if not synthetic:assert actual_meta['job']==job
        warm=path.with_name(path.stem+'_warm.npz');checks.append(check_raw(warm))
        with np.load(warm) as raw:
            assert len(raw['data'])==1 and raw['seed']==job['warmup_seed'] and raw['requested_reps']==1 and raw['reps_completed']==1
            for name,expected in dict(initial=case['initial'],a0=case['a0'],b0=case['b0'],omega=case['omega'],diff=case['diff'],src=g['src'],dst=g['dst'],weight=g['weight'],kappa=k,times=t,mask=m).items():np.testing.assert_array_equal(raw[name],expected)
    groups={key:np.concatenate(values) for key,values in groups.items()};counts={key:len(groups[key]) if key in groups else 0 for key in ((0,0),(1,0),(2,0),(3,0),(3,1))}
    group_full={key:all(complete[j['id']] for j in r['jobs'] if (j['cell'],j['engine'])==key) for key in counts}
    for key in counts:assert e['groups'][str(key)]==dict(n=counts[key],complete=group_full[key])
    analytic_full=all(group_full[c,0] for c in range(3));two_full=group_full[3,0] and group_full[3,1];assert analytic_full==e['analytic_family_complete'] and two_full==e['two_sample_family_complete']
    passed=[];epsilon=math.sqrt(math.log(24000)/(2*32768))
    for row in e['analytic']:
        c=row['cell'];data=groups.get((c,0));kind=('transport','A_source','B_contact')[c];case,g,k,t,m=cell(c);index=int(np.flatnonzero(t==row['time'])[0])
        if data is None or not len(data):assert row['distance'] is None and row['passed'] is None;continue
        labels=categorical(data,kind);raw=dict(initial=case['initial'],src=g['src'],dst=g['dst'],weight=g['weight'],kappa=k,omega=case['omega'],a0=case['a0'],b0=case['b0'])
        p=probability(raw,kind,row['time']);equal(p,row['reference']);distance=float(abs(cdf(labels[:,index],np.arange(len(p)))-np.cumsum(p)).max());equal(distance,row['distance'])
        if analytic_full:upper=distance+epsilon;equal(upper,row['upper']);assert row['passed']==bool(upper<=.03);passed.append(bool(upper<=.03))
        else:assert row['upper'] is None and row['passed'] is None
    eps=math.sqrt(math.log(224000)/(2*16384));x=groups.get((3,0));y=groups.get((3,1))
    for row in e['two_sample']:
        if x is None or y is None or not len(x) or not len(y):assert row['distance'] is None and row['passed'] is None;continue
        index=1 if row['time']==.05 else 2;s=row['species'];n=row['node'];a=x[:,index,s,n];b=y[:,index,s,n];grid=np.union1d(a,b);distance=float(abs(cdf(a,grid)-cdf(b,grid)).max());equal(distance,row['distance'])
        if two_full:upper=distance+2*eps;equal(upper,row['upper']);assert row['passed']==bool(upper<=.08);passed.append(bool(upper<=.08))
        else:assert row['upper'] is None and row['passed'] is None
    allcomplete=all(complete.values());expected='QUALIFIED WITHIN REGISTERED N1 SCOPE' if allcomplete and r['formal_passed'] and all(passed) else 'NOT QUALIFIED UNDER REGISTERED N1 CRITERIA' if allcomplete else 'INCOMPLETE N1 QUALIFICATION'
    assert expected==e['status'] and sum(complete.values())==e['jobs_complete'] and sum(counts.values())==e['paths_available']
    charge=None
    if not synthetic:
        ledger=json.loads((base/'production/ledger.json').read_text());assert ledger['ceiling']==3600 and ledger['prior_charge']==0;charge=ledger['charged_worker_seconds'];equal(charge,math.fsum(v['charged_worker_seconds'] for v in ledger['records']))
        if allcomplete:assert not ledger['stop_reason'] and len(ledger['records'])==32 and all(v['returncode']==0 and v['status']=='COMPLETE' for v in ledger['records'])
    atomic_json(base/'independent_check.json',dict(passed=True,status=expected,raw_checks=checks,charged_worker_seconds=charge,analytical_comparisons=6,two_sample_comparisons=28,
      dependencies='Same author/shared laws/constants/NumPy/SciPy and prior explicit raw checker. Separately assembled joint generator with expm_multiply, strict categorical state encoding and sorted/searchsorted CDFs. Not external replication or full temporal/large-lattice qualification.'))
if __name__=='__main__':main()

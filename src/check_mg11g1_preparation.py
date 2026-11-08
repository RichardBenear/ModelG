"""MG083 separate matching, state, generator-action, interval and cost audit."""
import os
for key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[key]='1'
from pathlib import Path
import json,hashlib,math,time,shutil
import numpy as np
from scipy.optimize import brentq
from scipy.special import betainc
from check_mg11n1_preparation import check_raw,probability
from mg11g1 import cell,CONTACT
ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/'results/MG11G1/preparation'
def read(p):return json.loads(Path(p).read_text(encoding='utf-8'))
def labels(data,resource):
    p=data[:,:,:7];n=p.shape[-1];expected=np.zeros_like(p)
    if resource==0:
        assert np.all(p[:,:,3].sum(axis=-1)+p[:,:,0].sum(axis=-1)==1);expected[:,:,3]=p[:,:,3];expected[:,:,0]=p[:,:,0]
        result=np.where(p[:,:,3].sum(axis=-1)==1,p[:,:,3].argmax(axis=-1),n+p[:,:,0].argmax(axis=-1))
        np.testing.assert_array_equal(data[:,:,7],p[:,:,0]);success=result==n+CONTACT
    else:
        c=data[:,:,10,CONTACT];assert np.all((c==0)|(c==1)) and np.all(p[:,:,4].sum(axis=-1)+c==1);expected[:,:,4]=p[:,:,4];expected[:,:,1,CONTACT]=1-c;expected[:,:,2,CONTACT]=c;expected[:,:,5,CONTACT]=c
        result=np.where(c==1,n,p[:,:,4].argmax(axis=-1));assert np.array_equal(data[:,:,10].sum(axis=-1),c);success=result==n
    np.testing.assert_array_equal(p,expected);np.testing.assert_array_equal(data[:,:,7 if resource==0 else 10,CONTACT],success.astype(np.int64));return result
def independent_interval(k,n):
    if n==0:return None
    lo=0. if k==0 else brentq(lambda p:betainc(k,n-k+1,p)-.0125,0,1,xtol=1e-14)
    hi=1. if k==n else brentq(lambda p:betainc(k+1,n-k,p)-.9875,0,1,xtol=1e-14)
    return [lo,hi]
def raw_case(case,g,k):return dict(initial=case['initial'],a0=case['a0'],b0=case['b0'],omega=case['omega'],src=g['src'],dst=g['dst'],weight=g['weight'],kappa=np.asarray(k))
def main():
    started=time.perf_counter();plan=read(BASE/'preparation_plan.json')
    for name,sha in plan['source_hashes'].items():assert hashlib.sha256((ROOT/name).read_text(encoding='utf-8').replace('\r\n','\n').encode()).hexdigest()==sha,name
    formal=read(BASE/'qualification/formal_references.json');q=read(BASE/'qualification/qualification_results.json');assert q['passed'] and all(x['passed'] for x in q['tests']);checks=[]
    for item in q['raw_checks']:
        path=ROOT/item['file'];v=check_raw(path);assert v['sha256']==item['sha256'];checks.append(v)
        if path.name.startswith('null_'):
            with np.load(path) as z:assert not np.any(z['data'][:,:,7:14,CONTACT]) and not np.any(z['environment'])
    for record in formal['geometry']:
        W=np.array(record['adjacency']);assert np.array_equal(W,W.T) and np.all(W.sum(axis=1)==3) and np.count_nonzero(W)==18
        reached={0};front={0};distance=None
        for depth in range(6):
            if CONTACT in front:distance=depth;break
            nextfront={int(j) for i in front for j in np.flatnonzero(W[i])}-reached;reached|=nextfront;front=nextfront
        assert distance==1;triangles=sum(W[i,j]*W[j,k]*W[k,i] for i in range(6) for j in range(i+1,6) for k in range(j+1,6));assert triangles==record['triangles']
    with np.load(BASE/'qualification/reference_inputs.npz') as stored:
        perm=stored['permutation']
        for r in range(2):
            for network in range(2):
                case,g,k,t,m=cell(r,network);kind='A_source' if r==0 else 'B_contact';p=probability(raw_case(case,g,k),kind,4.);key=f'r{r}_g{network}';np.testing.assert_allclose(p,stored[key+'_p'],atol=1e-12)
                mapping=np.r_[perm,6+perm] if r==0 else np.r_[perm,6];np.testing.assert_allclose(stored[key+'_Q_permuted'][np.ix_(mapping,mapping)],stored[key+'_Q'],atol=1e-15);np.testing.assert_allclose(stored[key+'_p_permuted'][mapping],p,atol=1e-12)
                fast=probability(raw_case(case,g,[1000.,1000.]),kind,4.);entry=next(x for x in formal['cases'] if x['resource']==r and x['network']==network);np.testing.assert_allclose(fast,entry['fast_probability'],atol=1e-10)
                assert abs((fast[6+CONTACT] if r==0 else fast[-1])-entry['fast_target_expected'])<.001
                if r==0:
                    for target in t:assert abs(probability(raw_case(case,g,k),kind,float(target))[6:].sum()-(1-math.exp(-.5*target)))<1e-12
    for x in read(BASE/'qualification/interval_cases.json'):np.testing.assert_allclose(independent_interval(x['k'],x['n']),x['interval'],atol=1e-11)
    design=read(BASE/'design_calculation.json');counts=design['primary_design_rounded_counts'];a=independent_interval(counts[0],65536);b=independent_interval(counts[1],65536);np.testing.assert_allclose([b[0]-a[1],b[1]-a[0]],design['primary_design_interval']['ci'],atol=1e-11)
    timings={};descriptive=[]
    for r in range(2):
        for network in range(2):
            meta=read(BASE/f'timing/timing_r{r}_g{network}.json');assert meta['status']=='COMPLETE';path=BASE/'timing'/meta['raw'];v=check_raw(path);assert v['sha256']==meta['sha256'];checks.append(v);checks.append(check_raw(path.with_name(path.stem+'_warm.npz')))
            case,g,k,t,m=cell(r,network)
            with np.load(path) as z:
                for name,expected in dict(initial=case['initial'],a0=case['a0'],b0=case['b0'],diff=case['diff'],src=g['src'],dst=g['dst'],weight=g['weight'],kappa=k,times=t,mask=m).items():np.testing.assert_array_equal(z[name],expected)
                assert z['seed']==meta['job']['seed'] and len(z['data'])==1024;tags=labels(z['data'],r)
                for i,target in enumerate(t[1:],1):
                    p=probability(raw_case(case,g,k),'A_source' if r==0 else 'B_contact',float(target));empirical=np.searchsorted(np.sort(tags[:,i]),np.arange(len(p)),side='right')/len(tags);descriptive.append(dict(resource=r,network=network,time=float(target),distance=float(abs(empirical-np.cumsum(p)).max()),scope='Development descriptive only, no production score'))
            timings[r,network]=meta
    forecasts={j['id']:timings[j['resource'],j['network']]['compute_seconds']*4*2+timings[j['resource'],j['network']]['io_seconds']*16*2+10 for j in plan['candidate_jobs']};total=math.fsum(forecasts.values())
    ledger=read(BASE/'timing/ledger.json');assert len(ledger['records'])==4 and not ledger['stop_reason'] and all(x['status']=='COMPLETE' and x['returncode']==0 for x in ledger['records'])
    allowance=ledger['charged_worker_seconds']+time.perf_counter()-started+5;assert allowance<900 and shutil.disk_usage(ROOT).free>=4*1024**3
    result=dict(passed=True,status='G1 DEVELOPMENT PASS; SCIENTIFIC ESTIMATION NOT REGISTERED',raw_checks=checks,formal_assertions=len(q['tests']),descriptive_cdf_distances=descriptive,
      forecast_seconds_per_job=forecasts,forecast_worker_seconds=total,forecast_core_hours=total/3600,eight_worker_wall_minutes=total/480,within_75_percent_reserve=total<=2700,candidate_cap_core_hours=1,
      raw_array_bound_bytes=262144*5*(20*6+2*18)*8,controlled_development_seconds=ledger['charged_worker_seconds'],total_preparation_allowance_seconds=allowance,
      dependencies='Same author/shared laws/NumPy/SciPy and frozen N1 exact sampler. Separate edge/counter/strict state reconstruction, BFS/triangle matching, joint-generator exponential action and beta-CDF inversion. No external replication or particle/native-geometry/quantum claim.')
    (BASE/'independent_check.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:v for k,v in result.items() if k not in ('raw_checks','descriptive_cdf_distances','forecast_seconds_per_job')},indent=2))
if __name__=='__main__':main()

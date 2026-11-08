"""MG084 separate raw accounting, strict labels, action generator and beta-CDF CI."""
import os
for key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[key]='1'
from pathlib import Path
import json,hashlib,math
import numpy as np
from mg11g1 import cell
from check_mg11g1_preparation import labels,raw_case,independent_interval
from check_mg11n1_preparation import check_raw,probability
from prepare_mg11g1 import candidate_jobs
from mg10d1_pool import atomic_json
ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/'results/MG11G1'
def read(p):return json.loads(Path(p).read_text(encoding='utf-8'))
def equal(a,b):np.testing.assert_allclose(a,b,rtol=1e-10,atol=1e-11)
def main(base=BASE,synthetic=False):
    base=Path(base);r=read(base/'registration.json');e=read(base/'evaluation.json');checks=[];groups={};complete={};counts={}
    assert r['jobs']==candidate_jobs() and len(e['analytic'])==16
    if not synthetic:
        for name,sha in r['source_hashes'].items():assert hashlib.sha256((ROOT/name).read_text(encoding='utf-8').replace('\r\n','\n').encode()).hexdigest()==sha
        for name,sha in r['formal_evidence'].items():assert hashlib.sha256((BASE/'preparation'/name).read_bytes()).hexdigest()==sha
    for job in r['jobs']:
        path=base/'production'/(job['id']+'.npz');key=(job['resource'],job['network']);case,g,k,t,m=cell(*key)
        if not path.exists():complete[job['id']]=False;continue
        v=check_raw(path);checks.append(v);assert v['sha256']==e['records'][job['id']]['sha256']
        with np.load(path) as raw:
            for name,expected in dict(**raw_case(case,g,k),diff=case['diff'],times=t,mask=m).items():np.testing.assert_array_equal(raw[name],expected)
            data=raw['data'];n=len(data);tag=labels(data,job['resource']);groups.setdefault(key,[]).append(tag)
            counter=data[:,:,7 if job['resource']==0 else 10,:];z=dict(n=n,successes=counter[:,:,3].sum(axis=0),contact=data[:,:,3+job['resource'],3].sum(axis=0),conversion=counter.sum(axis=(0,2)));counts.setdefault(key,[]).append(z)
            assert raw['seed']==job['seed'] and raw['requested_reps']==4096 and raw['reps_completed']==n and n<=4096
        meta=read(path.with_suffix('.json')) if path.with_suffix('.json').exists() else {};ok=n==4096 and meta.get('status')=='COMPLETE';complete[job['id']]=ok;assert ok==e['records'][job['id']]['complete']
        if ok:assert meta['sha256']==v['sha256'] and meta['job']==job
        warm=path.with_name(path.stem+'_warm.npz');checks.append(check_raw(warm))
        with np.load(warm) as raw:
            assert len(raw['data'])==1 and raw['seed']==job['warmup_seed'] and raw['requested_reps']==1 and raw['reps_completed']==1
            labels(raw['data'],job['resource'])
            for name,expected in dict(**raw_case(case,g,k),diff=case['diff'],times=t,mask=m).items():np.testing.assert_array_equal(raw[name],expected)
    tags={key:np.concatenate(v) for key,v in groups.items()};allfull=True;actual={}
    for key in ((0,0),(0,1),(1,0),(1,1)):
        n=sum(v['n'] for v in counts.get(key,[]));full=n==65536 and all(complete[j['id']] for j in r['jobs'] if (j['resource'],j['network'])==key);allfull &=full
        c=e['cells'][str(key)];assert c['n']==n and c['complete']==full
        actual[key]={name:sum((v[name] for v in counts.get(key,[])),np.zeros(5,dtype=np.int64)) for name in ('successes','contact','conversion')};equal(actual[key]['successes'],c['successes'])
        if n:
            for name,source in (('probability','successes'),('contact_mean','contact'),('total_conversion_mean','conversion')):equal(actual[key][source]/n,c[name])
            equal(actual[key]['contact']/n*(.5 if key[0]==0 else .8),c['expected_flux'])
    assert allfull==e['analytic_family_complete'];law=[];seen=set();eps=math.sqrt(math.log(32000)/(2*65536))
    for row in e['analytic']:
        key=(row['resource'],row['network']);assert (*key,row['time']) not in seen;seen.add((*key,row['time']));case,g,k,t,m=cell(*key);index=int(np.flatnonzero(t==row['time'])[0]);assert index>0
        p=probability(raw_case(case,g,k),'A_source' if key[0]==0 else 'B_contact',row['time']);equal(p,row['reference']);sample=tags.get(key)
        if sample is None or not len(sample):assert row['distance'] is None and row['passed'] is None;continue
        d=float(abs(np.searchsorted(np.sort(sample[:,index]),np.arange(len(p)),side='right')/len(sample)-np.cumsum(p)).max());equal(d,row['distance'])
        if allfull:equal(d+eps,row['upper']);assert row['passed']==bool(d+eps<=.03);law.append(d+eps<=.03)
        else:assert row['upper'] is None and row['epsilon'] is None and row['passed'] is None
    p=e['primary'];primary_full=all(e['cells'][str(key)]['complete'] for key in ((1,0),(1,1)));assert p['complete']==primary_full
    for name,key in (('successes_prism',(1,0)),('successes_bipartite',(1,1))):assert p[name]==actual[key]['successes'][-1]
    if primary_full:
        ka=int(actual[1,0]['successes'][-1]);kb=int(actual[1,1]['successes'][-1]);a=independent_interval(ka,65536);b=independent_interval(kb,65536);equal(p['ci'],[b[0]-a[1],b[1]-a[0]]);equal(p['mean'],(kb-ka)/65536)
    else:assert p['ci'] is None and p['mean'] is None
    full=allfull and all(complete.values());lawok=allfull and len(law)==16 and all(law);supported=bool(full and lawok and r['formal_passed'] and p['ci'][0]>0)
    assert e['complete']==full and e['law_supported']==lawok and e['connectivity_effect_supported']==supported
    assert e['status']==('COMPLETE' if full else 'INCOMPLETE')+' CONNECTIVITY-SUPPLY ESTIMATION' and e['jobs_complete']==sum(complete.values()) and e['paths_available']==sum(len(v) for v in tags.values())
    charge=None
    if not synthetic:
        ledger=read(base/'production/ledger.json');charge=ledger['charged_worker_seconds'];assert ledger['ceiling']==3600 and ledger['prior_charge']==0 and charge<=3600;equal(charge,math.fsum(v['charged_worker_seconds'] for v in ledger['records']))
        if full:assert not ledger['stop_reason'] and len(ledger['records'])==64 and all(v['returncode']==0 and v['status']=='COMPLETE' for v in ledger['records'])
    atomic_json(base/'independent_check.json',dict(passed=True,status=e['status'],raw_checks=checks,charged_worker_seconds=charge,analytical_comparisons=16,primary_interval_recomputed=primary_full,dependencies='Same author; shared constants/libraries and prior generic accounting checker. Strict state encoding, separately assembled generator action, sorted CDFs and beta-CDF root intervals. Not external replication.'))
if __name__=='__main__':main()

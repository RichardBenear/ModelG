"""MG084 chunkwise endpoint estimation; formal scores require full denominators."""
import os
for key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[key]='1'
from pathlib import Path
import json,hashlib,math
import numpy as np
from mg11g1 import cell,exact_probabilities,delivery,difference_interval,CONTACT
from mg10d1_pool import atomic_json
ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/'results/MG11G1'

def encode(data,resource):
    if resource==0:return np.where(data[:,:,3].sum(axis=2)>0,data[:,:,3].argmax(axis=2),6+data[:,:,0].argmax(axis=2))
    return np.where(data[:,:,10,CONTACT]>0,6,data[:,:,4].argmax(axis=2))

def summarize(r,records,groups):
    expected={str(key):dict(n=sum(len(v['labels']) for v in groups.get(key,[])),complete=all(records.get(j['id'],{}).get('complete',False) for j in r['jobs'] if (j['resource'],j['network'])==key)) for key in ((0,0),(0,1),(1,0),(1,1))}
    full=all(v['complete'] and v['n']==65536 for v in expected.values());eps=math.sqrt(math.log(32000)/(2*65536));rows=[];cells={}
    for resource in range(2):
        for network in range(2):
            key=(resource,network);items=groups.get(key,[]);n=expected[str(key)]['n'];case,g,k,t,m=cell(*key)
            labels=np.concatenate([v['labels'] for v in items]) if n else None
            sums={name:sum((v[name] for v in items),np.zeros(5,dtype=np.int64)) for name in ('successes','contact','conversion')}
            cells[str(key)]=dict(**expected[str(key)],successes=sums['successes'].tolist(),probability=(sums['successes']/n).tolist() if n else None,contact_mean=(sums['contact']/n).tolist() if n else None,expected_flux=(sums['contact']/n*(.5 if resource==0 else .8)).tolist() if n else None,total_conversion_mean=(sums['conversion']/n).tolist() if n else None)
            for index,target in enumerate(t[1:],1):
                p=exact_probabilities(case,g,k,resource,float(target));d=None if not n else float(abs(np.cumsum(np.bincount(labels[:,index],minlength=len(p)))/n-np.cumsum(p)).max())
                upper=d+eps if full else None
                rows.append(dict(resource=resource,network=network,time=float(target),n=n,reference=p.tolist(),distance=d,epsilon=eps if full else None,upper=upper,margin=.03,passed=bool(upper<=.03) if full else None))
    a=cells[str((1,0))];b=cells[str((1,1))];primary_full=all(v['complete'] and v['n']==65536 for v in (a,b))
    primary=difference_interval(a['successes'][-1],a['n'],b['successes'][-1],b['n']) if primary_full else dict(mean=None,ci=None)
    primary.update(complete=primary_full,successes_prism=a['successes'][-1],successes_bipartite=b['successes'][-1],n_prism=a['n'],n_bipartite=b['n'],analytic_prediction=.011529941280641665)
    law=full and all(v['passed'] for v in rows);complete=full and len(records)==64 and all(v['complete'] for v in records.values())
    supported=bool(complete and law and r['formal_passed'] and primary['ci'][0]>0)
    return dict(status=('COMPLETE' if complete else 'INCOMPLETE')+' CONNECTIVITY-SUPPLY ESTIMATION',complete=complete,jobs_complete=sum(v['complete'] for v in records.values()),paths_available=sum(v['n'] for v in records.values()),expected_paths=262144,records=records,cells=cells,analytic_family_complete=full,analytic=rows,law_supported=bool(law),primary=primary,connectivity_effect_supported=supported,formal_passed=r['formal_passed'])

def main(base=BASE):
    base=Path(base);r=json.loads((base/'registration.json').read_text(encoding='utf-8'));records={};groups={}
    for job in r['jobs']:
        path=base/'production'/(job['id']+'.npz')
        if not path.exists():continue
        with np.load(path) as raw:
            data=raw['data'];resource=job['resource'];labels=encode(data,resource)
            item=dict(labels=labels,successes=delivery(data,resource).sum(axis=0),contact=data[:,:,3+resource,CONTACT].sum(axis=0),conversion=data[:,:,7 if resource==0 else 10,:].sum(axis=(0,2)))
        meta=path.with_suffix('.json');ok=len(labels)==job['reps'] and meta.exists() and json.loads(meta.read_text(encoding='utf-8')).get('status')=='COMPLETE'
        records[job['id']]=dict(n=len(labels),complete=bool(ok),sha256=hashlib.sha256(path.read_bytes()).hexdigest());groups.setdefault((resource,job['network']),[]).append(item)
    result=summarize(r,records,groups);atomic_json(base/'evaluation.json',result)
    p=result['primary'];summary=('G1 supports a small connectivity-dependent change in resource delivery between two networks with matched degree, capacity, contact distance and hop rate. ' if result['connectivity_effect_supported'] else 'G1 has not established the registered connectivity-supply contrast with all required checks. ')
    summary+='This tests a mechanism by which an explicit environment can affect supply, reproducing a prediction of the chosen classical network model. It does not establish self-sustaining particles, a physical vacuum, energy or quantum entanglement. Separate raw-data audit must pass before acceptance.'
    text=f'# MG11-G1 results\n\n{summary}\n\nStatus: **{result["status"]}**; PI closure pending.\n\nCompleted {result["jobs_complete"]}/64 jobs; {result["paths_available"]}/262144 scientific paths. Warmups excluded.\n\nPrimary B capture at T=4, bipartite minus prism: {p["mean"]}; at least 95% joint-coverage interval {p["ci"]}. Predicted before sampling: {p["analytic_prediction"]}.\n\nEndpoint-law criteria: {sum(x["passed"] is True for x in result["analytic"])}/16 supported; 99.9% nominal family confidence for registered CDF bounds. Primary and law-family coverage are separate. Secondary A delivery, total conversion, occupancy, flux and curves are descriptive (evaluation.json), without simultaneous secondary confidence claims. A missing family withholds all its formal scores; failed laws withhold interpretation.\n'
    (base/'MG11G1_results.md').write_text(text,encoding='utf-8')
    return result
if __name__=='__main__':main()

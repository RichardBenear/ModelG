"""MG080 frozen categorical/two-sample endpoint qualification, no partial scoring."""
import os
for key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[key]='1'
from pathlib import Path
import json,hashlib,math
import numpy as np
from mg11n1_cases import cell
from mg11n1_oracle import linear_law,categories
from mg10d1_pool import atomic_json
ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/'results/MG11N1'
def main(base=BASE):
    base=Path(base);r=json.loads((base/'registration.json').read_text());groups={};records={}
    for job in r['jobs']:
        path=base/'production'/(job['id']+'.npz')
        if not path.exists():continue
        with np.load(path) as raw:data=raw['data'].copy()
        meta=path.with_suffix('.json');complete=len(data)==job['reps'] and meta.exists() and json.loads(meta.read_text()).get('status')=='COMPLETE'
        records[job['id']]=dict(n=len(data),complete=bool(complete),sha256=hashlib.sha256(path.read_bytes()).hexdigest())
        groups.setdefault((job['cell'],job['engine']),[]).append(data)
    pooled={key:np.concatenate(values) for key,values in groups.items()};expected={}
    for key in [(0,0),(1,0),(2,0),(3,0),(3,1)]:expected[str(key)]=dict(n=len(pooled[key]) if key in pooled else 0,complete=all(j['id'] in records and records[j['id']]['complete'] for j in r['jobs'] if (j['cell'],j['engine'])==key))
    analytic_full=all(expected[str((c,0))]['complete'] for c in range(3));two_full=all(expected[str((3,m))]['complete'] for m in range(2));analytic=[];two=[]
    epsilon=math.sqrt(math.log(2*6/.0005)/(2*32768))
    for c in range(3):
        case,g,k,t,m=cell(c);key=(c,0);data=pooled.get(key);labels=categories(data,case['name']) if data is not None and len(data) else None
        for index,target in enumerate(t[1:],1):
            p=linear_law(case,g,k,case['name'],float(target));distance=None;upper=None;passed=None
            if labels is not None:
                empirical=np.cumsum(np.bincount(labels[:,index],minlength=len(p)))/len(labels);distance=float(np.max(np.abs(empirical-np.cumsum(p))))
                if analytic_full:upper=distance+epsilon;passed=upper<=.03
            analytic.append(dict(cell=c,time=float(target),n=len(labels) if labels is not None else 0,distance=distance,epsilon=epsilon if analytic_full else None,upper=upper,margin=.03,passed=passed,reference=p.tolist()))
    eps=math.sqrt(math.log(4*28/.0005)/(2*16384));a=pooled.get((3,0));b=pooled.get((3,1));case,g,k,t,m=cell(3)
    for index,target in enumerate(t[1:],1):
        for species in range(7):
            for node in range(2):
                distance=None;upper=None;passed=None
                if a is not None and b is not None and len(a) and len(b):
                    x=a[:,index,species,node];y=b[:,index,species,node];size=max(int(x.max()),int(y.max()))+1
                    cx=np.cumsum(np.bincount(x,minlength=size))/len(x);cy=np.cumsum(np.bincount(y,minlength=size))/len(y);distance=float(abs(cx-cy).max())
                    if two_full:upper=distance+2*eps;passed=upper<=.08
                two.append(dict(time=float(target),species=species,node=node,distance=distance,epsilon_each=eps if two_full else None,upper=upper,margin=.08,passed=passed))
    complete=len(records)==32 and all(v['complete'] for v in records.values());allpass=complete and r['formal_passed'] and all(x['passed'] for x in analytic+two)
    status='QUALIFIED WITHIN REGISTERED N1 SCOPE' if allpass else 'NOT QUALIFIED UNDER REGISTERED N1 CRITERIA' if complete else 'INCOMPLETE N1 QUALIFICATION'
    result=dict(status=status,complete=complete,scientific_pass=bool(allpass),jobs_complete=sum(v['complete'] for v in records.values()),paths_available=sum(v['n'] for v in records.values()),expected_paths=131072,
       records=records,groups=expected,analytic_family_complete=analytic_full,two_sample_family_complete=two_full,analytic=analytic,two_sample=two,formal_passed=r['formal_passed'])
    atomic_json(base/'evaluation.json',result)
    paragraph=('N1 qualified the fixed environmental network on its registered small-system endpoint laws and declared formal/linear limits. Resources now move through an explicitly accounted network, providing tested machinery for later questions about geometry and supply. ' if allpass else 'N1 did not establish the complete registered environmental-network qualification; the missing or failed criteria are retained explicitly. ')
    paragraph+='This does not establish persistent particles on that network, full nonlinear fast-mixing recovery, emergent vacuum geometry, physical energy or quantum entanglement. Separate audit must pass before accepting the result.'
    text=f'# MG11-N1 results\n\n{paragraph}\n\nStatus: **{status}**; PI closure pending.\n\nCompleted{result["jobs_complete"]}/32jobs; available{result["paths_available"]}/131072scientific paths.\n\n'
    for title,items in (('Analytical endpoint family',analytic),('Full-tiny two-implementation family',two)):
        bounds=[x['upper'] for x in items if x['upper'] is not None];text+=f'- {title}: {sum(x["passed"] is True for x in items)}/{len(items)} criteria supported; largest upper={max(bounds) if bounds else None}; margin={items[0]["margin"]}.\n'
    text+='\nPrespecified DKW families have nominal>=99.9% combined confidence for the stated endpoint CDF-distance bounds under sampling assumptions. Warmups excluded; no incomplete-family score or full-joint-temporal/large-lattice claim. Raw/probability tables: evaluation.json; budget/progress/failures: production/ledger.json. Original MG10 studies unchanged.\n'
    (base/'MG11N1_results.md').write_text(text,encoding='utf-8')
if __name__=='__main__':main()

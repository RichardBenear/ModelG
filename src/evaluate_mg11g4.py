"""MG096 fixed finest-step and exact-oracle marginal CDF qualification."""
import os
for key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[key]='1'
from pathlib import Path
import json,math
import numpy as np
from mg11g4_cases import fixture,R
from mg11g4_io import digest
from mg10d1_pool import atomic_json
ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/'results/MG11G4'
def read(p):return json.loads(Path(p).read_text(encoding='utf-8'))
def cdf_distance(a,b):
    size=max(int(a.max()),int(b.max()))+1;return float(np.max(abs(np.cumsum(np.bincount(a,minlength=size))/len(a)-np.cumsum(np.bincount(b,minlength=size))/len(b))))
def main(base=BASE):
    base=Path(base);r=read(base/'registration.json');groups={};records={}
    for job in r['jobs']:
        p=base/'production'/(job['id']+'.npz')
        if not p.exists():continue
        with np.load(p) as raw:data=raw['data'][:,:,:7].copy()
        meta=p.with_suffix('.json');complete=len(data)==job['reps'] and meta.exists() and read(meta).get('status')=='COMPLETE';records[job['id']]=dict(n=len(data),complete=bool(complete),sha256=digest(p));groups.setdefault((job['cell'],job['engine']),[]).append(data)
    pooled={key:np.concatenate(values) for key,values in groups.items()};counts={str((c,m)):dict(n=len(pooled[c,m]) if (c,m) in pooled else 0,complete=all(records.get(j['id'],{}).get('complete',False) for j in r['jobs'] if (j['cell'],j['engine'])==(c,m))) for c in range(3) for m in range(5)}
    full=len(records)==120 and all(v['complete'] and v['n']==R for v in counts.values());eps=math.sqrt(math.log(4*196/.0005)/(2*R));formal=[];descriptive=[]
    for c in range(3):
        case,ga,gb,k,t,mask=fixture(c)
        for method,reference in ((3,0),(0,4),(1,0),(2,0)):
            x=pooled.get((c,method));y=pooled.get((c,reference))
            for index,target in enumerate(t[1:],1):
                for species in range(7):
                    for node in range(ga['n']):
                        distance=None if x is None or y is None or not len(x) or not len(y) else cdf_distance(x[:,index,species,node],y[:,index,species,node]);row=dict(cell=c,engine=method,reference_engine=reference,time=float(target),species=species,node=node,n=len(x) if x is not None else 0,n_reference=len(y) if y is not None else 0,distance=distance)
                        if method in (0,3):
                            upper=distance+2*eps if full else None;row.update(epsilon_each=eps if full else None,upper=upper,margin=.08,passed=bool(upper<=.08) if full else None);formal.append(row)
                        else:descriptive.append(row)
    qualified=full and r['formal_passed'] and all(x['passed'] for x in formal);status='QUALIFIED WITHIN REGISTERED SMALL-SYSTEM STOCHASTIC SCOPE' if qualified else 'NOT QUALIFIED UNDER REGISTERED G4 CRITERIA' if full else 'INCOMPLETE G4 STOCHASTIC QUALIFICATION'
    result=dict(status=status,complete=full,scientific_pass=bool(qualified),records=records,groups=counts,jobs_complete=sum(v['complete'] for v in records.values()),paths_available=sum(v['n'] for v in records.values()),expected_paths=245760,family_complete=full,formal=formal,coarse_medium_descriptive=descriptive)
    atomic_json(base/'evaluation.json',result);(base/'MG11G4_results.md').write_text('# MG11-G4 Step1 results\n\nThis qualification compares a bounded simultaneous stochastic network approximation with exact small-system simulations and a separate exact implementation. Its endpoint-law criteria do not establish full-lattice pattern accuracy, a fundamental tick, a physical vacuum or entanglement. Separate audit required.\n\nStatus: **'+status+'**; PI closure pending.\n',encoding='utf-8');return result
if __name__=='__main__':main()

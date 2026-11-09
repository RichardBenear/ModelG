"""MG100 separate incidence/CDF-based numerical recomputation; shared raw kernel."""
from mg11g5 import ROOT,OBS,OMEGA,jobs,factory,contexts,concentration_case
import json
from pathlib import Path
import numpy as np
from scipy.optimize import brentq
from scipy.special import stdtr
from mg11g3_oracle import accounting
from mg11g2_pattern_v2_oracle import dense_projection
from mg11g2_pattern_v2 import reference
from mg11g4_io import check_raw,digest
from mg10d1_pool import atomic_json
BASE=ROOT/'results/MG11G5'
def quantile(p,df):return brentq(lambda z:stdtr(df,z)-p,0,10000)
def audit(directory,evaluation,output):
    directory=Path(directory);e=json.loads(Path(evaluation).read_text(encoding='utf-8'));values={};checked=[];controls=0
    for j in jobs():
        p=directory/(j['id']+'.npz')
        if not p.exists():continue
        m=p.with_suffix('.json');meta=json.loads(m.read_text(encoding='utf-8')) if m.exists() else {};full=meta.get('status')=='COMPLETE' and meta.get('job')==j and meta.get('sha256')==digest(p)
        with np.load(p) as z:
            times=z['times'];full=full and np.array_equal(times,OBS);case,ga,gb,k=factory(j)
            if j['kind']=='stochastic':
                checked.append(check_raw(p));data=z['data'][0];env=z['environment'][0];state=np.concatenate([data.reshape(len(times),-1),env.reshape(len(times),-1),np.zeros((len(times),1000),np.int64)],axis=1)/OMEGA
                # Independent incidence reconstruction in concentration units, in addition to integer check.
                assert accounting(state,concentration_case(case),ga,gb)['passed'];pop=data[:,:7]/OMEGA
                if full:
                    projection=dense_projection(pop,reference()[0]);ix=times>=10;values[(j['setting'],j['block'],j['origin'],j['start'],j['arm'])]=float(np.sum((projection[ix,1][1:]+projection[ix,1][:-1])*.125)/5)
            else:
                assert accounting(z['state'],concentration_case(case),ga,gb)['passed'];controls+=int(full)
    independent=[];blocks=[]
    for v in range(2):
        ds=[]
        for b in range(4):
            keys=[(v,b,c['origin'],c['start'],c['arm']) for c in contexts()]
            if all(key in values for key in keys):
                origin=[]
                for o in range(2):origin.append(values[v,b,o,1,1]-values[v,b,o,0,1]-values[v,b,o,1,0]+values[v,b,o,0,0])
                ds.append(sum(origin)/2)
        blocks.append(ds);expected=e['settings'][v];assert len(ds)==expected['n']
        if len(ds)==4:
            mean=sum(ds)/4;sd=np.sqrt(sum((x-mean)**2 for x in ds)/3);half=quantile(.9875,3)*sd/2
            for key,val in (('mean',mean),('sd',sd),('lower',mean-half),('upper',mean+half)):np.testing.assert_allclose(expected[key],val,rtol=1e-10,atol=1e-10)
            independent.append(dict(mean=mean,lower=mean-half,upper=mean+half))
    agree=False
    if len(independent)==2:
        a=np.var(blocks[0],ddof=1)/4;b=np.var(blocks[1],ddof=1)/4;se=np.sqrt(a+b);df=(a+b)**2/(a*a/3+b*b/3) if se else 6.;d=np.mean(blocks[1])-np.mean(blocks[0]);half=quantile(.95,df)*se
        for key,val in (('mean',d),('lower',d-half),('upper',d+half),('df',df)):np.testing.assert_allclose(e['refinement'][key],val,rtol=1e-10,atol=1e-10)
        agree=d-half>=-.02 and d+half<=.02
    full=len(values)==64 and controls==16;assert e['complete']==full;assert e['direction_consistent']==bool(full and agree and all(x['lower']>0 for x in independent))
    result=dict(passed=True,integer_checks=checked,independent_intervals=independent,complete=full,shared_dependencies='same raw data, donors, template and accounting module; separate dense projection, block algebra and CDF-root intervals; same author, not external replication');atomic_json(output,result);return result
if __name__=='__main__':audit(BASE/'production',BASE/'evaluation.json',BASE/'independent_check.json')

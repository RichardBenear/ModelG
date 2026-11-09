"""MG100 balanced conditional pilot estimator; incomplete blocks never zero."""
from mg11g5 import *
import json
from mg10d1_pool import atomic_json
from mg11g4_io import check_raw,digest
from mg11g3_oracle import accounting
BASE=ROOT/'results/MG11G5'
def evaluate(directory,output):
    directory=Path(directory);records={};details=[];complete=0;controls=0
    for j in jobs():
        p=directory/(j['id']+'.npz');meta=p.with_suffix('.json')
        if not p.exists():continue
        status=json.loads(meta.read_text(encoding='utf-8')) if meta.exists() else {}
        with np.load(p) as raw:
            times=raw['times'];isfull=np.array_equal(times,OBS) and status.get('status')=='COMPLETE' and status.get('job')==j and status.get('sha256')==digest(p)
            if j['kind']=='stochastic':
                check_raw(p);case,ga,gb,k=factory(j)
                for key,value in (('initial',case['initial']),('a0',case['a0']),('b0',case['b0']),('src',ga['src']),('dst',ga['dst']),('srcB',gb['src']),('dstB',gb['dst']),('kappa',k)):
                    np.testing.assert_array_equal(raw[key],value)
                assert int(raw['seed'])==j['seed'] and float(raw['tau'])==TAUS[j['setting']] and float(raw['omega'])==OMEGA
                pop=raw['data'][0,:,:7]/OMEGA
                if isfull:records[(j['setting'],j['block'],j['origin'],j['start'],j['arm'])]=late_value(pop,times);complete+=1
            else:
                integer,ga,gb,k=factory(j);check=accounting(raw['state'],concentration_case(integer),ga,gb)
                if not check['passed']:raise ValueError('Bad deterministic ledger')
                if isfull:controls+=1
                pop=raw['state'][:,:3500].reshape(-1,7,500)
            details.append(dict(job=j,complete=bool(isfull),frames=len(times),diagnostics=metrics(pop,times,full=False),late_X=late_value(pop,times) if isfull else None))
    result=summarize(records);result.update(complete_branches=complete,complete_controls=controls,complete=complete==64 and controls==16,details=details)
    result['direction_consistent']=bool(result['direction_consistent'] and result['complete']);atomic_json(output,result);return result
if __name__=='__main__':evaluate(BASE/'production',BASE/'evaluation.json')

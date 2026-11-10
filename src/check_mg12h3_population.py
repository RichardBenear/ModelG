"""MG115 field/reference/balance/refinement and separately projected bridge audit."""
from mg12h3_population import *
import json
from mg10d1_pool import atomic_json
from mg11g2_pattern_v2_oracle import dense_projection
from mg11g2_pattern_v2 import reference as template
from mg10p1 import continuity
def evaluate(directory,output,full=True):
    directory=Path(directory);details=[];fields={};refs={};complete=0;maximum=0.;primary={};independent={}
    for job in jobs():
        path=directory/(job['id']+'.npz')
        if not path.exists():continue
        with np.load(path) as z:
            state=z['state'];times=z['times'];assert json.loads(str(z['job_json']))==job
        expected=np.arange(int(job['duration']*4)+1)/4;iscomplete=np.array_equal(times,expected);checked=check(job,state,times);assert checked['passed'];maximum=max(maximum,checked['reference_error']);pop=state[:,:3500].reshape(-1,7,500);details.append(dict(job=job,complete=bool(iscomplete),frames=len(times),check=checked));fields[job['id']]=pop
        if iscomplete:
            complete+=1
            if job['setting']==1:
                score=np.asarray(measure(job,state,times)['projection']);other=dense_projection(pop,template()[0]);np.testing.assert_allclose(score,other,atol=1e-10,rtol=1e-10);details[-1]['reference_continuity']=continuity(pop[:,1],times+(10 if job['kind']=='FINITE' else 0))
                if job['kind']=='FINITE':
                    keep=times>=10;key=(job['origin'],job['start'],job['arm']);primary[key]=float(np.trapezoid(score[keep,1],times[keep])/5);independent[key]=float(np.sum((other[keep,1][1:]+other[keep,1][:-1])*.125)/5)
                    ref=factory(job)[5]
                    with np.load(ref) as old:refpop=old['state'][:,:3500].reshape(-1,7,500)
                    rscore=dense_projection(refpop,template()[0]);refs[key]=float(np.sum((rscore[keep,1][1:]+rscore[keep,1][:-1])*.125)/5)
    refinement=[]
    for job in [j for j in jobs() if j['setting']==0]:
        name=job['id'];other=name[:-1]+'1'
        if name in fields and other in fields and fields[name].shape==fields[other].shape:
            initial=factory(job)[0]['initial'];scale=np.maximum(1,initial.max(axis=1));error=float(np.max(abs(fields[name]-fields[other])/scale[None,:,None]));assert error<=1e-5;refinement.append(dict(base=name,fine=other,error=error))
    contrast=None
    if len(primary)==8:
        def combine(values):return sum((1 if s else -1)*(1 if a else -1)*values[o,s,a] for o in range(2) for s in range(2) for a in range(2))/2
        value=combine(primary);audit=combine(independent);old=combine(refs);assert abs(value-audit)<=1e-10;assert abs(value-old)<=1e-5;contrast=dict(new_fine=value,closed_reference=old,separate_projection=audit,error=abs(value-old))
    done=complete==24;qualified=done and len(refinement)==12 and contrast is not None;result=dict(complete=done,qualified=bool(qualified),completed=complete,expected=24,max_reference_error=maximum,refinement=refinement,contrast=contrast,details=details,scope='deterministic population-limit trajectory/reference bridge; no finite-tick stochastic qualification');atomic_json(output,result);return result
if __name__=='__main__':evaluate(ROOT/'results/MG12H3/population/production',ROOT/'results/MG12H3/population/evaluation.json')

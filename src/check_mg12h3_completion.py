"""MG120 mixed-provenance matrix audit; original run is never rescored."""
from mg12h3_completion import *
from mg11g2_pattern_v2_oracle import dense_projection
from mg11g2_pattern_v2 import reference as template
from mg10p1 import continuity
def evaluate(entries,output,held=False):
    verify_reuse(entries);details=[];fields={};values={};independent={};reference_values={};complete=0;accepted=True
    for e in entries:
        job=e['job'];p=ROOT/e['path']
        if not p.exists():details.append(dict(id=job['id'],status='MISSING'));accepted=False;continue
        with np.load(p) as z:
            state=z['state'];times=z['times'];rawjob=json.loads(str(z['job_json']))
            if e['origin']=='reused_validated_MG116':assert rawjob==job
            else:
                assert json.loads(str(z['method_json']))==e['method']
                if held:assert rawjob['kind']==job['kind'] and rawjob['start']==job['start']
                else:assert rawjob==job and str(z['origin'])=='fresh_corrected_replay'
        case,ga,gb,k,mode,reference,parent=factory(job);expected=np.arange(int(job['duration']*4)+1)/4;full=np.array_equal(times,expected);checked=check(job,state,times);accepted=accepted and checked['passed'];pop=state[:,:3500].reshape(-1,7,500);fields[job['id']]=pop;complete+=int(full);row=dict(entry=e,complete=bool(full),check=checked)
        np.testing.assert_array_equal(pop[0],case['initial']);assert not np.any(state[0,3500:]);row['sha256']=hashlib.sha256(p.read_bytes()).hexdigest()
        if e['origin']=='fresh_corrected_replay':
            with np.load(ROOT/e['diagnostic_reference']) as ref:diagnostic=ref['state'];np.testing.assert_array_equal(times,ref['times'])
            scale=np.maximum(1,case['initial'].max(axis=1));replay_error=float(np.max(abs(pop-diagnostic[:,:3500].reshape(-1,7,500))/scale[None,:,None]));row['held_replay_error']=replay_error;accepted=accepted and replay_error<=1e-10
        if full and job['setting']==1:
            scores=np.asarray(measure(job,state,times)['projection']);audit=dense_projection(pop,template()[0]);np.testing.assert_allclose(scores,audit,rtol=1e-10,atol=1e-10);row['reference_continuity']=continuity(pop[:,1],times+(10 if job['kind']=='FINITE' else 0))
            if job['kind']=='FINITE':
                keep=times>=10;key=job['origin'],job['start'],job['arm'];values[key]=float(np.trapezoid(scores[keep,1],times[keep])/5);independent[key]=float(np.sum((audit[keep,1][1:]+audit[keep,1][:-1])*.125)/5)
                with np.load(reference) as ref:old=ref['state'][:,:3500].reshape(-1,7,500)
                score=dense_projection(old,template()[0]);reference_values[key]=float(np.sum((score[keep,1][1:]+score[keep,1][:-1])*.125)/5)
        details.append(row)
    pairs=[]
    for job in [j for j in jobs() if j['setting']==0]:
        other=job['id'][:-1]+'1'
        if job['id'] in fields and other in fields and fields[job['id']].shape==fields[other].shape:
            scale=np.maximum(1,factory(job)[0]['initial'].max(axis=1));error=float(np.max(abs(fields[job['id']]-fields[other])/scale[None,:,None]));pairs.append(dict(left=job['id'],right=other,error=error,passed=error<=1e-5));accepted=accepted and error<=1e-5
    interaction=None
    if len(values)==8:
        def combine(x):return sum((1 if s else -1)*(1 if a else -1)*x[o,s,a] for o in range(2) for s in range(2) for a in range(2))/2
        new=combine(values);old=combine(reference_values);audit=combine(independent);assert abs(new-audit)<=1e-10;interaction=dict(new=new,reference=old,audit=audit,error=abs(new-old));accepted=accepted and abs(new-old)<=1e-5
    done=complete==24 and len(entries)==24;result=dict(complete=done,qualified=bool(done and accepted and len(pairs)==12 and interaction is not None),completed=complete,expected=24,reused=sum(e['origin']=='reused_validated_MG116' for e in entries),fresh_slots=sum(e['origin']=='fresh_corrected_replay' for e in entries),held_reference_pipeline=held,details=details,pairs=pairs,interaction=interaction,original_study_status='MG117 failed/incomplete, unchanged',scope='case-specific corrected deterministic reference matrix, not stochastic/particle qualification');atomic_json(output,result);return result

"""MG092 conditional B-transport contrasts; no statistical replication claims."""
import os
for key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[key]='1'
from pathlib import Path
import json,hashlib
import numpy as np
from mg11g3 import jobs,factory,measure,donor_path,OBS,PARENT_INDEX,S0
from mg10d1_pool import atomic_json
ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/'results/MG11G3'
def read(p):return json.loads(Path(p).read_text(encoding='utf-8'))

def decide(records,refinement,shams):
    full=len(records)==24 and all(v['complete'] for v in records.values());by_setting=[];means=[];off=[]
    for v in (0,1):
        origin=[];off_origin=[]
        for o in (0,1):
            ids=[f'o{o}_s{s}_v{v}_a{a}' for s in (0,1) for a in (0,1)]
            if all(records.get(id,{}).get('complete',False) for id in ids):
                z=lambda s,a:np.array(records[f'o{o}_s{s}_v{v}_a{a}']['metrics']['late_projection'])
                origin.append((z(1,1)-z(0,1)-z(1,0)+z(0,0)).tolist())
            else:origin.append(None)
            ids_off=[f'o{o}_s{s}_v{v}_a{a}' for s in (0,1) for a in (o,2)]
            if all(records.get(id,{}).get('complete',False) for id in ids_off):off_origin.append((np.array(records[f'o{o}_s1_v{v}_a2']['metrics']['late_projection'])-np.array(records[f'o{o}_s0_v{v}_a2']['metrics']['late_projection'])-np.array(records[f'o{o}_s1_v{v}_a{o}']['metrics']['late_projection'])+np.array(records[f'o{o}_s0_v{v}_a{o}']['metrics']['late_projection'])).tolist())
            else:off_origin.append(None)
        by_setting.append(origin);means.append(None if any(x is None for x in origin) else np.mean(origin,axis=0).tolist());off.append(off_origin)
    numerical=full and len(refinement)==12 and all(x['active_error']<=1e-4 for x in refinement);sham_ok=full and len(shams)==8 and all(x['active_error']<=1e-4 for x in shams);delta=None if any(x is None for x in means) else abs(means[1][1]-means[0][1]);threshold=None if delta is None else max(.01,5*delta);resolved=bool(numerical and sham_ok and means[0][1]*means[1][1]>0 and abs(means[1][1])>threshold)
    verdict='INCOMPLETE' if not full else 'NUMERICAL/SHAM QUALIFICATION FAILED' if not numerical or not sham_ok else 'RESOLVED ONGOING B-TRANSPORT CONTRAST' if resolved else 'UNRESOLVED AT REGISTERED RESOLUTION'
    return dict(complete=full,origin_contrasts_GXY=by_setting,primary_GXY_by_setting=means,off_minus_sham_GXY=off,numerical_qualified=bool(numerical),sham_qualified=bool(sham_ok),primary_refinement_difference=delta,primary_resolution_threshold=threshold,primary_resolved=resolved,verdict=verdict)

def main(base=BASE):
    base=Path(base);r=read(base/'registration.json');records={};fields={};shams=[]
    for job in r['jobs']:
        p=base/'production'/(job['id']+'.npz')
        if not p.exists():continue
        case,ga,gb,k=factory(job)
        with np.load(p) as raw:t=raw['times'];state=raw['state'];pop=state[:,:3500].reshape(len(t),7,500);rx=state[:,3500:7000].reshape(len(t),7,500)
        meta=p.with_suffix('.json');full=np.array_equal(t,OBS) and meta.exists() and read(meta).get('status')=='COMPLETE';m=measure(job,pop,t,full);stock=pop[:,3:5]/np.array([case['a0'],case['b0']])[None]
        records[job['id']]=dict(complete=bool(full),checkpoints=len(t),sha256=hashlib.sha256(p.read_bytes()).hexdigest(),metrics=m,stock_mean=stock.mean(axis=-1).tolist(),stock_local=stock[:,:,(250+np.arange(-25,26))%500].mean(axis=-1).tolist(),total_post_A_B_uptake=rx[:,[0,3]].sum(axis=-1).tolist());fields[job['id']]=pop[:,:3].copy()
        if full and job['arm']==job['origin']:
            with np.load(donor_path(job)) as raw:target=raw['state'][PARENT_INDEX:,:1500].reshape(61,3,500)
            shams.append(dict(job=job['id'],active_error=float(np.max(abs(pop[:,:3]-target)/S0[None,:,None]))))
    refinements=[]
    for o in (0,1):
        for s in (0,1):
            for a in (0,1,2):
                ids=[f'o{o}_s{s}_v{v}_a{a}' for v in (0,1)]
                if all(records.get(id,{}).get('complete',False) for id in ids):refinements.append(dict(origin=o,start=s,arm=a,active_error=float(np.max(abs(fields[ids[1]]-fields[ids[0]])/S0[None,:,None]))))
    result=dict(status=('COMPLETE' if len(records)==24 and all(x['complete'] for x in records.values()) else 'INCOMPLETE')+' B-TRANSPORT SWITCH ESTIMATION',jobs_complete=sum(x['complete'] for x in records.values()),records=records,refinement=refinements,shams=shams,**decide(records,refinements,shams));atomic_json(base/'evaluation.json',result)
    (base/'MG11G3_results.md').write_text('# MG11-G3 B-transport switch\n\nThis study compares counterfactual changes to B transport from identical saved Model G states. The conditional deterministic contrasts address ongoing transport rather than a physical source for reservoirs or quantum memory. Separate raw audit is required.\n\nStatus: **'+result['status']+'**. Verdict: **'+result['verdict']+'**. PI closure pending.\n',encoding='utf-8');return result
if __name__=='__main__':main()

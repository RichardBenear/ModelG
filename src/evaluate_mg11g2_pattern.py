"""MG088 deterministic contrasts; refinements are not statistical replicates."""
import os
for key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[key]='1'
from pathlib import Path
import json,hashlib
import numpy as np
from mg11g2_pattern import factory,metrics,CONFIGS,OBS,S0
from mg10d1_pool import atomic_json
ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/'results/MG11G2/pattern_response'
def read(path):return json.loads(Path(path).read_text(encoding='utf-8'))

def decision(records,refinement):
    full=len(records)==28 and all(v['complete'] for v in records.values());d=[];secondary=[]
    def contrast(local,cross,version):
        ids=[f'c{c}_s{s}_v{version}' for c in (local,cross) for s in (0,1)]
        if not all(records.get(key,{}).get('complete',False) for key in ids):return None
        a=[np.array(records[f'c{c}_s1_v{version}']['metrics']['late_projection'])-np.array(records[f'c{c}_s0_v{version}']['metrics']['late_projection']) for c in (local,cross)];return (a[1]-a[0]).tolist()
    for version in (0,1):d.append(contrast(2,3,version));secondary.append(contrast(5,6,version))
    numerical=full and len(refinement)==14 and all(v['active_error']<=1e-4 for v in refinement)
    change=None if any(v is None for v in d) else abs(d[1][1]-d[0][1]);threshold=None if change is None else max(.01,5*change)
    resolved=bool(numerical and d[0][1]*d[1][1]>0 and abs(d[1][1])>threshold)
    maintained=all(records.get(f'c0_s1_v{v}',{}).get('complete',False) and records[f'c0_s1_v{v}']['metrics']['reference_continuity']['persistent_reference_signal'] for v in (0,1))
    verdict='INCOMPLETE' if not full else 'NUMERICAL RESOLUTION FAILED' if not numerical else 'RESOLVED CONNECTIVITY CONTRAST' if resolved else 'UNRESOLVED AT REGISTERED RESOLUTION'
    return dict(complete=full,numerical_qualified=bool(numerical),primary_GXY_by_setting=d,C100_secondary_GXY_by_setting=secondary,primary_refinement_difference=change,primary_resolution_threshold=threshold,primary_resolved=resolved,maintained_reference_control=bool(maintained),pattern_response_supported=bool(resolved and maintained),verdict=verdict)

def main(base=BASE):
    base=Path(base);r=read(base/'registration.json');records={};fields={}
    for job in r['jobs']:
        path=base/'production'/(job['id']+'.npz')
        if not path.exists():continue
        case,g,k,mode=factory(job['config'],job['start'])
        with np.load(path) as raw:times=raw['times'];state=raw['state'];pop=state[:,:7*g['n']].reshape(len(times),7,g['n']);rx=state[:,7*g['n']:14*g['n']].reshape(len(times),7,g['n'])
        meta=path.with_suffix('.json');full=np.array_equal(times,OBS) and meta.exists() and read(meta).get('status')=='COMPLETE';result=metrics(pop,times,full=full);stock=pop[:,3:5]/np.array([case['a0'],case['b0']])[None];local=(250+np.arange(-25,26))%500;remote=np.arange(-25,26)%500
        records[job['id']]=dict(complete=bool(full),checkpoints=len(times),sha256=hashlib.sha256(path.read_bytes()).hexdigest(),metrics=result,stock_mean=stock.mean(axis=-1).tolist(),stock_local=stock[:,:,local].mean(axis=-1).tolist(),stock_remote=stock[:,:,remote].mean(axis=-1).tolist(),total_A_B_uptake=rx[:,[0,3]].sum(axis=-1).tolist())
        fields[job['id']]=pop[:,:3].copy()
    refinement=[]
    for c in range(7):
        for s in (0,1):
            a=f'c{c}_s{s}_v0';b=f'c{c}_s{s}_v1'
            if all(records.get(key,{}).get('complete',False) for key in (a,b)):refinement.append(dict(condition=c,start=s,active_error=float(np.max(abs(fields[a]-fields[b])/S0[None,:,None]))))
    result=dict(status='COMPLETE DETERMINISTIC PATTERN-RESPONSE ESTIMATION' if len(records)==28 and all(v['complete'] for v in records.values()) else 'INCOMPLETE DETERMINISTIC PATTERN-RESPONSE ESTIMATION',jobs_complete=sum(v['complete'] for v in records.values()),records=records,refinement=refinement,**decision(records,refinement))
    atomic_json(base/'evaluation.json',result)
    paragraph='This study compares the same seeded Model G field under matched initial resource stocks and different environmental connectivity, with maintained and isolated-bath controls. Its deterministic contrasts and refinement checks address pattern response within the constructed model; no stochastic confidence interval, physical particle death, native vacuum or quantum claim follows. Separate raw-data audit is required before acceptance.'
    (base/'MG11G2_pattern_results.md').write_text(f'# MG11-G2 pattern response\n\n{paragraph}\n\nStatus: **{result["status"]}**. Verdict: **{result["verdict"]}**; PI closure pending.\n\nJobs {result["jobs_complete"]}/28. Primary G/X/Y contrast by base/fine setting: {result["primary_GXY_by_setting"]}. Registered numerical threshold: {result["primary_resolution_threshold"]}. Maintained reference control: {result["maintained_reference_control"]}. No statistical confidence interpretation of deterministic refinement.\n',encoding='utf-8');return result
if __name__=='__main__':main()

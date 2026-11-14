import os
for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS'):os.environ[k]='1'
import json,gzip,time,copy,traceback
import numpy as np
from run_mg180 import OUT,RAW,models,read,atomic_json,result
from mg180_background import make_tracker

def run():
 start=time.perf_counter();assert not(OUT/'nuisance_attempt.json').exists();atomic_json(OUT/'nuisance_attempt.json',dict(status='RUNNING'));rows=[]
 try:
  assert read(OUT/'attempt.json')['status']=='COMPLETE';model,L=next((m,L) for ri,m,L in models() if ri==2)
  for case in ('moving','displacedY'):
   tag='r2_rho0_scale1.25_'+case;d=json.loads(gzip.decompress((RAW/(tag+'_records.json.gz')).read_bytes()))
   with np.load(RAW/(tag+'.npz')) as z:fields=z['fields']
   for p,seq in zip(d['paths'],fields):
    for variable,ci in [('scale',p['calibration']['ci95']),('rho',p['calibration']['rho_ci95'])]:
     for index,value in enumerate(ci):
      cal=copy.deepcopy(p['calibration']);cal[variable]=float(np.clip(value,-.99,.99)) if variable=='rho' else value;obj=make_tracker(model,L,p['seed'],cal);records=[]
      for field in seq:
       r=obj.step(field);r['emitted']=r['label'] if cal['adequate'] else 'noise_unqualified';records.append(r)
      rows.append(dict(case=case,replicate=p['replicate'],variable=variable,endpoint='low' if index==0 else 'high',value=cal[variable],outcome=result(records,p['truth'],case,True),changed_raw_labels=sum(r['label']!=o['label'] for r,o in zip(records,p['records'])),deadline15=records[15]['emitted'],labels=[r['emitted'] for r in records],max_class_weight_difference=max(abs(r['class_weights'][k]-o['class_weights'][k]) for r,o in zip(records,p['records']) for k in r['class_weights'])))
      assert time.perf_counter()-start<300
  atomic_json(OUT/'nuisance_diagnostic.json',dict(rows=rows,elapsed_seconds=time.perf_counter()-start,no_gate=True,no_best_variant_selected=True,limits='One-parameter endpoint sensitivity only; not integration or a full joint confidence-region guarantee.'));atomic_json(OUT/'nuisance_attempt.json',dict(status='COMPLETE',elapsed_seconds=time.perf_counter()-start))
 except BaseException:atomic_json(OUT/'nuisance_attempt.json',dict(status='FAILED',error=traceback.format_exc(),elapsed_seconds=time.perf_counter()-start));raise
if __name__=='__main__':run()

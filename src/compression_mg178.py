import os
for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS'):os.environ[k]='1'
import json,gzip,time,traceback
import numpy as np
from run_mg178 import OUT,RAW,models,read,atomic_json,result
from mg178_background import make_tracker

def run():
 start=time.perf_counter();assert not(OUT/'compression_attempt.json').exists();atomic_json(OUT/'compression_attempt.json',dict(status='RUNNING'));rows=[]
 try:
  assert read(OUT/'attempt.json')['status']=='COMPLETE';model,L=next((m,L) for ri,m,L in models() if ri==2)
  for case in ('moving','split'):
   tag='r2_rho0_scale1.25_'+case;d=json.loads(gzip.decompress((OUT/(tag+'.json.gz')).read_bytes()))
   with np.load(RAW/(tag+'.npz')) as z:fields=z['fields']
   for p,seq in zip(d['paths'],fields):
    for keep in(4,8,16):
     obj=make_tracker(model,L,p['seed'],p['calibration'],keep);records=[]
     for field,old in zip(seq,p['records']):
      r=obj.step(field);r['emitted']=r['label'] if p['calibration']['adequate'] else 'noise_unqualified';records.append(r)
      if keep==8:assert r['label']==old['label'];np.testing.assert_allclose(list(r['class_weights'].values()),list(old['class_weights'].values()),rtol=1e-9,atol=1e-10)
     rows.append(dict(case=case,replicate=p['replicate'],keep=keep,outcome=result(records,p['truth'],case,True),changed_raw_labels=sum(r['label']!=o['label'] for r,o in zip(records,p['records'])),multiple_weights=[r['class_weights']['multiple'] for r in records],labels=[r['emitted'] for r in records],max_class_weight_difference=max(abs(r['class_weights'][k]-o['class_weights'][k]) for r,o in zip(records,p['records']) for k in r['class_weights'])))
     assert time.perf_counter()-start<300
  atomic_json(OUT/'compression_diagnostic.json',dict(rows=rows,elapsed_seconds=time.perf_counter()-start,primary_keep=8,no_new_gate=True,no_best_keep_selected=True));atomic_json(OUT/'compression_attempt.json',dict(status='COMPLETE',elapsed_seconds=time.perf_counter()-start))
 except BaseException:atomic_json(OUT/'compression_attempt.json',dict(status='FAILED',error=traceback.format_exc(),elapsed_seconds=time.perf_counter()-start));raise
if __name__=='__main__':run()

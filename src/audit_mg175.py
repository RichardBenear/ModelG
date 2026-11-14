import os
for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS'):os.environ[k]='1'
import json,gzip,hashlib,time,collections
from pathlib import Path
import numpy as np
from run_mg174 import ROOT,RAW,OUT as OLD,POSITIVES,atomic_json
from check_mg174 import inner,mu
OUT=ROOT/'results/MG12H3/initialization_noise_audit/v2'
def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def inputs():return sorted(OLD.glob('r*_rho*.json.gz'))+sorted(RAW.glob('r*_rho*_static_snap10.npz'))+sorted(RAW.glob('r*_rho*.npz'))
def register():
 OUT.mkdir(exist_ok=True,parents=True);assert not(OUT/'plan.json').exists(); files=inputs()+[ROOT/'src/audit_mg175.py',ROOT/'spec/MG175_initialization_noise_audit.md'];atomic_json(OUT/'plan.json',dict(files={str(p.relative_to(ROOT)):digest(p) for p in files},cap_seconds=120,no_draws=True,no_relabel=True))
def run():
 start=time.perf_counter();assert not(OUT/'attempt.json').exists();atomic_json(OUT/'attempt.json',dict(status='RUNNING'))
 try:
  plan=json.loads((OUT/'plan.json').read_text());assert all(digest(ROOT/p)==h for p,h in plan['files'].items());counts=collections.defaultdict(collections.Counter);uncertain=[];initial=[];oracle=[];checked=0;quadratics=0;maxerr=0.
  for p in sorted(OLD.glob('r*_rho*.json.gz')):
   d=json.loads(gzip.decompress(p.read_bytes()));case=d['case'];ri=0 if d['regime']==.1 else 2
   for path in d['paths']:
    for method,records in path['records'].items():
     for t,r in enumerate(records):
      fails=[]
      if r['class_weights'][r['winner']]<.95:fails.append('weight')
      if r['innovation_q']>r['qcut']:fails.append('innovation')
      if r['winner']=='single' and (r['single_center_spread']>4 or r['single_width_spread']>.25):fails.append('spread')
      expected='uncertain' if fails else r['winner'];assert expected==r['label'];checked+=1
      key=f"{d['regime']:g}/{d['rho']:g}/{case}/{method}";counts[key]['frames']+=1;counts[key]['label_'+r['label']]+=1
      for f in fails:counts[key]['failed_'+f]+=1
      if method=='modes' and d['primary'] and case in POSITIVES and r['label']!='single':uncertain.append(dict(regime=d['regime'],rho=d['rho'],case=case,replicate=path['replicate'],frame=t,fails=fails,**{k:r[k] for k in ['winner','class_weights','innovation_q','qcut','single_center_spread','single_width_spread']}))
      if method=='modes' and case in POSITIVES and t<4:initial.append(dict(regime=d['regime'],rho=d['rho'],case=case,replicate=path['replicate'],frame=t,label=r['label'],weights=r['class_weights'],motion=r['single_motion_weights']))
   if case in ('static','stop_start','noise_scaled'):
    with np.load(RAW/(p.name[:-8]+'.npz')) as z:fields=z['fields']
    with np.load(RAW/f'r{ri}_rho{d["rho"]:g}_static_snap10.npz') as z:ft=z['FT'];P=z['precision'];w=z['weights']
    for path,seq in zip(d['paths'],fields):
     residuals=[]
     for field,(k,width) in zip(seq,path['truth']):
      y=np.fft.rfft(field-field.mean(axis=1)[:,None],axis=1);residuals.append(y-mu('single',int(k),round((width-.5)/.125),None,ft))
     values=[]
     for t,e in enumerate(residuals):
      rho=d['rho'] if t else 0.;e=e-rho*residuals[t-1] if t else e;var=1-rho*rho;q=inner(e,e,P,w)/var
      # Independent real-space convolution with inverse spectral covariance.
      real=np.fft.irfft(e,n=500,axis=1);filtered=np.fft.irfft(np.einsum('fab,bf->af',P,e),n=500,axis=1);qr=float(np.sum(real*filtered)/var);np.testing.assert_allclose(q,qr,rtol=1e-8,atol=1e-6);maxerr=max(maxerr,abs(q-qr));quadratics+=1;values.append(q)
     cut=path['records']['modes'][0]['qcut'];factor=1.5625 if case=='noise_scaled' else 1.;oracle.append(dict(regime=d['regime'],rho=d['rho'],case=case,replicate=path['replicate'],qcut=cut,oracle_q=values,exceedances=sum(q>cut for q in values),known_scale_corrected_exceedances=sum(q/factor>cut for q in values),selected_q=[r['innovation_q'] for r in path['records']['modes']]))
   assert time.perf_counter()-start<120
  assert checked==9024;atomic_json(OUT/'audit.json',dict(guard_counts=dict(counts),uncertain_intact_primary=uncertain,initial_frames=initial,oracle=oracle,unchanged_labels_checked=checked,real_space_quadratics_checked=quadratics,max_absolute_quadratic_difference=maxerr,no_new_verdict=True,elapsed_seconds=time.perf_counter()-start));atomic_json(OUT/'attempt.json',dict(status='COMPLETE',elapsed_seconds=time.perf_counter()-start))
 except BaseException:
  import traceback
  atomic_json(OUT/'attempt.json',dict(status='FAILED',error=traceback.format_exc(),elapsed_seconds=time.perf_counter()-start));raise
if __name__=='__main__':register() if __import__('sys').argv[1]=='register' else run()

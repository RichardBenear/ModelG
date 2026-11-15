"""One bounded MG185 saved-history audit. No new dynamics or tuned instrument."""
import os
for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS'):os.environ[key]='1'
import json,time,hashlib
from pathlib import Path
from collections import Counter
import numpy as np
from mg184_phase2 import ref,initial
from run_mg185_phase2_G import ROOT,RAW,PUBLIC,read,binding
from mg10e2c1_pool import atomic_json
OUT=ROOT/'results/MG187_G_loss_audit'

def fade(times,amplitude):
 smooth=np.full(len(times),np.nan);smooth[20:]=np.convolve(amplitude,np.ones(21)/21,mode='valid')
 for i in range(20,len(times)-20):
  if np.all(smooth[i:i+21]<.5):return dict(onset=float(times[i]),confirmed=float(times[i+20]))
 return None

def ordering(g,x):
 if g is None or x is None:return 'unresolved/censored'
 return 'before' if g-x < -2 else 'after' if g-x>2 else 'together'

def tests():
 times=np.arange(201)*.5
 def drop(t):return np.where(times<t,1.,0.)
 a=fade(times,drop(40))['onset'];b=fade(times,drop(50))['onset'];assert ordering(a,b)=='before' and ordering(b,a)=='after' and ordering(a,a)=='together';assert fade(times,np.ones(201)) is None and ordering(None,a)=='unresolved/censored'

def metric(profile,unit,norm):
 amplitude=profile@unit/norm;orthogonal=profile-(profile@unit)[:,None]*unit
 residual=np.linalg.norm(orthogonal,axis=1)/norm
 np.testing.assert_allclose(np.sum(profile**2,axis=1)/norm**2,amplitude**2+residual**2,atol=1e-9,rtol=1e-10)
 return amplitude,residual

def run():
 tests();OUT.mkdir(exist_ok=False);start=time.perf_counter();atomic_json(OUT/'attempt.json',dict(status='RUNNING',new_simulations=False))
 try:
  assert binding()==read(PUBLIC/'plan.json')['source'];manifest=read(RAW/'analysis_segment_manifest.json');index={v['id']:v for v in manifest['runs']};export={v['id']:v for v in read(PUBLIC/'decision_exports.json')['runs']};plan=read(PUBLIC/'plan.json');reference=ref();template=reference['phi']-reference['phi'].mean(axis=1,keepdims=True);norm=np.linalg.norm(template,axis=1);unit=template/norm[:,None];times=np.arange(2001)*.5;k=2*np.pi*np.fft.fftfreq(500,d=.2);lossrows=[];high=[];curves=[];verified=0;projection_error=0.
  jobs=[j for j in plan['jobs'] if j['phase']=='main' and j['condition']=='prepared' and j['omega'] in (1000,3000,10000)]
  for ji,job in enumerate(jobs):
   assert time.perf_counter()-start<600;old=export[job['id']];directory=RAW/old['phase']/job['id'];saved=np.load(directory/'GXY_measurements.npz');valid=saved['unique_original'];positions=saved['positions'];loss=[v[1] for v in old['events'] if v[0]=='loss'];loss=min(loss) if loss else None
   anchor=np.zeros(2001);last=0.
   for i in range(2001):
    if valid[i] and (loss is None or times[i]<=loss):last=positions[i]
    anchor[i]=last
   counts=[initial(job['omega'],'prepared',job['seed'],reference)[None]]
   for segment in index[job['id']]['segments']:
    path=ROOT/segment['path'];assert hashlib.sha256(path.read_bytes()).hexdigest()==segment['sha256']
    with np.load(path) as a:counts.append(a['counts'].copy())
    verified+=1
   field=np.concatenate(counts)/job['omega'];field-=field.mean(axis=2,keepdims=True);shift=np.exp(1j*anchor[:,None]*k);aligned=np.fft.ifft(np.fft.fft(field,axis=2)*shift[:,None,:],axis=2).real
   amplitudes=np.sum(aligned*unit[None],axis=2)/norm[None];compare=valid&((times<=loss) if loss is not None else True);error=float(np.max(abs(amplitudes[compare]-saved['aligned_relative_amplitude'][compare])));projection_error=max(projection_error,error);np.testing.assert_allclose(amplitudes[compare],saved['aligned_relative_amplitude'][compare],atol=2e-10,rtol=1e-9)
   ga,gr=metric(aligned[:,0],unit[0],norm[0]);np.testing.assert_allclose(ga,amplitudes[:,0],atol=1e-12)
   if job['omega']==1000:
    gf=fade(times,amplitudes[:,0]);xf=fade(times,amplitudes[:,1]);row=dict(id=job['id'],X_detector_loss=loss,G_half_amplitude=gf,X_half_amplitude=xf,G_minus_detector_loss=None if gf is None or loss is None else gf['onset']-loss,G_minus_X_half=None if gf is None or xf is None else gf['onset']-xf['onset'],relative_to_detector=ordering(None if gf is None else gf['onset'],loss),relative_to_X_half=ordering(None if gf is None else gf['onset'],None if xf is None else xf['onset']),held_anchor_after_loss=True,windows=[])
    if loss is not None:
     for left,right in ((-20,-10),(-10,0),(0,10),(10,20)):
      mask=(times>=loss+left)&(times<loss+right);row['windows'].append(dict(left=left,right=right,n=int(mask.sum()),unique_original_frames=int(valid[mask].sum()),G_amplitude=float(ga[mask].mean()) if mask.any() else None,X_amplitude=float(amplitudes[mask,1].mean()) if mask.any() else None,G_residual=float(gr[mask].mean()) if mask.any() else None))
     relative=np.arange(-60,61)*.5;indices=np.rint((loss+relative)*2).astype(int);ok=(indices>=0)&(indices<len(times));curve=np.full((121,2),np.nan);curve[ok]=amplitudes[indices[ok],:2];curves.append(curve)
    lossrows.append(row)
   else:
    blocks=[]
    for name,mask in [('initial',(times>=0)&(times<100)),('final',(times>=900)&(times<=1000))]:
     keep=mask&valid;mean=aligned[keep,0].mean(axis=0) if keep.any() else None
     if mean is None:blocks.append(dict(name=name,valid=0,total=int(mask.sum()),amplitude=None,amplitude_SD=None,instantaneous_residual=None,averaged_profile_residual=None));continue
     ma,mr=metric(mean[None],unit[0],norm[0]);block_index=0 if name=='initial' else 9;np.testing.assert_allclose(ma[0],old['blocks'][block_index]['aligned'][0]/norm[0],atol=1e-10)
     blocks.append(dict(name=name,valid=int(keep.sum()),total=int(mask.sum()),amplitude=float(ma[0]),amplitude_SD=float(ga[keep].std(ddof=1)),instantaneous_residual=float(gr[keep].mean()),averaged_profile_residual=float(mr[0])))
    high.append(dict(id=job['id'],omega=job['omega'],X_persistent=old['X_persistent'],strict=old['strict'],blocks=blocks))
   if ji%10==0:print('audited',ji+1,'/90',flush=True)
  assert len(lossrows)==30 and sum(v['X_detector_loss'] is not None for v in lossrows)==29 and len(high)==60
  actual=[v for v in lossrows if v['X_detector_loss'] is not None];curve=np.array(curves);cohort=[]
  for omega in (3000,10000):
   for block in ('initial','final'):
    rows=[next(v for v in r['blocks'] if v['name']==block) for r in high if r['omega']==omega];values={}
    for name in ('amplitude','amplitude_SD','instantaneous_residual','averaged_profile_residual'):
     v=[r[name] for r in rows if r[name] is not None];values[name]=dict(n=len(v),median=float(np.median(v)),min=float(min(v)),max=float(max(v)))
    cohort.append(dict(omega=omega,block=block,metrics=values,total_valid_frames=sum(r['valid'] for r in rows),total_frames=sum(r['total'] for r in rows)))
  atomic_json(OUT/'audit.json',dict(loss_runs=lossrows,high_population_runs=high,cohort=cohort,counts_detector=dict(Counter(v['relative_to_detector'] for v in actual)),counts_half_amplitude=dict(Counter(v['relative_to_X_half'] for v in actual)),loss_curve=dict(relative_time=(np.arange(-60,61)*.5).tolist(),mean=np.nanmean(curve,axis=0).tolist(),contributors=np.isfinite(curve[:,:,0]).sum(axis=0).tolist()),posthoc_exploratory=True,no_new_G_gate=True))
  atomic_json(OUT/'verification.json',dict(passed=True,trajectories=90,segments=verified,maximum_saved_projection_error=projection_error,orthogonal_energy_identity=True,high_population_block_means_match_MG185=True,synthetic_timing_fixtures=True,external_replication=False,shared_dependencies=['reference','X track and loss event','NumPy FFT'],elapsed_seconds=time.perf_counter()-start))
  atomic_json(OUT/'attempt.json',dict(status='COMPLETE',elapsed_seconds=time.perf_counter()-start,new_simulations=False,bounded_passes=1))
 except BaseException:
  import traceback
  atomic_json(OUT/'attempt.json',dict(status='FAILED',error=traceback.format_exc(),elapsed_seconds=time.perf_counter()-start));raise
if __name__=='__main__':run()

"""MG191 post-run raw-field audit and frozen-gate reductions; no new integration."""
import os
for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[key]='1'
import time,hashlib
import numpy as np
from scipy.signal import find_peaks
from run_mg191 import ROOT,PUBLIC,RAW,read,binding
from mg10e2c1_pool import atomic_json

def independent_core(v,h,depth):
 n=len(v);padded=np.tile(v,3);peaks=find_peaks(-padded,height=.5*depth)[0];peaks=peaks[(peaks>=n)&(peaks<2*n)];peaks=peaks[padded[peaks]<-.5*depth];count=len(peaks)
 if count!=1:return np.array([np.nan,np.nan,np.nan,float(count)])
 center=int(peaks[0]);crossings=np.flatnonzero(np.signbit(padded[:-1])!=np.signbit(padded[1:]));left=crossings[crossings<center];right=crossings[crossings>=center]
 if not len(left) or not len(right):return np.array([np.nan,np.nan,np.nan,float(count)])
 zeros=[]
 for i in (int(left[-1]),int(right[0])):zeros.append(i-padded[i]/(padded[i+1]-padded[i]))
 coordinate=((np.mean(zeros)-n-n//2)*h+50)%100-50
 return np.array([coordinate,-padded[center],(zeros[1]-zeros[0])*h,float(count)])

def grade(rows,mirrors,level_flags):
 eligible=[];mesh_down=[];mesh_up=[]
 for eps in (.00025,.0005,.001):
  coarse=[r for r in rows if r['arm']=='A' and r['eps']==eps];reverse=[r for r in rows if r['arm']=='R' and abs(r['eps'])==eps];fine=[r for r in rows if r['arm']=='H' and r['eps']==eps]
  if len(coarse)==len(reverse)==8 and all(r['valid_late'] and r['resolved'] and r['mobility']>0 for r in coarse+reverse) and all(r['passed'] for r in mirrors if r['eps']==eps):eligible.append(eps)
  if len(coarse)==len(fine)==8 and all(r['valid_late'] and r['resolved'] and r['mobility'] is not None and r['mobility']!=0 for r in coarse+fine):
   agree=all(abs(a['mobility']-next(b for b in fine if b['phase']==a['phase'])['mobility'])/abs(next(b for b in fine if b['phase']==a['phase'])['mobility'])<=.1 for a in coarse)
   if agree and all(r['mobility']>0 for r in coarse+fine):mesh_down.append(eps)
   if agree and all(r['mobility']<0 for r in coarse+fine):mesh_up.append(eps)
 level_block=any(v['median_above_floor'] or v['maximum_block_above_floor'] for v in level_flags);p1=len(eligible)>=2;p4=len(mesh_down)>=2
 if p1 and p4 and not level_block:classification='SUPPORTED'
 elif len(mesh_up)>=2:classification='CONTRADICTED'
 elif len([r for r in rows if r['arm']=='H'])>=16 and all(r['valid_late'] for r in rows if r['arm'] in ('A','R','H')) and not any(r['resolved'] for r in rows if r['arm'] in ('A','R','H')):classification='NO RESOLVED RESPONSE'
 else:classification='UNRESOLVED'
 reasons=[]
 if not any(r['resolved'] for r in rows if r['arm'] in ('A','R')):reasons.append('No source/reversed cases exceeded the fixed ten-times-null floor')
 if not any(r['arm']=='H' for r in rows):reasons.append('Both refinement amplitudes excluded by uniform-source viability, so no fine gradient response exists for mesh comparison')
 if level_block:reasons.append('Level-control drift above matching null floor in some cases')
 return dict(classification=classification,P1_supported=p1,P4_supported=p4,P2='NOT_QUALIFIED_NO_RESOLVED_SOURCE_CASES' if not any(r['arm']=='A' and r['resolved'] for r in rows) else 'DESCRIPTIVE_ONLY',P3='NOT_QUALIFIED_NO_RESOLVED_SOURCE_OR_CRITICALITY_CASES' if not any(r['arm'] in ('A','C') and r['resolved'] for r in rows) else 'DESCRIPTIVE_ONLY',reasons=reasons,not_contradicted=classification!='CONTRADICTED',no_true_zero_response_claim=True)

def run():
 start=time.perf_counter();plan=read(PUBLIC/'plan.json');assert plan['source']==binding() and read(PUBLIC/'attempt.json')['status']=='COMPLETE';q=read(RAW/'qualification/ledger.json');p=read(RAW/'production/ledger.json');allowed=read(PUBLIC/'qualification.json');assert len(q['records'])==31 and len(p['records'])==len(allowed['admitted_jobs'])==104;assert all(r['status']=='COMPLETE' for r in q['records']+p['records']);cache=np.load(ROOT/'results/MG003/ref_profiles.npz');source_jobs={j['id']:j for j in plan['jobs']+plan['qualification_jobs']};measurements={};manifest=[];maxerror=0.;frames=0;segments=0
 for key,row in plan['inputs'].items():assert hashlib.sha256((RAW/'inputs'/(key+'.npz')).read_bytes()).hexdigest()==row['sha256']
 for phase,ledger in (('qualification',q),('production',p)):
  for ji,record in enumerate(ledger['records']):
   assert time.perf_counter()-start<plan['analysis_reserve'];job=record['job'];assert job==source_jobs[job['id']]
   if job['kind']=='background':continue
   directory=RAW/phase/job['id'];meta=read(directory/'metadata.json');ref=cache[f"h{job['h']}_b29.0"];refcore=independent_core(ref[2],job['h'],abs(ref[2].min()));np.testing.assert_allclose(refcore[1:3],[meta['reference']['depth'],meta['reference']['width']],atol=1e-10);inp=np.load(RAW/'inputs'/(job['background_key']+'.npz'));bg=inp['background'];phi=ref[1:];wave=2*np.pi*np.fft.fftfreq(phi.shape[1],d=job['h']);initial=bg+np.fft.ifft(np.fft.fft(phi,axis=1)*np.exp(-1j*wave*job.get('position',0.)),axis=1).real;previous_fields=initial;previous_position=job.get('position',0.);lost=False;observations=[];valid=[];hashes=[]
   for chunk in range(job['duration']//50):
    path=directory/f'chunk_{chunk:03d}.npz';hashes.append(hashlib.sha256(path.read_bytes()).hexdigest())
    with np.load(path) as archive:
     field=archive['fields'];assert field.shape==(101,3,bg.shape[1]) and np.isfinite(field).all() and field.min()>0;np.testing.assert_allclose(field[0],previous_fields,atol=1e-12,rtol=0);np.testing.assert_array_equal(archive['initial'],field[0]);np.testing.assert_array_equal(archive['final'],field[-1]);np.testing.assert_array_equal(archive['times'],chunk*50+np.arange(101)*.5);previous_fields=field[-1].copy()
     for frame in (field if chunk==0 else field[1:]):
      v=independent_core(frame[1]-bg[1],job['h'],meta['reference']['depth'])
      if not np.isfinite(v[:3]).all() or v[3]!=1:lost=True
      elif previous_position is not None and abs((v[0]-previous_position+50)%100-50)>meta['reference']['wavelength']/2:lost=True
      ok=not lost
      if ok:previous_position=float(v[0])
      observations.append(v);valid.append(ok);frames+=1
    segments+=1
   observations=np.array(observations);valid=np.array(valid);saved=np.load(directory/'observations.npz');np.testing.assert_array_equal(valid,saved['original_valid']);np.testing.assert_allclose(observations,saved['core'],atol=1e-10,rtol=0,equal_nan=True);finite=np.isfinite(observations)&np.isfinite(saved['core']);maxerror=max(maxerror,float(np.max(abs(observations[finite]-saved['core'][finite]))) if finite.any() else 0.);times=np.arange(job['duration']*2+1)*.5;np.testing.assert_array_equal(times,saved['times']);shape=valid&(observations[:,1]>=.6*refcore[1])&(observations[:,2]>=.75*refcore[2])&(observations[:,2]<=1.25*refcore[2])
   if phase=='qualification':assert record['passed']==bool(shape[times>=100].all() and record['uniform_stability']['passed'])
   else:
    assert record['valid_late']==bool(shape[times>=500].all());unwrapped=np.unwrap(observations[:,0]*np.pi/50)*50/np.pi;gradient=np.interp(observations[:,0],inp['x'],inp['gradient'],period=100);blocks=[]
    for b in record['blocks']:
     i=b['left']*2;j=b['right']*2;ok=bool(shape[i:j+1].all());assert ok==b['valid'];velocity=float((unwrapped[j]-unwrapped[i])/50) if ok else None;force=float(-gradient[i:j+1].mean()) if ok else None
     if ok:np.testing.assert_allclose(velocity,b['velocity'],atol=2e-12,rtol=1e-6);np.testing.assert_allclose(force,b['minus_G_gradient'],atol=1e-13)
     blocks.append(dict(left=b['left'],right=b['right'],velocity=velocity,force=force,valid=ok))
    good=[b for b in blocks if b['valid']];f=np.array([b['force'] for b in good]);v=np.array([b['velocity'] for b in good]);mobility=None if f@f<1e-24 else float(np.linalg.lstsq(f[:,None],v,rcond=None)[0][0]);speed=float(np.median(abs(v))) if len(v) else None
    if mobility is not None:np.testing.assert_allclose(mobility,record['mobility'],atol=2e-7,rtol=1e-6)
    if speed is not None:np.testing.assert_allclose(speed,record['median_abs_late_speed'],atol=2e-12,rtol=1e-6)
    measurements[job['id']]=dict(id=job['id'],arm=job['arm'],h=job['h'],eps=job['eps'],phase=job['phase'],valid_late=record['valid_late'],mobility=mobility,speed=speed,blocks=blocks,final_distance_to_minimum=record['final_distance_to_minimum'],first_near_minimum=record['first_near_minimum'],near_minimum_censored=record['near_minimum_censored'])
   manifest.append(dict(id=job['id'],chunk_sha256=hashes,observation_sha256=hashlib.sha256((directory/'observations.npz').read_bytes()).hexdigest()))
  print('audited',phase,flush=True)
 floors={(r['h'],r['phase']):max(1e-10,max(abs(b['velocity']) for b in r['blocks'] if b['valid'])) for r in measurements.values() if r['arm']=='N'};rows=[]
 for row in measurements.values():
  row['null_floor']=floors[row['h'],row['phase']];row['speed_to_floor']=None if row['speed'] is None else row['speed']/row['null_floor'];row['resolved']=bool(row['valid_late'] and row['speed'] is not None and row['speed']>10*row['null_floor']);rows.append(row)
 aggregate=[]
 for arm in ('A','R','C'):
  for eps in (.00025,.0005,.001):
   rr=[r for r in rows if r['arm']==arm and abs(r['eps'])==eps];mus=[r['mobility'] for r in rr if r['mobility'] is not None];aggregate.append(dict(arm=arm,eps=eps,n=len(rr),valid=sum(r['valid_late'] for r in rr),resolved=sum(r['resolved'] for r in rr),raw_positive_mu=sum(v>0 for v in mus),raw_negative_mu=sum(v<0 for v in mus),median_raw_mu=float(np.median(mus)),raw_mu_range=[float(min(mus)),float(max(mus))],median_speed=float(np.median([r['speed'] for r in rr])),max_speed_to_floor=float(max(r['speed_to_floor'] for r in rr))))
 mirrors=[]
 for row in rows:
  if row['arm']!='A':continue
  other=next(r for r in rows if r['arm']=='R' and r['phase']==row['phase'] and abs(r['eps'])==row['eps']);err=max(abs(b['velocity']+c['velocity']) for b,c in zip(row['blocks'],other['blocks']) if b['valid'] and c['valid']);mirrors.append(dict(eps=row['eps'],phase=row['phase'],max_velocity_sum=err,envelope=row['null_floor']+other['null_floor'],passed=err<=row['null_floor']+other['null_floor']))
 level_flags=[dict(id=r['id'],median_above_floor=r['speed']>r['null_floor'],maximum_block_above_floor=max(abs(b['velocity']) for b in r['blocks'] if b['valid'])>r['null_floor']) for r in rows if r['arm']=='L'];verdict=grade(rows,mirrors,level_flags)
 charged=read(PUBLIC/'attempt.json')['charged_worker_seconds'];seconds=time.perf_counter()-start;assert seconds<540 and charged+seconds<7200;manifest_path=RAW/'analysis_manifest.json';atomic_json(manifest_path,dict(trajectories=manifest));raw_size=sum(p.stat().st_size for p in RAW.rglob('*') if p.is_file())
 atomic_json(PUBLIC/'decision_exports.json',dict(verdict=verdict,runs=rows,aggregate=aggregate,mirror_checks=mirrors,level_flags=level_flags,coarse_null_floors=[dict(phase=k[1],floor=v) for k,v in floors.items() if k[0]==.2],qualified_production=104,excluded_production=16,charged_worker_seconds_before_analysis=charged,analysis_seconds=seconds,raw_bytes=raw_size,finite_phase_ranges_not_confidence_intervals=True))
 atomic_json(PUBLIC/'verification.json',dict(passed=True,raw_segments=segments,raw_snapshots=frames,trajectories=120,independent_zero_crossings_via_tiled_peak_and_sign_changes=True,maximum_core_measurement_error=maxerror,checkpoint_continuity=True,shape_and_identity_flags_recomputed=True,mobility_recomputed_by_lstsq=True,input_and_source_hashes_unchanged=True,manifest_sha256=hashlib.sha256(manifest_path.read_bytes()).hexdigest(),seconds=seconds,external_replication=False,shared_dependencies=['NumPy/SciPy','same cached reference and frozen detector thresholds'],limits='Nonuniform eigenvalue calculations are those of registered qualification, not independently re-solved here. This saved-field audit checks numerical records, not physical gravity.'))
 print(verdict,'seconds',seconds,'raw bytes',raw_size)
if __name__=='__main__':run()

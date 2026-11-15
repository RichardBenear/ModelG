"""Post-run joined MG185 analysis and separate raw projection recomputation; no draws."""
import os
for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS'):os.environ[key]='1'
import json,hashlib,time
from pathlib import Path
import numpy as np
from scipy.stats import t,beta
from mg184_phase2 import ref,initial
from run_mg185_phase2_G import ROOT,PUBLIC,RAW,read,binding
from mg10e2c1_pool import atomic_json

def interval(values):
 v=np.asarray(values,float);n=len(v)
 if not n:return dict(n=0,mean=None,CI95=None)
 mean=float(v.mean());half=float(t.ppf(.975,n-1)*v.std(ddof=1)/np.sqrt(n)) if n>1 else None
 return dict(n=n,mean=mean,CI95=None if half is None else [mean-half,mean+half])

def cp(k,n):return [0. if k==0 else float(beta.ppf(.025,k,n-k+1)),1. if k==n else float(beta.ppf(.975,k+1,n-k))]

def run():
 start=time.perf_counter();plan=read(PUBLIC/'plan.json');continuation=read(PUBLIC/'continuation_plan.json');assert binding()==plan['source'];assert hashlib.sha256((ROOT/'src/resume_mg185.py').read_bytes()).hexdigest()==continuation['resume_sha256']
 original=read(RAW/'main/ledger.json');resumed=read(RAW/'continuation/ledger.json');qual=read(RAW/'qualification/ledger.json');jobs=[j for j in plan['jobs'] if j['phase']=='main'];chosen={}
 for phase,ledger in [('main',original),('continuation',resumed)]:
  for record in ledger['records']:
   if record['status']=='COMPLETE':
    assert record['id'] not in chosen;chosen[record['id']]=(phase,record)
 assert len(chosen)==180 and set(chosen)=={j['id'] for j in jobs}
 reference=ref();template=reference['phi']-reference['phi'].mean(axis=1,keepdims=True);norm=np.sqrt(np.sum(template**2,axis=1));unit=template/norm[:,None];exports=[];manifests=[];checks=0;maxerror=0.;segments=0
 for index,job in enumerate(jobs):
  phase,record=chosen[job['id']];assert record['job']==job and record['clips']==0
  directory=RAW/phase/job['id'];summary=read(directory/'summary.json');measure=np.load(directory/'GXY_measurements.npz');assert summary['horizon']==1000 and summary['full_GXY_persistence']=='NOT_SCORED_FROM_X_ALONE'
  previous=initial(job['omega'],job['condition'],job['seed'],reference);fixed=[np.array([np.dot(previous[s]/job['omega']-np.mean(previous[s]/job['omega']),unit[s]) for s in range(3)])];paths=[]
  prefix=continuation['prefixes'].get(job['id'],[]) if phase=='continuation' else []
  for chunk in range(100):
   path=ROOT/prefix[chunk]['path'] if chunk<len(prefix) else directory/f'chunk_{chunk:03d}.npz';digest=hashlib.sha256(path.read_bytes()).hexdigest()
   if chunk<len(prefix):assert digest==prefix[chunk]['sha256']
   with np.load(path) as a:
    assert np.array_equal(previous,a['initial']) and int(a['seed'])==job['segment_seeds'][chunk] and int(a['clips'])==0
    counts=a['counts'];assert counts.shape==(20,3,500) and counts.min()>=0 and np.array_equal(a['final'],counts[-1]);previous=a['final'].copy();logical=chunk*10+np.arange(1,21)*.5
    np.testing.assert_allclose(a['times'],logical,atol=1e-7,rtol=0);np.testing.assert_allclose(a['kernel_times'],logical,atol=1e-7,rtol=0)
    # Separate implementation: direct dot product, no adapter einsum/summarize.
    for row in counts:
     fixed.append(np.array([np.dot(row[s]/job['omega']-np.mean(row[s]/job['omega']),unit[s]) for s in range(3)]))
   paths.append(dict(path=str(path.relative_to(ROOT)),sha256=digest));segments+=1
  fixed=np.array(fixed);times=np.arange(2001)*.5;np.testing.assert_array_equal(measure['times'],times);error=float(np.max(abs(fixed-measure['fixed_projection'])));maxerror=max(maxerror,error);np.testing.assert_allclose(fixed,measure['fixed_projection'],atol=2e-10,rtol=1e-10)
  legacy=summary['legacy'];blocks=[];valid=measure['unique_original'];aligned=measure['aligned_projection'];assert np.isnan(aligned[~valid]).all()
  if job['condition']=='prepared':
   directvalid=np.isfinite(np.array(legacy['positions']))&(np.array(legacy['n_dets'])==1);assert np.array_equal(valid,directvalid);assert bool(valid.all())==summary['strict_sampled_unique_X_track'];assert bool(legacy['eligible'])==summary['historical_X_persistence']
   # Real-space periodic interpolation via an independently constructed Fourier phase.
   for frame in (0,1000,2000):
    if not valid[frame]:continue
    chunk=max(0,(frame-1)//20);path=ROOT/prefix[chunk]['path'] if chunk<len(prefix) else directory/f'chunk_{chunk:03d}.npz'
    row=initial(job['omega'],job['condition'],job['seed'],reference) if frame==0 else np.load(path)['counts'][(frame-1)%20]
    centered=row/job['omega']-(row/job['omega']).mean(axis=1,keepdims=True);position=measure['positions'][frame];shifted=np.fft.ifft(np.fft.fft(centered,axis=1)*np.exp(2j*np.pi*np.fft.fftfreq(500,d=.2)*position),axis=1).real
    np.testing.assert_allclose(np.sum(shifted*unit,axis=1),aligned[frame],atol=2e-10,rtol=1e-9);checks+=1
  for block in range(10):
   mask=(times>=block*100)&((times<(block+1)*100) if block<9 else times<=1000);keep=mask&valid;fm=fixed[mask].mean(axis=0);am=aligned[keep].mean(axis=0) if keep.any() else None
   np.testing.assert_allclose(fm,summary['blocks'][block]['fixed_mean'],atol=2e-10);assert int(keep.sum())==summary['blocks'][block]['unique_original_snapshots']
   if am is not None:np.testing.assert_allclose(am,summary['blocks'][block]['aligned_mean'],atol=2e-10)
   blocks.append(dict(block=block,fixed=fm.tolist(),aligned=None if am is None else am.tolist(),relative_aligned=None if am is None else (am/norm).tolist(),valid_frames=int(keep.sum()),total_frames=int(mask.sum())))
  events=legacy.get('events',[]);births=legacy.get('births',[]);exports.append(dict(id=job['id'],omega=job['omega'],condition=job['condition'],phase=phase,clips=0,X_persistent=summary['historical_X_persistence'],strict=summary['strict_sampled_unique_X_track'],terminal_ambiguous=summary['end_unmatched_or_ambiguous'],events=events,births=births,blocks=blocks,final_fixed=fixed[-1].tolist(),final_aligned=None if not valid[-1] else aligned[-1].tolist(),initial_fixed=fixed[0].tolist()))
  manifests.append(dict(id=job['id'],summary_sha256=hashlib.sha256((directory/'summary.json').read_bytes()).hexdigest(),measurements_sha256=hashlib.sha256((directory/'GXY_measurements.npz').read_bytes()).hexdigest(),segments=paths))
  if index%20==0:print('verified',index+1,'/180',flush=True)
 population=[];series=[]
 for omega in (100,300,1000,3000,10000):
  prep=[e for e in exports if e['omega']==omega and e['condition']=='prepared'];controls=[e for e in exports if e['omega']==omega and e['condition']=='homogeneous'];k=sum(e['X_persistent'] for e in prep)
  population.append(dict(omega=omega,n=30,X_survived=k,CP95=cp(k,30),strict=sum(e['strict'] for e in prep),terminal_ambiguous=sum(e['terminal_ambiguous'] for e in prep),control_birth_runs=sum(bool(e['births']) for e in controls),control_n=6,loss_runs=sum(any(v[0]=='loss' for v in e['events']) for e in prep),replication_runs=sum(any(v[0]=='replication' for v in e['events']) for e in prep)))
  for block in range(10):
   for species in range(3):
    a=[e['blocks'][block]['fixed'][species] for e in prep];b=[e['blocks'][block]['fixed'][species] for e in controls];va=np.var(a,ddof=1)/30;vb=np.var(b,ddof=1)/6;df=(va+vb)**2/(va**2/29+vb**2/5);difference=float(np.mean(a)-np.mean(b));half=float(t.ppf(.975,df)*np.sqrt(va+vb));cond=[e['blocks'][block]['aligned'][species] for e in prep if e['blocks'][block]['aligned'] is not None]
    series.append(dict(omega=omega,block=block,species=('G','X','Y')[species],prepared_fixed=interval(a),control_fixed=interval(b),fixed_difference=dict(mean=difference,CI95=[difference-half,difference+half]),prepared_aligned=interval(cond),aligned_valid_frames=[e['blocks'][block]['valid_frames'] for e in prep],conditional_alignment=True))
 # Independently reconstruct the numerical qualification Welch comparisons from saved projections.
 qualification=[]
 for condition in ('prepared','homogeneous'):
  groups=[]
  for tau in (plan['jobs'][0]['tau'],plan['jobs'][24]['tau']):
   values=[]
   for record in qual['records']:
    job=record['job']
    if job['condition']==condition and job['tau']==tau:
     a=np.load(RAW/'qualification'/job['id']/'GXY_measurements.npz');values.append(float(a['fixed_projection'][(a['times']>=50)&(a['times']<=100),0].mean()));np.testing.assert_allclose(values[-1],record['numeric_G_block_mean'])
   assert len(values)==12;groups.append(np.array(values))
  a,b=groups;va=np.var(a,ddof=1)/12;vb=np.var(b,ddof=1)/12;df=(va+vb)**2/(va**2/11+vb**2/11);diff=float(b.mean()-a.mean());half=float(t.ppf(.975,df)*np.sqrt(va+vb));stored=next(x for x in read(PUBLIC/'qualification.json')['checks'] if x['condition']==condition);np.testing.assert_allclose(stored['CI95'],[diff-half,diff+half]);assert -norm[0]*.05<=diff-half and diff+half<=norm[0]*.05;qualification.append(stored)
 seconds=qual['charged_worker_seconds']+resumed['charged_worker_seconds'];assert seconds<360000;size=sum(p.stat().st_size for p in RAW.rglob('*') if p.is_file());assert size<10000000000
 atomic_json(PUBLIC/'decision_exports.json',dict(runs=exports,population=population,series=series,reference_norms=norm.tolist(),qualification=qualification,charged_worker_hours=seconds/3600,raw_bytes=size,full_GXY_persistence_claim=False))
 # Full segment hashes remain local; compact public manifest binds that file.
 manifest=RAW/'analysis_segment_manifest.json';atomic_json(manifest,dict(runs=manifests));atomic_json(PUBLIC/'verification.json',dict(passed=True,jobs=180,segments=segments,raw_snapshots=180*2001,maximum_fixed_projection_error=maxerror,aligned_sentinel_checks=checks,qualification_recomputed=True,manifest_sha256=hashlib.sha256(manifest.read_bytes()).hexdigest(),raw_bytes_before_analysis_manifest=size,elapsed_seconds=time.perf_counter()-start,external_replication=False,shared_dependencies=['NumPy/SciPy','reference','legacy X detector outputs'],limits='G projections and interval calculations independently recomputed; X tracking classifications checked against saved legacy outputs, not an independent detector. No new simulations.'))
 print(json.dumps(population));print('worker hours',seconds/3600,'raw GB',size/1e9)
if __name__=='__main__':run()

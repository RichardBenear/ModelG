import os
for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS'):os.environ[k]='1'
import time,hashlib,gzip,json,math,traceback
import numpy as np
from scipy.special import logsumexp,expit
from scipy.stats import beta,f
from run_mg176 import ROOT,OUT as PRIOR,RAW,read,atomic_json,models
from mg176_acquisition import tracker
from check_mg174 import inner,mu
OUT=ROOT/'results/MG12H3/split_reference_audit'
FILES=['src/audit_mg177.py','spec/MG177_split_reference_audit.md','src/run_mg176.py','src/mg176_acquisition.py','src/mg174_modes.py','src/mg172_corrected.py','src/analyze_mg170.py','src/check_mg174.py','src/mg167_xy.py','src/mg161_local_shape.py']
def register():
 OUT.mkdir(exist_ok=True);assert not(OUT/'plan.json').exists();paths=sorted(PRIOR.glob('r*.json.gz'))+sorted(RAW.glob('r*.npz'))+[PRIOR/'plan.json',PRIOR/'analysis.json'];atomic_json(OUT/'plan.json',dict(files={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths+[ROOT/x for x in FILES]},cap_seconds=600,no_draws=True,workers=1))
def directq(e,P):
 r=np.fft.irfft(e,n=500,axis=-1);v=np.fft.irfft(np.einsum('fab,...bf->...af',P,e),n=500,axis=-1);return np.sum(r*v,axis=(-1,-2))
def run():
 assert not(OUT/'attempt.json').exists();start=time.perf_counter();atomic_json(OUT/'attempt.json',dict(status='RUNNING'))
 try:
  assert all(hashlib.sha256((ROOT/k).read_bytes()).hexdigest()==v for k,v in read(OUT/'plan.json')['files'].items());cache={ri:(m,L) for ri,m,L in models()};refs=[];results=[];labels=gridchecks=0;maxerror=0.
  for file in sorted(PRIOR.glob('r*.json.gz')):
   d=json.loads(gzip.decompress(file.read_bytes()));ri=0 if d['regime']==.1 else 2;model,L=cache[ri]
   from mg174_modes import ModeTracker
   base=ModeTracker(model,L,d['rho'],100);P=base.P;w=base.freqweights;ft=base.ft;ev,V=np.linalg.eigh(P);sqrtP=np.einsum('fij,fj,fkj->fik',V,np.sqrt(np.maximum(ev,0)),np.conj(V))
   with np.load(RAW/(file.name[:-8]+'.npz')) as z:fields=z['fields'];calfields=z['calibration']
   for path,seq,calfield in zip(d['paths'],fields,calfields):
    y=np.fft.rfft(calfield-calfield.mean(axis=2)[:,:,None],axis=2);q=np.sum(directq(y,P));np.testing.assert_allclose(q,sum(path['calibration']['q']),rtol=1e-8,atol=1e-6);white=np.einsum('fab,rbf->raf',sqrtP,y);energies=np.sum(w[None,None,:]*abs(white)**2,axis=(0,2))/500;np.testing.assert_allclose(sum(energies),q,rtol=1e-8,atol=1e-6);fraction=float(energies[0]/sum(energies));pc=float(min(1.,2*min(beta.cdf(fraction,3992,3992),beta.sf(fraction,3992,3992))));mean=y.mean(axis=0);qm=float(16*directq(mean,P));qc=float(np.sum(directq(y-mean,P)));F=15*qm/qc;pm=float(f.sf(F,998,14970));refs.append(dict(regime=d['regime'],rho=d['rho'],scale=d['noise_amplitude'],case=d['case'],primary=d['primary'],replicate=path['replicate'],original_adequate=path['calibration']['adequate'],channel_fraction=fraction,channel_p=pc,mean_F=F,mean_p=pm,channel_tail=pc<.001/2,mean_tail=pm<.001/2))
    if ri!=2 or d['rho']!=0 or d['case'] not in ('split','moving'):continue
    obj=tracker(model,L,0,path['seed'],path['calibration']);scale=path['calibration']['scale'];Pe=P/scale;norms=np.array([inner(t,t,Pe,w) for t in ft]);K=np.einsum('fab,kbf->kaf',Pe,ft);offset=np.arange(500);distance=np.minimum(offset,500-offset)*.2;eligible=np.broadcast_to(distance[None,:]>6.4,(13,500));N=int(eligible.sum());assert N==5655;redsingle=redmulti=-math.log(2);perssingle=-math.log(2);persgrid=np.full((13,500),-np.inf);persgrid[eligible]=-math.log(2)-math.log(N);rows=[];previous={k:.2 for k in ('absent','single','Xonly','Yonly','multiple')}
    for t,(field,truth,stored) in enumerate(zip(seq,path['truth'],path['records'])):
     r=obj.step(field);np.testing.assert_allclose(list(r['class_weights'].values()),list(stored['class_weights'].values()),rtol=1e-9,atol=1e-10);assert r['label']==stored['label'];labels+=1;pred={k:.98*previous[k]+.005*(1-previous[k]) for k in previous};prior=math.log(pred['multiple']/pred['single']);post=math.log(r['class_weights']['multiple']/r['class_weights']['single']);previous=r['class_weights'];k,ss=truth;si=round((ss-.5)/.125);yf=np.fft.rfft(field-field.mean(axis=1)[:,None],axis=1);e=yf-mu('single',k,si,None,ft);ls=-.5*inner(e,e,Pe,w);corr=np.fft.irfft((e[None,:,:]*np.conj(K)).sum(axis=1),n=500,axis=1);grid=np.roll(corr-.5*norms[:,None],-k,axis=1);grid=np.where(eligible,grid,-np.inf);bf=float(logsumexp(grid)-math.log(N));known=float(grid[4,250]);truebf=known*scale/(d['noise_amplitude']**2)
     rpS=float(np.logaddexp(redsingle+math.log(.98),redmulti+math.log(.005)));rpM=float(np.logaddexp(redmulti+math.log(.98),redsingle+math.log(.005)));redsingle=rpS;redmulti=rpM+bf;z=float(np.logaddexp(redsingle,redmulti));redsingle-=z;redmulti-=z
     ps=float(np.logaddexp(perssingle+math.log(.98),logsumexp(persgrid)+math.log(.005)));pg=np.logaddexp(persgrid+math.log(.98),perssingle+math.log(.005)-math.log(N));pg=np.where(eligible,pg+grid,-np.inf);z=float(np.logaddexp(ps,logsumexp(pg)));perssingle=ps-z;persgrid=pg-z
     for off,sj in [(100,0),(250,4),(400,12)]:
      sec=mu('single',(k+off)%500,sj,None,ft);direct=float(-.5*(directq(e-sec,Pe)-directq(e,Pe)));np.testing.assert_allclose(direct,grid[sj,off],rtol=1e-8,atol=2e-6);maxerror=max(maxerror,abs(direct-grid[sj,off]));gridchecks+=1
     rows.append(dict(frame=t,label=r['label'],prior_log_multiple_single=prior,posterior_log_multiple_single=post,integrated_log_evidence_ratio=post-prior,oracle_primary_uniform_logBF=bf,oracle_pair_logBF=known,oracle_pair_true_scale_logBF=truebf,oracle_primary_redrawn_multiple_weight=float(np.exp(redmulti)),oracle_primary_persistent_multiple_weight=float(np.exp(logsumexp(persgrid))),selected_center=r['center'],truth_center=k,geometry_fraction=r['retained_pair_geometry_fraction'],class_weights=r['class_weights']))
    results.append(dict(scale=d['noise_amplitude'],case=d['case'],replicate=path['replicate'],estimated_variance=scale,rows=rows))
   assert time.perf_counter()-start<600
  assert labels==288 and len(refs)==274;atomic_json(OUT/'audit.json',dict(reference_diagnostics=refs,conditional_paths=results,unchanged_labels=labels,direct_grid_checks=gridchecks,max_grid_difference=maxerror,real_space_reference_checks=len(refs),elapsed_seconds=time.perf_counter()-start,conditional_truth_primary=True,no_new_verdict=True));atomic_json(OUT/'attempt.json',dict(status='COMPLETE',elapsed_seconds=time.perf_counter()-start))
 except BaseException:atomic_json(OUT/'attempt.json',dict(status='FAILED',error=traceback.format_exc(),elapsed_seconds=time.perf_counter()-start));raise
if __name__=='__main__':register() if __import__('sys').argv[1]=='register' else run()

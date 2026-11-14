"""Separate real-space calibration, label, acquisition and outcome recomputation."""
import time,hashlib,gzip,json,math,collections,traceback
import numpy as np
from scipy.stats import beta,chi2,f
from run_mg176 import ROOT,OUT,RAW,binding,models,read,atomic_json
from mg174_modes import ModeTracker
from check_mg174 import mu

def check():
 start=time.perf_counter();assert not(OUT/'check_attempt.json').exists();atomic_json(OUT/'check_attempt.json',dict(status='RUNNING'));labels=paths=calchecks=bchecks=qchecks=0;stats=[];maxq=0.
 try:
  assert read(OUT/'attempt.json')['status']=='COMPLETE' and binding()==read(OUT/'plan.json')['source'];allmodels={ri:(model,L) for ri,model,L in models()};rows=read(OUT/'analysis.json')['rows']
  for p in sorted(OUT.glob('r*_rho*.json.gz')):
   d=json.loads(gzip.decompress(p.read_bytes()));ri=0 if d['regime']==.1 else 2;model,L=allmodels[ri];base=ModeTracker(model,L,d['rho'],100);P=base.P;w=base.freqweights;ft=base.ft
   with np.load(RAW/(p.name[:-8]+'.npz')) as z:fields=z['fields'];cals=z['calibration']
   success=adequate_count=acquired=0
   for path,seq,reference in zip(d['paths'],fields,cals):
    ys=np.fft.rfft(reference-reference.mean(axis=2)[:,:,None],axis=2);real=np.fft.irfft(ys,n=500,axis=2);filtered=np.fft.irfft(np.einsum('fab,rbf->raf',P,ys),n=500,axis=2);q=np.sum(real*filtered,axis=(1,2));cal=path['calibration'];np.testing.assert_allclose(q,cal['q'],rtol=1e-8,atol=1e-6);calchecks+=16;total=math.fsum(q);scale=total/15968;np.testing.assert_allclose(scale,cal['scale'],rtol=1e-8,atol=1e-8);np.testing.assert_allclose(total/chi2.ppf([.975,.025],15968),cal['ci95'],rtol=1e-8);adequate=True
    for ix,stored in zip(np.array_split(np.arange(1,251),8),cal['bands']):
     mask=np.zeros(251);mask[ix]=1;filtered=np.fft.irfft(np.einsum('fab,rbf->raf',P,ys)*mask,n=500,axis=2);qb=float(np.sum(real*filtered));db=int(sum(2 if j==250 else 4 for j in ix)*16);lo,hi=beta.ppf([.001/16,1-.001/16],db/2,(15968-db)/2);okay=bool(lo<=qb/total<=hi);np.testing.assert_allclose(qb,stored['q'],rtol=1e-8,atol=1e-6);assert db==stored['df'] and okay==stored['passed'];adequate &= okay;bchecks+=1
    assert adequate==cal['adequate'];adequate_count+=adequate;records=path['records'];ysfield=np.fft.rfft(seq-seq.mean(axis=2)[:,:,None],axis=2);old=ModeTracker(model,L,d['rho'],path['seed']).beam
    for t,r in enumerate(records):
     pw=np.exp(np.array(r['full_pair_logmass'])-np.logaddexp.reduce(r['full_pair_logmass']));np.testing.assert_allclose(pw,[*r['pair_weights'].values()],rtol=1e-8,atol=1e-10);np.testing.assert_allclose(pw.reshape(5,2).sum(1),[*r['class_weights'].values()],rtol=1e-8,atol=1e-10)
     cutoff=998*f.ppf(.999,998,15968);np.testing.assert_allclose(cutoff,r['qcut']);winner=max(r['class_weights'],key=r['class_weights'].get);single=[b for b in r['beam'] if b['kind']=='single'];cs=max(min((x['center']-y['center'])%500,(y['center']-x['center'])%500) for x in single for y in single);ws=max(abs(x['width']-y['width']) for x in single for y in single);assert cs==r['single_center_spread'];np.testing.assert_allclose(ws,r['single_width_spread']);label=winner if r['class_weights'][winner]>=.95 and r['innovation_q']<=cutoff else 'uncertain'
     if label=='single' and(cs>4 or ws>.25):label='uncertain'
     assert r['label']==label and r['emitted']==(label if adequate else 'noise_unqualified');labels+=1
     if t in (0,10):
      best=max((b for b in r['beam'] if b['kind']==winner),key=lambda b:b['merged_log']);si=round((best['width']-.5)/.125);cur=mu(winner,best['center'],si,best['second'],ft);rho=d['rho'] if t else 0;innovation=ysfield[t]-rho*ysfield[t-1] if t else ysfield[t];incoming=[]
      for parent in old:
       ps=parent.get('si',round((parent.get('width',1)-.5)/.125));pk=parent['center'];pm=parent.get('mean')
       if pm is None:pm=mu(parent['kind'],pk,ps,parent['second'],ft)
       allowed=(best['center']==pk and si==ps) if winner=='absent' or best['motion']=='S' else min((best['center']-pk)%500,(pk-best['center'])%500)<=2 and abs(si-ps)<=1
       if not allowed:continue
       e=innovation-(cur-rho*pm);er=np.fft.irfft(e,n=500,axis=1);ef=np.fft.irfft(np.einsum('fab,bf->af',P,e),n=500,axis=1);qr=float(np.sum(er*ef)/cal['scale']/(1-rho*rho));score=parent.get('log',parent.get('log_weight'))+math.log(.98 if winner==parent['kind'] else .005)+math.log(.98 if best['motion']==parent['motion'] else .02)-.5*qr
       if winner!='absent' and best['motion']=='M':score-=math.log(5*(min(13,ps+2)-max(0,ps-1)))
       if winner=='multiple':score-=math.log(5655)
       incoming.append((score,qr))
      qr=max(incoming,key=lambda x:x[0])[1];np.testing.assert_allclose(qr,r['innovation_q'],rtol=1e-8,atol=1e-6);maxq=max(maxq,abs(qr-r['innovation_q']));qchecks+=1
     old=r['beam']
    times=[t for t in range(1,6) if records[t-1]['emitted']==records[t]['emitted']=='single'];acq=min(times) if times else None;assert acq==path['outcome']['acquisition'];acquired+=acq is not None;positive=d['case'] in ('static','moving','deforming','stop_start') or not d['primary']
    if positive:
     tail=range(acq+1,24) if acq is not None else [];ss=sum(records[t]['emitted']=='single' for t in tail);aa=sum(records[t]['emitted']=='single' and min((records[t]['center']-path['truth'][t][0])%500,(path['truth'][t][0]-records[t]['center'])%500)<=2 and abs(records[t]['width']-path['truth'][t][1])<=.12500001 for t in tail);okay=acq is not None and ss>=math.ceil(.9*len(tail)) and aa>=math.ceil(.9*len(tail));assert ss==path['outcome']['single'] and aa==path['outcome']['accurate']
    elif d['case']=='absent':okay=not any(records[t-1]['emitted']==records[t]['emitted']=='single' for t in range(1,24))
    else:okay=acq is not None and records[15]['emitted']!='single'
    assert okay==path['outcome']['success'];success+=okay;paths+=1;stats.append(dict(regime=d['regime'],rho=d['rho'],noise_amplitude=d['noise_amplitude'],case=d['case'],primary=d['primary'],replicate=path['replicate'],adequate=adequate,estimated_variance_scale=cal['scale'],true_variance_scale=d['noise_amplitude']**2,ci95=cal['ci95'],acquisition=acq,success=okay,labels=dict(collections.Counter(r['emitted'] for r in records))))
   row=next(r for r in rows if all(r[k]==d[k] for k in ('regime','rho','noise_amplitude','case','primary')));assert (success,adequate_count,acquired)==(row['success'],row['reference_adequate'],row['acquired']);assert time.perf_counter()-start<120
  assert paths==274 and labels==6576;atomic_json(OUT/'separate_check.json',dict(passed=True,labels=labels,paths=paths,real_space_reference_quadratics=calchecks,reference_band_checks=bchecks,selected_history_quadratics=qchecks,max_selected_quadratic_difference=maxq,path_summary=stats,elapsed_seconds=time.perf_counter()-start));atomic_json(OUT/'local_raw_manifest.json',dict(directory=str(RAW),files={p.name:dict(bytes=p.stat().st_size,sha256=hashlib.sha256(p.read_bytes()).hexdigest()) for p in RAW.iterdir()},total_bytes=sum(p.stat().st_size for p in RAW.iterdir()),full_fields_published=False));atomic_json(OUT/'check_attempt.json',dict(status='COMPLETE',elapsed_seconds=time.perf_counter()-start))
 except BaseException:atomic_json(OUT/'check_attempt.json',dict(status='FAILED',error=traceback.format_exc(),elapsed_seconds=time.perf_counter()-start));raise
if __name__=='__main__':check()

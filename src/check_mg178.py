import os
for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS'):os.environ[k]='1'
import time,hashlib,gzip,json,math,traceback,collections
import numpy as np
from scipy.stats import beta,f,norm,chi2
from run_mg178 import ROOT,OUT,RAW,read,atomic_json,models,binding
from mg174_modes import ModeTracker
from check_mg174 import mu

def rq(a,b,P):
 ar=np.fft.irfft(a,n=500,axis=-1);bf=np.fft.irfft(np.einsum('fab,...bf->...af',P,b),n=500,axis=-1);return np.sum(ar*bf,axis=(-1,-2))
def check():
 start=time.perf_counter();assert not(OUT/'check_attempt.json').exists();atomic_json(OUT/'check_attempt.json',dict(status='RUNNING'));labels=paths=qchecks=calqs=masschecks=0;summary=[];maxerr=0.;allmodels={ri:(m,L) for ri,m,L in models()}
 try:
  assert read(OUT/'attempt.json')['status']=='COMPLETE' and binding()==read(OUT/'plan.json')['source'];rows=read(OUT/'analysis.json')['rows']
  for file in sorted(OUT.glob('r*.json.gz')):
   d=json.loads(gzip.decompress(file.read_bytes()));ri=0 if d['regime']==.1 else 2;model,L=allmodels[ri];base=ModeTracker(model,L,0,100);P=base.P;w=base.freqweights;ft=base.ft;ev,V=np.linalg.eigh(P);W=np.einsum('fij,fj,fkj->fik',V,np.sqrt(np.maximum(ev,0)),np.conj(V));succeeded=adequate_count=acquired=0
   with np.load(RAW/(file.name[:-8]+'.npz')) as z:fields=z['fields'];spatials=z['spatial'];temporals=z['temporal']
   for path,seq,spatial,temporal in zip(d['paths'],fields,spatials,temporals):
    cal=path['calibration'];ys=np.fft.rfft(spatial-spatial.mean(axis=2)[:,:,None],axis=2);yt=np.fft.rfft(temporal-temporal.mean(axis=2)[:,:,None],axis=2);qs=rq(ys,ys,P);total=float(qs.sum());scale=total/15968;np.testing.assert_allclose(qs,cal['q'],rtol=1e-8,atol=2e-6);np.testing.assert_allclose(scale,cal['scale'],rtol=1e-8);np.testing.assert_allclose(total/chi2.ppf([.975,.025],15968),cal['ci95'],rtol=1e-8);rhat=float(rq(yt[1:],yt[:-1],P).sum()/rq(yt[:-1],yt[:-1],P).sum());np.testing.assert_allclose(rhat,cal['raw_rho'],rtol=1e-7,atol=1e-8);rho=float(np.clip(rhat,-.99,.99));innov=(yt[1:]-rho*yt[:-1])/np.sqrt(1-rho*rho);qt=float(rq(innov,innov,P).sum());df=63*998-1;se=float(np.sqrt(qt*(1-rho*rho)/df/rq(yt[:-1],yt[:-1],P).sum()));np.testing.assert_allclose([rhat-1.96*se,rhat+1.96*se],cal['rho_ci95'],rtol=1e-7,atol=1e-8);calqs+=79;alpha=.001/12;true_reference_rho=.8 if d['case']=='actual_rho' else d['rho'];true_innov=(yt[1:]-true_reference_rho*yt[:-1])/np.sqrt(1-true_reference_rho**2);true_ratio=float(rq(true_innov,true_innov,P).sum()/(63*998)/scale);true_pv=float(min(1.,2*min(f.cdf(true_ratio,63*998,15968),f.sf(true_ratio,63*998,15968))));bands=[]
    for ix,b in zip(np.array_split(np.arange(1,251),8),cal['bands']):
     mask=np.zeros(251);mask[ix]=1;qb=float(rq(ys,ys*mask,P).sum());db=int(sum(2 if j==250 else 4 for j in ix)*16);lo,hi=beta.ppf([alpha/2,1-alpha/2],db/2,(15968-db)/2);bands.append(bool(lo<=qb/total<=hi));np.testing.assert_allclose(qb,b['q'],rtol=1e-8,atol=2e-6)
    z=np.einsum('fab,rbf->raf',W,ys);zr=np.fft.irfft(z,n=500,axis=2);eng=(zr*zr).sum(axis=(0,2));frac=float(eng[0]/eng.sum());pc=float(min(1.,2*min(beta.cdf(frac,3992,3992),beta.sf(frac,3992,3992))));mean=ys.mean(axis=0);qm=float(16*rq(mean,mean,P));qc=float(rq(ys-mean,ys-mean,P).sum());pm=float(f.sf(15*qm/qc,998,14970));ratio=qt/df/scale;pv=float(min(1.,2*min(f.cdf(ratio,df,15968),f.sf(ratio,df,15968))));lag=float(rq(innov[2:],innov[:-2],P).sum()/np.sqrt(rq(innov[2:],innov[2:],P).sum()*rq(innov[:-2],innov[:-2],P).sum()));pl=float(2*norm.sf(abs(lag)*np.sqrt(61*998)));np.testing.assert_allclose([pc,pm,pv,pl],[cal['checks'][k] for k in ['channel_p','mean_p','temporal_variance_p','lag2_p']],rtol=1e-4,atol=2e-7);adequate=abs(rhat)<.99 and all(bands) and min(pc,pm,pv,pl)>=alpha;assert adequate==cal['adequate'];adequate_count+=adequate;records=path['records'];yf=np.fft.rfft(seq-seq.mean(axis=2)[:,:,None],axis=2);init=ModeTracker(model,L,0,path['seed']).beam;old=init;resets=0
    for t,r in enumerate(records):
     logs=np.array(r['full_pair_logmass']);pw=np.exp(logs-np.logaddexp.reduce(logs));np.testing.assert_allclose(pw,list(r['pair_weights'].values()),rtol=1e-8,atol=1e-12);np.testing.assert_allclose(pw.reshape(5,2).sum(1),list(r['class_weights'].values()),rtol=1e-8,atol=1e-12)
     for pair,pweight in zip(r['pair_weights'],pw):
      kind,motion=pair.split('/');mass=math.fsum(math.exp(b['log_weight']) for b in r['beam'] if b['kind']==kind and b['motion']==motion);np.testing.assert_allclose(mass,pweight,rtol=1e-8,atol=1e-12);masschecks+=1
     winner=max(r['class_weights'],key=r['class_weights'].get);single=[b for b in r['beam'] if b['kind']=='single'];cs=max(min((x['center']-y['center'])%500,(y['center']-x['center'])%500) for x in single for y in single);ws=max(abs(x['width']-y['width']) for x in single for y in single);cut=998*f.ppf(.999,998,15968);np.testing.assert_allclose(cut,r['qcut']);assert cs==r['single_center_spread'];np.testing.assert_allclose(ws,r['single_width_spread']);lab=winner if r['class_weights'][winner]>=.95 and r['innovation_q']<=cut else 'uncertain'
     if lab=='single' and(cs>4 or ws>.25):lab='uncertain'
     assert r['label']==lab and r['emitted']==(lab if adequate else 'noise_unqualified');labels+=1;resets+=r['secondary_kernel_resets']
     if t in (0,10):
      best=max((b for b in r['beam'] if b['kind']==winner),key=lambda b:b['merged_log']);si=round((best['width']-.5)/.125);current=mu(winner,best['center'],si,best['second'],ft);rr=rho if t else 0.;innovation=yf[t]-rr*yf[t-1] if t else yf[t];incoming=[];distance=np.minimum((np.arange(500)-best['center'])%500,(best['center']-np.arange(500))%500)*.2;eligible=np.broadcast_to(distance[None,:]>6.4,(13,500))
      for p in old:
       ps=p.get('si',round((p.get('width',1)-.5)/.125));pk=p['center'];allowed=(best['center']==pk and si==ps) if winner=='absent' or best['motion']=='S' else min((best['center']-pk)%500,(pk-best['center'])%500)<=2 and abs(si-ps)<=1
       if not allowed:continue
       lp=0.
       if winner=='multiple':
        lp=-math.log(5655)
        if p['kind']=='multiple' and p['second'] is not None:
         sk,ss=p['second'];options=[(sj,kj%500) for sj in range(max(0,ss-1),min(13,ss+2)) for kj in range(sk-2,sk+3) if eligible[sj,kj%500]]
         if options:
          pr={q:.5/len(options) for q in options}
          if eligible[ss,sk]:pr[(ss,sk)]=pr.get((ss,sk),0)+.5
          dest=(best['second'][1],best['second'][0])
          if dest not in pr:continue
          lp=math.log(pr[dest]/sum(pr.values()))
       pmu=p.get('mean')
       if pmu is None:pmu=mu(p['kind'],pk,ps,p['second'],ft)
       e=innovation-(current-rr*pmu);q=float(rq(e,e,P)/cal['scale']/(1-rr*rr));score=p.get('log',p.get('log_weight'))+math.log(.98 if winner==p['kind'] else .005)+math.log(.98 if best['motion']==p['motion'] else .02)+lp-.5*q
       if winner!='absent' and best['motion']=='M':score-=math.log(5*(min(13,ps+2)-max(0,ps-1)))
       incoming.append((score,q))
      q=max(incoming,key=lambda v:v[0])[1];np.testing.assert_allclose(q,r['innovation_q'],rtol=1e-8,atol=3e-6);maxerr=max(maxerr,abs(q-r['innovation_q']));qchecks+=1
     old=r['beam']
    case=d['case'];target='multiple' if case=='secondary_initial' else 'single';times=[t for t in range(1,6) if records[t-1]['emitted']==records[t]['emitted']==target];acq=min(times) if times else None;assert acq==path['outcome']['acquisition'];acquired+=acq is not None
    if case.startswith('secondary_'):
     checks={12:'multiple',20:'single'} if case=='secondary_vanish' else {15:'multiple'};okay=acq is not None and all(records[t]['emitted']==v for t,v in checks.items())
    elif case in ['static','moving','deforming','stop_start'] or not d['primary']:
     tail=range(acq+1,24) if acq is not None else [];ss=sum(records[t]['emitted']=='single' for t in tail);aa=sum(records[t]['emitted']=='single' and min((records[t]['center']-path['truth'][t][0])%500,(path['truth'][t][0]-records[t]['center'])%500)<=2 and abs(records[t]['width']-path['truth'][t][1])<=.12500001 for t in tail);okay=acq is not None and ss>=math.ceil(.9*len(tail)) and aa>=math.ceil(.9*len(tail));assert ss==path['outcome']['single'] and aa==path['outcome']['accurate']
    elif case=='absent':okay=not any(records[t-1]['emitted']==records[t]['emitted']=='single' for t in range(1,24))
    else:okay=acq is not None and records[15]['emitted']!='single'
    assert okay==path['outcome']['success'];succeeded+=okay;paths+=1;summary.append(dict(regime=d['regime'],rho=d['rho'],scale=d['noise_amplitude'],case=case,primary=d['primary'],replicate=path['replicate'],adequate=bool(adequate),rho_estimate=cal['rho'],rho_ci95=cal['rho_ci95'],scale_estimate=cal['scale'],scale_ci95=cal['ci95'],checks=cal['checks'],true_reference_rho_variance_p=true_pv,acquisition=acq,success=bool(okay),labels=dict(collections.Counter(r['emitted'] for r in records)),deadline15=records[15]['emitted'],secondary_resets=resets))
   row=next(r for r in rows if all(r[k]==d[k] for k in ['regime','rho','noise_amplitude','case','primary']));assert (succeeded,adequate_count,acquired)==(row['success'],row['reference_adequate'],row['acquired']);assert time.perf_counter()-start<600
  assert paths==292 and labels==7008;atomic_json(OUT/'separate_check.json',dict(passed=True,labels=labels,path_outcomes=paths,class_mode_mass_checks=masschecks,real_space_reference_quadratics=calqs,selected_history_quadratics=qchecks,max_selected_quadratic_difference=maxerr,path_summary=summary,elapsed_seconds=time.perf_counter()-start));atomic_json(OUT/'local_raw_manifest.json',dict(directory=str(RAW),files={p.name:dict(bytes=p.stat().st_size,sha256=hashlib.sha256(p.read_bytes()).hexdigest()) for p in RAW.iterdir()},total_bytes=sum(p.stat().st_size for p in RAW.iterdir()),full_fields_published=False));atomic_json(OUT/'check_attempt.json',dict(status='COMPLETE',elapsed_seconds=time.perf_counter()-start))
 except BaseException:atomic_json(OUT/'check_attempt.json',dict(status='FAILED',error=traceback.format_exc(),elapsed_seconds=time.perf_counter()-start));raise
if __name__=='__main__':check()

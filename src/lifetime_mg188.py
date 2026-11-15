"""Frozen MG188 survival analysis. Deterministic fits, no new simulations."""
import os
for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS'):os.environ[key]='1'
import json,hashlib,csv
from pathlib import Path
import numpy as np
from scipy.optimize import minimize,minimize_scalar,brentq
from scipy.special import logsumexp
from scipy.stats import chi2,norm,beta as betadist
from run_mg185_phase2_G import ROOT,PUBLIC,RAW,read
from mg10e2c1_pool import atomic_json
OUT=ROOT/'results/MG188_lifetime';HORIZON=980.

def exponential(times,event,x):
 D=event.sum();sx=x[event].sum()
 def profile(b):
  alpha=np.log(D)-logsumexp(np.log(times)+b*x)
  return float(D*alpha+b*sx-D)
 fit=minimize_scalar(lambda b:-profile(b),bounds=(-50,50),method='bounded',options={'xatol':1e-12});b=float(fit.x);alpha=float(np.log(D)-logsumexp(np.log(times)+b*x));maximum=profile(b);cut=chi2.ppf(.95,1)/2
 def root(direction):
  width=.1
  while maximum-profile(b+direction*width)<cut:width*=2
  return float(brentq(lambda v:maximum-profile(v)-cut,*sorted((b,b+direction*width))))
 ci=[root(-1),root(1)];general=minimize(lambda p:float(np.sum(np.exp(p[0]+p[1]*x)*times)-np.sum((p[0]+p[1]*x)[event])),[alpha,b],method='BFGS');np.testing.assert_allclose(-general.fun,maximum,atol=1e-7)
 return dict(alpha=alpha,beta=b,beta_CI95=ci,loglik=maximum,parameters=2)

def weibull(times,event,x,start):
 logt=np.log(times)
 def objective(p):
  k=np.exp(p[2]);eta=p[0]+p[1]*x;z=k*(logt-eta)
  if np.max(z)>600:return 1e260
  return float(np.sum(np.exp(z))-np.sum((np.log(k)-k*eta+(k-1)*logt)[event]))
 fit=minimize(objective,[-start['alpha'],-start['beta'],0.],method='Nelder-Mead',options={'maxiter':10000,'xatol':1e-9,'fatol':1e-9});assert fit.success
 p=fit.x;step=1e-4;eye=np.eye(3)*step;hess=np.empty((3,3))
 for i in range(3):
  for j in range(3):hess[i,j]=(objective(p+eye[i]+eye[j])-objective(p+eye[i]-eye[j])-objective(p-eye[i]+eye[j])+objective(p-eye[i]-eye[j]))/(4*step**2)
 covariance=np.linalg.inv(hess);assert np.all(np.linalg.eigvalsh(hess)>0);se=np.sqrt(covariance[2,2]);shape=np.exp(p[2]);return dict(alpha=float(p[0]),beta=float(p[1]),shape=float(shape),shape_CI95=np.exp([p[2]-1.96*se,p[2]+1.96*se]).tolist(),loglik=float(-fit.fun),parameters=3)

def km(times,event):
 points=[dict(time=0.,survival=1.,lower=1.,upper=1.,at_risk=len(times),deaths=0)];s=1.;g=0.;last=0.;area=0.;median=None
 for time in sorted(set(times)):
  risk=int(np.sum(times>=time));deaths=int(np.sum((times==time)&event));area+=s*(time-last)
  if deaths:
   s*=1-deaths/risk
   if deaths<risk:g+=deaths/(risk*(risk-deaths))
   if median is None and s<=.5:median=float(time)
  if s==0:lo=hi=0.
  elif s==1:lo=hi=1.
  else:
   c=np.log(-np.log(s));se=np.sqrt(g)/abs(np.log(s));lo=float(np.exp(-np.exp(c+1.96*se)));hi=float(np.exp(-np.exp(c-1.96*se)))
  points.append(dict(time=float(time),survival=float(s),lower=lo,upper=hi,at_risk=risk,deaths=deaths));last=time
 area+=s*(HORIZON-last);np.testing.assert_allclose(area,times.mean(),atol=1e-9)
 return points,float(area),median

def run():
 assert not OUT.exists();OUT.mkdir();atomic_json(OUT/'attempt.json',dict(status='RUNNING',new_simulations=False))
 try:
  exports=read(PUBLIC/'decision_exports.json');rows=[]
  for r in exports['runs']:
   if r['condition']!='prepared':continue
   summary=read(RAW/r['phase']/r['id']/'summary.json');assert summary['legacy']['events']==r['events'];events=r['events'];loss=[float(v[1]) for v in events if v[0]=='loss' and v[1]<=HORIZON];solitary=[float(v[1]) for v in events if v[0] in ('loss','replication','merge') and v[1]<=HORIZON]
   rows.append(dict(id=r['id'],omega=r['omega'],time=min(loss) if loss else HORIZON,event=bool(loss),secondary_time=min(solitary) if solitary else HORIZON,secondary_event=bool(solitary),replication_or_merge=any(v[0] in ('replication','merge') and v[1]<=HORIZON for v in events)))
  assert len(rows)==150;times=np.array([r['time'] for r in rows]);assert times.min()>0;event=np.array([r['event'] for r in rows]);omega=np.array([r['omega'] for r in rows]);x=(omega-1000)/1000;primary=exponential(times,event,x);power=exponential(times,event,np.log(omega/1000));weib=weibull(times,event,x,primary);secondary=exponential(np.array([r['secondary_time'] for r in rows]),np.array([r['secondary_event'] for r in rows]),x);primary['factor_per_1000']=float(np.exp(-primary['beta']));primary['factor_CI95']=np.exp(-np.array(primary['beta_CI95'][::-1])).tolist();secondary['factor_per_1000']=float(np.exp(-secondary['beta']))
  models={'exponential_in_Omega':primary,'power_in_Omega':power,'Weibull_in_Omega':weib}
  for v in models.values():k=v['parameters'];v['AICc']=float(2*k-2*v['loglik']+2*k*(k+1)/(150-k-1))
  cohorts=[]
  for om in (100,300,1000,3000,10000):
   mask=omega==om;tt=times[mask];ee=event[mask];curve,rmst,median=km(tt,ee);d=int(ee.sum());exposure=float(tt.sum());rate=d/exposure;probCI=[0. if d==0 else float(betadist.ppf(.025,d,31-d)),1. if d==30 else float(betadist.ppf(.975,d+1,30-d))];predictions={}
   for name,model in models.items():
    xx=np.log(om/1000) if name=='power_in_Omega' else (om-1000)/1000
    if name=='Weibull_in_Omega':scale=np.exp(model['alpha']+model['beta']*xx);hazard=(HORIZON/scale)**model['shape'];mean=float(scale*np.exp(__import__('scipy').special.gammaln(1+1/model['shape'])))
    else:rr=np.exp(model['alpha']+model['beta']*xx);hazard=rr*HORIZON;mean=float(1/rr)
    probability=float(-np.expm1(-hazard));predictions[name]=dict(expected_events=30*probability,event_probability=probability,extrapolated_mean=mean,outside_observed_CP95=not(probCI[0]<=probability<=probCI[1]))
   cohorts.append(dict(omega=om,events=d,censored=30-d,exposure=exposure,constant_hazard_rate=rate,constant_hazard_mean=None if d==0 else 1/rate,zero_event_mean_lower95=exposure/-np.log(.05) if d==0 else None,restricted_mean_980=rmst,KM_median=median,KM=curve,observed_event_fraction_CP95=probCI,replication_merge_runs=sum(r['replication_or_merge'] for r in rows if r['omega']==om),secondary_events=sum(r['secondary_event'] for r in rows if r['omega']==om),predictions=predictions))
  flags=[]
  for name in ('power_in_Omega','Weibull_in_Omega'):
   if primary['AICc']-models[name]['AICc']>=6:flags.append(name+' favored by AICc >= 6')
  if not weib['shape_CI95'][0]<=1<=weib['shape_CI95'][1]:flags.append('Weibull shape interval excludes constant hazard')
  for c in cohorts:
   if c['predictions']['exponential_in_Omega']['outside_observed_CP95']:flags.append('primary event fraction outside observed interval at Omega '+str(c['omega']))
  atomic_json(OUT/'fits.json',dict(models=models,secondary_solitary_exponential=secondary,cohorts=cohorts,primary_diagnostic_flags=flags,rows=rows,horizon=HORIZON,primary_quotable_without_qualification=not flags))
  # Separate rate identities, KM/RMST check, and censoring fixture, without random draws.
  fixture=np.ones(30)*HORIZON;curve,area,med=km(fixture,np.zeros(30,bool));assert area==HORIZON and med is None and curve[-1]['survival']==1
  assert sum(c['events'] for c in cohorts)==int(event.sum());score=np.sum(event)-np.sum(times*np.exp(primary['alpha']+primary['beta']*x));assert abs(score)<1e-7
  atomic_json(OUT/'verification.json',dict(passed=True,summary_event_times_match_exports=True,KM_integrals_match_capped_time_means=True,analytic_profile_matches_general_optimizer=True,zero_event_fixture=True,likelihood_intercept_score=float(score),input_sha256=hashlib.sha256((PUBLIC/'decision_exports.json').read_bytes()).hexdigest(),external_replication=False,shared_dependencies=['legacy X event detector','NumPy/SciPy'],new_simulations=False))
  atomic_json(OUT/'attempt.json',dict(status='COMPLETE',new_simulations=False));print(json.dumps(dict(models=models,flags=flags,secondary=secondary),indent=2))
 except BaseException:
  import traceback
  atomic_json(OUT/'attempt.json',dict(status='FAILED',traceback=traceback.format_exc()));raise
if __name__=='__main__':run()

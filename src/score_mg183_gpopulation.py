"""Direct G population and paired-tick analysis; no fitted particle observer."""
import os
for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS'):os.environ[k]='1'
import json,math
import numpy as np
from scipy.stats import t
from mg182_gscale import ROOT,setup
from mg10e2c1_pool import atomic_json
PUBLIC=ROOT/'results/MG12H3/G_population_pilot';RAW=ROOT/'local_raw/MG183_G_population'
def interval(values):
 v=np.asarray(values,float);mean=float(v.mean());se=float(v.std(ddof=1)/np.sqrt(len(v)));half=float(t.ppf(.975,len(v)-1)*se);return dict(n=len(v),mean=mean,sd=float(v.std(ddof=1)),ci95=[mean-half,mean+half],assumption='t interval, small n; not exact distribution-free coverage')
def welch(left,right):
 a=np.asarray(left,float);b=np.asarray(right,float);estimate=float(a.mean()-b.mean());va=float(a.var(ddof=1)/len(a));vb=float(b.var(ddof=1)/len(b));se=math.sqrt(va+vb);df=(va+vb)**2/(va*va/(len(a)-1)+vb*vb/(len(b)-1)) if se else None;half=float(t.ppf(.975,df)*se) if se else 0.;return dict(estimate=estimate,ci95=[estimate-half,estimate+half],df=df,assumption='Welch t interval, independent conditions, small n')
def run():
 if not(PUBLIC/'attempt.json').exists():raise RuntimeError('No scientific pilot has been run.')
 status=json.loads((PUBLIC/'attempt.json').read_text());assert status['status']=='COMPLETE';plan=json.loads((PUBLIC/'plan.json').read_text());from run_mg183_gpopulation import binding;assert binding()==plan['source'];scores={};profiles={};audits=0
 for job in plan['jobs']:
  c=setup(job['omega']);path=RAW/'production'/job['id'];prior=None
  for step in range(100):
   with np.load(path/f'checkpoint_{step:03d}.npz') as z:
    previous=z['initial'];final=z['final'];hist=z['histogram'];assert prior is None or np.array_equal(previous,prior)
    for side in range(len(previous)):
     counts=previous[side].ravel().copy();weights=hist[0]+hist[side+1] if job['paired'] else hist
     for ch,n in zip(c['model'].channels,weights):
      for pos,v in ch.change:counts[pos]+=int(n)*v
     np.testing.assert_array_equal(counts,final[side].ravel());assert counts.min()>=0 and counts.sum()==previous[side].sum();audits+=1
    assert z['available'].min()>=0;prior=final.copy()
  with np.load(path/'observations.npz') as z:counts=z['counts'];times=z['times'];omega=float(z['omega'])
  assert counts.shape[0]==11 and len(times)==11 and omega==job['omega'];np.testing.assert_array_equal(counts[-1],prior);u=c['reference'][0]-c['reference'][0].mean();u-=u.mean();u/=np.linalg.norm(u);q=np.einsum('tsn,n->ts',counts[:,:,0]/omega,u);key=(job['omega'],job['condition']);scores.setdefault(key,[]).append(q);profiles.setdefault(key,[]).append(counts[:,:,0]/omega)
 results=[];tick=[];arrays={}
 for omega in(10,40,160):
  prepared=np.array(scores[(omega,'prepared')]);control=np.array(scores[(omega,'homogeneous')]);assert len(prepared)==len(control)==4
  for j,when in enumerate(times):results.append(dict(omega=omega,time=float(when),prepared=interval(prepared[:,j,0]),control=interval(control[:,j,0]),difference=welch(prepared[:,j,0],control[:,j,0]),control_sd_snr=(float((prepared[:,j,0].mean()-control[:,j,0].mean())/control[:,j,0].std(ddof=1)) if control[:,j,0].std(ddof=1)>0 else None)))
  for condition in('prepared','homogeneous'):
   v=np.array(profiles[(omega,condition)]);arrays[f'Omega{omega}_{condition}_Gmean']=v.mean(axis=0);arrays[f'Omega{omega}_{condition}_Gsd']=v.std(axis=0,ddof=1)
 margin=.05*np.linalg.norm(setup(160)['reference'][0]-setup(160)['reference'][0].mean())
 for condition in('prepared','homogeneous'):
  v=np.array(scores[(160,condition)])
  for j in(5,10):
   d=interval(v[:,j,1]-v[:,j,0]);lo,hi=d['ci95'];tick.append(dict(condition=condition,time=float(times[j]),difference=d,margin=float(margin),classification='EQUIVALENT_IN_SCOPE' if lo>=-margin and hi<=margin else 'OUTSIDE_MARGIN' if lo>margin or hi<-margin else 'INCONCLUSIVE'))
 np.savez_compressed(PUBLIC/'G_profiles_summary.npz',times=times,**arrays);atomic_json(PUBLIC/'G_direct_results.json',dict(population_estimates=results,paired_tick_checks=tick,marginal_endpoints_checked=audits,full_particle_persistence_claim=False,variance_scaling_qualified=False,limits='Direct fixed-geometry G contrasts; n4 per condition. T.01 only, conditional model, no long-horizon G conclusion.'))
if __name__=='__main__':run()

import os
for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS'):os.environ[k]='1'
import time,json,gzip,hashlib,math,statistics,traceback,gc
from fractions import Fraction
import numpy as np
from scipy.sparse import coo_matrix
from scipy.stats import t
from run_mg183_gpopulation import ROOT,PUBLIC,RAW,binding,read
from mg182_gscale import setup
from mg10e2c1_pool import atomic_json

def interval(v):
 m=statistics.mean(map(float,v));sd=statistics.stdev(map(float,v));half=float(t.ppf(.975,len(v)-1))*sd/math.sqrt(len(v));return m,sd,[m-half,m+half]
def check():
 start=time.perf_counter();attempt=PUBLIC/'verification_attempt.json';assert not attempt.exists();atomic_json(attempt,dict(status='RUNNING'));endpoints=available_checks=rng_checks=trace_events=traces=0;groups={};max_projection_error=0.;max_ci_error=0.
 try:
  plan=read(PUBLIC/'plan.json');assert plan['source']==binding();ledger=read(RAW/'production'/'ledger.json');assert len(ledger['records'])==24 and all(r['status']=='COMPLETE' and r['returncode']==0 for r in ledger['records']);models={o:setup(o) for o in (10,40,160)};matrices={}
  for omega,c in models.items():
   sr=[];sc=[];sv=[];rr=[];rc=[];rv=[]
   for i,ch in enumerate(c['model'].channels):
    for pos,v in ch.change:sr.append(pos);sc.append(i);sv.append(v)
    for pos,v in ch.reserve:rr.append(pos);rc.append(i);rv.append(v)
   matrices[omega]=(coo_matrix((sv,(sr,sc)),shape=(3500,len(c['model'].channels))).tocsr(),coo_matrix((rv,(rr,rc)),shape=(3500,len(c['model'].channels))).tocsr())
  for job in plan['jobs']:
   c=models[job['omega']];model=c['model'];S,R=matrices[job['omega']];path=RAW/'production'/job['id'];meta=read(path/'metadata.json');assert meta['source']==plan['source'] and meta['job']==job;previous=None;previous_rng=None;observations=[]
   for step in range(100):
    with np.load(path/f'checkpoint_{step:03d}.npz') as z:
     initial=z['initial'];final=z['final'];av=z['available'];h=z['histogram'];before=json.loads(str(z['rng_before']));after=json.loads(str(z['rng_after']));assert before==(previous_rng if previous_rng is not None else np.random.default_rng(job['seed']).bit_generator.state);rng_checks+=1
     if previous is not None:np.testing.assert_array_equal(initial,previous)
     if step==0:observations.append(initial.copy())
     for side in range(len(initial)):
      hist=h[0]+h[side+1] if job['paired'] else h;np.testing.assert_array_equal(initial[side].ravel()+S@hist,final[side].ravel());assert final[side].min()>=0 and int(final[side].sum())==int(initial[side].sum());endpoints+=1
      if not job['paired'] or side==0:np.testing.assert_array_equal(initial[side].ravel()-R@hist,av[side].ravel());available_checks+=1
     assert av.min()>=0;previous=final.copy();previous_rng=after
     if (step+1)%10==0:observations.append(final.copy())
   with np.load(path/'observations.npz') as z:counts=z['counts'];times=z['times'];np.testing.assert_array_equal(counts,np.array(observations));np.testing.assert_allclose(times,np.linspace(0,.01,11),rtol=0,atol=1e-16)
   g=c['reference'][0];u=(g-g.mean());u/=np.linalg.norm(u);v=counts[:,:,0].astype(float)/job['omega'];centered=v-v.mean(axis=2,keepdims=True);proj=np.sum(centered*u,axis=2);groups.setdefault((job['omega'],job['condition']),[]).append(proj)
   if job['sentinel']:
    with np.load(path/'initial.npz') as z:initial=z['counts']
    if not job['paired']:
     with np.load(path/'sentinel.npz') as z:ids=z['ids'];marks=z['marks']
     av=initial.ravel().copy();hist=np.zeros(len(model.channels),np.int64);last=0.
     for index,mark in zip(ids,marks):
      assert last<mark<.0001;last=mark;ch=model.channels[index];assert ch.rate(av)>0
      for pos,n in ch.reserve:assert av[pos]>=n;av[pos]-=n
      hist[index]+=1
     with np.load(path/'checkpoint_000.npz') as z:np.testing.assert_array_equal(hist,z['histogram']);np.testing.assert_array_equal(av,z['available'][0].ravel());np.testing.assert_array_equal(initial.ravel()+S@hist,z['final'][0].ravel())
     trace_events+=len(ids);traces+=1
    else:
     out=json.loads(gzip.decompress((path/'sentinel.json.gz').read_bytes()));av=np.stack([initial.ravel(),initial.ravel()]);output=av.copy();hist=np.zeros((3,len(model.channels)),np.int64);index=0;last=0.;calendar=[Fraction('0.0001'),Fraction('0.00002')];ds=calendar.copy()
     for boundary in out['boundaries']:
      limit=boundary['time'];expected=min(calendar);assert limit==float(expected);which=[s for s in (0,1) if calendar[s]==expected];assert boundary['mask']==sum(1<<s for s in which)
      while index<len(out['trace']) and out['trace'][index]['time']<limit:
       event=out['trace'][index];mark=event['time'];assert last<mark;last=mark;ch=model.channels[event['channel']];mask=event['mask'];assert mask in(1,2,3);hist[{3:0,1:1,2:2}[mask],event['channel']]+=1
       for side in(0,1):
        if mask&(1<<side):
         assert ch.rate(av[side])>0
         for pos,n in ch.reserve:assert av[side,pos]>=n;av[side,pos]-=n
         for pos,vv in ch.change:output[side,pos]+=vv
       index+=1
      np.testing.assert_array_equal(av.reshape(2,7,500),boundary['available']);np.testing.assert_array_equal(output.reshape(2,7,500),boundary['output'])
      for side in which:av[side]=output[side].copy();calendar[side]+=ds[side]
     assert index==len(out['trace']);np.testing.assert_array_equal(hist,out['histogram'])
     with np.load(path/'checkpoint_000.npz') as z:np.testing.assert_array_equal(hist,z['histogram']);np.testing.assert_array_equal(output.reshape(2,7,500),z['final'])
     trace_events+=index;traces+=1;del out;gc.collect()
   assert time.perf_counter()-start<500
  result=read(PUBLIC/'G_direct_results.json');tickchecks=[]
  for row in result['population_estimates']:
   j=int(round(row['time']/.001));p=np.array(groups[(row['omega'],'prepared')])[:,j,0];h=np.array(groups[(row['omega'],'homogeneous')])[:,j,0]
   for v,key in[(p,'prepared'),(h,'control')]:
    m,sd,ci=interval(v);np.testing.assert_allclose([m,sd,*ci],[row[key]['mean'],row[key]['sd'],*row[key]['ci95']],rtol=1e-9,atol=1e-11);max_projection_error=max(max_projection_error,abs(m-row[key]['mean']))
   m=statistics.mean(p)-statistics.mean(h);va=statistics.variance(p)/4;vb=statistics.variance(h)/4;df=(va+vb)**2/(va*va/3+vb*vb/3) if va+vb>0 else None;half=float(t.ppf(.975,df)*math.sqrt(va+vb)) if df else 0.;ci=[m-half,m+half];np.testing.assert_allclose(ci,row['difference']['ci95'],rtol=1e-9,atol=2e-11);max_ci_error=max(max_ci_error,max(abs(x-y) for x,y in zip(ci,row['difference']['ci95'])))
  margin=.05*np.linalg.norm(models[160]['reference'][0]-models[160]['reference'][0].mean())
  for row in result['paired_tick_checks']:
   j=int(round(row['time']/.001));values=np.array(groups[(160,row['condition'])]);m,sd,ci=interval(values[:,j,1]-values[:,j,0]);np.testing.assert_allclose(ci,row['difference']['ci95'],rtol=1e-9,atol=2e-11);classification='EQUIVALENT_IN_SCOPE' if ci[0]>=-margin and ci[1]<=margin else 'OUTSIDE_MARGIN' if ci[0]>margin or ci[1]<-margin else 'INCONCLUSIVE';assert classification==row['classification'];tickchecks.append(classification)
  assert endpoints==3200 and available_checks==2400 and rng_checks==2400 and traces==6;atomic_json(PUBLIC/'separate_check.json',dict(passed=True,completed_jobs=24,marginal_endpoints=endpoints,coarse_or_single_reservation_checks=available_checks,RNG_continuity_checks=rng_checks,ordered_sentinel_traces=traces,ordered_sentinel_events=trace_events,projection_estimates=33,paired_tick_intervals=len(tickchecks),max_projection_mean_difference=max_projection_error,max_CI_difference=max_ci_error,elapsed_seconds=time.perf_counter()-start,limits='Independent sparse ledger, centered projection/statistics and chronological sentinel replay; shared input fields/model/Student reference and same author, not external replication. Fine reservations replayed in sentinel segments, not all fine ticks.'));atomic_json(attempt,dict(status='COMPLETE',elapsed_seconds=time.perf_counter()-start));files={str(p.relative_to(RAW)):dict(bytes=p.stat().st_size,sha256=hashlib.sha256(p.read_bytes()).hexdigest()) for p in RAW.rglob('*') if p.is_file()};atomic_json(PUBLIC/'local_raw_manifest.json',dict(directory=str(RAW),files=files,total_bytes=sum(v['bytes'] for v in files.values()),full_raw_published=False))
 except BaseException:atomic_json(attempt,dict(status='FAILED',elapsed_seconds=time.perf_counter()-start,error=traceback.format_exc()));raise
if __name__=='__main__':check()

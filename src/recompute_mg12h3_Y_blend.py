"""MG154 independent compact archive and fixed-blend variance recomputation."""
import os
for key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[key]='1'
from pathlib import Path
import numpy as np,json,math,hashlib
from scipy.sparse import coo_matrix
from scipy.stats import skew,kurtosis,t
from mg12h3_Y_blend import ROOT,BASE,hashes,ALPHA,core
setup=core.setup
power=core.power

def run():
 r=json.loads((BASE/'registration.json').read_text());official=json.loads((BASE/'evaluation.json').read_text());assert hashes()==r['source_hashes'];initial,contrast=setup();model=contrast.model;n=initial.size;m=len(model.channels);rr=[];rc=[];rv=[];sr=[];sc=[];sv=[]
 for j,c in enumerate(model.channels):
  for i,v in c.reserve:rr.append(i);rc.append(j);rv.append(v)
  for i,v in c.change:sr.append(i);sc.append(j);sv.append(v)
 R=coo_matrix((np.asarray(rv,np.int64),(rr,rc)),shape=(n,m)).tocsr();D=coo_matrix((np.asarray(sv,np.int64),(sr,sc)),shape=(n,m)).tocsr();groups={level:[] for level in range(3)};records=[];eventsum=0;ticksum=0;maxdrift=0.
 for job in r['jobs']:
  path=BASE/'production'/job['id'];meta=json.loads((BASE/'production'/(job['id']+'.json')).read_text());savedmeta=json.loads((path/'metadata.json').read_text());assert meta['status']=='COMPLETE' and meta['job']==job and meta['seed']==r['seeds'][job['id']] and savedmeta['source_hashes']==core.hashes();assert savedmeta['environment_sha256']==hashlib.sha256(json.dumps(savedmeta['environment'],sort_keys=True).encode()).hexdigest();coeff=np.array(savedmeta['contrast']['coefficients']);np.testing.assert_array_equal(coeff,contrast.coefficients);w=np.array(savedmeta['contrast']['weights']);np.testing.assert_array_equal(w,contrast.weights);pop=initial.copy();cumulative=contrast.direct(initial);frames=[pop.copy()];direct=[cumulative.copy()];est=[cumulative.copy()];powers=[np.zeros(3)];step=0;events=0;prior_rng=None
  for chunk in sorted(path.glob('chunk_*.npz')):
   with np.load(chunk) as z:
    steps=z['steps'];hist=z['histogram'].astype(np.int64);final=z['final'];available=z['available'];exposure=z['exposures'];assert hist.shape==(len(steps),m) and hist.min()>=0 and np.isfinite(exposure).all() and exposure.min()>=0;assert steps.tolist()==list(range(step,step+len(steps)));before=np.concatenate((pop[None],final[:-1]),axis=0).reshape(len(steps),n);np.testing.assert_array_equal(before-(R@hist.T).T,available.reshape(len(steps),n));np.testing.assert_array_equal(before+(D@hist.T).T,final.reshape(len(steps),n));assert final.min()>=0 and available.min()>=0 and np.all(final.sum(axis=(1,2))==initial.sum());integral=np.einsum('ij,jk->ik',exposure,coeff);np.testing.assert_allclose(integral,z['integral'],rtol=1e-11,atol=1e-12);np.testing.assert_allclose(final.reshape(len(steps),n)@w.T,z['direct'],rtol=1e-11,atol=1e-12);assert z['events'].tolist()==hist.sum(axis=1).tolist();offsets=z['trace_offsets'];assert len(offsets)==len(steps)+1 and offsets[0]==0 and offsets[-1]==len(z['trace_ids'])==len(z['trace_marks'])
    for k in range(len(steps)):
     assert 0<=z['last_accepted_time'][k]<job['delta'] and (np.isnan(z['next_mark'][k]) or z['next_mark'][k]>=job['delta']);before_rng=str(z['rng_before'][k]);assert (json.loads(before_rng)==np.random.default_rng(meta['seed']).bit_generator.state if prior_rng is None else before_rng==prior_rng);prior_rng=str(z['rng_after'][k]);cumulative+=integral[k];np.testing.assert_allclose(cumulative,z['cumulative'][k],rtol=1e-11,atol=1e-12)
     checks=json.loads(str(z['checks'][k]));assert checks and checks[-1]['events']==int(z['events'][k]);maxdrift=max(maxdrift,max(x['max_drift_difference'] for x in checks))
     if job['sentinel']:
      ids=z['trace_ids'][offsets[k]:offsets[k+1]];marks=z['trace_marks'][offsets[k]:offsets[k+1]];assert np.all(np.diff(marks)>0) and np.all(marks>0) and np.all(marks<job['delta']);np.testing.assert_array_equal(np.bincount(ids,minlength=m),hist[k])
     else:assert offsets[k]==offsets[k+1]
     step+=1;events+=int(z['events'][k])
     if step%job['stride']==0:frames.append(final[k].copy());direct.append(contrast.direct(final[k]));est.append(cumulative.copy());powers.append(power(final[k],initial))
    pop=final[-1].copy()
  assert step==round(job['duration']/job['delta']) and events==meta['events']
  with np.load(path/'observations.npz') as z:
   assert int(z['seed'])==meta['seed'];np.testing.assert_array_equal(frames,z['counts']);np.testing.assert_allclose(z['direct'],direct,rtol=1e-11,atol=1e-12);np.testing.assert_allclose(z['estimates'],est,rtol=1e-11,atol=1e-12);np.testing.assert_allclose(z['power'],powers,rtol=1e-11,atol=1e-12);np.testing.assert_allclose(z['times'],np.arange(11)*.001,rtol=0,atol=1e-14)
  groups[job['level']].append((np.asarray(direct)[[5,10]],np.asarray(est)[[5,10]],np.asarray(powers)[[5,10]]));records.append(dict(id=job['id'],ticks=step,events=events,body_seconds=meta['elapsed_seconds'],sentinel=job['sentinel'],passed=True));eventsum+=events;ticksum+=step
 direct=np.array([[v[0][:,2] for v in groups[l]] for l in range(3)]);integrated=np.array([[v[1][:,2] for v in groups[l]] for l in range(3)]);blended=ALPHA[None,None,:]*direct+(1-ALPHA)[None,None,:]*integrated
 with np.load(BASE/'per_path_Y.npz') as z:np.testing.assert_allclose(z['direct'],direct,rtol=1e-12,atol=1e-12);np.testing.assert_allclose(z['integrated'],integrated,rtol=1e-12,atol=1e-12);np.testing.assert_allclose(z['blend'],blended,rtol=1e-12,atol=1e-12);np.testing.assert_array_equal(z['alpha'],ALPHA)
 index_path=BASE/'preparation_v2/bootstrap_indices.npz';assert hashlib.sha256(index_path.read_bytes()).hexdigest()==r['bootstrap_sha256']
 with np.load(index_path) as z:indices=z['indices']
 def variance(x):
  mean=math.fsum(map(float,x))/len(x);return math.fsum((float(v)-mean)**2 for v in x)/(len(x)-1)
 def quantile(x,p):
  a=np.sort(x);at=(len(a)-1)*p;k=int(math.floor(at));return float(a[k]+(at-k)*(a[min(k+1,len(a)-1)]-a[k]))
 statistics=[];bootstraps=[]
 for ti,key in enumerate(('secondary_early','primary_late')):
  a=direct[:,:,ti];b=blended[:,:,ti];v1=np.array([variance(row) for row in a]);v2=np.array([variance(row) for row in b]);pooled1=math.fsum(v1)/3;pooled2=math.fsum(v2)/3;per=v2/v1;ratios=[]
  for sample in indices:
   va=math.fsum(variance(a[g,sample[g]]) for g in range(3))/3;vb=math.fsum(variance(b[g,sample[g]]) for g in range(3))/3;ratios.append(vb/va if va>0 else math.nan)
  ratios=np.asarray(ratios);assert np.isfinite(ratios).all();lo,hi=quantile(ratios,.025),quantile(ratios,.975)
  if hi<1-1e-12 and np.all(per<1-1e-12):verdict='VARIANCE BENEFIT SUPPORTED IN THIS VALIDATION'
  elif lo>1+1e-12:verdict='EVIDENCE OF HIGHER POOLED VARIANCE'
  elif hi<1-1e-12:verdict='POOLED BENEFIT; TICK-SPECIFIC VALIDATION UNRESOLVED'
  else:verdict='UNRESOLVED VARIANCE BENEFIT'
  old=official[key];np.testing.assert_allclose([pooled1,pooled2,pooled2/pooled1,lo,hi],[old['pooled_direct_variance'],old['pooled_blend_variance'],old['blend_over_direct'],*old['bootstrap_95_interval']],rtol=1e-11,atol=1e-12);np.testing.assert_allclose(per,old['per_tick_blend_over_direct'],rtol=1e-11,atol=1e-12);assert verdict==old['verdict'];statistics.append(dict(time=(.005,.01)[ti],alpha_direct=float(ALPHA[ti]),pooled_direct_variance=pooled1,pooled_blend_variance=pooled2,blend_over_direct=pooled2/pooled1,variance_gain=pooled1/pooled2,per_tick_blend_over_direct=per.tolist(),bootstrap_95_interval=[lo,hi],verdict=verdict));bootstraps.append(ratios)
 with np.load(BASE/'bootstrap_ratios.npz') as z:np.testing.assert_allclose(bootstraps[0],z['early'],rtol=1e-11,atol=1e-12);np.testing.assert_allclose(bootstraps[1],z['late'],rtol=1e-11,atol=1e-12)
 result=dict(passed=True,trajectories=24,ticks=ticksum,events=eventsum,records=records,statistics=statistics,maximum_saved_drift_difference=maxdrift,source_hashes_unchanged=True,weights_refitted=False,equivalence_verdict=None,shared='Frozen compact model/coefficients and registered bootstrap indexes; independent sparse-matrix count and exposure arithmetic, scalar variances and resampling loops, sorted linear quantiles. Same author, not external replication; ordered exposure replay remains the registered sentinel audit.')
 (BASE/'separate_recomputation.json').write_text(json.dumps(result,indent=2)+'\n');files=[dict(path=p.relative_to(ROOT).as_posix(),bytes=p.stat().st_size,sha256=hashlib.sha256(p.read_bytes()).hexdigest()) for p in sorted((BASE/'production').rglob('*')) if p.is_file()];(BASE/'raw_manifest.json').write_text(json.dumps(dict(files=files,total_bytes=sum(x['bytes'] for x in files)),indent=2)+'\n');print(json.dumps(dict(passed=True,trajectories=24,ticks=ticksum,events=eventsum,statistics=statistics,files=len(files),bytes=sum(x['bytes'] for x in files))))
if __name__=='__main__':run()

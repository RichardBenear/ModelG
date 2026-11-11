"""MG149 saved compact-pilot recomputation, no sampling or qualification verdict."""
import os
for key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[key]='1'
from pathlib import Path
import numpy as np,json,math,hashlib
from scipy.sparse import coo_matrix
from scipy.stats import skew,kurtosis,t
from mg12h3_finer_precision import ROOT,BASE,setup,hashes,power

def run():
 r=json.loads((BASE/'registration.json').read_text());official=json.loads((BASE/'evaluation.json').read_text());assert hashes()==r['source_hashes'];initial,contrast=setup();model=contrast.model;n=initial.size;m=len(model.channels);rr=[];rc=[];rv=[];sr=[];sc=[];sv=[]
 for j,c in enumerate(model.channels):
  for i,v in c.reserve:rr.append(i);rc.append(j);rv.append(v)
  for i,v in c.change:sr.append(i);sc.append(j);sv.append(v)
 R=coo_matrix((np.asarray(rv,np.int64),(rr,rc)),shape=(n,m)).tocsr();D=coo_matrix((np.asarray(sv,np.int64),(sr,sc)),shape=(n,m)).tocsr();groups={level:[] for level in range(3)};records=[];eventsum=0;ticksum=0;maxdrift=0.
 for job in r['jobs']:
  path=BASE/'production'/job['id'];meta=json.loads((BASE/'production'/(job['id']+'.json')).read_text());savedmeta=json.loads((path/'metadata.json').read_text());assert meta['status']=='COMPLETE' and meta['job']==job and meta['seed']==r['seeds'][job['id']] and savedmeta['source_hashes']==hashes();assert savedmeta['environment_sha256']==hashlib.sha256(json.dumps(savedmeta['environment'],sort_keys=True).encode()).hexdigest();coeff=np.array(savedmeta['contrast']['coefficients']);np.testing.assert_array_equal(coeff,contrast.coefficients);w=np.array(savedmeta['contrast']['weights']);np.testing.assert_array_equal(w,contrast.weights);pop=initial.copy();cumulative=contrast.direct(initial);frames=[pop.copy()];direct=[cumulative.copy()];est=[cumulative.copy()];powers=[np.zeros(3)];step=0;events=0;prior_rng=None
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
 precision=[]
 for level,values in groups.items():
  a=np.array([x[0] for x in values]);b=np.array([x[1] for x in values]);p=np.array([x[2] for x in values]);v1=a.var(0,ddof=1);v2=b.var(0,ddof=1);row=dict(level=level,replicates=len(values),direct_variance=v1.tolist(),integral_variance=v2.tolist(),variance_ratio=(v1/v2).tolist(),direct_mean=a.mean(0).tolist(),integral_mean=b.mean(0).tolist(),power_mean=p.mean(0).tolist(),power_variance=p.var(0,ddof=1).tolist());prior=official['precision'][level]
  for key in ('direct_variance','integral_variance','variance_ratio','direct_mean','integral_mean','power_mean','power_variance'):np.testing.assert_allclose(row[key],prior[key],rtol=1e-10,atol=1e-12)
  row.update(direct_skew=skew(a,axis=0,bias=False).tolist(),integral_skew=skew(b,axis=0,bias=False).tolist(),power_skew=skew(p,axis=0,bias=False).tolist());precision.append(row)
 # Optional hybrid diagnostics are explicitly exploratory, fitted on pilot data.
 hybrid=[]
 for time_index in range(2):
  for s in range(3):
   direct=np.concatenate([np.asarray([x[0] for x in values])[:,time_index,s]-np.mean([x[0][time_index,s] for x in values]) for values in groups.values()]);integral=np.concatenate([np.asarray([x[1] for x in values])[:,time_index,s]-np.mean([x[1][time_index,s] for x in values]) for values in groups.values()]);va=float(direct@direct/21);vb=float(integral@integral/21);cov=float(direct@integral/21);alpha=(vb-cov)/(va+vb-2*cov);variance=alpha*alpha*va+(1-alpha)**2*vb+2*alpha*(1-alpha)*cov;hybrid.append(dict(time=(.005,.01)[time_index],species=('G','X','Y')[s],alpha_direct=alpha,pooled_direct_variance=va,pooled_integral_variance=vb,pooled_hybrid_in_sample_variance=variance,scope='Exploratory same-pilot fitted coefficient; optimistic in-sample variance, not validated or deployed'))
 result=dict(passed=True,trajectories=24,ticks=ticksum,events=eventsum,records=records,precision=precision,exploratory_hybrid=hybrid,maximum_saved_drift_difference=maxdrift,source_hashes_unchanged=True,equivalence_verdict=None,shared='Frozen compact channel compilation and saved coefficients; sparse-matrix integer aggregation and independent exposure/statistic arithmetic. Same author, not external replication. Ordered exposure replay remains the registered sentinel audit.')
 (BASE/'separate_recomputation.json').write_text(json.dumps(result,indent=2)+'\n');files=[dict(path=p.relative_to(ROOT).as_posix(),bytes=p.stat().st_size,sha256=hashlib.sha256(p.read_bytes()).hexdigest()) for p in sorted((BASE/'production').rglob('*')) if p.is_file()];(BASE/'raw_manifest.json').write_text(json.dumps(dict(files=files,total_bytes=sum(x['bytes'] for x in files)),indent=2)+'\n');print(json.dumps(dict(passed=True,trajectories=24,ticks=ticksum,events=eventsum,files=len(files),bytes=sum(x['bytes'] for x in files),maximum_saved_drift_difference=maxdrift)))
if __name__=='__main__':run()

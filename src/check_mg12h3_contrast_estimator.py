"""MG146 bounded estimator/compact writer validation, no production pilot."""
import os
for key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[key]='1'
from pathlib import Path
import numpy as np,json,time,math,hashlib,sys,traceback
from mg12h3_contrast_estimator import Contrast,tick,aggregate_audit,trace_audit,save_tick
from mg12h3_sparse import Model,tick as frozen_tick
from prepare_mg12h3_bridge import fixture,config,ROOT,seed_for
from mg10e2c1_pool import atomic_json
BASE=ROOT/'results/MG12H3/contrast_estimator'
FILES=['src/mg12h3_contrast_estimator.py','src/check_mg12h3_contrast_estimator.py','src/mg12h3_sparse.py','src/prepare_mg12h3_bridge.py','src/mg12h1.py','src/mg12h2.py','src/runner.py','spec/MG12H3_contrast_estimator_development.md']
def hashes():return {p:hashlib.sha256((ROOT/p).read_text(encoding='utf-8').replace('\r\n','\n').encode()).hexdigest() for p in FILES}
def plan():
 assert not (BASE/'plan.json').exists();used=set()
 def collect(x):
  if isinstance(x,dict):
   for v in x.values():collect(v)
  elif isinstance(x,list):
   for v in x:collect(v)
  elif type(x) is int and 0<=x<2**32:used.add(x)
 for p in (ROOT/'results/MG12H3').rglob('registration.json'):collect(json.loads(p.read_text()).get('seeds',{}))
 def seeds(c,stage,n):
  out=[]
  for r in range(n):
   nonce=0;v=seed_for(12,c,0,stage,r)
   while v in used:nonce+=1;v=seed_for(12,c,0,stage,r,nonce)
   used.add(v);out.append(v)
  return out
 atomic_json(BASE/'plan.json',dict(source_hashes=hashes(),single=seeds(0,14600,4096),paired=[seeds(c,14601,32) for c in range(4)],mixed=[seeds(c,14602,512) for c in range(4)],ceiling=60,workers=1,formal=False))

def weights(n):
 w=np.zeros((3,7*n));w[0,:n]=np.arange(1,n+1)/n;w[1,n:2*n]=1/n;w[2,2*n:3*n]=np.where(np.arange(n)%2,1.,-1.);return w

def one_model():
 cfg=(1.,[(14.,4.)],[],(0.,0.,0.),(0.,0.),(True,False,False,False,False,False));model=Model(1,cfg);a=np.zeros((7,1),np.int64);a[3,0]=1;w=np.zeros((1,7));w[0,0]=1;return a,Contrast(model,w)

def static():
 class BG:state={'scripted':True}
 class Fixed:
  bit_generator=BG()
  def __init__(self,wait):self.wait=wait
  def exponential(self,scale):return self.wait
  def random(self):return .5
 a,c=one_model();out=tick(c,a,.05,Fixed(.02),sentinel=True,check_every=1);aggregate_audit(c,a,out);trace_audit(c,a,out,.05);np.testing.assert_allclose(out['integral'],[.02],rtol=0,atol=1e-15);assert out['counts'][0,0]==1
 boundary=tick(c,a,.05,Fixed(.05),sentinel=True);assert boundary['events']==0;np.testing.assert_allclose(boundary['integral'],[.05]);trace_audit(c,a,boundary,.05)
 zero=np.zeros_like(a);empty=tick(c,zero,.05,Fixed(.02),sentinel=True);assert empty['events']==0 and empty['integral'][0]==0
 zero_weights=Contrast(c.model,np.zeros((1,7)));z=tick(zero_weights,a,.05,Fixed(.02),sentinel=True);assert z['integral'][0]==0 and len(z['exposures'])==0
 bad=dict(out,exposures=out['exposures']+.1)
 try:aggregate_audit(c,a,bad);raise RuntimeError('Corrupt exposure accepted')
 except AssertionError:pass
 bad=dict(out,histogram=out['histogram']*0)
 try:aggregate_audit(c,a,bad);raise RuntimeError('Corrupt histogram accepted')
 except AssertionError:pass
 return dict(passed=True,scripted_checks=6)

def run():
 p=json.loads((BASE/'plan.json').read_text());assert p['source_hashes']==hashes() and not (BASE/'attempt.json').exists();clock=time.perf_counter();deadline=clock+60;atomic_json(BASE/'attempt.json',dict(status='STARTED'));scripted=static();paired=[];saved=[]
 for c in range(4):
  state=fixture(c);model=Model(len(state.sites),config(state));contrast=Contrast(model,weights(model.n))
  for replicate,seed in enumerate(p['paired'][c]):
   a=state.counts();rng=np.random.default_rng(seed);old_rng=np.random.default_rng(seed)
   for step in range(3):
    out=tick(contrast,a,.02,rng,deadline=deadline,sentinel=True,check_every=1);old=frozen_tick(model,a,.02,old_rng,deadline=deadline);np.testing.assert_array_equal(out['counts'],old['counts']);np.testing.assert_array_equal(out['available'],old['available']);assert [(x['channel'],x['time']) for x in out['trace']]==[(x['channel'],x['time']) for x in old['trace']];assert json.dumps(rng.bit_generator.state)==json.dumps(old_rng.bit_generator.state);aggregate_audit(contrast,a,out);trace_audit(contrast,a,out,.02)
    if replicate==0:
     file=BASE/f'paired_c{c}_tick{step}.npz';save_tick(file,contrast,a,out);saved.append(file.name)
    a=out['counts']
   paired.append(dict(fixture=c,seed=seed,replicate=replicate,ticks=3,passed=True))
 a,contrast=one_model();initial=float(contrast.direct(a)[0]);single=[];hist=[];exposures=[];final=[]
 for seed in p['single']:
  out=tick(contrast,a,.05,np.random.default_rng(seed),deadline=deadline);aggregate_audit(contrast,a,out);single.append([float(out['direct'][0]),initial+float(out['integral'][0])]);hist.append(out['histogram']);exposures.append(out['exposures']);final.append(out['counts'])
 values=np.asarray(single);expected=-math.expm1(-.05);checks=[]
 for k,name in enumerate(('direct','integral')):
  mean=float(values[:,k].mean());se=float(values[:,k].std(ddof=1)/math.sqrt(len(values)));assert abs(mean-expected)<=6*se+1e-6;checks.append(dict(estimator=name,mean=mean,exact_mean=expected,MCSE=se,passed=True))
 np.savez_compressed(BASE/'single_channel.npz',values=values,histogram=np.asarray(hist),exposures=np.asarray(exposures),final=np.asarray(final),seeds=p['single']);ratio=float(values[:,0].var(ddof=1)/values[:,1].var(ddof=1));mixed=[]
 for c in range(4):
  state=fixture(c);a=state.counts();model=Model(len(state.sites),config(state));contrast=Contrast(model,weights(model.n));initial=contrast.direct(a);values=[];hist=[];exposures=[];final=[];available=[]
  for replicate,seed in enumerate(p['mixed'][c]):
   out=tick(contrast,a,.1,np.random.default_rng(seed),deadline=deadline,sentinel=replicate<2,check_every=1);aggregate_audit(contrast,a,out)
   if replicate<2:trace_audit(contrast,a,out,.1);save_tick(BASE/f'mixed_c{c}_r{replicate}.npz',contrast,a,out)
   values.append(np.stack([out['direct'],initial+out['integral']]));hist.append(out['histogram']);exposures.append(out['exposures']);final.append(out['counts']);available.append(out['available'])
  values=np.asarray(values);difference=values[:,0]-values[:,1];mean=difference.mean(axis=0);se=difference.std(axis=0,ddof=1)/math.sqrt(len(values));assert np.all(abs(mean)<=6*se+1e-6);np.savez_compressed(BASE/f'mixed_c{c}.npz',values=values,histogram=np.asarray(hist),exposures=np.asarray(exposures),final=np.asarray(final),available=np.asarray(available),seeds=p['mixed'][c],initial=a,metadata=json.dumps(contrast.metadata()));mixed.append(dict(fixture=c,paths=512,paired_mean_difference=mean.tolist(),paired_MCSE=se.tolist(),passed=True))
 assert time.perf_counter()-clock<60 and hashes()==p['source_hashes'];atomic_json(BASE/'results.json',dict(passed=True,elapsed_seconds=time.perf_counter()-clock,static=scripted,paired_paths=128,paired_ticks=384,paired=paired,single_channel=checks,single_channel_variance_reduction_ratio=ratio,mixed=mixed,scope='Bounded implementation diagnostics, not formal stochastic qualification or actual-ring precision',source_hashes=hashes(),saved_sentinels=saved))
if __name__=='__main__':
 try:plan() if '--plan' in sys.argv else run()
 except BaseException:atomic_json(BASE/'failure.json',dict(error=traceback.format_exc()));raise

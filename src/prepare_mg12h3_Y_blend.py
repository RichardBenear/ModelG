"""MG151 saved validation/MG152 prospective seed/bootstrap plan."""
from mg12h3_Y_blend import *
from runner import seed_for
PREP=BASE/'preparation'
def plan():
 assert not (BASE/'plan.json').exists();used=set()
 def collect(x):
  if isinstance(x,dict):
   for v in x.values():collect(v)
  elif isinstance(x,list):
   for v in x:collect(v)
  elif type(x) is int and 0<=x<2**32:used.add(x)
 for p in (ROOT/'results/MG12H3').rglob('registration.json'):collect(json.loads(p.read_text()).get('seeds',{}))
 collect(json.loads((ROOT/'results/MG12H3/contrast_estimator/plan.json').read_text()));seeds={};collisions=[]
 for j in core.jobs():
  nonce=0;value=seed_for(12,j['level'],0,15201,j['replicate'])
  while value in used:nonce+=1;value=seed_for(12,j['level'],0,15201,j['replicate'],nonce);collisions.append(dict(id=j['id'],nonce=nonce))
  used.add(value);seeds[j['id']]=value
 nonce=0;analysis=seed_for(12,0,0,15202,0)
 while analysis in used:nonce+=1;analysis=seed_for(12,0,0,15202,0,nonce)
 old=json.loads((ROOT/'results/MG12H3/finer_tick_precision_pilot/production/ledger.json').read_text());records={x['id']:x for x in old['records']};assert len(records)==24 and all(x['status']=='COMPLETE' for x in records.values());forecast={j['id']:1.5*records[j['id']]['charged_worker_seconds']+3 for j in core.jobs()};atomic_json(BASE/'plan.json',dict(source_hashes=hashes(),jobs=core.jobs(),seeds=seeds,bootstrap_seed=analysis,alpha_direct=ALPHA.tolist(),forecast=forecast,collisions=collisions))

def run():
 p=json.loads((BASE/'plan.json').read_text());assert p['source_hashes']==hashes() and p['alpha_direct']==ALPHA.tolist() and not (PREP/'attempt.json').exists();clock=time.perf_counter();atomic_json(PREP/'attempt.json',dict(status='STARTED',new_token_paths=0));indices=np.random.default_rng(p['bootstrap_seed']).integers(0,8,size=(10000,3,8),dtype=np.uint8);np.savez_compressed(PREP/'bootstrap_indices.npz',indices=indices,seed=p['bootstrap_seed']);known=[]
 with np.load(ROOT/'results/MG12H3/contrast_estimator/single_channel.npz') as z:values=z['values']
 expected=-math.expm1(-.05);mixed_values=[]
 for time_index in range(2):
  value=blend(values[:,0],values[:,1],time_index);mean=float(value.mean());se=float(value.std(ddof=1)/math.sqrt(len(value)));assert abs(mean-expected)<=6*se+1e-6;known.append(dict(weight_time=(.005,.01)[time_index],known_case_horizon=.05,mean=mean,known_mean=expected,MCSE=se,passed=True));mixed_values.append(value)
 np.savez_compressed(PREP/'known_case_blends.npz',values=np.asarray(mixed_values).T,expected=expected,alpha=ALPHA)
 small=[]
 for c in range(4):
  with np.load(ROOT/f'results/MG12H3/contrast_estimator/mixed_c{c}.npz') as z:a=z['values'][:,0,2];b=z['values'][:,1,2]
  rows=[]
  for time_index in range(2):
   w=float(ALPHA[time_index]);value=blend(a,b,time_index);np.testing.assert_allclose(value.mean(),w*a.mean()+(1-w)*b.mean(),rtol=1e-12,atol=1e-12);cov=float(np.cov(a,b,ddof=1)[0,1]);variance=w*w*a.var(ddof=1)+(1-w)**2*b.var(ddof=1)+2*w*(1-w)*cov;np.testing.assert_allclose(value.var(ddof=1),variance,rtol=1e-12,atol=1e-12);diff=value-a;mean=float(diff.mean());se=float(diff.std(ddof=1)/math.sqrt(len(diff)));assert abs(mean)<=6*se+1e-6;rows.append(dict(weight_time=(.005,.01)[time_index],paired_mean_difference=mean,MCSE=se,passed=True))
  small.append(dict(fixture=c,paths=len(a),checks=rows))
 # Synthetic statistics with paired resampling, no simulated token trajectories.
 x=np.arange(8,dtype=float)-3.5;a=np.array([x+10*c for c in range(3)]);b=np.array([-x+10*c for c in range(3)]);good,_=validate_variance(a,b,indices);assert good['verdict']=='VARIANCE BENEFIT SUPPORTED IN THIS VALIDATION';bad,_=validate_variance(a,np.array([3*x+10*c for c in range(3)]),indices);assert bad['verdict']=='EVIDENCE OF HIGHER POOLED VARIANCE';equal,_=validate_variance(a,a,indices);assert equal['verdict']=='UNRESOLVED VARIANCE BENEFIT';zero,_=validate_variance(np.ones((3,8)),np.ones((3,8)),indices);assert zero['verdict'].startswith('UNRESOLVED');guards=dict(synthetic_supported=True,synthetic_worse=True,synthetic_unresolved=True,zero_variance_unresolved=True)
 try:validate_variance(a[:,:7],b[:,:7],indices);raise RuntimeError('Incomplete accepted')
 except ValueError:guards['incomplete_rejected']=True
 mutated=ALPHA.copy();mutated[1]+=.1;guards['weight_mutation_detected']=mutated.tolist()!=p['alpha_direct'];assert all(guards.values())
 atomic_json(PREP/'synthetic_statistics.json',dict(supported=good,worse=bad,unresolved=equal,zero_variance=zero));forecast=sum(p['forecast'].values());assert forecast<10350 and max(p['forecast'].values())<600 and hashes()==p['source_hashes'] and time.perf_counter()-clock<60;atomic_json(BASE/'readiness.json',dict(passed=True,elapsed_seconds=time.perf_counter()-clock,known_mean_checks=known,small_case_checks=small,guards=guards,source_hashes=hashes(),alpha_direct=ALPHA.tolist(),plan_sha256=hashlib.sha256((BASE/'plan.json').read_bytes()).hexdigest(),bootstrap_sha256=hashlib.sha256((PREP/'bootstrap_indices.npz').read_bytes()).hexdigest(),forecast_worker_seconds=forecast,new_token_paths=0,production_registered=False))
if __name__=='__main__':
 try:plan() if '--plan' in sys.argv else run()
 except BaseException:atomic_json(PREP/'failure.json',dict(error=traceback.format_exc()));raise

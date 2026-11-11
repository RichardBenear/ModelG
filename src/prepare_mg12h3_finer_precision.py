"""MG147 actual compact writer preparation, no precision-pilot draws."""
from mg12h3_finer_precision import *
import shutil
PREP=BASE/'preparation_v2'
PLAN=BASE/'plan_v2.json'
def plan():
 assert not (PLAN).exists();used=set()
 def collect(x):
  if isinstance(x,dict):
   for v in x.values():collect(v)
  elif isinstance(x,list):
   for v in x:collect(v)
  elif type(x) is int and 0<=x<2**32:used.add(x)
 for p in (ROOT/'results/MG12H3').rglob('registration.json'):collect(json.loads(p.read_text()).get('seeds',{}))
 collect(json.loads((ROOT/'results/MG12H3/contrast_estimator/plan.json').read_text()));formal={};prep={}
 for j in jobs():
  nonce=0;v=seed_for(12,j['level'],0,14701,j['replicate'])
  while v in used:nonce+=1;v=seed_for(12,j['level'],0,14701,j['replicate'],nonce)
  used.add(v);formal[j['id']]=v
 for level in range(3):
  nonce=0;v=seed_for(12,level,0,14702,0)
  while v in used:nonce+=1;v=seed_for(12,level,0,14702,0,nonce)
  used.add(v);prep[str(level)]=v
 atomic_json(PLAN,dict(jobs=jobs(),seeds=formal,preparation_seeds=prep,source_hashes=hashes()))
def run():
 p=json.loads((PLAN).read_text());assert hashes()==p['source_hashes'] and not (PREP/'attempt.json').exists();clock=time.perf_counter();deadline=clock+60;atomic_json(PREP/'attempt.json',dict(status='STARTED'));timings={};storage={};guards={};sample=None
 for j in jobs():
  if j['replicate']!=0:continue
  short=dict(j,duration=.0001,observation=.0001,stride=round(.0001/j['delta']));start=time.perf_counter();write(short,p['preparation_seeds'][str(j['level'])],PREP,deadline);audit(short,p['preparation_seeds'][str(j['level'])],PREP,deadline);timings[str(j['level'])]=time.perf_counter()-start;storage[str(j['level'])]=sum(x.stat().st_size for x in (PREP/j['id']).rglob('*') if x.is_file());sample=short if sample is None else sample
 seed=p['preparation_seeds'][str(sample['level'])]
 def reject(name,fn):
  try:fn()
  except (AssertionError,FileNotFoundError):guards[name]=True;return
  raise ValueError('Guard accepted '+name)
 reject('wrong_seed',lambda:audit(sample,seed+1,PREP,deadline));reject('missing',lambda:audit(dict(sample,id='missing'),seed,PREP,deadline))
 for name in ('corrupt_exposure','incomplete_chunk'):
  target=PREP/name;target.mkdir();shutil.copyfile(PREP/(sample['id']+'.json'),target/(sample['id']+'.json'));shutil.copytree(PREP/sample['id'],target/sample['id']);file=next((target/sample['id']).glob('chunk_*.npz'))
  with np.load(file) as z:values={k:z[k].copy() for k in z.files}
  if name=='corrupt_exposure':values['exposures'][0,0]+=1.
  else:values['steps']=values['steps'][:0]
  np.savez_compressed(file,**values);reject(name,lambda:audit(sample,seed,target,deadline))
 forecast={j['id']:2*timings[str(j['level'])]*100+3 for j in jobs()};assert time.perf_counter()-clock<60 and hashes()==p['source_hashes'];atomic_json(BASE/'readiness.json',dict(passed=True,elapsed_seconds=time.perf_counter()-clock,timings=timings,storage=storage,guards=guards,forecast=forecast,forecast_worker_seconds=sum(forecast.values()),source_hashes=hashes(),plan_sha256=hashlib.sha256((PLAN).read_bytes()).hexdigest(),production_registered=False))
if __name__=='__main__':
 try:plan() if '--plan' in sys.argv else run()
 except BaseException:atomic_json(PREP/'failure.json',dict(error=traceback.format_exc()));raise

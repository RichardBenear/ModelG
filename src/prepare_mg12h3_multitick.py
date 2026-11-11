"""MG140 measured short actual writer preparation."""
from mg12h3_multitick import *
import shutil
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
 seeds={};prep={};collisions=[]
 for index,job in enumerate(jobs()):
  for stage,target in [(14001,seeds),(14000,prep)]:
   nonce=0;v=seed_for(12,index,0,stage,0)
   while v in used:nonce+=1;v=seed_for(12,index,0,stage,0,nonce);collisions.append(dict(id=job['id'],stage=stage,nonce=nonce))
   used.add(v);target[job['id']]=v
 atomic_json(BASE/'plan.json',dict(jobs=jobs(),seeds=seeds,preparation_seeds=prep,source_hashes=hashes(),collisions=collisions))
def run():
 p=json.loads((BASE/'plan.json').read_text());assert hashes()==p['source_hashes'] and not (PREP/'attempt.json').exists();clock=time.perf_counter();deadline=clock+60;atomic_json(PREP/'attempt.json',dict(status='STARTED'));timings={};guards={}
 chosen=[j for j in p['jobs'] if j['case'].startswith('o0_a0_')]
 for job in chosen:
  short=dict(job,duration=.001);start=time.perf_counter();run_job(short,p['preparation_seeds'][job['id']],PREP,deadline);audit_job(short,p['preparation_seeds'][job['id']],PREP);timings[job['id']]=time.perf_counter()-start
 sample=dict(chosen[0],duration=.001);seed=p['preparation_seeds'][sample['id']]
 def reject(name,fn):
  try:fn()
  except (AssertionError,FileNotFoundError):guards[name]=True;return
  raise ValueError('Guard accepted '+name)
 reject('wrong_seed',lambda:audit_job(sample,seed+1,PREP));reject('missing',lambda:audit_job(dict(sample,id='missing'),seed,PREP))
 for label in ('corrupt_marks','incomplete_frames'):
  target=PREP/label;target.mkdir();shutil.copyfile(PREP/(sample['id']+'.json'),target/(sample['id']+'.json'));shutil.copytree(PREP/sample['id'],target/sample['id'])
  file=target/sample['id']/('tick_0000.npz' if label=='corrupt_marks' else 'frames.npz')
  with np.load(file) as z:values={k:z[k].copy() for k in z.files}
  if label=='corrupt_marks':values['marks'][0]=-1.
  else:values['counts']=values['counts'][:-1]
  np.savez_compressed(file,**values);reject(label,lambda:audit_job(sample,seed,target))
 altered=dict(p['source_hashes']);altered[next(iter(altered))]='0'*64;guards['source_mutation']=altered!=hashes()
 forecast={}
 for job in jobs():
  representative='o0_a0_Omega'+('10' if job['case'].endswith('Omega10') else '1')+'_d'+str(job['level']);forecast[job['id']]=2*timings[representative]*10+3
 assert hashes()==p['source_hashes'] and time.perf_counter()-clock<60;atomic_json(BASE/'readiness.json',dict(passed=True,guards=guards,timings=timings,forecast=forecast,forecast_worker_seconds=sum(forecast.values()),elapsed_seconds=time.perf_counter()-clock,source_hashes=hashes(),plan_sha256=hashlib.sha256((BASE/'plan.json').read_bytes()).hexdigest(),production_registered=False))
if __name__=='__main__':
 try:plan() if '--plan' in sys.argv else run()
 except BaseException:atomic_json(PREP/'failure.json',dict(error=traceback.format_exc()));raise

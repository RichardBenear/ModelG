"""MG137 frozen actual-writer preparation and seed plan."""
from mg12h3_sparse_qualification import *
import shutil
PREP=BASE/'preparation_v2'
PLAN=BASE/'plan_v2.json'
def plan():
 assert not (PLAN).exists();used=set();collisions=[]
 def collect(x):
  if isinstance(x,dict):
   for v in x.values():collect(v)
  elif isinstance(x,list):
   for v in x:collect(v)
  elif type(x) is int and 0<=x<2**32:used.add(x)
 for p in (ROOT/'results/MG12H3').rglob('*.json'):
  if p.name not in ('registration.json','plan.json','readiness.json'):continue
  obj=json.loads(p.read_text())
  for key in ('seeds','tiny_seeds','actual_seeds','tiny','formal_seeds','preparation_seeds'):
   if key in obj:collect(obj[key])
 prior=len(used);formal={};prep={}
 def seed(c,stage,r):
  nonce=0;v=seed_for(12,c,0,stage,r)
  while v in used:collisions.append(dict(cell=c,stage=stage,rep=r,nonce=nonce));nonce+=1;v=seed_for(12,c,0,stage,r,nonce)
  used.add(v);return v
 for j in jobs():
  if j['kind']=='endpoint':formal[j['id']]=[seed(j['cell'],13800,j['chunk']*4096+r) for r in range(4096)]
 for c in range(4):prep[str(c)]=[seed(c,13700,r) for r in range(32)]
 atomic_json(PLAN,dict(source_hashes=hashes(),jobs=jobs(),formal_seeds=formal,preparation_seeds=prep,prior_seed_values=prior,collisions_resolved=collisions))

def run():
 p=json.loads((PLAN).read_text());assert p['source_hashes']==hashes() and not (PREP/'attempt.json').exists();clock=time.perf_counter();deadline=clock+60;atomic_json(PREP/'attempt.json',dict(status='STARTED'));control=controls();timings={};held={};guards={}
 for c in range(4):
  job=dict(id=f'prep_c{c}',kind='endpoint',cell=c,chunk=0,reps=32);start=time.perf_counter();meta=write_endpoint(job,p['preparation_seeds'][str(c)],PREP,deadline);timings[str(c)]=time.perf_counter()-start;held[c]=[check(job,p['preparation_seeds'][str(c)],PREP,deadline)]
  if c==0:
   def reject(label,fn,expected=(AssertionError,FileNotFoundError)):
    try:fn()
    except expected:guards[label]=True;return
    raise ValueError('Guard failed '+label)
   reject('wrong_seed',lambda:check(job,[v+1 for v in p['preparation_seeds']['0']],PREP,deadline))
   reject('missing',lambda:check(dict(job,id='missing'),[],PREP,deadline))
   reject('wrong_job',lambda:check(dict(job,reps=33),p['preparation_seeds']['0'],PREP,deadline))
   bad=PREP/'corrupt_trace';bad.mkdir();[shutil.copyfile(PREP/(job['id']+suffix),bad/(job['id']+suffix)) for suffix in ('.json','.jsonl','.npz')];path=bad/(job['id']+'.jsonl');lines=path.read_text().splitlines();changed=False
   for n,line in enumerate(lines):
    row=json.loads(line)
    if row['trace']:row['trace'][0]['increment'][0][0]+=1;lines[n]=json.dumps(row);changed=True;break
   assert changed;path.write_text('\n'.join(lines)+'\n');reject('corrupt_trace',lambda:check(job,p['preparation_seeds']['0'],bad,deadline))
   prefix=PREP/'prefix';prefix.mkdir();[shutil.copyfile(PREP/(job['id']+suffix),prefix/(job['id']+suffix)) for suffix in ('.json','.jsonl','.npz')];path=prefix/(job['id']+'.npz')
   with np.load(path) as z:counts=z['counts'][:16];seeds=z['seeds'][:16]
   np.savez_compressed(path,counts=counts,seeds=seeds,job=json.dumps(job));meta['sha256']=hashlib.sha256(path.read_bytes()).hexdigest();atomic_json(prefix/(job['id']+'.json'),meta);reject('prefix',lambda:check(job,p['preparation_seeds']['0'],prefix,deadline))
 held_score=score(held,required=32);atomic_json(PREP/'held_score.json',dict(pipeline_only=True,not_formal_result=True,result=held_score))
 try:score(held);raise ValueError('Incomplete accepted')
 except AssertionError:guards['incomplete']=True
 stress_timing={}
 for job in p['jobs']:
  if job['kind']=='stress':
   short=dict(job,updates=10000);start=time.perf_counter();stress(short,PREP,deadline);stress_check(short,PREP);stress_timing[job['case']]=time.perf_counter()-start
 forecast={j['id']:2*timings[str(j['cell'])]*128+3 if j['kind']=='endpoint' else 2*stress_timing[j['case']]*100+3 for j in p['jobs']};elapsed=time.perf_counter()-clock;assert elapsed<60 and hashes()==p['source_hashes'];atomic_json(BASE/'readiness.json',dict(passed=True,elapsed_seconds=elapsed,source_hashes=hashes(),plan_sha256=hashlib.sha256((PLAN).read_bytes()).hexdigest(),controls=control,guards=guards,endpoint_timing=timings,stress_timing=stress_timing,forecast=forecast,forecast_worker_seconds=sum(forecast.values()),formal_registered=False,formal_paths=0))
if __name__=='__main__':
 try:plan() if '--plan' in sys.argv else run()
 except BaseException:atomic_json(PREP/'failure.json',dict(error=traceback.format_exc()));raise

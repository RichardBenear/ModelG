"""MG147 approval-bound compact precision pilot pool."""
from mg12h3_finer_precision import *
from mg10e2c1_pool import run as pool
import argparse,subprocess,shutil

def register(approval):
 p=json.loads((BASE/'plan_v2.json').read_text());q=json.loads((BASE/'readiness.json').read_text());assert approval and q['passed'] and hashes()==q['source_hashes']==p['source_hashes'] and p['jobs']==jobs() and q['plan_sha256']==hashlib.sha256((BASE/'plan_v2.json').read_bytes()).hexdigest() and sum(q['forecast'].values())<15750 and max(q['forecast'].values())<600 and not (BASE/'registration.json').exists() and shutil.disk_usage(ROOT).free>3*1024**3;atomic_json(BASE/'registration.json',dict(approval=approval,jobs=jobs(),seeds=p['seeds'],source_hashes=hashes(),forecast=q['forecast'],plan_sha256=q['plan_sha256'],readiness_sha256=hashlib.sha256((BASE/'readiness.json').read_bytes()).hexdigest(),worker_cap=21000,parent_cap=600,total_cap=21600,wall_cap=3600,workers=8,job_timeout=600,scope='precision/cost pilot, no equivalence qualification',compact_archive=True,sentinel_replicate=0))
def worker(job,directory):
 r=json.loads((BASE/'registration.json').read_text());assert hashes()==r['source_hashes'] and job in r['jobs'];write(job,r['seeds'][job['id']],directory,time.perf_counter()+590)
def run():
 r=json.loads((BASE/'registration.json').read_text());assert hashes()==r['source_hashes'] and r['jobs']==jobs() and r['readiness_sha256']==hashlib.sha256((BASE/'readiness.json').read_bytes()).hexdigest() and r['plan_sha256']==hashlib.sha256((BASE/'plan_v2.json').read_bytes()).hexdigest() and not (BASE/'run_state.json').exists();subprocess.run(['git','ls-files','--error-unmatch',(BASE/'registration.json').relative_to(ROOT).as_posix()],check=True,cwd=ROOT,stdout=subprocess.DEVNULL);clock=time.perf_counter();atomic_json(BASE/'run_state.json',dict(status='RUNNING',started=datetime.now(timezone.utc).isoformat()));directory=BASE/'production';ledger=pool(jobs(),directory,Path(__file__).resolve(),ROOT,21000,workers=8,timeout=600,forecast=r['forecast']);auditstart=time.perf_counter();checks=[]
 for j in jobs():
  if time.perf_counter()-auditstart>600 or time.perf_counter()-clock>3600:raise TimeoutError('Pilot audit/wall cap')
  checks.append(audit(j,r['seeds'][j['id']],directory,auditstart+600))
 precision=[]
 for level in range(3):
  direct=[];est=[];power_values=[]
  for j in jobs():
   if j['level']==level:
    with np.load(directory/j['id']/'observations.npz') as z:direct.append(z['direct'][[5,10]]);est.append(z['estimates'][[5,10]]);power_values.append(z['power'][[5,10]])
  a=np.asarray(direct);b=np.asarray(est);v1=a.var(axis=0,ddof=1);v2=b.var(axis=0,ddof=1);precision.append(dict(level=level,replicates=8,direct_variance=v1.tolist(),integral_variance=v2.tolist(),variance_ratio=np.divide(v1,v2,out=np.full_like(v1,np.nan),where=v2>0).tolist(),direct_mean=a.mean(0).tolist(),integral_mean=b.mean(0).tolist(),power_mean=np.mean(power_values,axis=0).tolist(),power_variance=np.var(power_values,axis=0,ddof=1).tolist()))
 elapsed=time.perf_counter()-auditstart;assert elapsed<=600 and time.perf_counter()-clock<=3600 and elapsed+ledger['charged_worker_seconds']<=21600 and hashes()==r['source_hashes'];atomic_json(BASE/'evaluation.json',dict(integrity_passed=all(x['passed'] for x in checks),precision=precision,equivalence_verdict=None));atomic_json(BASE/'run_state.json',dict(status='RESULTS READY; PI closure pending',wall_seconds=time.perf_counter()-clock,worker_seconds=ledger['charged_worker_seconds'],parent_seconds=elapsed,total_seconds=elapsed+ledger['charged_worker_seconds'],finished=datetime.now(timezone.utc).isoformat()))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--register',action='store_true');p.add_argument('--approval');p.add_argument('--run',action='store_true');p.add_argument('--worker',nargs=2);a=p.parse_args()
 try:
  if a.worker:worker(json.loads(a.worker[0]),Path(a.worker[1]))
  elif a.register:register(a.approval)
  elif a.run:run()
 except BaseException:
  if a.worker:atomic_json(Path(a.worker[1])/(json.loads(a.worker[0])['id']+'.json'),dict(status='FAILED',job=json.loads(a.worker[0]),error=traceback.format_exc()))
  else:atomic_json(BASE/'failure.json',dict(error=traceback.format_exc()));atomic_json(BASE/'run_state.json',dict(status='STOPPED; no retry, partial data retained'))
  raise

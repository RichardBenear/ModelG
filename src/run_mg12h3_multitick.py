"""MG140 prospective actual-ring pilot controller."""
from mg12h3_multitick import *
from mg10e2c1_pool import run as pool
import argparse,subprocess,shutil

def register(approval):
 p=json.loads((BASE/'plan.json').read_text());q=json.loads((BASE/'readiness.json').read_text());assert approval and q['passed'] and not (BASE/'registration.json').exists() and hashes()==p['source_hashes']==q['source_hashes'] and p['jobs']==jobs() and sum(q['forecast'].values())<4950 and q['plan_sha256']==hashlib.sha256((BASE/'plan.json').read_bytes()).hexdigest() and shutil.disk_usage(ROOT).free>3*1024**3;atomic_json(BASE/'registration.json',dict(approval=approval,jobs=jobs(),seeds=p['seeds'],source_hashes=hashes(),forecast=q['forecast'],plan_sha256=q['plan_sha256'],readiness_sha256=hashlib.sha256((BASE/'readiness.json').read_bytes()).hexdigest(),worker_cap=6600,parent_cap=600,total_cap=7200,workers=8,threads_per_library=1,job_timeout=600,public_raw=True))
def worker(job,directory):
 r=json.loads((BASE/'registration.json').read_text());assert hashes()==r['source_hashes'] and job in r['jobs'];run_job(job,r['seeds'][job['id']],directory,time.perf_counter()+590)
def run():
 r=json.loads((BASE/'registration.json').read_text());assert r['jobs']==jobs() and hashes()==r['source_hashes'] and r['readiness_sha256']==hashlib.sha256((BASE/'readiness.json').read_bytes()).hexdigest() and r['plan_sha256']==hashlib.sha256((BASE/'plan.json').read_bytes()).hexdigest() and not (BASE/'run_state.json').exists();subprocess.run(['git','ls-files','--error-unmatch',(BASE/'registration.json').relative_to(ROOT).as_posix()],cwd=ROOT,check=True,stdout=subprocess.DEVNULL);clock=time.perf_counter();started=datetime.now(timezone.utc).isoformat();atomic_json(BASE/'run_state.json',dict(status='RUNNING',started=started));directory=BASE/'production';ledger=pool(jobs(),directory,Path(__file__).resolve(),ROOT,6600,workers=8,timeout=600,forecast=r['forecast']);auditstart=time.perf_counter();results=[]
 for job in jobs():
  if time.perf_counter()-auditstart>600:raise TimeoutError('Parent audit cap')
  results.append(dict(id=job['id'],audit=audit_job(job,r['seeds'][job['id']],directory)))
 parent=time.perf_counter()-auditstart;assert hashes()==r['source_hashes'] and parent<=600 and parent+ledger['charged_worker_seconds']<=7200;atomic_json(BASE/'evaluation.json',dict(feasibility_integrity_passed=all(x['audit']['passed'] for x in results),records=results,particle_verdict=None,refinement_verdict=None));atomic_json(BASE/'run_state.json',dict(status='RESULTS READY; PI closure pending',started=started,finished=datetime.now(timezone.utc).isoformat(),wall_seconds=time.perf_counter()-clock,worker_seconds=ledger['charged_worker_seconds'],parent_seconds=parent,total_seconds=parent+ledger['charged_worker_seconds'],stop_reason=ledger['stop_reason']))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--register',action='store_true');p.add_argument('--approval');p.add_argument('--run',action='store_true');p.add_argument('--worker',nargs=2);a=p.parse_args()
 try:
  if a.worker:worker(json.loads(a.worker[0]),Path(a.worker[1]))
  elif a.register:register(a.approval)
  elif a.run:run()
 except BaseException:
  if a.worker:atomic_json(Path(a.worker[1])/(json.loads(a.worker[0])['id']+'.json'),dict(status='FAILED',job=json.loads(a.worker[0]),error=traceback.format_exc()))
  else:atomic_json(BASE/'failure.json',dict(error=traceback.format_exc()));atomic_json(BASE/'run_state.json',dict(status='STOPPED; no retry, partial records retained'))
  raise

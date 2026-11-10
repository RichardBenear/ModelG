"""MG137 approval-bound sparse endpoint/stress pool driver."""
from mg12h3_sparse_qualification import *
from mg10e2c1_pool import run as pool
import argparse,subprocess,shutil

def register(approval):
 p=json.loads((BASE/'plan.json').read_text());q=json.loads((BASE/'readiness.json').read_text());assert approval and not (BASE/'registration.json').exists() and p['jobs']==jobs() and q['passed'] and hashes()==q['source_hashes']==p['source_hashes'] and q['plan_sha256']==hashlib.sha256((BASE/'plan.json').read_bytes()).hexdigest() and sum(q['forecast'].values())<4950 and shutil.disk_usage(ROOT).free>2*1024**3
 atomic_json(BASE/'registration.json',dict(approval=approval,source_hashes=hashes(),jobs=jobs(),seeds=p['formal_seeds'],forecast=q['forecast'],plan_sha256=q['plan_sha256'],readiness_sha256=hashlib.sha256((BASE/'readiness.json').read_bytes()).hexdigest(),worker_cap=6600,parent_cap=600,total_cap=7200,workers=8,timeout=600,endpoint_margin=.1,endpoint_alpha=.001,public_raw=True))

def worker(job,directory):
 r=json.loads((BASE/'registration.json').read_text());assert hashes()==r['source_hashes'] and job in r['jobs'];deadline=time.perf_counter()+590
 if job['kind']=='endpoint':write_endpoint(job,r['seeds'][job['id']],directory,deadline)
 else:stress(job,directory,deadline)

def run():
 r=json.loads((BASE/'registration.json').read_text());assert hashes()==r['source_hashes'] and r['jobs']==jobs() and not (BASE/'run_state.json').exists() and r['readiness_sha256']==hashlib.sha256((BASE/'readiness.json').read_bytes()).hexdigest() and r['plan_sha256']==hashlib.sha256((BASE/'plan.json').read_bytes()).hexdigest();subprocess.run(['git','ls-files','--error-unmatch',(BASE/'registration.json').relative_to(ROOT).as_posix()],check=True,cwd=ROOT,stdout=subprocess.DEVNULL);atomic_json(BASE/'run_state.json',dict(status='RUNNING'));directory=BASE/'production';ledger=pool(r['jobs'],directory,Path(__file__).resolve(),ROOT,6600,workers=8,timeout=600,forecast=r['forecast']);clock=time.perf_counter();arrays={c:[] for c in range(4)};stress_results=[]
 for job in r['jobs']:
  if time.perf_counter()-clock>600:raise TimeoutError('Parent audit cap')
  if job['kind']=='endpoint':arrays[job['cell']].append(check(job,r['seeds'][job['id']],directory,clock+600))
  else:stress_results.append(stress_check(job,directory))
 result=score(arrays);result.update(stress_passed=all(x['passed'] for x in stress_results),stress_results=stress_results);result['combined_qualified']=result['qualified'] and result['stress_passed'];atomic_json(BASE/'evaluation.json',result);elapsed=time.perf_counter()-clock;assert hashes()==r['source_hashes'] and elapsed<=600 and elapsed+ledger['charged_worker_seconds']<=7200;atomic_json(BASE/'run_state.json',dict(status='RESULTS READY; PI closure pending',combined_qualified=result['combined_qualified'],worker_seconds=ledger['charged_worker_seconds'],parent_seconds=elapsed,total_seconds=elapsed+ledger['charged_worker_seconds'],stop_reason=ledger['stop_reason']))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--register',action='store_true');p.add_argument('--approval');p.add_argument('--run',action='store_true');p.add_argument('--worker',nargs=2);a=p.parse_args()
 try:
  if a.worker:worker(json.loads(a.worker[0]),Path(a.worker[1]))
  elif a.register:register(a.approval)
  elif a.run:run()
 except BaseException:
  if a.worker:atomic_json(Path(a.worker[1])/(json.loads(a.worker[0])['id']+'.json'),dict(status='FAILED',error=traceback.format_exc()))
  else:atomic_json(BASE/'failure.json',dict(error=traceback.format_exc()));atomic_json(BASE/'run_state.json',dict(status='STOPPED; no retry, partial data preserved'))
  raise

"""MG152 prospective fresh fixed-blend validation using frozen raw writer."""
from mg12h3_Y_blend import *
from mg10e2c1_pool import run as pool
import argparse,subprocess,shutil

def register(approval):
 p=json.loads((BASE/'plan.json').read_text());q=json.loads((BASE/'readiness.json').read_text());assert approval and q['passed'] and hashes()==q['source_hashes']==p['source_hashes'] and q['alpha_direct']==ALPHA.tolist()==p['alpha_direct'] and p['jobs']==core.jobs() and q['plan_sha256']==hashlib.sha256((BASE/'plan.json').read_bytes()).hexdigest() and q['bootstrap_sha256']==hashlib.sha256((BASE/'preparation/bootstrap_indices.npz').read_bytes()).hexdigest() and sum(p['forecast'].values())<10350 and max(p['forecast'].values())<600 and shutil.disk_usage(ROOT).free>3*1024**3 and not (BASE/'registration.json').exists();atomic_json(BASE/'registration.json',dict(approval=approval,jobs=core.jobs(),seeds=p['seeds'],source_hashes=hashes(),alpha_direct=ALPHA.tolist(),forecast=p['forecast'],plan_sha256=q['plan_sha256'],bootstrap_sha256=q['bootstrap_sha256'],readiness_sha256=hashlib.sha256((BASE/'readiness.json').read_bytes()).hexdigest(),worker_cap=13800,parent_cap=600,total_cap=14400,workers=8,timeout=600,wall_checks=3600,scope='Fresh fixed-estimator variance validation, no tick equivalence verdict',compact_archive=True,sentinel_replicate=0))
def worker(job,directory):
 r=json.loads((BASE/'registration.json').read_text());assert hashes()==r['source_hashes'] and r['alpha_direct']==ALPHA.tolist() and job in r['jobs'];core.write(job,r['seeds'][job['id']],directory,time.perf_counter()+590)
def run():
 r=json.loads((BASE/'registration.json').read_text());assert r['jobs']==core.jobs() and hashes()==r['source_hashes'] and r['readiness_sha256']==hashlib.sha256((BASE/'readiness.json').read_bytes()).hexdigest() and r['plan_sha256']==hashlib.sha256((BASE/'plan.json').read_bytes()).hexdigest() and r['bootstrap_sha256']==hashlib.sha256((BASE/'preparation/bootstrap_indices.npz').read_bytes()).hexdigest() and not (BASE/'run_state.json').exists();subprocess.run(['git','ls-files','--error-unmatch',(BASE/'registration.json').relative_to(ROOT).as_posix()],check=True,cwd=ROOT,stdout=subprocess.DEVNULL);clock=time.perf_counter();atomic_json(BASE/'run_state.json',dict(status='RUNNING',started=core.datetime.now(core.timezone.utc).isoformat()));directory=BASE/'production';ledger=pool(core.jobs(),directory,Path(__file__).resolve(),ROOT,13800,workers=8,timeout=600,forecast=r['forecast']);auditstart=time.perf_counter();checks=[]
 for job in core.jobs():
  if time.perf_counter()-auditstart>600 or time.perf_counter()-clock>3600:raise TimeoutError('Audit/wall budget')
  checks.append(core.audit(job,r['seeds'][job['id']],directory,auditstart+600))
 direct=np.empty((3,8,2));integrated=direct.copy()
 for job in core.jobs():
  with np.load(directory/job['id']/'observations.npz') as z:direct[job['level'],job['replicate']]=z['direct'][[5,10],2];integrated[job['level'],job['replicate']]=z['estimates'][[5,10],2]
 with np.load(BASE/'preparation/bootstrap_indices.npz') as z:indices=z['indices']
 early,early_boot=validate_variance(direct[:,:,0],integrated[:,:,0],indices,0);late,late_boot=validate_variance(direct[:,:,1],integrated[:,:,1],indices,1);mixed=np.stack([blend(direct[:,:,i],integrated[:,:,i],i) for i in range(2)],axis=-1);np.savez_compressed(BASE/'per_path_Y.npz',direct=direct,integrated=integrated,blend=mixed,times=[.005,.01],alpha=ALPHA);np.savez_compressed(BASE/'bootstrap_ratios.npz',early=early_boot,late=late_boot);elapsed=time.perf_counter()-auditstart;assert hashes()==r['source_hashes'] and elapsed<=600 and elapsed+ledger['charged_worker_seconds']<=14400 and time.perf_counter()-clock<=3600;atomic_json(BASE/'evaluation.json',dict(integrity_passed=all(x['passed'] for x in checks),primary_late=late,secondary_early=early,equivalence_verdict=None,weights_refitted=False));atomic_json(BASE/'run_state.json',dict(status='RESULTS READY; PI closure pending',wall_seconds=time.perf_counter()-clock,worker_seconds=ledger['charged_worker_seconds'],parent_seconds=elapsed,total_seconds=elapsed+ledger['charged_worker_seconds'],finished=core.datetime.now(core.timezone.utc).isoformat()))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--register',action='store_true');p.add_argument('--approval');p.add_argument('--run',action='store_true');p.add_argument('--worker',nargs=2);a=p.parse_args()
 try:
  if a.worker:worker(json.loads(a.worker[0]),Path(a.worker[1]))
  elif a.register:register(a.approval)
  elif a.run:run()
 except BaseException:
  if a.worker:atomic_json(Path(a.worker[1])/(json.loads(a.worker[0])['id']+'.json'),dict(status='FAILED',job=json.loads(a.worker[0]),error=traceback.format_exc()))
  else:atomic_json(BASE/'failure.json',dict(error=traceback.format_exc()));atomic_json(BASE/'run_state.json',dict(status='STOPPED; partial outputs preserved, no retry'))
  raise

"""MG130 prospective qualification controller."""
from mg12h3_efficiency_qualification import *
from prepare_mg12h3_efficiency import hashes as timing_hashes
from mg10e2c1_pool import run as pool
import os,sys,argparse,subprocess,traceback
FILES=['src/mg12h3_efficiency_qualification.py','src/run_mg12h3_efficiency_qualification.py','spec/MG12H3_efficiency_qualification.md']
def hashes():
    h=timing_hashes();h.update({p:hashlib.sha256((ROOT/p).read_text(encoding='utf-8').replace('\r\n','\n').encode()).hexdigest() for p in FILES})
    for c in range(4):
        p=f'results/MG12H3/mixed_qualification/references/c{c}_reference.npz';h[p]=hashlib.sha256((ROOT/p).read_bytes()).hexdigest()
    return h

def register(approval):
    q=json.loads((BASE/'readiness.json').read_text());assert approval and q['passed'] and q['source_hashes']==hashes() and not (BASE/'registration.json').exists()
    atomic_json(BASE/'registration.json',dict(approval=approval,jobs=jobs(),seeds=q['seeds'],source_hashes=hashes(),worker_cap=4200,parent_cap=300,total_cap=4500,workers=8,timeout=600,forecast=q['forecast'],margin=.1,alpha=.001))
def worker(job,directory):
    r=json.loads((BASE/'registration.json').read_text());assert hashes()==r['source_hashes'] and job in r['jobs'];write(job,r['seeds'][job['id']],directory,time.perf_counter()+590)
def run():
    r=json.loads((BASE/'registration.json').read_text());assert hashes()==r['source_hashes'] and not (BASE/'run_state.json').exists();subprocess.run(['git','ls-files','--error-unmatch',str((BASE/'registration.json').relative_to(ROOT))],check=True,cwd=ROOT,stdout=subprocess.DEVNULL);atomic_json(BASE/'run_state.json',dict(status='RUNNING'));directory=BASE/'production';ledger=pool(r['jobs'],directory,Path(__file__).resolve(),ROOT,4200,workers=8,timeout=600,forecast=r['forecast']);clock=time.perf_counter();arrays={c:[] for c in range(4)}
    for job in r['jobs']:arrays[job['cell']].append(check(job,r['seeds'][job['id']],directory,clock+300))
    result=score(arrays);assert hashes()==r['source_hashes'];atomic_json(BASE/'evaluation.json',result);elapsed=time.perf_counter()-clock;assert elapsed<=300 and ledger['charged_worker_seconds']+elapsed<=4500;atomic_json(BASE/'run_state.json',dict(status='RESULTS READY; PI closure pending',qualified=result['qualified'],worker_seconds=ledger['charged_worker_seconds'],parent_seconds=elapsed,total_seconds=ledger['charged_worker_seconds']+elapsed))
if __name__=='__main__':
    from pathlib import Path
    p=argparse.ArgumentParser();p.add_argument('--register',action='store_true');p.add_argument('--approval');p.add_argument('--run',action='store_true');p.add_argument('--worker',nargs=2);a=p.parse_args()
    try:
        if a.worker:worker(json.loads(a.worker[0]),Path(a.worker[1]))
        elif a.register:register(a.approval)
        elif a.run:run()
    except BaseException:
        if a.worker:atomic_json(Path(a.worker[1])/(json.loads(a.worker[0])['id']+'.json'),dict(status='FAILED',error=traceback.format_exc()))
        else:atomic_json(BASE/'failure.json',dict(error=traceback.format_exc()))
        raise

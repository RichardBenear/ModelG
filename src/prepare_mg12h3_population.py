"""MG115 prospective T1 preparation; full pattern reference run approval gated."""
from mg12h3_population import *
import json,hashlib,argparse,subprocess,traceback
from mg11g2 import initial as old_initial,rhs as common_rhs
from mg11g3 import rhs as separate_rhs
from mg10d1_pool import atomic_json
BASE=ROOT/'results/MG12H3/population';PREP=BASE/'preparation'
FILES=('src/mg12h3_population.py','src/prepare_mg12h3_population.py','src/run_mg12h3_population.py','src/check_mg12h3_population.py','spec/MG12H3_pattern_bridge.md','src/mg11g2.py','src/mg11g3.py','src/mg11g2_pattern_v2.py','src/mg11g2_pattern_v2_oracle.py','src/mg10p1.py','src/mg11n1.py')
def read(p):return json.loads(Path(p).read_text(encoding='utf-8'))
def hashes():
    names=set(FILES)
    for j in jobs():
        case,ga,gb,k,mode,reference,parent=factory(j);names.add(reference.relative_to(ROOT).as_posix())
        if parent:names.add(parent.relative_to(ROOT).as_posix())
    names.add('results/MG003/ref_profiles.npz')
    return {name:hashlib.sha256((ROOT/name).read_bytes() if name.endswith('.npz') else (ROOT/name).read_text(encoding='utf-8').replace('\r\n','\n').encode()).hexdigest() for name in sorted(names)}
def save(path,job,times,state):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True);tmp=path.with_suffix('.npz.tmp')
    with tmp.open('wb') as f:np.savez_compressed(f,state=state,times=times,job_json=json.dumps(job))
    os.replace(tmp,path)
def plan():
    if (PREP/'plan.json').exists():raise ValueError('Already planned')
    atomic_json(PREP/'plan.json',dict(stage='MG115',scope='T1preparation only; full deterministic reference run awaits concrete approval',source_hashes=hashes(),jobs=jobs(),times=(np.arange(5)/4).tolist(),ceiling_elapsed_seconds=300,workers=1,threads_per_library=1,candidate_full_cap_seconds=3600))
def run():
    p=read(PREP/'plan.json')
    if p['source_hashes']!=hashes() or (PREP/'attempt.json').exists():raise ValueError('Changed source or existing attempt')
    subprocess.run(['git','ls-files','--error-unmatch','results/MG12H3/population/preparation/plan.json'],cwd=ROOT,check=True,stdout=subprocess.DEVNULL);start=time.perf_counter();deadline=start+300;atomic_json(PREP/'attempt.json',dict(status='RUNNING'));records=[];forecast={};times=np.array(p['times']);static=[]
    for job in p['jobs']:
        clock=time.perf_counter();case,ga,gb,k,mode,reference,parent=factory(job);n=ga['n'];y=np.r_[case['initial'].ravel(),np.zeros(2*n)];new=rhs(y,case,matrices(case,ga,gb),k,mode);old=common_rhs(0,old_initial(case,ga),case,ga,k,mode) if job['kind']!='FINITE' else separate_rhs(0,old_initial(case,ga),case,ga,gb,k);np.testing.assert_allclose(new,np.r_[old[:7*n],old[-2*n:]],rtol=1e-12,atol=1e-9);static.append(dict(id=job['id'],max_initial_RHS_difference=float(max(abs(new-np.r_[old[:7*n],old[-2*n:]])))));io=[0.]
        def checkpoint(t,z):
            tick=time.perf_counter();save(PREP/(job['id']+'.npz'),job,t,z);io[0]+=time.perf_counter()-tick
        result=solve(job,times,deadline,checkpoint);elapsed=time.perf_counter()-clock;checked=check(job,result,times);assert checked['passed'];fullframes=int(job['duration']*4)+1;compute=max(0,elapsed-io[0]);estimate=1.5*(compute*job['duration']+io[0]*fullframes*(fullframes+1)/30)+2;forecast[job['id']]=estimate;records.append(dict(job=job,elapsed_seconds=elapsed,io_seconds=io[0],check=checked,forecast_seconds=estimate))
    # Test invalid balance/negative fields with held preparation data, no new integrations.
    j=jobs()[0]
    with np.load(PREP/(j['id']+'.npz')) as z:held=z['state'].copy()
    bad=held.copy();bad[-1,5*500]+=1e8;assert not check(j,bad,times)['passed'];bad=held.copy();bad[-1,6*500]=-1;assert not check(j,bad,times)['passed']
    from check_mg12h3_population import evaluate
    incomplete=evaluate(PREP,PREP/'prefix_evaluation.json',full=False);assert not incomplete['complete'] and not incomplete['qualified']
    empty=PREP/'empty';empty.mkdir(exist_ok=True);missing=evaluate(empty,PREP/'missing_evaluation.json',full=False);assert not missing['complete'] and missing['completed']==0
    elapsed=time.perf_counter()-start;assert elapsed<300 and hashes()==p['source_hashes'];atomic_json(PREP/'results.json',dict(passed=True,source_hashes=hashes(),jobs=records,static=static,forecast=forecast,forecast_seconds=sum(forecast.values())+30,elapsed_seconds=elapsed,invalid_balance_rejected=True,negative_rejected=True,prefix_suppressed=True,missing_suppressed=True,scope='T1only; full late trajectories not run'));atomic_json(PREP/'attempt.json',dict(status='COMPLETE'))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--plan',action='store_true');p.add_argument('--run',action='store_true');a=p.parse_args()
    if a.plan:plan()
    elif a.run:
        try:run()
        except BaseException:atomic_json(PREP/'failure.json',dict(error=traceback.format_exc()));raise
    else:p.error('Choose action')

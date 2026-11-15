import os
for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS'):os.environ[k]='1'
import time,json,hashlib,gzip,collections,traceback
import numpy as np
from scipy.integrate import solve_ivp
from mg182_gscale import setup,projection,rhs,sparsity,PATH,ROOT
from mg12h3_sparse import tick
from mg12h3_tick_coupling import pair,replay
from run_mg176 import read,atomic_json
OUT=ROOT/'results/MG12H3/G_population_development';RAW=ROOT/'local_raw/MG182_development'
FILES=['src/mg182_gscale.py','src/prepare_mg182.py','src/mg12h3_sparse.py','src/mg12h3_tick_coupling.py','src/mg11n1.py','spec/MG182_G_population_development.md']
def binding():return {x:hashlib.sha256((ROOT/x).read_bytes()).hexdigest() for x in FILES}
def register():
 OUT.mkdir(exist_ok=True);assert not(OUT/'plan.json').exists();atomic_json(OUT/'plan.json',dict(source=binding(),inputs={str(PATH.relative_to(ROOT)):hashlib.sha256(PATH.read_bytes()).hexdigest()},cap_seconds=180,workers=1,development_only=True,production_approved=False,omegas=[10,40,160,640],benchmark_seed=18200))
def run():
 plan=read(OUT/'plan.json');assert plan['source']==binding() and not(OUT/'attempt.json').exists();assert hashlib.sha256(PATH.read_bytes()).hexdigest()==plan['inputs'][str(PATH.relative_to(ROOT))];RAW.mkdir(exist_ok=False);start=time.perf_counter();deadline=start+180;atomic_json(OUT/'attempt.json',dict(status='RUNNING'));rows=[];bench=[]
 try:
  cases={};references=[]
  for omega in plan['omegas']:
   t=time.perf_counter();c=setup(omega);setupseconds=time.perf_counter()-t;cases[omega]=c;pop=c['pop'];z=pop/omega;drift=np.zeros(3500);byrule=collections.defaultdict(float)
   for ch in c['model'].channels:
    rate=ch.rate(pop.ravel());byrule[ch.rule]+=rate
    for i,v in ch.change:drift[i]+=rate*v/omega
   finite=rhs(z.ravel(),c,omega);np.testing.assert_allclose(drift,finite,rtol=1e-8,atol=3e-8);continuum=rhs(z.ravel(),c);correction=np.zeros((7,500));correction[1]=-z[1]*z[2]/omega;correction[2]=-correction[1];np.testing.assert_allclose(drift,continuum+correction.ravel(),rtol=1e-8,atol=3e-8);np.testing.assert_allclose(drift.reshape(7,500)[0],(14*z[3]/c['a']-z[0]+.1*z[1]+c['diff'][0]*(np.roll(z[0],1)+np.roll(z[0],-1)-2*z[0])),rtol=1e-9,atol=1e-9)
   if omega==10:
    with np.load(PATH) as original:np.testing.assert_array_equal(pop,original['counts'])
   initial=projection(z,c['reference']);ref=projection(c['reference'],c['reference']);rounderror=float(np.linalg.norm((z[0]-z[0].mean())-(c['reference'][0]-c['reference'][0].mean())));rows.append(dict(omega=omega,setup_seconds=setupseconds,initial_projection=initial,unrounded_projection=ref,centered_G_quantization_error=rounderror,relative_G_quantization_error=rounderror/ref['G_reference_norm'],max_finite_drift_error=float(np.max(abs(drift-finite))),hazards=dict(byrule),total_hazard=float(sum(byrule.values())),G_count_range=[int(pop[0].min()),int(pop[0].max())],total_tokens=int(pop.sum())));np.savez_compressed(RAW/f'inputs_Omega{omega}.npz',counts=pop,reference=c['reference'],omega=omega,a=c['a'],b=c['b'])
  c=cases[10];t=time.perf_counter()
  def bounded(t,y):
   if time.perf_counter()>deadline:raise TimeoutError('Development ODEdeadline')
   return rhs(y,c)
  ode=solve_ivp(bounded,(0,.1),c['reference'].ravel(),t_eval=[0,.005,.01,.1],method='BDF',jac_sparsity=sparsity(c['model']),rtol=1e-9,atol=1e-11,max_step=.005);assert ode.success and ode.y.shape==(3500,4);np.savez_compressed(RAW/'continuum_reference.npz',times=ode.t,fields=ode.y.T.reshape(-1,7,500));refs=[dict(time=float(t),**projection(y.reshape(7,500),c['reference'])) for t,y in zip(ode.t,ode.y.T)];atomic_json(OUT/'scaling_checks.json',dict(passed=True,rows=rows,continuum_reference=refs,ode_seconds=time.perf_counter()-t,ode_nfev=ode.nfev,no_stationarity_claim=True))
  jobs=[(omega,1e-6,False) for omega in plan['omegas']]+[(10,1e-4,False),(160,1e-6,True)]
  for index,(omega,duration,paired) in enumerate(jobs):
   c=cases[omega];seed=18200+index;rng=np.random.default_rng(seed);before=json.dumps(rng.bit_generator.state);t=time.perf_counter()
   if paired:
    out=pair(c['model'],c['pop'],c['pop'],[duration,duration/5],[0,duration],seed,deadline=deadline,rng=rng,sentinel=True);replay(c['model'],c['pop'],c['pop'],[duration,duration/5],[0,duration],out);final=out['counts'][:,-1];atomic_json(RAW/f'benchmark_{index}_paired.json',out)
   else:
    out=tick(c['model'],c['pop'],duration,rng,deadline=deadline,trace=True);a=c['pop'].copy().ravel();end=a.copy();hist=np.zeros(len(c['model'].channels),np.int64)
    for event in out['trace']:
     ch=c['model'].channels[event['channel']];assert ch.rate(a)>0
     for i,n in ch.reserve:assert a[i]>=n;a[i]-=n
     for i,v in ch.change:end[i]+=v
     hist[event['channel']]+=1
    np.testing.assert_array_equal(a,out['available'].ravel());np.testing.assert_array_equal(end,out['counts'].ravel());assert end.sum()==c['pop'].sum();final=out['counts'];np.savez_compressed(RAW/f'benchmark_{index}.npz',initial=c['pop'],final=final,available=out['available'],histogram=hist,ids=[e['channel'] for e in out['trace']],marks=[e['time'] for e in out['trace']],rng_before=before,rng_after=json.dumps(rng.bit_generator.state))
   elapsed=time.perf_counter()-t;bench.append(dict(index=index,omega=omega,duration=duration,paired=paired,seed=seed,events=out['events'],seconds_including_replay=elapsed,events_per_second=out['events']/elapsed,model_time_per_second=duration/elapsed,not_scientific_sample=True));atomic_json(OUT/'benchmark.json',dict(rows=bench))
  atomic_json(OUT/'attempt.json',dict(status='COMPLETE',elapsed_seconds=time.perf_counter()-start,production_started=False));atomic_json(OUT/'local_raw_manifest.json',dict(directory=str(RAW),files={p.name:dict(bytes=p.stat().st_size,sha256=hashlib.sha256(p.read_bytes()).hexdigest()) for p in RAW.iterdir()},full_raw_published=False))
 except BaseException:atomic_json(OUT/'attempt.json',dict(status='FAILED',error=traceback.format_exc(),elapsed_seconds=time.perf_counter()-start,production_started=False));raise
if __name__=='__main__':register() if __import__('sys').argv[1]=='register' else run()

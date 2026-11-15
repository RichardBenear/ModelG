import os
for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS'):os.environ[k]='1'
import time,json,hashlib,traceback
import numpy as np
from pathlib import Path
from lattice import P17,tau_leap_euler_1d
from measure import analyze_run
from mg184_phase2 import ref,initial,advance,summarize,OMEGAS,TAU,H
from mg184_open_bridge import OpenBridge,generator
from mg10e2c1_pool import atomic_json
from run_mg176 import ROOT,read
OUT=ROOT/'results/MG12H3/phase2_G_development';RAW=ROOT/'local_raw/MG184_development';REFERENCE=ROOT/'results/MG003/ref_profiles.npz'
FILES=['src/mg184_phase2.py','src/mg184_open_bridge.py','src/prepare_mg184.py','src/lattice.py','src/measure.py','src/runner.py','spec/MG184_phase2_development.md']
def binding():return {p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in FILES}
def register():
 OUT.mkdir(exist_ok=True);assert not(OUT/'plan.json').exists();atomic_json(OUT/'plan.json',dict(source=binding(),reference_sha256=hashlib.sha256(REFERENCE.read_bytes()).hexdigest(),cap_seconds=240,workers=1,benchmark_seeds={str(o):18400+i for i,o in enumerate(OMEGAS)},production_started=False))
def run():
 start=time.perf_counter();plan=read(OUT/'plan.json');assert plan['source']==binding() and not(OUT/'attempt.json').exists();assert hashlib.sha256(REFERENCE.read_bytes()).hexdigest()==plan['reference_sha256'];RAW.mkdir(exist_ok=False);atomic_json(OUT/'attempt.json',dict(status='RUNNING'));checks=[]
 def test(name,condition):
  assert condition,name;checks.append(name)
 try:
  r=ref();pop=initial(10000,'prepared',0,r);times=np.arange(2001)*.5;full=np.broadcast_to(pop,(2001,3,500)).copy();base,scores=summarize(full,times,10000,'prepared',1000,r);test('constant reference passes full horizon',base['historical_X_persistence']);test('constant sampled track unique',base['strict_sampled_unique_X_track']);homogeneous=np.rint(10000*r['S0'][:,None])*np.ones((3,500));late=full.copy();late[times>=500]=homogeneous;legacy=analyze_run(late[:,0],late[:,1],late[:,2],times,10000,H,r['S0'],r['A_s'],r['lam'],T=100);correct,late_scores=summarize(late,times,10000,'prepared',1000,r);test('old T100 would miss late loss',legacy['eligible']);test('explicit T1000 detects late loss',not correct['historical_X_persistence']);test('unmatched G is NaN',np.isnan(late_scores['aligned_projection'][times>=530]).all());missing=full.copy();missing[:,0]=homogeneous[0];m,ms=summarize(missing,times,10000,'prepared',1000,r);test('X persistence not full G evidence',m['historical_X_persistence'] and np.max(abs(ms['aligned_projection'][:,0]))<1e-10 and m['full_GXY_persistence']=='NOT_SCORED_FROM_X_ALONE');shifted=np.roll(full,1,axis=2);sm,ss=summarize(shifted,times,10000,'prepared',1000,r);np.testing.assert_allclose(ss['aligned_projection'],scores['aligned_projection'],rtol=1e-8,atol=1e-9);checks.append('common-alignment translation invariant');tail=full.copy();tail[times>=995]=homogeneous;tm,ts=summarize(tail,times,10000,'prepared',1000,r);test('right-tail ambiguity explicit',tm['end_unmatched_or_ambiguous'] and not tm['strict_sampled_unique_X_track']);bridge=[]
  for omega in (100,160,3000):
   n=8;p=np.zeros((7,n),np.int64);p[:3]=np.rint(omega*np.array([np.linspace(14,16,n),np.linspace(10,20,n),np.linspace(1,3,n)]));model=OpenBridge(n,omega);drift=np.zeros_like(p,float)
   for ch in model.channels:
    rate=ch.rate(p.ravel())
    for pos,v in ch.change:drift.ravel()[pos]+=v*rate/omega
   expected=generator(p,omega);np.testing.assert_allclose(drift[:3],expected,rtol=1e-10,atol=1e-10);balance=P17['a']*n-(P17['q']-P17['p'])*np.sum(p[0]/omega)-(1-P17['g'])*np.sum(p[1]/omega);np.testing.assert_allclose(drift[:3].sum(),balance,atol=1e-10);bridge.append(dict(omega=omega,max_generator_error=float(np.max(abs(drift[:3]-expected))),external_active_balance=float(balance),finite_step_equivalence_claim=False))
  pg=TAU*(P17['q']+2/H**2);px=TAU*(1+P17['b']+2*P17['dx']/H**2);test('fixed G/X exit probabilities valid',pg<1 and px<1);atomic_json(OUT/'static_checks.json',dict(passed=True,checks=checks,bridge=bridge,reference_G_norm=float(np.linalg.norm(r['phi'][0]-r['phi'][0].mean())),G_exit_probability=pg,X_exit_probability=px,late_loss_events=correct['legacy']['events'],tail_legacy_eligible=tm['historical_X_persistence'],full_GXY_claim=False));del full,late,missing,shifted,tail
  warm=time.perf_counter();a=np.array([[16]*8,[16]*8,[2]*8],np.int64);tau_leap_euler_1d(*a,100.,14.,29.,.1,1.,1.,1.,12.,.2,TAU,1,0.,0.,np.arange(8)*.2,1,18499);warm_seconds=time.perf_counter()-warm;rows=[]
  for omega in OMEGAS:
   p=initial(omega,'prepared',18400,r);t=time.perf_counter();o=advance(p,omega,2.,0.,plan['benchmark_seeds'][str(omega)]);elapsed=time.perf_counter()-t;assert o['counts'].min()>=0 and np.array_equal(o['final'],o['counts'][-1]);np.savez_compressed(RAW/f'Omega{omega}_benchmark.npz',initial=p,counts=o['counts'],times=o['times'],final=o['final']);rows.append(dict(omega=omega,seed=plan['benchmark_seeds'][str(omega)],duration=2.,seconds=elapsed,seconds_per_model_unit=elapsed/2,clips=o['clips'],T1000_projected_seconds=elapsed*500,not_a_production_replicate=True));atomic_json(OUT/'benchmark.json',dict(rows=rows,JIT_warmup_seconds=warm_seconds));assert time.perf_counter()-start<240
  assert hashlib.sha256(REFERENCE.read_bytes()).hexdigest()==plan['reference_sha256'];atomic_json(OUT/'attempt.json',dict(status='COMPLETE',elapsed_seconds=time.perf_counter()-start,production_started=False))
 except BaseException:atomic_json(OUT/'attempt.json',dict(status='FAILED',elapsed_seconds=time.perf_counter()-start,error=traceback.format_exc()));raise
if __name__=='__main__':register() if __import__('sys').argv[1]=='register' else run()

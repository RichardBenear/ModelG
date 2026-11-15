import os
for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS'):os.environ[k]='1'
import time,json,hashlib,traceback
import numpy as np
from mg182_gscale import setup,ROOT
from prepare_mg182 import OUT,RAW,read,atomic_json
from mg12h3_tick_coupling import pair

def register():
 assert not(OUT/'full_tick_plan.json').exists();source=ROOT/'src/benchmark_mg182_fulltick.py';atomic_json(OUT/'full_tick_plan.json',dict(source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),omega=160,deltas=[.0001,.00002],duration=.0001,seed=18206,cap_seconds=150,max_events=2000000,development_timing_only=True,production_started=False,previous_microbenchmark_complete=read(OUT/'attempt.json')['status']=='COMPLETE'))
def run():
 p=read(OUT/'full_tick_plan.json');assert not(OUT/'full_tick_attempt.json').exists();assert hashlib.sha256((ROOT/'src/benchmark_mg182_fulltick.py').read_bytes()).hexdigest()==p['source_sha256'];start=time.perf_counter();atomic_json(OUT/'full_tick_attempt.json',dict(status='RUNNING'))
 try:
  c=setup(160);o=pair(c['model'],c['pop'],c['pop'],p['deltas'],[0,.0001],p['seed'],deadline=start+150,sentinel=False,max_events=2000000);hist=o['histogram'];end=o['counts'][:,-1]
  for side in range(2):
   direct=c['pop'].ravel().copy()
   for i,ch in enumerate(c['model'].channels):
    for coordinate,v in ch.change:direct[coordinate]+=int(hist[0,i]+hist[side+1,i])*v
   np.testing.assert_array_equal(direct,end[side].ravel());assert direct.min()>=0 and direct.sum()==c['pop'].sum()
  np.savez_compressed(RAW/'full_tick_probe.npz',initial=c['pop'],final=end,histogram=hist,available=o['final_available_before_reset'],rng_state=o['rng_state']);elapsed=time.perf_counter()-start;atomic_json(OUT/'full_tick_result.json',dict(passed=True,omega=160,deltas=p['deltas'],duration=.0001,events=o['events'],seconds=elapsed,projected_seconds_T0_01=elapsed*100,projected_seconds_T0_1=elapsed*1000,not_a_scientific_sample=True,aggregate_ledger_checked=True,no_full_event_trace=True,event_cap_raised_only_as_computational_bound=True));atomic_json(OUT/'full_tick_attempt.json',dict(status='COMPLETE',elapsed_seconds=elapsed))
 except BaseException:atomic_json(OUT/'full_tick_attempt.json',dict(status='FAILED',elapsed_seconds=time.perf_counter()-start,error=traceback.format_exc()));raise
if __name__=='__main__':register() if __import__('sys').argv[1]=='register' else run()

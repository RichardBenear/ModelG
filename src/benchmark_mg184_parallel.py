import os
for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS'):os.environ[k]='1'
import time,json,hashlib,traceback,multiprocessing
from concurrent.futures import ProcessPoolExecutor,as_completed
from pathlib import Path
import numpy as np
from mg184_phase2 import initial,advance
from prepare_mg184 import ROOT,OUT,RAW,read,atomic_json

def worker(rep):
 start=time.perf_counter();p=initial(3000,'prepared',18420+rep);o=advance(p,3000,5.,0.,18420+rep);np.savez_compressed(RAW/f'parallel_probe_r{rep}.npz',initial=p,counts=o['counts'],times=o['times']);return dict(replicate=rep,seed=18420+rep,omega=3000,duration=5.,seconds=time.perf_counter()-start,clips=o['clips'],not_production=True)
def register():
 assert not(OUT/'parallel_probe_plan.json').exists();p=ROOT/'src/benchmark_mg184_parallel.py';atomic_json(OUT/'parallel_probe_plan.json',dict(source_sha256=hashlib.sha256(p.read_bytes()).hexdigest(),workers=8,duration=5.,omega=3000,seeds=list(range(18420,18428)),wall_cap_seconds=90,development_only=True))
def run():
 p=read(OUT/'parallel_probe_plan.json');assert hashlib.sha256((ROOT/'src/benchmark_mg184_parallel.py').read_bytes()).hexdigest()==p['source_sha256'] and not(OUT/'parallel_probe_attempt.json').exists();start=time.perf_counter();atomic_json(OUT/'parallel_probe_attempt.json',dict(status='RUNNING'))
 try:
  with ProcessPoolExecutor(max_workers=8,mp_context=multiprocessing.get_context('spawn')) as pool:rows=list(pool.map(worker,range(8)))
  elapsed=time.perf_counter()-start;assert elapsed<90;atomic_json(OUT/'parallel_probe.json',dict(rows=rows,wall_seconds=elapsed,worker_seconds=sum(r['seconds'] for r in rows),average_worker_seconds_per_model_unit=sum(r['seconds'] for r in rows)/40,includes_process_startup_in_wall=True,production_started=False));atomic_json(OUT/'parallel_probe_attempt.json',dict(status='COMPLETE',elapsed_seconds=elapsed))
 except BaseException:atomic_json(OUT/'parallel_probe_attempt.json',dict(status='FAILED',elapsed_seconds=time.perf_counter()-start,error=traceback.format_exc()));raise
if __name__=='__main__':register() if __import__('sys').argv[1]=='register' else run()

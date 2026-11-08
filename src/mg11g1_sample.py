"""MG083 generic atomic ensemble writer, reused N1 stream law, no production CLI."""
import os,time,json,hashlib
import numpy as np
from pathlib import Path
from mg10d1_pool import atomic_json
from mg11g1 import cell
from mg11n1 import validate
from mg11n1_stream import stream_chunk
from check_mg11n1_preparation import check_raw
def raw_save(path,case,g,k,t,m,seed,data,env,stats,requested):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True);temp=path.with_suffix('.npz.tmp')
    with temp.open('wb') as f:np.savez_compressed(f,data=data,environment=env,stats=stats,initial=case['initial'],a0=case['a0'],b0=case['b0'],omega=case['omega'],diff=case['diff'],src=g['src'],dst=g['dst'],weight=g['weight'],kappa=k,times=t,mask=m,seed=seed,requested_reps=requested,reps_completed=len(data))
    os.replace(temp,path)
def generate(job,directory,block=128):
    directory=Path(directory);directory.mkdir(parents=True,exist_ok=True);case,g,k,t,m=cell(job['resource'],job['network']);validate(case,g,k,m);total=job['reps'];seed=job['seed']
    warm=stream_chunk(case['initial'],case['a0'],case['b0'],case['omega'],case['diff'],m,t,1,job['warmup_seed'],g['src'],g['dst'],g['weight'],k,1000000)
    raw_save(directory/(job['id']+'_warm.npz'),case,g,k,t,m,job['warmup_seed'],*warm,1)
    data=np.zeros((total,len(t),20,6),dtype=np.int64);env=np.zeros((total,len(t),2,len(g['src'])),dtype=np.int64);stats=np.zeros((total,3));path=directory/(job['id']+'.npz');compute=0.;io=0.;started=time.perf_counter()
    raw_save(path,case,g,k,t,m,seed,data[:0],env[:0],stats[:0],total)
    for left in range(0,total,block):
        right=min(total,left+block);clock=time.perf_counter();result=stream_chunk(case['initial'],case['a0'],case['b0'],case['omega'],case['diff'],m,t,right-left,seed if left==0 else -1,g['src'],g['dst'],g['weight'],k,1000000);compute+=time.perf_counter()-clock
        data[left:right],env[left:right],stats[left:right]=result;clock=time.perf_counter();raw_save(path,case,g,k,t,m,seed,data[:right],env[:right],stats[:right],total)
        atomic_json(directory/(job['id']+'_progress.json'),dict(job=job,replicas_completed=right,requested=total));io+=time.perf_counter()-clock
    checked=check_raw(path);meta=dict(status='COMPLETE',job=job,reps=total,compute_seconds=compute,io_seconds=io,elapsed_seconds=time.perf_counter()-started,sha256=hashlib.sha256(path.read_bytes()).hexdigest(),raw=path.name,raw_check=checked)
    atomic_json(path.with_suffix('.json'),meta);return meta

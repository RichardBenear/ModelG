"""MG096 atomic integer path writer, private continued RNG and raw reconstruction."""
import os,time,json,hashlib
from pathlib import Path
import numpy as np
from mg11g4 import simulate,validate,arrays,leap_chunk,exact_chunk,MAX_STEPS
from mg11g4_cases import fixture,TAUS,R,CHUNK
from mg11g4_oracle import direct_ssa
from mg10d1_pool import atomic_json
ROOT=Path(__file__).resolve().parents[1]
def digest(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def jobs():
    from runner import seed_for
    return [dict(id=f'c{c}_e{e}_r{r}',cell=c,engine=e,chunk=r,reps=CHUNK,seed=seed_for(11,c,e,9603,r),warmup_seed=seed_for(11,c,e,9603,r,1)) for c in range(3) for e in range(5) for r in range(R//CHUNK)]

def sample(inputs,engine,reps,seed,rng=None):
    case,ga,gb,k,t,mask=inputs;validate(case,ga,gb,k,mask)
    if engine==4:return direct_ssa(case,ga,gb,k,t,reps,seed,mask,_rng=rng)
    src,dst,w,degree,offsets,edge_ids=arrays(ga,gb,k)
    if engine==0:return exact_chunk(case['initial'],case['a0'],case['b0'],case['omega'],case['diff'],mask,t,reps,seed,src,dst,w,k,1000000)
    return leap_chunk(case['initial'],case['a0'],case['b0'],case['omega'],case['diff'],mask,TAUS[engine],t,reps,seed,src,dst,w,k,degree,offsets,edge_ids,MAX_STEPS)

def save(path,inputs,seed,result,requested,engine,tau_override=None):
    case,ga,gb,k,t,mask=inputs;data,environment,stats=result;path=Path(path);path.parent.mkdir(parents=True,exist_ok=True);tmp=path.with_suffix('.npz.tmp')
    with tmp.open('wb') as f:np.savez_compressed(f,data=data,environment=environment,stats=stats,initial=case['initial'],a0=case['a0'],b0=case['b0'],omega=case['omega'],diff=case['diff'],src=ga['src'],dst=ga['dst'],weight=ga['weight'],srcB=gb['src'],dstB=gb['dst'],weightB=gb['weight'],kappa=k,times=t,mask=mask,seed=seed,requested_reps=requested,reps_completed=len(data),engine=engine,tau=TAUS[engine] if tau_override is None else tau_override)
    os.replace(tmp,path)

def check_raw(path):
    path=Path(path)
    with np.load(path) as z:
        data=z['data'];ef=z['environment'];stats=z['stats'];initial=z['initial'];r,t,rows,n=data.shape;e=len(z['src']);assert rows==20 and ef.shape==(r,t,2,e) and stats.shape==(r,3);assert data.dtype==ef.dtype==initial.dtype==np.int64 and np.all(data>=0) and np.all(ef>=0) and np.all(np.isfinite(stats));assert z['reps_completed']==r and z['requested_reps']>=r
        if r:
            np.testing.assert_array_equal(data[:,0,:7],np.broadcast_to(initial,(r,7,n)));assert not np.any(data[:,0,7:]) and not np.any(ef[:,0]);assert np.all(np.diff(data[:,:,7:],axis=1)>=0) and np.all(np.diff(ef,axis=1)>=0)
            rx=data[:,:,7:14];f=rx;pred=initial[None,None]+np.stack((f[:,:,0]-f[:,:,1]+f[:,:,2]-f[:,:,6],f[:,:,1]-f[:,:,2]-f[:,:,3]+f[:,:,4]-f[:,:,5],f[:,:,3]-f[:,:,4],-f[:,:,0]+f[:,:,6],-f[:,:,3],f[:,:,3],f[:,:,5]),axis=2)
            hop=data[:,:,14:20]
            for s in range(3):
                left,right=hop[:,:,2*s],hop[:,:,2*s+1];pred[:,:,s]+=np.roll(left,-1,axis=-1)+np.roll(right,1,axis=-1)-left-right
            for s,(src,dst) in enumerate(((z['src'],z['dst']),(z['srcB'],z['dstB']))):
                for edge,(i,j) in enumerate(zip(src,dst)):pred[:,:,3+s,i]-=ef[:,:,s,edge];pred[:,:,3+s,j]+=ef[:,:,s,edge]
            np.testing.assert_array_equal(pred,data[:,:,:7]);np.testing.assert_array_equal(data[:,:,:7].sum(axis=(2,3)),np.full((r,t),initial.sum(),np.int64));assert not np.any(data[:,:,13])
            engine=int(z['engine'])
            if engine in (0,4):np.testing.assert_array_equal(stats[:,0],data[:,-1,7:].sum(axis=(1,2))+ef[:,-1].sum(axis=(1,2)))
            else:assert np.all(stats[:,0]<=MAX_STEPS) and np.all(stats[:,1]<=stats[:,0]) and np.all(stats[:,2]<=.25+1e-14)
    return dict(file=path.relative_to(ROOT).as_posix(),sha256=digest(path),replicas=r,checkpoints=t,integer_ledger=True,closed_total=True,nonnegative=True)

def generate(job,directory,block=128):
    directory=Path(directory);directory.mkdir(parents=True,exist_ok=True);inputs=fixture(job['cell']);case,ga,gb,k,t,mask=inputs;engine=job['engine'];start=time.perf_counter();warm=sample(inputs,engine,1,job['warmup_seed']);save(directory/(job['id']+'_warm.npz'),inputs,job['warmup_seed'],warm,1,engine)
    total=job['reps'];data=np.zeros((total,len(t),20,ga['n']),np.int64);env=np.zeros((total,len(t),2,len(ga['src'])),np.int64);stats=np.zeros((total,3));p=directory/(job['id']+'.npz');save(p,inputs,job['seed'],(data[:0],env[:0],stats[:0]),total,engine);compute=0.;io=0.;rng=np.random.default_rng(job['seed']) if engine==4 else None
    for left in range(0,total,block):
        right=min(total,left+block);clock=time.perf_counter();result=sample(inputs,engine,right-left,job['seed'] if left==0 or engine==4 else -1,rng=rng);compute+=time.perf_counter()-clock;data[left:right],env[left:right],stats[left:right]=result;clock=time.perf_counter();save(p,inputs,job['seed'],(data[:right],env[:right],stats[:right]),total,engine);atomic_json(p.with_name(p.stem+'_progress.json'),dict(job=job,completed=right,requested=total));io+=time.perf_counter()-clock
    checked=check_raw(p);meta=dict(status='COMPLETE',job=job,reps=total,compute_seconds=compute,io_seconds=io,elapsed_seconds=time.perf_counter()-start,sha256=digest(p),raw_check=checked);atomic_json(p.with_suffix('.json'),meta);return meta

"""MG157 full-horizon paired pilot; frozen event samplers, compact checkpoints."""
import os
for key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[key]='1'
from pathlib import Path
from collections import defaultdict
import json,math,time,hashlib,subprocess,sys,platform,traceback,argparse
import numpy as np
from mg12h3_tick_coupling import pair,replay
from mg12h3_sparse import Model,tick
from mg12h3_finer_precision import setup
from run_mg12h3_tick_coupling import exact
from prepare_mg12h3_bridge import instantiate,config,ROOT,seed_for
from mg12h2 import HyperState,Relation
from mg10e2c1_pool import atomic_json,run as pool_run
BASE=ROOT/'results/MG12H3/coupling_long_pilot'
DELTAS=((.0001,.00002),(.00005,.00002));R=8
FILES=['src/check_mg12h3_coupling_long.py','src/mg12h3_coupling_long.py','src/mg12h3_tick_coupling.py','src/mg12h3_sparse.py','src/mg12h3_finer_precision.py','src/mg12h3_contrast_estimator.py','src/run_mg12h3_tick_coupling.py','src/prepare_mg12h3_bridge.py','src/mg12h1.py','src/mg12h2.py','src/mg10e2c1_pool.py','src/runner.py','spec/MG12H3_coupling_long_pilot.md']
def read(p):return json.loads(Path(p).read_text(encoding='utf-8'))
def hashes():
    h={p:hashlib.sha256((ROOT/p).read_text(encoding='utf-8').replace('\r\n','\n').encode()).hexdigest() for p in FILES}
    p=ROOT/'results/MG12H3/ring_feasibility/o0_a0_Omega10.npz';h[p.relative_to(ROOT).as_posix()]=hashlib.sha256(p.read_bytes()).hexdigest();return h
def plan():
    assert not (BASE/'plan.json').exists();prior=set()
    def collect(x):
        if isinstance(x,dict):
            for v in x.values():collect(v)
        elif isinstance(x,list):
            for v in x:collect(v)
        elif type(x) is int and 0<=x<2**32:prior.add(x)
    for p in (ROOT/'results').rglob('*.json'):
        if 'plan' in p.name or 'registration' in p.name:collect(read(p))
    used=set(prior);collisions=[]
    def seed(stage,c,r):
        nonce=0;v=seed_for(12,c,0,stage,r,nonce)
        while v in used:collisions.append([stage,c,r,nonce,v]);nonce+=1;v=seed_for(12,c,0,stage,r,nonce)
        used.add(v);return v
    jobs=[]
    for c in range(2):
        for kind in ('paired','independent'):
            for r in range(R):jobs.append(dict(id=f'c{c}_{kind}_r{r:02d}',comparison=c,kind=kind,replicate=r,deltas=DELTAS[c],sentinel=kind=='paired' and r==0,seeds=[seed(15701 if kind=='paired' else 15702,c,2*r+s) for s in range(1 if kind=='paired' else 2)]))
    seeds=dict(control=[seed(15700,0,r) for r in range(8192)],segmentation=[seed(15704,0,r) for r in range(64)],bootstrap=seed(15703,0,0))
    atomic_json(BASE/'plan.json',dict(stage='MG157',authorization='go ahead and do the longer test.',scope='full-horizon paired pilot; no final tick equivalence',source_hashes=hashes(),jobs=jobs,seeds=seeds,prior_seed_values=len(prior),collisions_resolved=collisions,workers=8,threads=1,control_seconds=120,worker_ceiling=18000,parent_allowance=600,job_seconds=1200,forecast_seconds={j['id']:(360 if j['kind']=='paired' else 260) for j in jobs},duration=.01,checkpoint=.001,environment=dict(python=sys.version,numpy=np.__version__,platform=platform.platform())))
def state():
    return HyperState((0,),tuple(Relation('Occupies',(f't{i}',0,sp)) for i,sp in enumerate(('A','B','X','X','Y'))),2.,((0,10.,4.),),(0.,0.,0.),(0.,0.),(True,)*6).validate()
def controls():
    p=read(BASE/'plan.json');assert p['source_hashes']==hashes() and not (BASE/'control_attempt.json').exists();start=time.perf_counter();deadline=start+120;atomic_json(BASE/'control_attempt.json',dict(status='RUNNING'))
    try:
        s=state();pop=s.counts();model=Model(1,config(s));cache={};laws=[]
        for side,d in enumerate((.01,.004)):
            distribution={tuple(pop.ravel()):1.};snapshots=[]
            for step in range(round(.04/d)):
                nxt=defaultdict(float)
                for key,mass in distribution.items():
                    ck=(d,key)
                    if ck not in cache:cache[ck]=exact(instantiate(s,np.array(key,np.int64).reshape(7,1)),d,deadline)
                    ref=cache[ck]
                    for dest,q in zip(ref['support'],ref['probabilities']):nxt[tuple(dest)]+=mass*float(q)
                distribution=dict(nxt)
                if step+1 in (round(.02/d),round(.04/d)):
                    assert abs(sum(distribution.values())-1)<1e-10;snapshots.append(distribution.copy())
            laws.append(snapshots)
        for i,(ck,ref) in enumerate(cache.items()):np.savez_compressed(BASE/f'control_kernel_{i:03d}.npz',delta=ck[0],initial=ck[1],**ref)
        data=[];segmented=[]
        for r,seed in enumerate(p['seeds']['segmentation']):
            whole=pair(model,pop,pop,(.01,.004),[0,.02,.04],seed,deadline,sentinel=True);rng=np.random.default_rng(seed);first=pair(model,pop,pop,(.01,.004),[0,.02],seed,deadline,rng=rng,sentinel=True);second=pair(model,*first['counts'][:,-1],(.01,.004),[0,.02],seed,deadline,rng=rng,sentinel=True);np.testing.assert_array_equal(whole['counts'][:,-1],second['counts'][:,-1]);assert whole['rng_state']==second['rng_state'];replay(model,pop,pop,(.01,.004),[0,.02,.04],whole);segmented.append(whole['counts'])
        np.savez_compressed(BASE/'segmentation.npz',counts=segmented,seeds=p['seeds']['segmentation'])
        for seed in p['seeds']['control']:
            result=pair(model,pop,pop,(.01,.004),[0,.02,.04],seed,deadline);data.append(result['counts'])
        data=np.asarray(data);np.savez_compressed(BASE/'control_samples.npz',counts=data,seeds=p['seeds']['control']);rows=[]
        for side in range(2):
            for f,law in enumerate(laws[side],1):
                support=np.asarray(list(law),np.int64);prob=np.asarray(list(law.values()));np.savez_compressed(BASE/f'control_law_s{side}_f{f}.npz',support=support,probabilities=prob);observed=defaultdict(int)
                for pop in data[:,side,f]:observed[tuple(pop.ravel())]+=1
                assert set(observed)<=set(law);tv=sum(abs(observed.get(k,0)/len(data)-v) for k,v in law.items())/2;bound=math.sqrt((len(law)*math.log(2)+math.log(4/.001))/(2*len(data)));rows.append(dict(side=side,frame=f,TV=tv,bound=bound,support=len(law)));assert tv<=bound
        elapsed=time.perf_counter()-start;assert elapsed<120;atomic_json(BASE/'control_results.json',dict(passed=True,elapsed_seconds=elapsed,rows=rows,segmentation_paths=64,pairs=8192,kernels=len(cache),source_hashes=hashes()));atomic_json(BASE/'control_attempt.json',dict(status='COMPLETE'));print(json.dumps(rows))
    except BaseException:
        atomic_json(BASE/'control_failure.json',dict(error=traceback.format_exc(),elapsed_seconds=time.perf_counter()-start));atomic_json(BASE/'control_attempt.json',dict(status='FAILED'));raise
def save_trace(path,result):
    bs=result['boundaries'];es=result['trace'];np.savez_compressed(path,times=[e['time'] for e in es],masks=np.asarray([e['mask'] for e in es],np.uint8),channels=np.asarray([e['channel'] for e in es],np.uint16),boundary_times=[b['time'] for b in bs],boundary_masks=np.asarray([b['mask'] for b in bs],np.uint8),boundary_available=np.asarray([b['available'] for b in bs],np.int64),boundary_output=np.asarray([b['output'] for b in bs],np.int64))
def worker(job,directory):
    directory=Path(directory);p=read(BASE/'plan.json');assert p['source_hashes']==hashes();assert job in p['jobs'];path=directory/job['id'];path.mkdir(exist_ok=False);start=time.perf_counter();deadline=start+1190
    try:
        initial,contrast=setup();model=contrast.model;pop=np.stack([initial,initial]);rngs=[np.random.default_rng(seed) for seed in job['seeds']];observations=[pop.copy()];contrasts=[[contrast.direct(x) for x in pop]];total_events=0;segment_seconds=[]
        atomic_json(path/'metadata.json',dict(job=job,source_hashes=hashes(),weights=contrast.weights.tolist(),initial_total=sum(map(int,initial.ravel()))))
        for segment in range(10):
            clock=time.perf_counter();before=[json.dumps(r.bit_generator.state) for r in rngs];previous=pop.copy()
            if job['kind']=='paired':
                result=pair(model,*pop,job['deltas'],[0,.001],job['seeds'][0],deadline,rng=rngs[0],sentinel=job['sentinel']);pop=result['counts'][:,-1];hist=np.stack([result['histogram'][0]+result['histogram'][s+1] for s in range(2)]);available=result['final_available_before_reset'];events=result['events']
                if job['sentinel']:save_trace(path/f'trace_{segment:02d}.npz',result)
                del result
            else:
                hist=np.zeros((2,len(model.channels)),np.int64);available=np.zeros_like(pop);events=0
                for side,d in enumerate(job['deltas']):
                    for step in range(round(.001/d)):
                        result=tick(model,pop[side],d,rngs[side],deadline,trace=True);pop[side]=result['counts'];available[side]=result['available'];ids=np.fromiter((e['channel'] for e in result['trace']),dtype=np.int32,count=result['events']);hist[side]+=np.bincount(ids,minlength=len(model.channels));events+=result['events']
                    del result
            for side in range(2):
                reconstructed=previous[side].ravel().copy()
                for i,ch in enumerate(model.channels):
                    for coordinate,change in ch.change:reconstructed[coordinate]+=int(hist[side,i])*change
                np.testing.assert_array_equal(reconstructed,pop[side].ravel());assert sum(map(int,reconstructed))==sum(map(int,initial.ravel()))
            after=[json.dumps(r.bit_generator.state) for r in rngs];np.savez_compressed(path/f'checkpoint_{segment:02d}.npz',initial=previous,final=pop,available=available,histogram=hist,rng_before=before,rng_after=after,events=events);observations.append(pop.copy());contrasts.append([contrast.direct(x) for x in pop]);total_events+=events;segment_seconds.append(time.perf_counter()-clock);atomic_json(path/'progress.json',dict(segments_complete=segment+1,segments_required=10,events=total_events))
        np.savez_compressed(path/'observations.npz',counts=observations,contrasts=contrasts,times=np.arange(11)*.001,seeds=job['seeds']);atomic_json(directory/(job['id']+'.json'),dict(status='COMPLETE',job=job,elapsed_seconds=time.perf_counter()-start,events=total_events,segment_seconds=segment_seconds,source_hashes=hashes()))
    except BaseException:
        atomic_json(directory/(job['id']+'.json'),dict(status='FAILED',job=job,elapsed_seconds=time.perf_counter()-start,error=traceback.format_exc()));raise
def production():
    p=read(BASE/'plan.json');assert p['source_hashes']==hashes();assert read(BASE/'control_results.json')['passed'];assert not (BASE/'production_attempt.json').exists()
    subprocess.run(['git','ls-files','--error-unmatch',str((BASE/'plan.json').relative_to(ROOT))],cwd=ROOT,check=True,stdout=subprocess.DEVNULL);atomic_json(BASE/'production_attempt.json',dict(status='RUNNING'));start=time.perf_counter()
    ledger=pool_run(p['jobs'],BASE/'production',Path(__file__).resolve(),ROOT,p['worker_ceiling'],workers=8,timeout=p['job_seconds'],forecast=p['forecast_seconds']);complete=len(ledger['records'])==32 and all(r['status']=='COMPLETE' and r['returncode']==0 for r in ledger['records']);atomic_json(BASE/'production_attempt.json',dict(status='COMPLETE' if complete else 'STOPPED',wall_seconds=time.perf_counter()-start,charged_worker_seconds=ledger['charged_worker_seconds'],stop_reason=ledger['stop_reason']));print(json.dumps(read(BASE/'production_attempt.json')))
if __name__=='__main__':
    if len(sys.argv)>1 and sys.argv[1]=='--worker':worker(json.loads(sys.argv[2]),sys.argv[3])
    else:
        parser=argparse.ArgumentParser();parser.add_argument('action',choices=['plan','controls','run']);a=parser.parse_args();dict(plan=plan,controls=controls,run=production)[a.action]()

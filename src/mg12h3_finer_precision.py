"""MG147 compact actual-ring precision-pilot writer, frozen sampler wrapper."""
import os
for key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[key]='1'
from pathlib import Path
from datetime import datetime,timezone
import numpy as np,json,time,hashlib,math,sys,traceback,platform
from mg12h3_contrast_estimator import Contrast,tick,aggregate_audit
from run_mg12h3_sparse_combined import ROOT,load
from mg10e2c1_pool import atomic_json
from runner import seed_for
BASE=ROOT/'results/MG12H3/finer_tick_precision_pilot';CASE='o0_a0_Omega10';S=np.array([140/9,140/9,29/(140/9)])
FILES=['src/mg12h3_finer_precision.py','src/prepare_mg12h3_finer_precision.py','src/run_mg12h3_finer_precision.py','spec/MG12H3_finer_precision_pilot.md','src/mg12h3_contrast_estimator.py','src/mg12h3_sparse.py','src/run_mg12h3_sparse_combined.py','src/mg10e2c1_pool.py']
def hashes():
 h={p:hashlib.sha256((ROOT/p).read_text(encoding='utf-8').replace('\r\n','\n').encode()).hexdigest() for p in FILES}
 p=ROOT/f'results/MG12H3/ring_feasibility/{CASE}.npz';h[p.relative_to(ROOT).as_posix()]=hashlib.sha256(p.read_bytes()).hexdigest();return h

def environment():
 return dict(python=sys.version,numpy=np.__version__,platform=platform.platform(),byteorder=sys.byteorder,default_rng='PCG64')

def setup():
 p=ROOT/f'results/MG12H3/ring_feasibility/{CASE}.npz';initial,model=load(p)
 with np.load(p) as z:reference=z['concentrations']
 weights=np.zeros((3,initial.size))
 for s in range(3):
  template=reference[s]/S[s];template-=template.mean();template/=np.sqrt(np.mean(template**2));weights[s,s*500:(s+1)*500]=template/(500*10*S[s])
 return initial,Contrast(model,weights)

def jobs():return [dict(id=f'd{level}_r{r:02d}',level=level,replicate=r,delta=delta,stride=stride,duration=.01,observation=.001,sentinel=r==0) for level,(delta,stride) in enumerate(((.0001,10),(.00005,20),(.00002,50))) for r in range(8)]

def power(pop,initial):
 difference=(pop[:3]-initial[:3])/10.;mask=np.abs(np.fft.fftfreq(500)*500)>=125;filtered=np.fft.ifft(np.fft.fft(difference,axis=1)*mask,axis=1).real;return 10/S*np.mean(filtered**2,axis=1)

def save_chunk(path,rows):
 keys=('final','available','histogram','exposures','integral','direct','cumulative','events','last_accepted_time')
 values={key:np.asarray([r[key] for r in rows]) for key in keys};values['steps']=np.asarray([r['step'] for r in rows]);values['next_mark']=np.asarray([np.nan if r['next_mark'] is None else r['next_mark'] for r in rows]);values['rng_before']=np.asarray([r['rng_before'] for r in rows]);values['rng_after']=np.asarray([r['rng_after'] for r in rows]);values['checks']=np.asarray([json.dumps(r['checks']) for r in rows]);lengths=[len(r['trace']) for r in rows];values['trace_offsets']=np.r_[0,np.cumsum(lengths)];values['trace_ids']=np.asarray([e['channel'] for r in rows for e in r['trace']],np.int32);values['trace_marks']=np.asarray([e['time'] for r in rows for e in r['trace']],float);np.savez_compressed(path,**values)

def ordered_audit(contrast,initial,result,delta):
 # Replay only recorded events, independently reconstruct lazy factor exposures.
 available=initial.copy().ravel();groups=contrast.factors;coeff=contrast.coefficients;exposure=np.zeros(len(groups));last=np.zeros(len(groups));values=np.array([contrast.factor(g,available) for g in range(len(groups))]);hist=np.zeros(len(contrast.model.channels),np.uint32);previous=0.
 for event in result['trace']:
  mark=event['time'];assert previous<mark<delta;channel=contrast.model.channels[event['channel']];assert channel.rate(available)>0;touched=set()
  for i,n in channel.reserve:
   assert available[i]>=n;touched.update(contrast.dependencies[i])
  for g in touched:exposure[g]+=values[g]*(mark-last[g]);last[g]=mark
  for i,n in channel.reserve:available[i]-=n
  for g in touched:values[g]=contrast.factor(g,available)
  hist[event['channel']]+=1;previous=mark
 exposure+=values*(delta-last);np.testing.assert_array_equal(hist,result['histogram']);np.testing.assert_array_equal(available,result['available'].ravel());np.testing.assert_allclose(exposure,result['exposures'],rtol=1e-10,atol=1e-10)

def write(job,seed,directory,deadline):
 directory=Path(directory);directory.mkdir(parents=True,exist_ok=True);path=directory/job['id'];path.mkdir(exist_ok=False);clock=time.perf_counter();started=datetime.now(timezone.utc).isoformat();initial,contrast=setup();atomic_json(path/'metadata.json',dict(job=job,seed=seed,contrast=contrast.metadata(),source_hashes=hashes(),environment=environment(),environment_sha256=hashlib.sha256(json.dumps(environment(),sort_keys=True).encode()).hexdigest()));pop=initial.copy();nsteps=round(job['duration']/job['delta']);assert abs(nsteps*job['delta']-job['duration'])<1e-14 and abs(job['stride']*job['delta']-job['observation'])<1e-14;rng=np.random.default_rng(seed);cumulative=contrast.direct(initial);rows=[];chunk=0;frames=[initial.copy()];direct=[cumulative.copy()];estimates=[cumulative.copy()];powers=[np.zeros(3)];events=0
 try:
  for step in range(nsteps):
   out=tick(contrast,pop,job['delta'],rng,deadline=deadline,sentinel=job['sentinel']);aggregate_audit(contrast,pop,out);cumulative=cumulative+out['integral'];rows.append(dict(out,step=step,cumulative=cumulative.copy(),final=out['counts']));pop=out['counts'];events+=out['events']
   if (step+1)%job['stride']==0:frames.append(pop.copy());direct.append(out['direct']);estimates.append(cumulative.copy());powers.append(power(pop,initial))
   if len(rows)==20:save_chunk(path/f'chunk_{chunk:04d}.npz',rows);chunk+=1;rows=[];atomic_json(path/'progress.json',dict(ticks_complete=step+1,ticks_required=nsteps,events=events))
 finally:
  if rows:save_chunk(path/f'chunk_{chunk:04d}.npz',rows)
 np.savez_compressed(path/'observations.npz',initial=initial,counts=np.asarray(frames),direct=np.asarray(direct),estimates=np.asarray(estimates),power=np.asarray(powers),times=np.arange(len(frames))*job['observation'],job=json.dumps(job),seed=seed);meta=dict(status='COMPLETE',job=job,seed=seed,ticks=nsteps,events=events,elapsed_seconds=time.perf_counter()-clock,started=started,finished=datetime.now(timezone.utc).isoformat());atomic_json(directory/(job['id']+'.json'),meta);return meta

def audit(job,seed,directory,deadline=math.inf):
 directory=Path(directory);path=directory/job['id'];meta=json.loads((directory/(job['id']+'.json')).read_text());assert meta['status']=='COMPLETE' and meta['job']==job and meta['seed']==seed;initial,contrast=setup();savedmeta=json.loads((path/'metadata.json').read_text());assert savedmeta['job']==job and savedmeta['seed']==seed and savedmeta['contrast']==json.loads(json.dumps(contrast.metadata())) and savedmeta['source_hashes']==hashes() and savedmeta['environment']==environment() and savedmeta['environment_sha256']==hashlib.sha256(json.dumps(environment(),sort_keys=True).encode()).hexdigest();pop=initial.copy();frames=[pop.copy()];direct=[contrast.direct(pop)];cumulative=direct[0].copy();estimates=[cumulative.copy()];powers=[np.zeros(3)];step=0;events=0;prior_rng=None
 for file in sorted(path.glob('chunk_*.npz')):
  with np.load(file) as z:
   assert z['steps'].tolist()==list(range(step,step+len(z['steps'])));offsets=z['trace_offsets'];assert offsets[0]==0 and offsets[-1]==len(z['trace_ids'])==len(z['trace_marks'])
   for k in range(len(z['steps'])):
    if time.perf_counter()>deadline:raise TimeoutError('Compact audit cap')
    out=dict(counts=z['final'][k],available=z['available'][k],histogram=z['histogram'][k],exposures=z['exposures'][k],integral=z['integral'][k]);aggregate_audit(contrast,pop,out);assert int(z['events'][k])==int(out['histogram'].sum());np.testing.assert_allclose(z['direct'][k],contrast.direct(out['counts']),rtol=1e-12,atol=1e-12);assert 0<=z['last_accepted_time'][k]<job['delta'] and (np.isnan(z['next_mark'][k]) or z['next_mark'][k]>=job['delta'])
    before=str(z['rng_before'][k]);after=str(z['rng_after'][k]);assert (json.loads(before)==np.random.default_rng(seed).bit_generator.state if prior_rng is None else prior_rng==before);prior_rng=after
    checks=json.loads(str(z['checks'][k]));assert checks and checks[-1]['events']==int(z['events'][k])
    if job['sentinel']:
     out['trace']=[dict(channel=int(i),time=float(t)) for i,t in zip(z['trace_ids'][offsets[k]:offsets[k+1]],z['trace_marks'][offsets[k]:offsets[k+1]])];ordered_audit(contrast,pop,out,job['delta'])
    else:assert offsets[k]==offsets[k+1]
    pop=out['counts'].copy();cumulative+=out['integral'];np.testing.assert_allclose(cumulative,z['cumulative'][k],rtol=1e-12,atol=1e-12);events+=int(z['events'][k]);step+=1
    if step%job['stride']==0:frames.append(pop.copy());direct.append(contrast.direct(pop));estimates.append(cumulative.copy());powers.append(power(pop,initial))
 assert step==round(job['duration']/job['delta']) and events==meta['events']
 with np.load(path/'observations.npz') as z:
  assert int(z['seed'])==seed and json.loads(str(z['job']))==job;np.testing.assert_array_equal(z['counts'],frames);np.testing.assert_allclose(z['direct'],direct,rtol=1e-12,atol=1e-12);np.testing.assert_allclose(z['estimates'],estimates,rtol=1e-12,atol=1e-12);np.testing.assert_allclose(z['power'],powers,rtol=1e-12,atol=1e-12);np.testing.assert_allclose(z['times'],np.arange(len(frames))*job['observation'],rtol=0,atol=1e-14)
 return dict(passed=True,events=events,ticks=step)

"""MG142 independent saved channel-ID aggregation and descriptive metrics."""
from mg12h3_multitick import ROOT,BASE,hashes,load,np,json,hashlib,math

def run():
 r=json.loads((BASE/'registration.json').read_text());assert hashes()==r['source_hashes'];records=[];finals={};manifest=[]
 for job in r['jobs']:
  directory=BASE/'production'/job['id'];meta=json.loads((BASE/'production'/(job['id']+'.json')).read_text());assert meta['status']=='COMPLETE' and meta['job']==job and meta['seed']==r['seeds'][job['id']];initial,model=load(ROOT/f'results/MG12H3/ring_feasibility/{job["case"]}.npz');pop=initial.copy();frames=[pop.copy()];events=0;drift=0.;leaf=0.
  for step in range(round(job['duration']/job['delta'])):
   with np.load(directory/f'tick_{step:04d}.npz') as z:
    assert int(z['step'])==step and int(z['seed'])==meta['seed'] and float(z['delta'])==job['delta'];np.testing.assert_array_equal(pop,z['initial']);ids=z['channel_ids'];marks=z['marks'];assert len(ids)==len(marks) and np.all(np.diff(marks)>0) and np.all(marks>0) and np.all(marks<job['delta']);assert ids.min()>=0 and ids.max()<len(model.channels)
    ordered=np.sort(ids);freq=np.searchsorted(ordered,np.arange(len(model.channels)),side='right')-np.searchsorted(ordered,np.arange(len(model.channels)),side='left');reserve=np.zeros(pop.size,np.int64);change=reserve.copy()
    for c,n in zip(model.channels,freq):
     if n:
      np.add.at(reserve,[i for i,v in c.reserve],[int(n)*v for i,v in c.reserve]);np.add.at(change,[i for i,v in c.change],[int(n)*v for i,v in c.change])
    np.testing.assert_array_equal(pop.ravel()-reserve,z['available'].ravel());np.testing.assert_array_equal(pop.ravel()+change,z['final'].ravel());assert z['available'].min()>=0 and z['final'].min()>=0 and z['final'].sum()==initial.sum();pop=z['final'].copy();events+=len(ids)
   checks=json.loads((directory/f'tick_{step:04d}_checks.json').read_text())['records'];assert checks and checks[-1]['events']==len(ids)
   for x in checks:
    assert abs(x['tree_total']-x['total'])<=1e-6+1e-10*abs(x['total']);drift=max(drift,abs(x['tree_total']-x['total'])/max(1,abs(x['total'])));leaf=max(leaf,x['max_leaf_error'])
   if (step+1)%(job['level']+1)==0:frames.append(pop.copy())
  with np.load(directory/'frames.npz') as z:np.testing.assert_array_equal(z['counts'],frames);np.testing.assert_allclose(z['times'],np.arange(11)*.001,rtol=0,atol=1e-14)
  assert events==meta['events'];omega=10 if job['case'].endswith('Omega10') else 1;scale=np.maximum(1,np.sqrt(np.mean((initial[:3]/omega)**2,axis=1)));movement=float(np.sqrt(np.mean(((pop[:3]-initial[:3])/omega/scale[:,None])**2)));finals[job['id']]=pop.copy();records.append(dict(id=job['id'],case=job['case'],level=job['level'],ticks=len(frames)-1 if job['level']==0 else (len(frames)-1)*2,events=events,worker_body_seconds=meta['elapsed_seconds'],maximum_relative_tree_drift=drift,maximum_leaf_error=leaf,normalized_active_RMS_change=movement,passed=True))
 pairs=[]
 for case in sorted({j['case'] for j in r['jobs']}):
  initial,_=load(ROOT/f'results/MG12H3/ring_feasibility/{case}.npz');omega=10 if case.endswith('Omega10') else 1;scale=np.maximum(1,np.sqrt(np.mean((initial[:3]/omega)**2,axis=1)));difference=float(np.sqrt(np.mean(((finals[case+'_d0'][:3]-finals[case+'_d1'][:3])/omega/scale[:,None])**2)));pairs.append(dict(case=case,normalized_coarse_fine_RMS=difference,scope='Descriptive independent trajectories, not convergence or bias estimate'))
 out=dict(passed=True,trajectories=len(records),ticks=sum(x['ticks'] for x in records),events=sum(x['events'] for x in records),records=records,coarse_fine=pairs,source_hashes_unchanged=True,shared='Saved compact channel compilation and inputs; sorted-ID and np.add.at aggregation separate from registered bincount/scalar loops. Same author, not external replication. Does not rederive worker leaf/selection checks.')
 (BASE/'separate_recomputation.json').write_text(json.dumps(out,indent=2)+'\n')
 for p in sorted((BASE/'production').rglob('*')):
  if p.is_file():manifest.append(dict(path=p.relative_to(ROOT).as_posix(),bytes=p.stat().st_size,sha256=hashlib.sha256(p.read_bytes()).hexdigest()))
 (BASE/'raw_manifest.json').write_text(json.dumps(dict(files=manifest,total_bytes=sum(x['bytes'] for x in manifest)),indent=2)+'\n');print(json.dumps({k:v for k,v in out.items() if k not in ('records','coarse_fine')}));print('rawbytes',sum(x['bytes'] for x in manifest))
if __name__=='__main__':run()

"""Same-author separate recomputation; no new stochastic trajectories."""
import hashlib,json,tempfile
import numpy as np
from pathlib import Path
from mg184_phase2 import ref,initial,summarize
from mg10e2c1_pool import atomic_json
from run_mg176 import ROOT

def run():
 out=ROOT/'results/MG12H3/phase2_G_development'
 read=lambda p:json.loads(p.read_text())
 plan=read(out/'plan.json');checks=[]
 for name,digest in plan['source'].items():
  assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==digest
 assert hashlib.sha256((ROOT/'results/MG003/ref_profiles.npz').read_bytes()).hexdigest()==plan['reference_sha256']
 checks.append('frozen kernel, adapter, detector and reference hashes unchanged')
 r=ref();pop=initial(10000,'prepared',0,r);times=np.arange(201)*5.
 fields=np.broadcast_to(pop,(201,3,500)).copy();summary,a=summarize(fields,times,10000,'prepared',1000,r)
 # Independent explicit dot products, rather than the adapter einsum.
 for species in range(3):
  template=r['phi'][species]-np.mean(r['phi'][species]);unit=template/np.sqrt(sum(template**2))
  direct=sum((pop[species]/10000-np.mean(pop[species]/10000))*unit)
  np.testing.assert_allclose(a['fixed_projection'][:,species],direct,atol=1e-12)
 checks.append('G/X/Y projection recomputed with direct sums')
 fields[times>=500]=np.rint(10000*r['S0'][:,None])
 late,b=summarize(fields,times,10000,'prepared',1000,r)
 assert not late['historical_X_persistence'] and np.isnan(b['aligned_projection'][times>=530]).all()
 checks.append('full-horizon late-loss and missing-alignment regression')
 fields[:]=pop;fields[:,0]=np.rint(10000*r['S0'][0])
 missing,c=summarize(fields,times,10000,'prepared',1000,r)
 assert missing['historical_X_persistence'] and np.max(abs(c['fixed_projection'][:,0]))<1e-12
 assert missing['full_GXY_persistence']=='NOT_SCORED_FROM_X_ALONE'
 checks.append('X-only persistence cannot certify absent G')
 for row in read(out/'benchmark.json')['rows']:
  raw=np.load(ROOT/f"local_raw/MG184_development/Omega{row['omega']}_benchmark.npz")
  assert raw['counts'].shape==(4,3,500) and np.array_equal(raw['final'],raw['counts'][-1])
  np.testing.assert_allclose(raw['times'],[.5,1,1.5,2],atol=1e-10)
  assert raw['counts'].min()>=0 and row['clips']==0
 checks.append('all five saved benchmark endpoints, sampling grids and zero-clip records checked')
 from run_mg185_phase2_G import qualification_gate,TAU
 records=[dict(status='COMPLETE',clips=0,job=dict(condition=c,tau=tau),numeric_G_block_mean=0.) for c in ('prepared','homogeneous') for tau in (TAU,TAU/2) for i in range(12)]
 with tempfile.TemporaryDirectory() as d:
  path=Path(d)/'ledger.json';atomic_json(path,dict(records=records));assert qualification_gate(d)['passed']
  for record in records:
   if record['job']['tau']==TAU/2:record['numeric_G_block_mean']=1.
  atomic_json(path,dict(records=records));assert not qualification_gate(d)['passed']
  records[0]['status']='FAILED';atomic_json(path,dict(records=records))
  try:qualification_gate(d)
  except AssertionError:pass
  else:raise AssertionError('incomplete qualification admitted')
 checks.append('qualification accepts equality, rejects discrepancy and incomplete jobs')
 atomic_json(out/'separate_check.json',dict(passed=True,checks=checks,new_stochastic_draws=False,external_replication=False,shared_dependencies=['summarize','reference','numpy'],limits='Generator fixture checked by preparation; no finite-step bridge equivalence or long-horizon accuracy demonstrated.'))
if __name__=='__main__':run()

"""Independent direct-stencil background residual checks, no new integrations."""
import json
import numpy as np
from gravity_mg190 import PUBLIC,RAW,binding,bc
from mg10e2c1_pool import atomic_json

def run():
 read=lambda p:json.loads(p.read_text())
 assert read(PUBLIC/'development_plan.json')['source']==binding();a=read(PUBLIC/'preflight.json');maximum=0.
 for row in a['backgrounds']:
  with np.load(RAW/f"background_{row['arm']}_{row['eps']}.npz") as f:
   G,X,Y=f['fields'];aa=f['a'];bb=f['b'];lap=lambda v:(np.roll(v,1)+np.roll(v,-1)-2*v)/.2**2
   rhs=np.array([lap(G)-G+.1*X+aa,lap(X)+G-(1+bb)*X+X*X*Y,12*lap(Y)+bb*X-X*X*Y]);maximum=max(maximum,float(abs(rhs).max()));assert f['fields'].min()>0 and abs(rhs).max()<1e-9
   for name,v in [('G',G),('X',X)]:np.testing.assert_allclose(row[name+'_half_range'],(max(v)-min(v))/2)
 np.testing.assert_allclose(bc(14),30.04407748,atol=1e-7);ratios=[]
 for eps in (.00025,.0005,.001):
  A=next(r for r in a['backgrounds'] if r['arm']=='A' and r['eps']==eps);C=next(r for r in a['backgrounds'] if r['arm']=='C' and r['eps']==eps);ratios.append(dict(eps=eps,C_to_A_G_gradient=C['maximum_G_gradient']/A['maximum_G_gradient']))
 atomic_json(PUBLIC/'separate_check.json',dict(passed=True,backgrounds=16,maximum_direct_stencil_residual=maximum,threshold_matches_repository=True,source_and_reference_hashes_unchanged=True,criticality_gradient_leakage=ratios,external_replication=False,shared_dependencies=['NumPy','same background solutions and equations'],no_gradient_particle_runs=True))
if __name__=='__main__':run()

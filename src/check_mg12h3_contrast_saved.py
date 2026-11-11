"""MG146 independent saved exposure/contrast arithmetic; no simulation."""
from pathlib import Path
import numpy as np,json,math,hashlib
from prepare_mg12h3_bridge import ROOT
BASE=ROOT/'results/MG12H3/contrast_estimator'
def run():
 p=json.loads((BASE/'plan.json').read_text());r=json.loads((BASE/'results.json').read_text())
 for rel,sha in r['source_hashes'].items():assert hashlib.sha256((ROOT/rel).read_text(encoding='utf-8').replace('\r\n','\n').encode()).hexdigest()==sha
 with np.load(BASE/'single_channel.npz') as z:values=z['values'];exposure=z['exposures'];final=z['final'];assert z['seeds'].tolist()==p['single'];np.testing.assert_array_equal(values[:,0],final[:,0,0]);np.testing.assert_allclose(values[:,1],exposure[:,0],rtol=0,atol=1e-15)
 checks=[]
 for c in range(4):
  with np.load(BASE/f'mixed_c{c}.npz') as z:values=z['values'];final=z['final'];available=z['available'];hist=z['histogram'];exposure=z['exposures'];initial=z['initial'];meta=json.loads(str(z['metadata']));assert z['seeds'].tolist()==p['mixed'][c]
  w=np.array(meta['weights']);coef=np.array(meta['coefficients']);direct=final.reshape(len(final),-1)@w.T;initialvalue=initial.ravel()@w.T;integral=np.einsum('ij,jk->ik',exposure,coef);np.testing.assert_allclose(direct,values[:,0],rtol=1e-12,atol=1e-12);np.testing.assert_allclose(initialvalue+integral,values[:,1],rtol=1e-12,atol=1e-12);difference=direct-(initialvalue+integral);mean=difference.mean(0);se=difference.std(0,ddof=1)/math.sqrt(len(values));assert np.all(abs(mean)<=6*se+1e-6);assert final.min()>=0 and available.min()>=0;checks.append(dict(fixture=c,paths=len(values),passed=True))
 (BASE/'separate_recomputation.json').write_text(json.dumps(dict(passed=True,checks=checks,source_hashes_unchanged=True,shared='Saved factor coefficients and exposure arrays; independent matrix/einsum arithmetic, same author. Trace integration and sampler comparisons are the registered execution-time checks.'),indent=2)+'\n')
if __name__=='__main__':run()

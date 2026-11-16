"""MG191 static/synthetic implementation validation; no gradient trajectories."""
import tempfile
from pathlib import Path
from unittest.mock import patch
import numpy as np
import run_mg191 as driver
from mg191_gravity import reference_data,model_for,initialize,core,phase,velocity_blocks,wrap,uniform_stability
from mg10e2c1_pool import atomic_json

def run():
 old=driver.read(driver.PUBLIC/'implementation_checks.json');checks=[]
 for h in (.2,.1):
  r,ref=reference_data(h);n=len(r['x']);np.testing.assert_allclose(r['phi'],r['phi'][:,(-np.arange(n))%n],atol=1e-10,rtol=0);m=model_for(dict(arm='N',h=h));bg,_,_=m.steady();assert uniform_stability(m)['passed']
  for j in range(8):
   p=phase(j);positive=initialize(dict(h=h,position=p),bg);negative=initialize(dict(h=h,position=-p),bg);np.testing.assert_allclose(positive,negative[:,(-np.arange(n))%n],atol=1e-10,rtol=0);one=core(positive[1],bg.reshape(3,-1)[1],h,ref['depth']);two=core(negative[1],bg.reshape(3,-1)[1],h,ref['depth']);assert abs(wrap(one[0]+two[0]))<1e-10
 checks.append('16 reflected phase pairs preserve core measurements and template symmetry')
 times=np.arange(4001)*.5;centers=wrap(40+.01*times);obs=np.column_stack((centers,np.ones(4001)*ref['depth'],np.ones(4001)*ref['width'],np.ones(4001)));v=velocity_blocks(obs,times,np.ones(4001,bool),ref);assert v['valid_late'];np.testing.assert_allclose([b['velocity'] for b in v['blocks']],.01,atol=1e-14)
 lost=np.ones(4001,bool);lost[100:]=False;v=velocity_blocks(obs,times,lost,ref);assert not v['valid_late'] and all(b['velocity'] is None for b in v['blocks']);checks.append('periodic unwrap velocity exact; lost original never supplies valid late velocities')
 with tempfile.TemporaryDirectory() as directory:
  original_raw=driver.RAW;driver.RAW=Path(directory);(driver.RAW/'stable').write_bytes(b'1234');(driver.RAW/'ledger.json.tmp').write_bytes(b'ignore');assert driver.storage_bytes()==4;(driver.RAW/'vanished').write_bytes(b'x');original_stat=driver.os.stat
  def race(path,*args,**kwargs):
   if str(path).endswith('vanished'):raise FileNotFoundError(path)
   return original_stat(path,*args,**kwargs)
  try:
   with patch.object(driver.os,'stat',side_effect=race):assert driver.storage_bytes()==4
  finally:driver.RAW=original_raw
 checks.append('temporary-ledger rename race does not crash or inflate storage count')
 old.update(passed=True,synthetic_checks=checks,no_gradient_particle_runs=True);atomic_json(driver.PUBLIC/'implementation_checks.json',old)
if __name__=='__main__':run()

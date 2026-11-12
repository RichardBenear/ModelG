"""Deterministic inactive-model writer fixture: no random marks or production draws."""
from pathlib import Path
from types import SimpleNamespace
import json
import numpy as np
import run_mg160_persistence as driver
from mg12h3_sparse import Model
from mg10e2c1_pool import atomic_json
def check():
    base=driver.BASE/'writer_fixture';base.mkdir(exist_ok=False);initial=np.zeros((7,4),np.int64);initial[1]=[1,2,4,2];model=Model(4,(1.,[(10.,4.)]*4,[],(0.,0.,0.),(0.,0.),(False,)*6));job=dict(id='inactive_fixture',condition='prepared',seed=0,sentinel=True);binding={'fixture_only':'inactive/no_random_marks'};np.savez_compressed(base/'prepared_inputs.npz',initial=initial,homogeneous=initial,reference=initial.astype(float),omega=1.);atomic_json(base/'P1_plan.json',dict(jobs=[job],source_hashes=binding));old_base,old_setup,old_hash=driver.BASE,driver.setup,driver.hashes
    try:
        driver.BASE=base;driver.setup=lambda:(initial,SimpleNamespace(model=model));driver.hashes=lambda:binding;driver.worker(job,base)
    finally:driver.BASE=old_base;driver.setup=old_setup;driver.hashes=old_hash
    path=base/job['id'];meta=json.loads((base/(job['id']+'.json')).read_text());assert meta['status']=='COMPLETE' and meta['events']==0
    with np.load(path/'observations.npz') as z:assert z['counts'].shape==(101,2,7,4);assert np.all(z['counts']==initial);np.testing.assert_allclose(z['times'],np.arange(101)*.001)
    with np.load(path/'checkpoint_000.npz') as a,np.load(path/'checkpoint_099.npz') as b:assert str(a['rng_before'])==str(b['rng_after'])
    assert len(list(path.glob('checkpoint_*.npz')))==100 and len(list(path.glob('trace_*.npz')))==100
    atomic_json(driver.BASE/'writer_fixture_results.json',dict(passed=True,checkpoints=100,observation_frames=101,events=0,RNG_state_unchanged=True,scope='Inactive small-model file/framing fixture with mocked directory/setup/source binding. Not stochastic validation or P1 evidence.'))
if __name__=='__main__':check()

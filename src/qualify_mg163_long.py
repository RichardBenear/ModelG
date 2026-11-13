"""No-active-model framing, regression and resource-guard qualification."""
import time,json,tempfile,hashlib,traceback
from pathlib import Path
from types import SimpleNamespace
import numpy as np
import run_mg163_long as new
import run_mg160_persistence as old
from mg12h3_sparse import Model
from mg10e2c1_pool import atomic_json

class Detector:
    def __init__(self,*args):pass
    def measure(self,pop):return dict(total=int(pop.sum()))

def qualify():
    start=time.perf_counter();attempt=new.PUBLIC/'qualification_attempt.json';assert not attempt.exists();atomic_json(attempt,dict(status='RUNNING'))
    try:
        model=Model(1,(10,[(1.,1.)],[],[0.,0.,0.],[0.,0.],(False,)*6));initial=np.arange(7,dtype=np.int64).reshape(7,1);job=dict(id='inactive',condition='prepared',replicate=0,seed=163,sentinel=True)
        original={mod:{k:getattr(mod,k) for k in ('BASE','setup','hashes','Coherence')} for mod in (new,old)}
        with tempfile.TemporaryDirectory(prefix='mg163_qualification_') as temp:
            temp=Path(temp)
            try:
                for mod in (new,old):mod.setup=lambda:(None,SimpleNamespace(model=model));mod.hashes=lambda:{};mod.Coherence=Detector
                for name,segments,mod in [('long',1000,new),('short',100,new),('reference',100,old)]:
                    base=temp/name;base.mkdir();mod.BASE=base
                    np.savez_compressed(base/'prepared_inputs.npz',initial=initial,homogeneous=initial,reference=initial,omega=10.)
                    plan=dict(source_hashes={},jobs=[job],segments=segments,job_cap=100000)
                    atomic_json(base/('P1_plan.json' if mod is old else 'plan.json'),plan)
                    out=base/'output';out.mkdir();mod.worker(job,out)
                    with np.load(out/'inactive/observations.npz') as z:
                        np.testing.assert_array_equal(z['counts'],np.broadcast_to(initial,(segments+1,2,7,1)));np.testing.assert_array_equal(z['times'],np.arange(segments+1)*.001)
                    assert len(list((out/'inactive').glob('trace_*.npz')))==(segments if mod is old else segments//10)
                    expected=json.dumps(np.random.default_rng(job['seed']).bit_generator.state)
                    for i in range(segments):
                        with np.load(out/f'inactive/checkpoint_{i:03d}.npz') as z:
                            assert str(z['rng_before'])==str(z['rng_after'])==expected and int(z['events'])==0;np.testing.assert_array_equal(z['initial'],z['final']);np.testing.assert_array_equal(z['available'],z['final'])
                for i in range(100):
                    with np.load(temp/f'short/output/inactive/checkpoint_{i:03d}.npz') as a,np.load(temp/f'reference/output/inactive/checkpoint_{i:03d}.npz') as b:
                        for k in a.files:np.testing.assert_array_equal(a[k],b[k])
                for i in range(0,100,10):
                    with np.load(temp/f'short/output/inactive/trace_{i:03d}.npz') as a,np.load(temp/f'reference/output/inactive/trace_{i:03d}.npz') as b:
                        for k in a.files:np.testing.assert_array_equal(a[k],b[k])
                for call in (lambda:new.resource_guard(7999999999),lambda:new.storage_guard(temp,cap=1)):
                    try:call()
                    except RuntimeError:pass
                    else:raise AssertionError('Guard did not reject')
                new.resource_guard(8000000000)
                assert 563540.728079<.75*800000 and 1000*.001==1 and len([i for i in range(1000) if new.trace_segment(job,i)])==100
            finally:
                for mod,attrs in original.items():
                    for key,value in attrs.items():setattr(mod,key,value)
        atomic_json(new.PUBLIC/'qualification.json',dict(passed=True,elapsed_seconds=time.perf_counter()-start,inactive_segments=1200,active_model_draws=0,checks=['1000segment actual-engine inactive framing/count conservation/availability/RNG continuity','100segment checkpoint regression vs frozen old writer','selected full traces array-identical vs old writer','fixed100of1000trace schedule','local storage/free-space rejection boundaries','budget admission arithmetic'],limits='No active-model microbenchmark; frozen engine unchanged. Inactive checks do not establish long-horizon outcome.',source_hashes=new.hashes()))
        atomic_json(attempt,dict(status='COMPLETE',elapsed_seconds=time.perf_counter()-start))
    except BaseException:atomic_json(attempt,dict(status='FAILED',elapsed_seconds=time.perf_counter()-start,error=traceback.format_exc()));raise
if __name__=='__main__':qualify()

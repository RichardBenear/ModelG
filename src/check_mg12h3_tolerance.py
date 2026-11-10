"""MG118 saved-raw tolerance/reference/solver comparisons, no new integration."""
import json,hashlib
from itertools import combinations
import numpy as np
from prepare_mg12h3_population import read,factory,ROOT
from mg10d1_pool import atomic_json
BASE=ROOT/'results/MG12H3/tolerance_check'
def evaluate():
    from run_mg12h3_tolerance import diagnostic_hashes,JOB,METHODS
    registration=read(BASE/'registration.json');assert registration['source_original_hashes']==diagnostic_hashes();case,ga,gb,k,mode,reference,parent=factory(JOB);n=ga['n'];scale=np.maximum(1,case['initial'].max(axis=1));times=np.arange(101)/4
    with np.load(reference) as z:old=z['state'][:,:7*n].reshape(-1,7,n);np.testing.assert_array_equal(z['times'],times)
    data={};details=[];manifest=[]
    for method in METHODS:
        path=BASE/(method['id']+'.npz');meta=read(path.with_suffix('.json'));assert meta['method']==method and meta['status']=='COMPLETE' and meta['sha256']==hashlib.sha256(path.read_bytes()).hexdigest()
        with np.load(path) as z:state=z['state'];np.testing.assert_array_equal(z['times'],times);assert json.loads(str(z['job_json']))==JOB and json.loads(str(z['method_json']))==method
        assert state.shape==(101,9*n) and np.all(np.isfinite(state));pop=state[:,:7*n].reshape(-1,7,n);assert not np.any(state[:,7*n:]);np.testing.assert_array_equal(pop[0],case['initial']);error=np.abs(pop-old)/scale[None,:,None];maximum=float(error.max());balance=max(abs(sum(float(frame[s].sum()) for s in range(7))-float(case['initial'].sum())) for frame in pop)/float(case['initial'].sum());minimum=float(state.min());location=np.unravel_index(np.argmax(error),error.shape);passed=maximum<=1e-5 and balance<=1e-10 and minimum>=-1e-9
        np.testing.assert_allclose(maximum,meta['check']['reference_error'],atol=1e-14,rtol=1e-12);details.append(dict(method=method,reference_error=maximum,balance=balance,minimum=minimum,passed=bool(passed),species_errors=error.max(axis=(0,2)).tolist(),max_location=dict(time=float(times[location[0]]),species=('G','X','Y','A','B','Z','W')[location[1]],site=int(location[2]))));data[method['id']]=pop;manifest.append(dict(file=path.relative_to(ROOT).as_posix(),bytes=path.stat().st_size,sha256=hashlib.sha256(path.read_bytes()).hexdigest()))
    pairs=[]
    for left,right in combinations(data,2):
        error=float(np.max(abs(data[left]-data[right])/scale[None,:,None]));pairs.append(dict(left=left,right=right,error=error,passed=error<=1e-5))
    with np.load(ROOT/'results/MG12H3/population/production/ISOLATED_C1_s1_v0.npz') as z:failed=z['state'][:,:7*n].reshape(-1,7,n)
    original=float(np.max(abs(failed-old)/scale[None,:,None]));supported=all(r['passed'] for r in details+pairs);result=dict(complete=True,supported=bool(supported),verdict='TOLERANCE EXPLANATION SUPPORTED IN TESTED ISOLATED CONTROL' if supported else 'TOLERANCE DIAGNOSTIC NOT SUPPORTED UNDER FIXED CRITERIA',methods=details,pairs=pairs,original_failed_error=original,margin=1e-5,source_original_hashes_unchanged=True,original_study_status='failed/incomplete, unchanged; original fine isolated job remains unstarted',shared_dependencies='same physical initial fields/model and SciPy; independent Radau equations/sparsity, separate saved-array error/balance arithmetic, same author not external replication');atomic_json(BASE/'evaluation.json',result);atomic_json(BASE/'raw_manifest.json',dict(files=3,bytes=sum(r['bytes'] for r in manifest),records=manifest));return result
if __name__=='__main__':evaluate()

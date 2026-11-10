"""MG115 full analysis pipeline using archived reference fields only.

No ODE integration or RNG: these are held reference replays, not new trajectories.
"""
from prepare_mg12h3_population import *
from check_mg12h3_population import evaluate
def main():
    q=read(PREP/'results.json');assert q['passed'] and q['source_hashes']==hashes();start=time.perf_counter();directory=PREP/'held_reference_pipeline'
    if (PREP/'held_pipeline_results.json').exists():raise ValueError('Held pipeline already attempted')
    for job in jobs():
        if time.perf_counter()-start+q['elapsed_seconds']>=300:raise TimeoutError('Preparation aggregate cap')
        reference=factory(job)[5]
        with np.load(reference) as z:state=np.c_[z['state'][:,:3500],z['state'][:,-1000:]];times=z['times'].copy()
        save(directory/(job['id']+'.npz'),job,times,state)
    result=evaluate(directory,PREP/'held_full_evaluation.json');assert result['complete'] and result['qualified'];elapsed=time.perf_counter()-start;assert elapsed+q['elapsed_seconds']<300;assert hashes()==q['source_hashes'];atomic_json(PREP/'held_pipeline_results.json',dict(passed=True,scope='HELD CLOSED REFERENCE FIELDS ONLY; NOT NEW INTEGRATIONS OR SCIENTIFIC TRAJECTORIES',complete_jobs=24,elapsed_seconds=elapsed,aggregate_preparation_seconds=elapsed+q['elapsed_seconds'],source_hashes_unchanged=True,contrast=result['contrast'],max_reference_error=result['max_reference_error']))
if __name__=='__main__':main()

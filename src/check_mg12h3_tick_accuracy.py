"""MG158 complete archive checks and frozen paired accuracy intervals."""
import os
for k in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):os.environ[k]='1'
import time,json,traceback
import numpy as np
from scipy.stats import skew
from mg12h3_tick_accuracy import BASE,read,hashes
from mg12h3_accuracy_readiness import intervals
from mg12h3_finer_precision import setup,power
from mg12h3_tick_coupling import replay
from mg10e2c1_pool import atomic_json
def check():
    start=time.perf_counter();attempt=BASE/'analysis_attempt.json';atomic_json(attempt,dict(status='RUNNING'))
    try:
        plan=read(BASE/'plan.json');assert plan['source_hashes']==hashes();initial,contrast=setup();margins=np.r_[.05*abs(contrast.direct(initial)),[.05]*3];samples=[[],[]];endpoints=0;events=0
        for session in (1,2,3):assert read(BASE/f'session_{session}_attempt.json')['status']=='COMPLETE'
        for job in plan['jobs']:
            directory=BASE/f"session_{job['session']}";path=directory/job['id'];meta=read(directory/(job['id']+'.json'));assert meta['status']=='COMPLETE' and meta['source_hashes']==plan['source_hashes'];previous=np.stack([initial,initial]);frames=[previous.copy()];last_rng=None
            for segment in range(10):
                with np.load(path/f'checkpoint_{segment:02d}.npz') as z:
                    np.testing.assert_array_equal(z['initial'],previous);final=z['final'];hist=z['histogram']
                    if last_rng is not None:np.testing.assert_array_equal(last_rng,z['rng_before'])
                    last_rng=z['rng_after'].copy()
                for side in range(2):
                    flat=previous[side].ravel().copy()
                    for i,ch in enumerate(contrast.model.channels):
                        for coordinate,v in ch.change:flat[coordinate]+=int(hist[side,i])*v
                    np.testing.assert_array_equal(flat,final[side].ravel());assert flat.min()>=0 and sum(map(int,flat))==sum(map(int,initial.ravel()));endpoints+=1
                if job['sentinel']:
                    with np.load(path/f'trace_{segment:02d}.npz') as archive:t={k:archive[k] for k in archive.files}
                    branch=np.stack([np.bincount(t['channels'][t['masks']==mask],minlength=len(contrast.model.channels)) for mask in (3,1,2)]);np.testing.assert_array_equal(np.stack([branch[0]+branch[1],branch[0]+branch[2]]),hist)
                    result=dict(histogram=branch,counts=np.stack([np.stack([previous[s],final[s]]) for s in range(2)]),trace=[dict(time=float(mark),mask=int(mask),channel=int(channel)) for mark,mask,channel in zip(t['times'],t['masks'],t['channels'])],boundaries=[dict(time=float(mark),mask=int(mask),available=available,output=out) for mark,mask,available,out in zip(t['boundary_times'],t['boundary_masks'],t['boundary_available'],t['boundary_output'])]);replay(contrast.model,*previous,job['deltas'],[0,.001],result);events+=len(t['times']);del result,t
                previous=final.copy();frames.append(previous.copy())
            with np.load(path/'observations.npz') as z:np.testing.assert_array_equal(z['counts'],frames);np.testing.assert_array_equal(z['seeds'],job['seeds'])
            samples[job['comparison']].append([[*(contrast.direct(frames[f][0])-contrast.direct(frames[f][1])),*(power(frames[f][0],initial)-power(frames[f][1],initial))] for f in (5,10)])
        rows=[]
        for c,values in enumerate(samples):
            a=np.asarray(values);assert len(a)==plan['pairs'][c];row=intervals(a,margins);row.update(comparison=c,N=len(a),skewness=skew(a,axis=0,bias=False).tolist());np.savez_compressed(BASE/f'paired_observables_c{c}.npz',values=a,margins=margins,times=[.005,.01]);rows.append(row)
        atomic_json(BASE/'accuracy_results.json',dict(comparisons=rows,margins=margins.tolist(),source_hashes=hashes(),checkpoint_chain_endpoints=endpoints,sentinel_events=events,coverage='Nominal simultaneous95% paired-t/Bonferroni,24comparisons; approximate for nonnormal observations, no full-field equivalence.'));atomic_json(attempt,dict(status='COMPLETE',elapsed_seconds=time.perf_counter()-start))
    except BaseException:
        atomic_json(attempt,dict(status='FAILED',elapsed_seconds=time.perf_counter()-start,error=traceback.format_exc()));raise
if __name__=='__main__':check()

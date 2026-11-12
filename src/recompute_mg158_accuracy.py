"""MG158 separate raw-field/Parseval/interval analysis; no new simulation draws."""
import os
for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):os.environ[k]='1'
import json,time,math,hashlib
import numpy as np
from scipy.stats import t
from mg12h3_tick_accuracy import BASE,read,hashes
from mg12h3_finer_precision import setup,S
from mg10e2c1_pool import atomic_json
def recompute():
    started=time.perf_counter();attempt=BASE/'separate_analysis_attempt.json';atomic_json(attempt,dict(status='RUNNING'))
    try:
        plan=read(BASE/'plan.json');assert hashes()==plan['source_hashes'];initial,contrast=setup();margins=np.r_[.05*abs(contrast.weights@initial.ravel()),[.05]*3];values=[[],[]];means=np.zeros((2,11,2,7,500));inputs={};band=abs(np.fft.fftfreq(500)*500)>=125
        for job in plan['jobs']:
            path=BASE/f"session_{job['session']}"/job['id']/'observations.npz';inputs[path.relative_to(BASE).as_posix()]=hashlib.sha256(path.read_bytes()).hexdigest()
            with np.load(path) as z:
                counts=z['counts'];np.testing.assert_array_equal(z['seeds'],job['seeds']);assert counts.shape==(11,2,7,500)
            means[job['comparison']]+=counts;measure=[]
            for frame in (5,10):
                cs=[contrast.weights@counts[frame,s].ravel() for s in (0,1)];ps=[]
                for side in (0,1):
                    change=(counts[frame,side,:3]-initial[:3])/10.;coeff=np.fft.fft(change,axis=1);ps.append(10/S*np.sum(abs(coeff[:,band])**2,axis=1)/500**2)
                measure.append(np.r_[cs[0]-cs[1],ps[0]-ps[1]])
            values[job['comparison']].append(measure)
        rows=[]
        for c,group in enumerate(values):
            a=np.asarray(group);n=len(a);assert n==plan['pairs'][c];mean=np.sum(a,axis=0)/n;centered=a-mean;variance=np.sum(centered**2,axis=0)/(n-1);critical=float(t.isf(.05/48,n-1));hw=critical*np.sqrt(variance/n);lower=mean-hw;upper=mean+hw;labels=[]
            for f in range(2):
                row=[]
                for j,margin in enumerate(margins):
                    if variance[f,j]==0:label='unresolved'
                    elif lower[f,j]>margin or upper[f,j]<-margin:label='incompatible'
                    elif lower[f,j]>=-margin and upper[f,j]<=margin:label='equivalent'
                    else:label='unresolved'
                    row.append(label)
                labels.append(row)
            flat=sum(labels,[]);overall='equivalent' if all(x=='equivalent' for x in flat) else ('incompatible' if 'incompatible' in flat else 'unresolved');denom=np.sum(centered**2,axis=0);share=np.divide(np.max(centered**2,axis=0),denom,out=np.zeros_like(denom),where=denom>0);fourth=np.mean(centered**4,axis=0);second=np.mean(centered**2,axis=0);kurt=np.divide(fourth,second**2,out=np.zeros_like(second),where=second>0)-3
            rows.append(dict(comparison=c,N=n,mean=mean.tolist(),variance=variance.tolist(),critical=critical,halfwidth=hw.tolist(),lower=lower.tolist(),upper=upper.tolist(),labels=labels,overall=overall,max_variance_share=share.tolist(),population_moment_excess_kurtosis=kurt.tolist()));means[c]/=n
        np.savez_compressed(BASE/'mean_observation_fields.npz',mean_counts=means,times=np.arange(11)*.001,comparison_ticks=[.0001,.00005],reference_tick=.00002)
        atomic_json(BASE/'separate_accuracy_recomputation.json',dict(passed=True,comparisons=rows,margins=margins.tolist(),input_sha256=inputs,elapsed_seconds=time.perf_counter()-started,method='Raw endpoint projections and Fourier Parseval powers; centered sums of squares and inverse-survival t criticals. Same author/shared frozen input, weights and SciPy library; not external replication.'));atomic_json(attempt,dict(status='COMPLETE',elapsed_seconds=time.perf_counter()-started));print(json.dumps([dict(comparison=r['comparison'],overall=r['overall'],labels=r['labels']) for r in rows]))
    except BaseException:
        import traceback;atomic_json(attempt,dict(status='FAILED',elapsed_seconds=time.perf_counter()-started,error=traceback.format_exc()));raise
if __name__=='__main__':recompute()

"""MG160 exploratory saved-data diagnostics; never replaces frozen per-path scoring."""
import os
for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):os.environ[k]='1'
import time,traceback
import numpy as np
from run_mg160_persistence import BASE,read
from check_mg160_pilot import recomputer
from mg10e2c1_pool import atomic_json
from mg160_local_archive import local_archive
BASE=local_archive(BASE)
def diagnose():
    started=time.perf_counter();attempt=BASE/'P1_diagnostic_attempt.json';atomic_json(attempt,dict(status='RUNNING'))
    try:
        analysis=read(BASE/'P1_analysis_results.json');assert analysis['passed'];plan=read(BASE/'P1_plan.json')
        with np.load(BASE/'prepared_inputs.npz') as z:initial=z['initial'];reference=z['reference'];omega=float(z['omega'])
        measure=recomputer(reference,initial,omega);mean_fields=[];rows=[];stock=[]
        for condition in ('prepared','homogeneous'):
            group=[j for j in plan['jobs'] if j['condition']==condition];arrays=[]
            for job in group:
                with np.load(BASE/'P1_production'/job['id']/'observations.npz') as z:arrays.append(z['counts'])
            mean=np.mean(arrays,axis=0);mean_fields.append(mean);metrics=[[measure(pop) for pop in frame] for frame in mean];rows.append(dict(condition=condition,metrics=metrics,final=metrics[-1],independent_pairs=4))
            for job,a in zip(group,arrays):
                totals=a.sum(axis=3)
                for side in range(2):stock.append(dict(id=job['id'],side=side,A_depleted_fraction=float((totals[0,side,3]-totals[-1,side,3])/totals[0,side,3]),B_depleted_fraction=float((totals[0,side,4]-totals[-1,side,4])/totals[0,side,4])))
        np.savez_compressed(BASE/'P1_exploratory_mean_fields.npz',mean_counts=mean_fields,times=np.arange(101)*.001,conditions=['prepared','homogeneous']);atomic_json(BASE/'P1_exploratory_diagnostics.json',dict(label='Exploratory post-result diagnostics, not a new qualification or rescore.',conditions=rows,stocks=stock,method='Average raw unaligned fields across four independent pairs within each condition and tick, then apply the same measurement. This changes the estimand to an ensemble-mean field and cannot establish individual survival.',limits='Visible ensemble structure is consistent with noise obscuring per-trajectory shape but does not identify a unique mechanism; genuine deformation and local stock redistribution remain possible.',elapsed_seconds=time.perf_counter()-started));atomic_json(attempt,dict(status='COMPLETE',elapsed_seconds=time.perf_counter()-started));print([dict(condition=r['condition'],final=r['final']) for r in rows])
    except BaseException:
        atomic_json(attempt,dict(status='FAILED',elapsed_seconds=time.perf_counter()-started,error=traceback.format_exc()));raise
if __name__=='__main__':diagnose()

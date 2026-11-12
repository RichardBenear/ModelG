"""MG161 exploratory P1 application only after synthetic benchmark checks pass."""
import os
for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):os.environ[k]='1'
import time
import numpy as np
from run_mg161_shape import BASE,ARCHIVE,read,setup
from mg10e2c1_pool import atomic_json
def apply():
    start=time.perf_counter();assert read(BASE/'benchmark_results.json')['qualified_within_synthetic_scope'] and read(BASE/'separate_recomputation.json')['passed'];model,_,_=setup();cal=read(BASE/'calibration.json');plan=read(ARCHIVE/'P1_plan.json');features=[];labels=[];rows=[]
    for job in plan['jobs']:
        with np.load(ARCHIVE/'P1_production'/job['id']/'observations.npz') as z:counts=z['counts'];times=z['times']
        job_features=[];job_labels=[]
        for frame in counts:
            f=[];lab=[]
            for pop in frame:
                m=model.fit(pop[1]/model.omega);f.append([m[k] for k in ('amplitude','scale','center_site','shift_sites','score','residual','secondary_score','secondary_ratio','bank_edge')]);lab.append(model.classify(m,cal['detection_threshold'],cal['residual_threshold']))
            job_features.append(f);job_labels.append(lab)
        features.append(job_features);labels.append(job_labels);outcomes=[]
        for s in (0,1):
            path=[f[s] for f in job_labels];bad=[i for i,l in enumerate(path) if l!='core_compatible'];outcomes.append(dict(side=s,compatible_frames=path.count('core_compatible'),all_frames_compatible=not bad,first_noncompatible_time=float(times[bad[0]]) if bad else None,final_label=path[-1],final_amplitude=job_features[-1][s][0]/model.initial_amplitude,final_scale=job_features[-1][s][1]))
        rows.append(dict(id=job['id'],condition=job['condition'],outcomes=outcomes))
    np.savez_compressed(BASE/'P1_exploratory_features.npz',features=features,labels=labels,job_ids=[j['id'] for j in plan['jobs']],times=times);atomic_json(BASE/'P1_exploratory_application.json',dict(exploratory=True,jobs=rows,noise_scale_fixed=1.,scope='Post-result new-instrument compatibility on previously seen P1 data. Not held-out real trajectories, not a persistence rescore; original0/4continuous-retention result stays unchanged.',elapsed_seconds=time.perf_counter()-start));print(rows)
if __name__=='__main__':apply()

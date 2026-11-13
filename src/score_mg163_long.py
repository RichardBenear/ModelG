import time,traceback
import numpy as np
from run_mg163_long import PUBLIC as BASE,BASE as RAW,ROOT,read,hashes,atomic_json
from mg161_local_shape import LocalShape
def score():
    started=time.perf_counter();attempt=BASE/'fresh_score_attempt.json';assert not attempt.exists();atomic_json(attempt,dict(status='RUNNING'))
    try:
        plan=read(BASE/'plan.json');assert hashes()==plan['source_hashes'] and read(BASE/'attempt.json')['status']=='COMPLETE';assert read(RAW/'P1_analysis_results.json')['passed'];cal=read(ROOT/'results/MG12H3/noise_aware_shape/calibration.json')
        with np.load(RAW/'prepared_inputs.npz') as z:model=LocalShape(z['reference'][1],z['initial'][1],float(z['omega']))
        rows=[]
        for j in plan['jobs']:
            with np.load(RAW/'production'/j['id']/'observations.npz') as z:counts=z['counts'];times=z['times'];assert counts.shape==(1001,2,7,500) and int(z['seed'])==j['seed']
            sides=[]
            for s in (0,1):
                fits=[model.fit(frame[s,1]/model.omega) for frame in counts];labels=[model.classify(f,cal['detection_threshold'],cal['residual_threshold']) for f in fits];bad=[i for i,l in enumerate(labels) if l!='core_compatible'];edge=[i for i,f in enumerate(fits) if labels[i]=='core_compatible' and (f['amplitude']/model.initial_amplitude<=.6 or min(abs(f['scale']-.5),abs(f['scale']-2))<=.125)];sides.append(dict(side=s,compatible_frames=labels.count('core_compatible'),all_compatible=not bad,any_compatible='core_compatible' in labels,first_noncompatible_bracket=[float(times[max(0,bad[0]-1)]),float(times[bad[0]])] if bad else None,boundary_advisory_frames=edge,labels=labels,fits=fits))
            rows.append(dict(id=j['id'],condition=j['condition'],outcomes=sides))
        supported=all(s['all_compatible'] for r in rows if r['condition']=='prepared' for s in r['outcomes']) and all(not s['any_compatible'] for r in rows if r['condition']=='homogeneous' for s in r['outcomes']) and all(not s['boundary_advisory_frames'] for r in rows if r['condition']=='prepared' for s in r['outcomes']);atomic_json(BASE/'fresh_candidate_scores.json',dict(pilot_separation_supported=supported,jobs=rows,source_hashes=hashes(),qualification_status='Candidate scores only; causal/raw audit and separate metric recomputation required before final report.',boundary_advisory='Engineering sensitivity flags only; no fitted-parameter confidence intervals and no hidden score changes.',elapsed_seconds=time.perf_counter()-started));atomic_json(attempt,dict(status='COMPLETE',elapsed_seconds=time.perf_counter()-started))
    except BaseException:atomic_json(attempt,dict(status='FAILED',elapsed_seconds=time.perf_counter()-started,error=traceback.format_exc()));raise

if __name__=='__main__':score()

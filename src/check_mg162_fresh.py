"""Separate direct-sum recomputation of saved MG162B fitted metrics and frozen verdict."""
import time, traceback
import numpy as np
from run_mg162_fresh import BASE, RAW, ROOT, read, atomic_json, LocalShape

def check():
    start=time.perf_counter();attempt=BASE/'fresh_recompute_attempt.json'
    assert not attempt.exists();atomic_json(attempt,dict(status='RUNNING'))
    try:
        scores=read(BASE/'fresh_candidate_scores.json');cal=read(ROOT/'results/MG12H3/noise_aware_shape/calibration.json')
        with np.load(RAW/'prepared_inputs.npz') as z:m=LocalShape(z['reference'][1],z['initial'][1],float(z['omega']))
        kernels=[];norms=[]
        for template in m.bank:
            for lag in range(m.n):
                w=np.roll(m.window,lag);g=np.roll(template,lag);g=g-w@g/m.ws;kernels.append(w*g);norms.append(w@(g*g))
        kernels=np.asarray(kernels);norms=np.asarray(norms)
        rows=[];checked=0
        for job in scores['jobs']:
            with np.load(RAW/'P1_production'/job['id']/'observations.npz') as z:counts=z['counts'];times=z['times']
            sides=[]
            for side,stored in enumerate(job['outcomes']):
                labels=[];edges=[]
                for t,fit in enumerate(stored['fits']):
                    f=np.fft.irfft(np.fft.rfft(counts[t,side,1]/m.omega-(counts[t,side,1]/m.omega).mean())*m.smoother,n=m.n)
                    k=fit['center_site'];si=int(np.flatnonzero(m.scales==fit['scale'])[0]);w=np.roll(m.window,k);g=np.roll(m.bank[si],k);fm=w@f/m.ws;gm=w@g/m.ws
                    dot=w@((f-fm)*(g-gm));norm=w@((g-gm)**2);amplitude=dot/norm;gain=max(0,dot)**2/norm;residual=w@((f-fm-amplitude*(g-gm))**2)
                    np.testing.assert_allclose([amplitude,gain,residual],[fit['amplitude'],fit['score'],fit['residual']],rtol=1e-8,atol=1e-8)
                    # Direct weighted sums at every scale/position verify global and secondary search too.
                    dots=kernels@f;gains=(np.maximum(dots,0)**2/norms).reshape(len(m.scales),m.n)
                    np.testing.assert_allclose(gains.max(),fit['score'],rtol=1e-8,atol=1e-8)
                    dist=np.minimum((np.arange(m.n)-k)%m.n,(k-np.arange(m.n))%m.n)*m.h;second=gains[:,dist>m.background_exclusion].max()
                    np.testing.assert_allclose(second,fit['secondary_score'],rtol=1e-8,atol=1e-8)
                    if gain<=cal['detection_threshold'] or amplitude/m.initial_amplitude<.5:label='absent_or_weak'
                    elif second>cal['detection_threshold'] and second/max(gain,1e-30)>=.5:label='ambiguous_multiple_cores'
                    elif si in (0,len(m.scales)-1) or not .5<=fit['scale']<=2:label='deformed_scale'
                    elif residual>cal['residual_threshold']:label='shape_mismatch'
                    else:label='core_compatible'
                    labels.append(label)
                    if label=='core_compatible' and (amplitude/m.initial_amplitude<=.6 or min(abs(fit['scale']-.5),abs(fit['scale']-2))<=.125):edges.append(t)
                    checked+=1
                assert labels==stored['labels'] and edges==stored['boundary_advisory_frames']
                bad=[i for i,v in enumerate(labels) if v!='core_compatible'];bracket=[float(times[max(0,bad[0]-1)]),float(times[bad[0]])] if bad else None
                assert bracket==stored['first_noncompatible_bracket'] and (not bad)==stored['all_compatible'] and ('core_compatible' in labels)==stored['any_compatible']
                sides.append(dict(all_compatible=not bad,any_compatible='core_compatible' in labels,boundary_frames=len(edges)))
            rows.append(dict(id=job['id'],condition=job['condition'],outcomes=sides))
        supported=all(s['all_compatible'] and not s['boundary_frames'] for r in rows if r['condition']=='prepared' for s in r['outcomes']) and all(not s['any_compatible'] for r in rows if r['condition']=='homogeneous' for s in r['outcomes'])
        assert supported==scores['pilot_separation_supported']
        atomic_json(BASE/'fresh_separate_recomputation.json',dict(passed=True,fields_checked=checked,pilot_separation_supported=supported,jobs=rows,elapsed_seconds=time.perf_counter()-start,method='Direct spatial weighted sums over all scale/position fits; independent classification, advisory and outcome aggregation. Shared template bank, smoothing and input fields; same author, not external replication.'))
        atomic_json(attempt,dict(status='COMPLETE',elapsed_seconds=time.perf_counter()-start))
    except BaseException:
        atomic_json(attempt,dict(status='FAILED',elapsed_seconds=time.perf_counter()-start,error=traceback.format_exc()));raise
if __name__=='__main__':check()

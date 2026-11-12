"""MG161 saved-field scalar-fit and independent binomial-bound verification."""
import os
for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):os.environ[k]='1'
import time,math,hashlib
import numpy as np
from scipy.stats import binom
from scipy.optimize import brentq
from run_mg161_shape import BASE,RAW,ARCHIVE,read,setup,hashes
from mg10e2c1_pool import atomic_json
def check():
    start=time.perf_counter();model,psd,_=setup();plan=read(BASE/'plan.json');actual=read(BASE/'benchmark_results.json');cal=read(BASE/'calibration.json');assert hashes()==plan['source_hashes'];scalar_checks=0;bounds=[];raw_manifest={}
    for p,digest in plan['inputs'].items():assert hashlib.sha256((ARCHIVE/p).read_bytes()).hexdigest()==digest
    cal_features=[]
    for case in plan['cases']:
        path=RAW/(case['id']+'.npz');raw_manifest[path.name]=dict(bytes=path.stat().st_size,sha256=hashlib.sha256(path.read_bytes()).hexdigest())
        with np.load(path) as z:fields=z['fields'];features=z['features'];assert int(z['seed'])==case['seed']
        assert len(fields)==case['N'];sigma=case.get('noise_scale',1.)
        with np.load(BASE/(case['id']+'_features.npz')) as z:np.testing.assert_array_equal(z['features'],features)
        for r in (0,len(fields)//2,len(fields)-1):
            v=features[r];scale=v[1];k=int(v[2]);s=int(np.flatnonzero(model.scales==scale)[0]);f=np.fft.irfft(np.fft.rfft(fields[r]-fields[r].mean())*model.smoother,n=500);w=np.roll(model.window,k);g=np.roll(model.bank[s],k);fm=float(w@f/model.ws);gm=float(w@g/model.ws);norm=float(w@((g-gm)**2));dot=float(w@((f-fm)*(g-gm)));amp=dot/norm;gain=max(0,dot)**2/norm;residual=float(w@((f-fm-amp*(g-gm))**2));np.testing.assert_allclose([amp,gain/sigma**2,residual/sigma**2],[v[0],v[4],v[5]],rtol=1e-9,atol=1e-9);scalar_checks+=1
        if case['id'].startswith('cal_'):cal_features.append(features);continue
        labels=[]
        for v in features:
            if v[4]<=cal['detection_threshold'] or v[0]/model.initial_amplitude<.5:label='absent_or_weak'
            elif v[6]>cal['detection_threshold'] and v[7]>=.5:label='ambiguous_multiple_cores'
            elif v[8] or not .5<=v[1]<=2:label='deformed_scale'
            elif v[5]>cal['residual_threshold']:label='shape_mismatch'
            else:label='core_compatible'
            labels.append(label)
        row=next(r for r in actual['cases'] if r['case']['id']==case['id']);k=labels.count('core_compatible');n=len(labels);assert k==row['compatible'];assert {s:labels.count(s) for s in set(labels)}==row['label_counts']
        if case.get('primary'):
            alpha=plan['bound_alpha'];lower=brentq(lambda p:binom.sf(k-1,n,p)-alpha,1e-14,1-1e-14,xtol=1e-13) if k>0 else 0.;upper=brentq(lambda p:binom.cdf(k,n,p)-alpha,1e-14,1-1e-14,xtol=1e-13) if k<n else 1.;np.testing.assert_allclose([lower,upper],[row['lower'],row['upper']],rtol=1e-9,atol=1e-11);passed=lower>=.9 if case['truth']=='compatible' else upper<=(.01 if case['id']=='null' else .05);assert passed==row['passed'];bounds.append(dict(case=case['id'],lower=lower,upper=upper,passed=bool(passed)))
    assert cal['detection_threshold']==float(cal_features[0][:,4].max()) and cal['residual_threshold']==float(cal_features[1][:,5].max());assert all(b['passed'] for b in bounds)==actual['qualified_within_synthetic_scope'];atomic_json(BASE/'local_raw_manifest.json',dict(local_directory='local_raw/MG161_noise_shape',files=raw_manifest,total_bytes=sum(p['bytes'] for p in raw_manifest.values()),retention='Actual full synthetic fields retained locally, not public.'));atomic_json(BASE/'separate_recomputation.json',dict(passed=True,scalar_fit_checks=scalar_checks,bounds=bounds,source_hashes=hashes(),elapsed_seconds=time.perf_counter()-start,method='Direct weighted local means/dot products/residual sums; separate threshold logic; binomial tail inversion instead of beta quantiles. Same author/shared fields and template bank.'));print('Verified',scalar_checks,'scalar fits and',len(bounds),'exact-binomial case bounds.')
if __name__=='__main__':check()

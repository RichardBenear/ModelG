"""MG162A saved-field direct sums and exact finite-cohort recomputation."""
import os
for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):os.environ[k]='1'
import time,hashlib
import numpy as np
from run_mg162_challenge import BASE,RAW,ARCHIVE,read,hashes
from run_mg161_shape import setup,ROOT
from mg10e2c1_pool import atomic_json
def check():
    start=time.perf_counter();p=read(BASE/'plan.json');a=read(BASE/'challenge_results.json');assert p['source_hashes']==hashes();m,_,_=setup();cal=read(ROOT/'results/MG12H3/noise_aware_shape/calibration.json');rows=[];manifest={};fit_checks=0
    for name,digest in p['inputs'].items():assert hashlib.sha256((ARCHIVE/name).read_bytes()).hexdigest()==digest
    for c in p['cases']:
        path=RAW/(c['id']+'.npz');manifest[path.name]=dict(bytes=path.stat().st_size,sha256=hashlib.sha256(path.read_bytes()).hexdigest())
        with np.load(path) as z:fields=z['fields'];features=z['features'];saved=z['labels'];assert fields.shape==(2000,500)
        with np.load(BASE/(c['id']+'_features.npz')) as z:np.testing.assert_array_equal(z['features'],features);np.testing.assert_array_equal(z['labels'],saved)
        for r in (0,999,1999):
            v=features[r];i=int(np.flatnonzero(m.scales==v[1])[0]);lag=int(v[2]);f=np.fft.irfft(np.fft.rfft(fields[r]-fields[r].mean())*m.smoother,n=500);w=np.roll(m.window,lag);g=np.roll(m.bank[i],lag);fm=w@f/m.ws;gm=w@g/m.ws;norm=w@((g-gm)**2);dot=w@((f-fm)*(g-gm));amp=dot/norm;gain=max(0,dot)**2/norm;residual=w@((f-fm-amp*(g-gm))**2);np.testing.assert_allclose([amp,gain,residual],[v[0],v[4],v[5]],rtol=1e-9,atol=1e-9);fit_checks+=1
        labels=[]
        for v in features:
            if v[4]<=cal['detection_threshold'] or v[0]/m.initial_amplitude<.5:s='absent_or_weak'
            elif v[6]>cal['detection_threshold'] and v[7]>=.5:s='ambiguous_multiple_cores'
            elif v[8] or not .5<=v[1]<=2:s='deformed_scale'
            elif v[5]>cal['residual_threshold']:s='shape_mismatch'
            else:s='core_compatible'
            labels.append(s)
        np.testing.assert_array_equal(labels,saved);record=next(r for r in a['cases'] if r['case']['id']==c['id']);fractions=[labels[i:i+500].count('core_compatible')/500 for i in (0,500,1000,1500)];np.testing.assert_array_equal(fractions,[g['fraction'] for g in record['groups']]);assert record['compatible']==labels.count('core_compatible')
        if c['primary']:
            passed=min(fractions)>=.9 if c['truth']=='compatible' else max(fractions)<=(.01 if c['id']=='null' else .05);assert passed==record['passed'];rows.append(dict(case=c['id'],passed=passed))
    assert all(r['passed'] for r in rows)==a['supported_within_finite_challenge'];atomic_json(BASE/'separate_recomputation.json',dict(passed=True,scalar_fit_checks=fit_checks,primary_cases=rows,elapsed_seconds=time.perf_counter()-start,scope='Separate weighted scalar sums, classification branches and exhaustive finite-group fractions; same author/shared fields/bank, no independent-trajectory population inference.'));atomic_json(BASE/'local_raw_manifest.json',dict(local_directory='local_raw/MG162_shape_challenge',files=manifest,total_bytes=sum(v['bytes'] for v in manifest.values()),full_fields_published=False));print('Verified',fit_checks,'direct fits and15finite-case verdicts.')
if __name__=='__main__':check()

"""MG160 deterministic/saved-only measurement and cost preparation."""
import os
for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):os.environ[k]='1'
import json,hashlib,time
from pathlib import Path
import numpy as np
from mg12h3_pattern_persistence import Coherence,homogenize,continuity
from mg12h3_finer_precision import setup,ROOT
from mg10e2c1_pool import atomic_json
BASE=ROOT/'results/MG12H3/pattern_persistence'
def prepare():
    started=time.perf_counter();n=500;grid=np.arange(n);well=np.exp(-((grid-n//2)/9.)**2);reference=np.zeros((7,n));reference[1]=10-well;detector=Coherence(reference,reference,1.);checks=[]
    def test(name,condition):
        if not condition:raise AssertionError(name)
        checks.append(name)
    test('reference retained',detector.measure(reference)['coherent'])
    shifted=np.roll(reference,73,axis=1);m=detector.measure(shifted);test('translation retained with correct lag',m['coherent'] and m['shift_sites']==73)
    low=reference.copy();low[1]=10-.4*well;test('amplitude loss rejected',not detector.measure(low)['coherent'])
    flat=reference.copy();flat[1]=10;test('flat rejected',not detector.measure(flat)['coherent'])
    inverse=reference.copy();inverse[1]=10+well;test('inverted localized feature rejected',not detector.measure(inverse)['coherent'])
    doubled=reference.copy();doubled[1]=10-well-np.roll(well,150);test('split profile rejected',not detector.measure(doubled)['coherent'])
    broad=reference.copy();broad[1]=10-np.exp(-((grid-n//2)/22.5)**2);test('broadening distinguishes shape from coherence',detector.measure(broad)['coherent'] and not detector.measure(broad)['pattern_retained'])
    offset=reference.copy();offset[1]+=2;test('uniform background shift invariant',abs(detector.measure(offset)['relative_amplitude']-1)<1e-12)
    periodic=reference.copy();periodic[1]=10+np.cos(2*np.pi*grid/100);periodic_detector=Coherence(periodic,periodic,1.);test('periodic phase ambiguity flagged',periodic_detector.measure(periodic)['phase_ambiguous'])
    test('sampled retention right censored',continuity([0,1],[detector.measure(reference)]*2)['status']=='sampled coherence retained');test('loss interval bracket',continuity([0,1],[detector.measure(reference),detector.measure(flat)])['first_loss_interval']==[0.,1.])
    try:continuity([0,np.nan],[detector.measure(reference)]*2)
    except ValueError:checks.append('nonfinite time rejected')
    else:raise AssertionError('Missing time guard')
    initial,contrast=setup();path=ROOT/'results/MG12H3/ring_feasibility/o0_a0_Omega10.npz'
    with np.load(path) as z:ref=z['concentrations'];omega=float(z['omega']);a0=z['a0'];b0=z['b0']
    real=Coherence(ref,initial,omega);control=homogenize(initial);test('mass-matched homogeneous control',np.array_equal(initial.sum(1),control.sum(1)) and max(np.ptp(control,axis=1))<=1);test('homogeneous control initially rejected',not real.measure(control)['coherent'])
    np.savez_compressed(BASE/'prepared_inputs.npz',initial=initial,homogeneous=control,reference=ref,omega=omega,template=real.template)
    # Saved snapshots are development checks only, never new survival observations.
    saved=[];inputs={path.relative_to(ROOT).as_posix():hashlib.sha256(path.read_bytes()).hexdigest()}
    for c in (0,1):
        p=ROOT/f'results/MG12H3/tick_accuracy/session_1/c{c}_paired_r000/observations.npz';inputs[p.relative_to(ROOT).as_posix()]=hashlib.sha256(p.read_bytes()).hexdigest()
        with np.load(p) as z:counts=z['counts']
        saved.append(dict(comparison=c,metrics=[[real.measure(counts[f,s]) for s in (0,1)] for f in (0,5,10)]))
    old=ROOT/'results/MG11G3/production/o0_s1_v1_a0.npz';inputs[old.relative_to(ROOT).as_posix()]=hashlib.sha256(old.read_bytes()).hexdigest()
    with np.load(old) as z:
        old_times=z['times'];old_fields=z['state'][:,:3500].reshape(-1,7,500);old_metrics=[real.measure(field*omega) for field in old_fields]
    reference_context=dict(kind='archived deterministic population-limit fields; descriptive new diagnostic, no rescore',times=old_times.tolist(),metrics=old_metrics,continuity=continuity(old_times,old_metrics),shape_continuity=continuity(old_times,old_metrics,key='pattern_retained'))
    measured=json.loads((ROOT/'results/MG12H3/tick_accuracy/budget_summary.json').read_text())['worker_seconds']/320
    forecast=8*10*measured*1.15;projection15=8*1500*measured/8/86400
    atomic_json(BASE/'preparation_results.json',dict(passed=True,tests=checks,initial_metrics=real.measure(initial),homogeneous_metrics=real.measure(control),saved_prefix_checks=saved,archived_reference_context=reference_context,stock_reference_ranges=dict(A=[float(a0.min()),float(a0.max())],B=[float(b0.min()),float(b0.max())]),measured_pair_seconds_per_0_01=measured,candidate_P1=dict(duration=.1,prepared_pairs=4,homogeneous_pairs=4,ticks=[.0001,.00002],observation=.001,forecast_worker_seconds=forecast,forecast_wall_hours_at8=forecast/8/3600,proposed_worker_cap=100000,proposed_parent_cap=1800,proposed_job_cap=12000),linear_15unit_8pair_wall_days=projection15,input_hashes=inputs,random_draws=0,integrations=0,elapsed_seconds=time.perf_counter()-started,limits='Deterministic fixture/saved-prefix checks only; not empirical false-positive calibration or a persistence outcome.'))
    print(json.dumps(dict(tests=len(checks),prepared=real.measure(initial),homogeneous=real.measure(control),P1_forecast_hours=forecast/8/3600,full_15unit_days=projection15)))
if __name__=='__main__':BASE.mkdir(parents=True,exist_ok=True);prepare()

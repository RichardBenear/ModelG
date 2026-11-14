"""Saved-only detectability and averaging; competing-core candidate assessment."""
import os
for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):os.environ[k]='1'
from pathlib import Path
import time,json,hashlib,traceback
import numpy as np
from scipy.stats import norm
from mg165_joint import Joint
from run_mg163_long import ROOT,BASE as OLD,read,atomic_json
OUT=ROOT/'results/MG12H3/detectability_and_veto';SYN=ROOT/'local_raw/MG165_multivariate_v2';PRIOR=ROOT/'results/MG12H3/multivariate_instrument/v2'
def sources():return {p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in ['spec/MG166_detectability_and_veto.md','src/analyze_mg166.py','src/mg165_joint.py']}
def register():
    OUT.mkdir(exist_ok=True);assert not(OUT/'plan.json').exists();paths=[OLD/'prepared_inputs.npz',PRIOR/'benchmark.json']+[OLD/'production'/f'homogeneous_r{r}'/'observations.npz' for r in range(4)]+[SYN/f'r{i}_{case}.npz' for i in range(3) for case in ['null','split','intact_0','intact_1','intact_2']];atomic_json(OUT/'plan.json',dict(source=sources(),inputs={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},cap=600,no_new_draws=True))
def run():
    start=time.perf_counter();p=read(OUT/'plan.json');assert p['source']==sources() and not(OUT/'attempt.json').exists();atomic_json(OUT/'attempt.json',dict(status='RUNNING'))
    try:
        for name,h in p['inputs'].items():assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==h
        with np.load(OLD/'prepared_inputs.npz') as z:ref=z['reference'];initial=z['initial'];omega=float(z['omega'])
        benchmark=read(PRIOR/'benchmark.json');detect=[];averaging=[];veto=[];projectionchecks=0;artifacts={};target=float(norm.ppf(.95)+norm.ppf(.9));assert abs(norm.cdf(target-norm.ppf(.95))-.9)<1e-12
        for regime,rr in enumerate(benchmark['regime_inputs']):
            model=Joint(ref,initial,omega,np.array(rr['covariance']));m=model.x;g=model.g[int(np.flatnonzero(m.scales==1)[0])];contrast=np.zeros_like(g);contrast[0]=g[0];kernels=[contrast*m.window, (model.inv@contrast)*m.window];gap=[float((contrast*k).sum()) for k in kernels]
            with np.load(SYN/f'r{regime}_null.npz') as z:null=z['fields']
            for kind,k,delta in zip(['G_only','G_conditional_GXY'],kernels,gap):
                scores=np.einsum('rsn,sn->r',null,k);std=float(np.std(scores[:512],ddof=1));checkstd=float(np.std(scores[512:],ddof=1));snr=delta/std;count=int(np.ceil((target/snr)**2));detect.append(dict(regime=rr['time'],projection=kind,signal_difference=delta,calibration_std=std,validation_std=checkstd,SNR=snr,ideal_independent_count=count,gaussian_power_one_frame=float(norm.cdf(snr-norm.ppf(.95)))))
                np.testing.assert_allclose(scores[[0,511,1023]],[sum(float(null[r,s]@k[s]) for s in range(3)) for r in(0,511,1023)],rtol=1e-12,atol=1e-12);projectionchecks+=3
            k=kernels[1];positions=np.floor(np.arange(16)*500/16).astype(int);timeslices=[(50,151),(450,551),(900,1001)];a,b=timeslices[regime];allvalues=[]
            for r in range(4):
                with np.load(OLD/'production'/f'homogeneous_r{r}'/'observations.npz') as z:f=z['counts'][a:b,0,:3]/omega
                f=f-f.mean(axis=2)[:,:,None];smooth=np.fft.irfft(np.fft.rfft(f,axis=2)*m.smoother,n=500,axis=2);projection=np.fft.irfft((np.fft.rfft(smooth,axis=2)*np.conj(np.fft.rfft(k,axis=1))[None,:,:]).sum(axis=1),n=500,axis=1)[:,positions].T;allvalues.append(projection)
                for center in(0,7,15):
                    val=float(np.sum(smooth[0]*np.roll(k,int(positions[center]),axis=1)));np.testing.assert_allclose(val,projection[center,0],rtol=1e-10,atol=1e-10);projectionchecks+=1
            allvalues=np.array(allvalues);artifacts[f'r{regime}_projections']=allvalues
            for group,values in [('development_r0r1',allvalues[:2]),('seen_check_r2r3',allvalues[2:])]:
                baseline=float(np.var(values,ddof=1));drift=float(np.sqrt(np.mean((values[:,:,51:].mean(axis=2)-values[:,:,:50].mean(axis=2))**2)))
                for length in(1,2,5,10,20,50):
                    nblocks=101//length;means=values[:,:,:nblocks*length].reshape(2,16,nblocks,length).mean(axis=-1);variance=float(np.var(means,ddof=1));averaging.append(dict(regime=rr['time'],group=group,M=length,blocks_per_center=nblocks,variance=variance,variance_reduction=baseline/variance,standard_deviation_gain=float(np.sqrt(baseline/variance)),independent_reference_gain=float(np.sqrt(length)),half_window_RMS_drift=drift,baseline_std=float(np.sqrt(baseline))))
            for case in('intact_0','intact_1','intact_2','split'):
                with np.load(SYN/f'r{regime}_{case}.npz') as z:fields=z['fields']
                flags=[];secondrat=[]
                for field in fields:
                    f=field-field.mean(axis=1)[:,None];dots=np.fft.irfft((np.fft.rfft(f,axis=1)[None,:,:]*model.kfft).sum(axis=1),n=500,axis=1);gains=np.maximum(dots,0)**2/model.norm[:,None];si,center=np.unravel_index(gains.argmax(),gains.shape);distance=np.minimum((np.arange(500)-center)%500,(center-np.arange(500))%500)*m.h;second=float(gains[:,distance>m.background_exclusion].max());first=float(gains[si,center]);ratio=second/max(first,1e-30);flags.append(second>benchmark['calibrations'][regime]['score'] and ratio>=.5);secondrat.append(ratio)
                    if time.perf_counter()-start>600:raise TimeoutError('MG166 active600s cap')
                fraction=sum(flags)/len(flags);passed=fraction>=.9 if case=='split' else fraction<=.05;veto.append(dict(regime=rr['time'],case=case,n=len(flags),vetoes=int(sum(flags)),fraction=fraction,screen_passed=passed));artifacts[f'r{regime}_{case}_veto']=np.array(flags);artifacts[f'r{regime}_{case}_ratios']=np.array(secondrat)
        artifacts['positions']=positions;np.savez_compressed(OUT/'compact_projections_and_veto.npz',**artifacts);atomic_json(OUT/'analysis.json',dict(detectability=detect,averaging=averaging,veto=veto,veto_development_screen_passed=all(r['screen_passed'] for r in veto),gaussian_target_separation=target,projection_checks=projectionchecks,no_new_draws=True,elapsed_seconds=time.perf_counter()-start));atomic_json(OUT/'attempt.json',dict(status='COMPLETE',elapsed_seconds=time.perf_counter()-start))
    except BaseException:atomic_json(OUT/'attempt.json',dict(status='FAILED',elapsed_seconds=time.perf_counter()-start,error=traceback.format_exc()));raise
if __name__=='__main__':
    if __import__('sys').argv[1]=='register':register()
    else:run()

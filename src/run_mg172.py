"""Independent synthetic validation of exact-current-prior, compressed-history tracker."""
import os
for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):os.environ[k]='1'
from pathlib import Path
import time,hashlib,traceback,json
import numpy as np
from scipy.special import logsumexp
from mg172_corrected import Corrected
from analyze_mg170 import Beam,ROOT,OLD,BENCH,CLASSES,read,atomic_json
from mg167_xy import XYCore
from run_mg169_tracking import truth
OUT=ROOT/'results/MG12H3/corrected_model_selection';RAW=ROOT/'local_raw/MG172_corrected'
FILES=['src/mg172_corrected.py','src/run_mg172.py','src/analyze_mg170.py','src/mg167_xy.py','src/run_mg169_tracking.py','src/mg161_local_shape.py','spec/MG172_corrected_model_selection.md']
def binding():return {p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in FILES}
def register():
    OUT.mkdir(exist_ok=True);assert not(OUT/'plan.json').exists();paths=[OLD/'prepared_inputs.npz',BENCH]+[OLD/'production'/f'homogeneous_r{r}'/'observations.npz' for r in(0,1)];atomic_json(OUT/'plan.json',dict(source=binding(),inputs={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},seed=17200,cap=600,workers=1,threads=1,frames=24,replicates=4,regimes=[.1,1.],rho=[0.,.95],cases=['static','moving','deforming','absent','vanish','missingX','missingY','displacedY','split','jump'],no_ModelG_draws=True))
def fixture_check(tracker):
    s=tracker.audit_snapshot;key=next(iter(s['ends']['multiple']));k,si=key;eligible=np.broadcast_to(np.minimum((np.arange(500)-k)%500,(k-np.arange(500))%500)[None,:]*.2>tracker.model.x.background_exclusion,(13,500));N=eligible.sum();checks=0;rho=0 if getattr(tracker,'fixture_first',False) else tracker.rho
    for secondary,ss in [((k+100)%500,0),((k+250)%500,6),((k+400)%500,12)]:
        mu=tracker.mean('multiple',k,si,(secondary,ss));values=[]
        for parent,(_,_,nloc) in zip(s['old'],s['parentinfo']):
            if min((k-parent['center'])%500,(parent['center']-k)%500)>2 or abs(si-parent['si'])>1:continue
            err=s['innovation']-(mu-rho*parent['mean']);q=float(tracker.inner(err,err)/s['variance']);values.append(parent['log']+np.log(.98 if parent['kind']=='multiple' else .005)-np.log(nloc)-np.log(N)-.5*q)
        expected=float(logsumexp(values));np.testing.assert_allclose(expected,s['ends']['multiple'][key][ss,secondary],rtol=1e-9,atol=1e-7);checks+=1
    return checks

def save_snapshot(tracker,path,rho):
    s=tracker.audit_snapshot;keys=list(s['ends']['multiple']);data=dict(multiple_keys=np.array(keys),multiple_logs=np.array([s['ends']['multiple'][k] for k in keys]),old_kind=np.array([CLASSES.index(p['kind']) for p in s['old']]),old_center=np.array([p['center'] for p in s['old']]),old_si=np.array([p['si'] for p in s['old']]),old_mean=np.array([p['mean'] for p in s['old']]),old_log=np.array([p['log'] for p in s['old']]),innovation=s['innovation'],variance=s['variance'],rho=rho,precision=tracker.P,FT=tracker.ft,norms=tracker.norms,weights=tracker.freqweights,classlogs=s['classlogs'])
    for kind in CLASSES[:-1]:data[kind+'_keys']=np.array(list(s['ends'][kind]));data[kind+'_logs']=np.array(list(s['ends'][kind].values()))
    np.savez_compressed(path,**data)
def run():
    p=read(OUT/'plan.json');assert p['source']==binding() and not(OUT/'attempt.json').exists();RAW.mkdir(exist_ok=False);start=time.perf_counter();atomic_json(OUT/'attempt.json',dict(status='RUNNING'));checks=0;rows=[];preflight_seconds=[]
    try:
        for name,h in p['inputs'].items():assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==h
        with np.load(OLD/'prepared_inputs.npz') as z:ref=z['reference'];initial=z['initial'];omega=float(z['omega'])
        data=[]
        for r in(0,1):
            with np.load(OLD/'production'/f'homogeneous_r{r}'/'observations.npz') as z:data.append(z['counts'][:,0,:3]/omega)
        cov=read(BENCH)['regime_inputs'];rng=np.random.default_rng(p['seed'])
        for ri in(0,2):
            when=[.1,.5,1.][ri];model=XYCore(ref,initial,omega,np.array(cov[ri]['covariance'])[1:,1:]);samples=np.array([x[round(when*1000)] for x in data]);samples-=samples.mean(axis=2)[:,:,None];ft=np.fft.rfft(samples,axis=2);S=np.einsum('rsf,rtf->fst',ft,np.conj(ft))/(2*500);S=np.stack([S[max(0,i-2):min(251,i+3)].mean(axis=0) for i in range(251)]);diag=np.zeros_like(S)
            for sp in range(3):diag[:,sp,sp]=S[:,sp,sp].real
            S=.9*S+.1*diag;ev,V=np.linalg.eigh(S);ev=np.maximum(ev,1e-8*ev.max(axis=1)[:,None]);L=V*np.sqrt(ev)[:,None,:];L[0]=0
            # No-draw fixtures precede validation in each regime.
            pre=time.perf_counter();test=Corrected(model,L,.95,100);test.fixture_first=True;test.step(model.signal(1,100));checks+=fixture_check(test);test.fixture_first=False;test.step(model.signal(1,101));checks+=fixture_check(test);preflight_seconds.append(time.perf_counter()-pre)
            def noise():return np.fft.irfft(np.einsum('fij,jf->if',L,np.fft.rfft(rng.normal(size=(3,500)),axis=1))*model.x.smoother,n=500,axis=1)[1:]
            for rho in(0.,.95):
                for case in p['cases']:
                    paths=[];fields=[]
                    for r in range(4):
                        seed=(100+71*r)%500;trackers=dict(corrected=Corrected(model,L,rho,seed),original=Beam(model,L,rho,seed));records={k:[] for k in trackers};n=noise();sequence=[];truths=[]
                        for t in range(24):
                            if t:n=rho*n+np.sqrt(1-rho*rho)*noise()
                            center,scale=truth(case,t,seed);signal=model.signal(scale,center)
                            if case=='absent' or(case=='vanish' and t>=10):signal[:]=0
                            if case=='missingX' and t>=10:signal[0]=0
                            if case=='missingY' and t>=10:signal[1]=0
                            if case=='displacedY' and t>=10:signal[1]=np.roll(signal[1],16)
                            if case=='split' and t>=10:signal+=model.signal(scale,(center+250)%500)
                            field=signal+n;sequence.append(field);truths.append([center,scale])
                            for name,tracker in trackers.items():records[name].append(tracker.step(field))
                            corrected=trackers['corrected'];corrected.fixture_first=t==0
                            if r==0 and t in(0,10):checks+=fixture_check(corrected)
                            if r==0 and case in('static','split') and t in(0,10,15):save_snapshot(corrected,RAW/f'r{ri}_rho{rho:g}_{case}_snap{t}.npz',0 if t==0 else rho)
                            if time.perf_counter()-start>600:raise TimeoutError('MG172 active600s cap')
                        outcome={};positive=case in('static','moving','deforming');onset=0 if case=='absent' else 10
                        for name,rec in records.items():
                            singles=sum(v['label']=='single' for v in rec);accurate=sum(v['label']=='single' and min((v['center']-tr[0])%500,(tr[0]-v['center'])%500)<=2 and abs(v['width']-tr[1])<=.12500001 for v,tr in zip(rec,truths));f=rec[onset+5]['label']=='single' if not positive else None;success=singles>=22 and accurate>=22 if positive else not f;outcome[name]=dict(single_frames=singles,accurate_frames=accurate,false_single_deadline=f,success=bool(success))
                        paths.append(dict(replicate=r,seed=seed,truth=truths,records=records,outcomes=outcome));fields.append(sequence)
                    tag=f'r{ri}_rho{rho:g}_{case}';np.savez_compressed(RAW/(tag+'.npz'),fields=fields);atomic_json(OUT/(tag+'.json'),dict(regime=when,rho=rho,case=case,paths=paths));row=dict(regime=when,rho=rho,case=case,n=4)
                    for name in trackers:
                        successes=sum(v['outcomes'][name]['success'] for v in paths);row[name]=dict(successful_paths=successes,false_single_deadline=sum(bool(v['outcomes'][name]['false_single_deadline']) for v in paths) if not positive else None,screen_passed=successes==4)
                    rows.append(row)
                    if sum(p.stat().st_size for p in RAW.iterdir())>100000000:raise RuntimeError('MG172 local0.1GBreserve')
        atomic_json(OUT/'analysis.json',dict(rows=rows,corrected_screen_passed=all(r['corrected']['screen_passed'] for r in rows),original_screen_passed=all(r['original']['screen_passed'] for r in rows),explicit_endpoint_history_sum_checks=checks,preflight_seconds=preflight_seconds,elapsed_seconds=time.perf_counter()-start,G_status='UNRESOLVED',new_persistence_verdict=None));atomic_json(OUT/'attempt.json',dict(status='COMPLETE',elapsed_seconds=time.perf_counter()-start))
    except BaseException:atomic_json(OUT/'attempt.json',dict(status='FAILED',elapsed_seconds=time.perf_counter()-start,error=traceback.format_exc()));raise
if __name__=='__main__':
    if __import__('sys').argv[1]=='register':register()
    else:run()

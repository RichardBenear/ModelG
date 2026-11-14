"""Independent full-suite validation of stationary/moving-mode tracker."""
import os
for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):os.environ[k]='1'
from pathlib import Path
import time,math,hashlib,traceback
import numpy as np
from scipy.special import logsumexp
from mg174_modes import ModeTracker,PAIRS,MOTIONS
from mg172_corrected import Corrected
from run_mg172 import ROOT,OLD,BENCH,read,atomic_json
from mg167_xy import XYCore
from run_mg169_tracking import truth
from analyze_mg170 import CLASSES
OUT=ROOT/'results/MG12H3/stationary_moving_tracker';RAW=ROOT/'local_raw/MG174_modes'
POSITIVES=['static','moving','deforming','stop_start'];NEGATIVES=['absent','vanish','missingX','missingY','displacedY','split','jump'];SENS=['weak_static','weak_moving','noise_scaled']
FILES=['src/mg174_modes.py','src/run_mg174.py','src/mg172_corrected.py','src/analyze_mg170.py','src/mg167_xy.py','src/mg161_local_shape.py','src/run_mg169_tracking.py','spec/MG174_stationary_moving_tracker.md']
def binding():return {p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in FILES}
def register():
    OUT.mkdir(exist_ok=True);assert not(OUT/'plan.json').exists();paths=[OLD/'prepared_inputs.npz',BENCH]+[OLD/'production'/f'homogeneous_r{r}'/'observations.npz' for r in(0,1)];atomic_json(OUT/'plan.json',dict(source=binding(),inputs={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},seed=17400,cap=600,workers=1,threads=1,primary_cases=POSITIVES+NEGATIVES,sensitivity_cases=SENS,primary_replicates=4,primary_frames=24,sensitivity_replicates=2,sensitivity_frames=12,mode_stay=.98,no_ModelG_draws=True))
def geometry(case,t,seed):
    if case=='stop_start':return (seed+min(max(t-7,0),8))%500,1.
    if case in('weak_static','noise_scaled'):case='static'
    if case=='weak_moving':case='moving'
    return truth(case,t,seed)
def verify(tracker,first):
    s=tracker.audit_snapshot;rho=0 if first else tracker.rho;checks=0
    for motion in MOTIONS:
        pair=('multiple',motion);key=next(iter(s['ends'][pair]));k,si=key;N=5655
        for secondary,ss in [((k+100)%500,0),((k+250)%500,6),((k+400)%500,12)]:
            mu=tracker.mean('multiple',k,si,(secondary,ss));vals=[]
            for parent in s['old']:
                allowed=(k==parent['center'] and si==parent['si']) if motion=='S' else min((k-parent['center'])%500,(parent['center']-k)%500)<=2 and abs(si-parent['si'])<=1
                if not allowed:continue
                err=s['innovation']-(mu-rho*parent['mean']);q=float(tracker.inner(err,err)/s['variance']);nloc=1 if motion=='S' else 5*(min(13,parent['si']+2)-max(0,parent['si']-1));vals.append(parent['log']+np.log(.98 if parent['kind']=='multiple' else .005)+np.log(.98 if parent['motion']==motion else .02)-np.log(nloc)-np.log(N)-.5*q)
            np.testing.assert_allclose(float(logsumexp(vals)),s['ends'][pair][key][ss,secondary],rtol=1e-9,atol=2e-6);checks+=1
    return checks

def snapshot(tracker,path,rho):
    s=tracker.audit_snapshot;data=dict(old_kind=np.array([CLASSES.index(x['kind']) for x in s['old']]),old_motion=np.array([MOTIONS.index(x['motion']) for x in s['old']]),old_center=np.array([x['center'] for x in s['old']]),old_si=np.array([x['si'] for x in s['old']]),old_mean=np.array([x['mean'] for x in s['old']]),old_log=np.array([x['log'] for x in s['old']]),innovation=s['innovation'],variance=s['variance'],rho=rho,precision=tracker.P,FT=tracker.ft,norms=tracker.norms,weights=tracker.freqweights,pairlogs=s['pairlogs'])
    for kind,motion in PAIRS:
        key=kind+'_'+motion;data[key+'_keys']=np.array(list(s['ends'][(kind,motion)]));data[key+'_logs']=np.array(list(s['ends'][(kind,motion)].values()))
    np.savez_compressed(path,**data)
def run():
    p=read(OUT/'plan.json');assert p['source']==binding() and not(OUT/'attempt.json').exists();RAW.mkdir(exist_ok=False);start=time.perf_counter();atomic_json(OUT/'attempt.json',dict(status='RUNNING'));rows=[];checks=0
    try:
        for name,h in p['inputs'].items():assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==h
        with np.load(OLD/'prepared_inputs.npz') as z:reference=z['reference'];initial=z['initial'];omega=float(z['omega'])
        inputs=[]
        for r in(0,1):
            with np.load(OLD/'production'/f'homogeneous_r{r}'/'observations.npz') as z:inputs.append(z['counts'][:,0,:3]/omega)
        cov=read(BENCH)['regime_inputs'];rng=np.random.default_rng(p['seed'])
        assert abs(.98+.02-1)<1e-12 and abs(.98+4*.005-1)<1e-12
        for ri in(0,2):
            when=[.1,.5,1.][ri];model=XYCore(reference,initial,omega,np.array(cov[ri]['covariance'])[1:,1:]);samples=np.array([x[round(when*1000)] for x in inputs]);samples-=samples.mean(axis=2)[:,:,None];ft=np.fft.rfft(samples,axis=2);S=np.einsum('rsf,rtf->fst',ft,np.conj(ft))/(2*500);S=np.stack([S[max(0,i-2):min(251,i+3)].mean(axis=0) for i in range(251)]);diag=np.zeros_like(S)
            for sp in range(3):diag[:,sp,sp]=S[:,sp,sp].real
            S=.9*S+.1*diag;ev,V=np.linalg.eigh(S);ev=np.maximum(ev,1e-8*ev.max(axis=1)[:,None]);L=V*np.sqrt(ev)[:,None,:];L[0]=0
            test=ModeTracker(model,L,.95,100);test.step(model.signal(1,100));checks+=verify(test,True);test.step(model.signal(1,101));checks+=verify(test,False)
            def noise():return np.fft.irfft(np.einsum('fij,jf->if',L,np.fft.rfft(rng.normal(size=(3,500)),axis=1))*model.x.smoother,n=500,axis=1)[1:]
            for rho in(0.,.95):
                for case in p['primary_cases']+p['sensitivity_cases']:
                    primary=case not in SENS;N=p['primary_replicates'] if primary else p['sensitivity_replicates'];T=p['primary_frames'] if primary else p['sensitivity_frames'];positive=case in POSITIVES or not primary;onset=0 if case=='absent' else 10;paths=[];fields=[]
                    for rep in range(N):
                        seed=(100+71*rep)%500;trackers=dict(modes=ModeTracker(model,L,rho,seed),MG172=Corrected(model,L,rho,seed));records={k:[] for k in trackers};n=noise();seq=[];truths=[]
                        for t in range(T):
                            if t:n=rho*n+np.sqrt(1-rho*rho)*noise()
                            k,s=geometry(case,t,seed);signal=model.signal(s,k)*(.8 if case.startswith('weak_') else 1.)
                            if case=='absent' or(case=='vanish' and t>=10):signal[:]=0
                            if case=='missingX' and t>=10:signal[0]=0
                            if case=='missingY' and t>=10:signal[1]=0
                            if case=='displacedY' and t>=10:signal[1]=np.roll(signal[1],16)
                            if case=='split' and t>=10:signal+=model.signal(s,(k+250)%500)
                            field=signal+(1.25 if case=='noise_scaled' else 1.)*n;seq.append(field);truths.append([k,s])
                            for name,tracker in trackers.items():records[name].append(tracker.step(field))
                            if rep==0 and t in(0,10):checks+=verify(trackers['modes'],t==0)
                            if primary and rep==0 and case in('static','moving','split') and t==10:snapshot(trackers['modes'],RAW/f'r{ri}_rho{rho:g}_{case}_snap{t}.npz',rho)
                            if time.perf_counter()-start>600:raise TimeoutError('MG174 active600s cap')
                        outcomes={}
                        for name,rec in records.items():
                            single=sum(x['label']=='single' for x in rec);acc=sum(x['label']=='single' and min((x['center']-tr[0])%500,(tr[0]-x['center'])%500)<=2 and abs(x['width']-tr[1])<=.12500001 for x,tr in zip(rec,truths));false=rec[onset+5]['label']=='single' if not positive else None;success=single>=math.ceil(.9*T) and acc>=math.ceil(.9*T) if positive else not false;outcomes[name]=dict(single_frames=single,accurate_frames=acc,false_single_deadline=false,success=bool(success))
                        paths.append(dict(replicate=rep,seed=seed,truth=truths,records=records,outcomes=outcomes));fields.append(seq)
                    tag=f'r{ri}_rho{rho:g}_{case}';np.savez_compressed(RAW/(tag+'.npz'),fields=fields);atomic_json(OUT/(tag+'.json'),dict(regime=when,rho=rho,case=case,primary=primary,frames=T,n=N,paths=paths));row=dict(regime=when,rho=rho,case=case,primary=primary,frames=T,n=N)
                    for name in trackers:
                        successes=sum(x['outcomes'][name]['success'] for x in paths);row[name]=dict(successful_paths=successes,false_single_deadline=sum(bool(x['outcomes'][name]['false_single_deadline']) for x in paths) if not positive else None,screen_passed=successes==N)
                    rows.append(row)
                    if sum(p.stat().st_size for p in RAW.iterdir())>200000000:raise RuntimeError('MG174 local0.2GBreserve')
        atomic_json(OUT/'analysis.json',dict(rows=rows,primary_modes_screen_passed=all(r['modes']['screen_passed'] for r in rows if r['primary']),primary_MG172_screen_passed=all(r['MG172']['screen_passed'] for r in rows if r['primary']),sensitivity_not_qualification=True,explicit_history_endpoint_checks=checks,G_status='UNRESOLVED',ready_for_real_persistence=False,new_persistence_verdict=None,elapsed_seconds=time.perf_counter()-start));atomic_json(OUT/'attempt.json',dict(status='COMPLETE',elapsed_seconds=time.perf_counter()-start))
    except BaseException:atomic_json(OUT/'attempt.json',dict(status='FAILED',elapsed_seconds=time.perf_counter()-start,error=traceback.format_exc()));raise
if __name__=='__main__':
    if __import__('sys').argv[1]=='register':register()
    else:run()

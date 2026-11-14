"""Bounded conditional tracking prototype and complete finite validation denominator."""
import os
for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):os.environ[k]='1'
from pathlib import Path
import time,json,hashlib,traceback
import numpy as np
from mg167_xy import XYCore
from analyze_mg168 import oracle,guard_flags
from run_mg163_long import ROOT,BASE as OLD,read,atomic_json
OUT=ROOT/'results/MG12H3/tracking_feasibility';RAW=ROOT/'local_raw/MG169_tracking';PRIOR=ROOT/'results/MG12H3/XY_core_instrument'
FILES=['src/run_mg169_tracking.py','src/mg167_xy.py','src/analyze_mg168.py','spec/MG169_tracking_feasibility.md']
def binding():return {p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in FILES}
def register():
    OUT.mkdir(exist_ok=True);assert not(OUT/'plan.json').exists();inputs=[OLD/'prepared_inputs.npz',PRIOR/'benchmark.json']+[OLD/'production'/f'homogeneous_r{r}'/'observations.npz' for r in(0,1)];atomic_json(OUT/'plan.json',dict(source=binding(),inputs={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs},seed=16900,cap=600,workers=1,threads=1,frames=40,replicates=16,regimes=[.1,.5,1.],rho=[0.,.95],cases=['static','moving','deforming','absent','vanish','missingX','missingY','displacedY','split','jump'],no_ModelG_draws=True))
def truth(case,t,seed):
    moving=case!='static';center=(seed+(t if moving else 0))%500;scale=1.
    if case=='deforming':center=(seed+t//2)%500;scale=[1.,1.125,1.25,1.125][t//10]
    if case=='jump' and t>=10:center=(center+20)%500
    return center,scale

def localfit(model,field,center,scale):
    m=model.x;f=field-field.mean(axis=1)[:,None];dot=np.fft.irfft((np.fft.rfft(f,axis=1)[None,:,:]*model.kfft).sum(axis=1),n=500,axis=1);gain=np.maximum(dot,0)**2/model.norm[:,None];distance=np.minimum((np.arange(500)-center)%500,(center-np.arange(500))%500);allowed=(abs(m.scales-scale)<=.12500001)&(m.scales>=.5)&(m.scales<=2);masked=np.where(allowed[:,None]&(distance[None,:]<=2),gain,-np.inf);si,k=np.unravel_index(masked.argmax(),masked.shape);return oracle(model,field,int(k),float(m.scales[si]))

def preflight():
    for seed in(0,499):
        for case in('static','moving','deforming'):
            prev=truth(case,0,seed)
            for t in range(1,40):
                now=truth(case,t,seed);assert min((now[0]-prev[0])%500,(prev[0]-now[0])%500)<=2 and abs(now[1]-prev[1])<=.125;prev=now
    assert abs(np.sqrt(1-.95**2)**2+.95**2-1)<1e-12
    # three failures stop; accepted observations reset the gap count
    gap=0;alive=True
    for accepted in(False,False,True,False,False,False):
        if accepted:gap=0
        else:gap+=1
        if gap>=3:alive=False
    assert not alive and gap==3

def run():
    p=read(OUT/'plan.json');assert p['source']==binding() and not(OUT/'attempt.json').exists();preflight();RAW.mkdir(exist_ok=False);start=time.perf_counter();atomic_json(OUT/'attempt.json',dict(status='RUNNING'));checks=0
    try:
        for name,h in p['inputs'].items():assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==h
        with np.load(OLD/'prepared_inputs.npz') as z:reference=z['reference'];initial=z['initial'];omega=float(z['omega'])
        saved=[]
        for r in(0,1):
            with np.load(OLD/'production'/f'homogeneous_r{r}'/'observations.npz') as z:saved.append(z['counts'][:,0,:3]/omega)
        old=read(PRIOR/'benchmark.json');rng=np.random.default_rng(p['seed']);summaries=[];calibrations=[]
        for ri,when in enumerate(p['regimes']):
            model=XYCore(reference,initial,omega,np.array(old['regime_inputs'][ri]['covariance'])[1:,1:]);m=model.x;samples=np.array([x[round(when*1000)] for x in saved]);samples-=samples.mean(axis=2)[:,:,None];ft=np.fft.rfft(samples,axis=2);spectra=np.einsum('rsf,rtf->fst',ft,np.conj(ft))/(2*500);spectra=np.stack([spectra[max(0,i-2):min(251,i+3)].mean(axis=0) for i in range(251)]);diag=np.zeros_like(spectra)
            for sp in range(3):diag[:,sp,sp]=spectra[:,sp,sp].real
            spectra=.9*spectra+.1*diag;ev,V=np.linalg.eigh(spectra);ev=np.maximum(ev,1e-8*ev.max(axis=1)[:,None]);L=V*np.sqrt(ev)[:,None,:];L[0]=0
            def noise():return np.fft.irfft(np.einsum('fij,jf->if',L,np.fft.rfft(rng.normal(size=(3,500)),axis=1))*m.smoother,n=500,axis=1)[1:]
            def deadline():
                if time.perf_counter()-start>p['cap']:raise TimeoutError('MG169 active600s cap')
            null=[];good=[];calfields=[]
            for i in range(256):
                c=(37*i)%500;s=.5+.125*(i%13);f=noise();null.append(localfit(model,f,c,s));calfields.append(f)
                s=[1.,1.125,1.25][i%3];f=noise()+model.signal(s,c);good.append(localfit(model,f,(c+(i%3)-1)%500,s));calfields.append(f);deadline()
            cal=dict(score=max(f['score'] for f in null),residual=max(f['residual'] for f in good),component=max(f['component_stat'] for f in good),secondary=max(f['secondary_score'] for f in good));atomic_json(OUT/f'r{ri}_calibration.json',dict(thresholds=cal,null_fits=null,shape_fits=good));np.savez_compressed(RAW/f'r{ri}_calibration.npz',fields=calfields);calibrations.append(cal)
            for rho in p['rho']:
                for case in p['cases']:
                    paths=[];fields=[]
                    for replicate in range(16):
                        seed=(100+71*replicate)%500;states={name:dict(center=seed,scale=1.,gap=0,alive=True,stop=None) for name in('local','free')};records={name:[] for name in states};truths=[];sequence=[];n=noise()
                        for t in range(40):
                            if t:n=rho*n+np.sqrt(1-rho*rho)*noise()
                            center,scale=truth(case,t,seed);signal=model.signal(scale,center)
                            if case=='absent' or(case=='vanish' and t>=10):signal[:]=0
                            if case=='missingX' and t>=10:signal[0]=0
                            if case=='missingY' and t>=10:signal[1]=0
                            if case=='displacedY' and t>=10:signal[1]=np.roll(signal[1],16)
                            if case=='split' and t>=10:signal+=model.signal(scale,(center+250)%500)
                            field=signal+n;sequence.append(field);truths.append([center,scale])
                            for name,state in states.items():
                                if not state['alive']:records[name].append(dict(active=False,accepted=False,center=None,scale=None,fit=None));continue
                                if name=='local':f=localfit(model,field,state['center'],state['scale']);flags=guard_flags(model,f,cal);accepted=not any(flags.values())
                                else:f=model.fit(field,already_smooth=True);accepted=model.label(f,old['calibrations'][ri])=='candidate_compatible';flags=None
                                if accepted:state.update(center=f.get('center',f.get('center_site')),scale=f['scale'],gap=0)
                                else:state['gap']+=1
                                if state['gap']>=3:state['alive']=False;state['stop']=t
                                records[name].append(dict(active=True,accepted=accepted,center=state['center'],scale=state['scale'],fit=f,flags=flags))
                                if t in(0,10,39) and replicate==0:
                                    w=np.roll(m.window,f['center']);template=model.signal(f['scale'],f['center']);A=np.zeros((2,500,3));A[:,:,0]=template
                                    for sp in range(2):A[sp,:,sp+1]=1
                                    normal=np.einsum('ina,ij,jnb,n->ab',A,model.inv,A,w);rhs=np.einsum('ina,ij,jn,n->a',A,model.inv,field,w);coef=np.linalg.solve(normal,rhs);np.testing.assert_allclose(coef[0],f['amplitude'],rtol=1e-8,atol=1e-8);checks+=1
                            deadline()
                        outcomes={}
                        for name,state in states.items():
                            rec=records[name];compatible=sum(r['accepted'] for r in rec);accurate=sum(r['center'] is not None and min((r['center']-tr[0])%500,(tr[0]-r['center'])%500)<=2 and abs(r['scale']-tr[1])<=.12500001 for r,tr in zip(rec,truths));positive=case in('static','moving','deforming');onset=0 if case=='absent' else 10;continued=state['stop'] is None or state['stop']>onset+5;success=state['alive'] and compatible>=36 and accurate>=36 if positive else not continued;outcomes[name]=dict(compatible_frames=compatible,accurate_frames=accurate,stop=state['stop'],alive_final=state['alive'],success=bool(success),false_continuation=bool(continued) if not positive else None)
                        paths.append(dict(replicate=replicate,initial_seed=seed,truth=truths,records=records,outcomes=outcomes));fields.append(sequence)
                    tag=f'r{ri}_rho{rho:g}_{case}';np.savez_compressed(RAW/(tag+'.npz'),fields=fields);atomic_json(OUT/(tag+'.json'),dict(regime=when,rho=rho,case=case,paths=paths));positive=case in('static','moving','deforming');summary=dict(regime=when,rho=rho,case=case,n=16)
                    for name in('local','free'):
                        successes=sum(x['outcomes'][name]['success'] for x in paths);false=sum(bool(x['outcomes'][name]['false_continuation']) for x in paths) if not positive else None;summary[name]=dict(successful_paths=successes,false_continuations=false,screen_passed=successes/16>=.9 if positive else false/16<=.05)
                    summaries.append(summary)
                    if sum(p.stat().st_size for p in RAW.iterdir() if p.is_file())>500000000:raise RuntimeError('MG169 local0.5GB reserve exceeded')
        atomic_json(OUT/'analysis.json',dict(rows=summaries,calibrations=calibrations,independent_primary_checks=checks,local_feasibility_screen=all(r['local']['screen_passed'] for r in summaries),free_screen=all(r['free']['screen_passed'] for r in summaries),G_status='UNRESOLVED',new_persistence_verdict=None,elapsed_seconds=time.perf_counter()-start));atomic_json(OUT/'attempt.json',dict(status='COMPLETE',elapsed_seconds=time.perf_counter()-start))
    except BaseException:atomic_json(OUT/'attempt.json',dict(status='FAILED',elapsed_seconds=time.perf_counter()-start,error=traceback.format_exc()));raise
if __name__=='__main__':
    if __import__('sys').argv[1]=='register':register()
    else:run()

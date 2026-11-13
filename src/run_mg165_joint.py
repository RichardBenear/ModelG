"""Frozen conditional benchmark and exploratory saved-data application."""
import os
for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):os.environ[k]='1'
from pathlib import Path
import numpy as np,json,time,hashlib,traceback
from scipy.stats import beta
from mg165_joint import Joint
from run_mg163_long import ROOT,BASE as OLD,PUBLIC as PRIOR,read,atomic_json
OUT=ROOT/'results/MG12H3/multivariate_instrument';RAW=ROOT/'local_raw/MG165_multivariate'
FILES=['src/mg165_joint.py','src/run_mg165_joint.py','spec/MG165_multivariate_instrument.md']
def binding():return {p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in FILES}
def register():
    OUT.mkdir(exist_ok=True);assert not(OUT/'plan.json').exists();inputs={}
    for p in [OLD/'prepared_inputs.npz']+[OLD/'production'/j['id']/'observations.npz' for j in read(PRIOR/'plan.json')['jobs']]:inputs[str(p.relative_to(ROOT))]=hashlib.sha256(p.read_bytes()).hexdigest()
    atomic_json(OUT/'plan.json',dict(source=binding(),inputs=inputs,seed=16500,cap=600,regimes=[.1,.5,1.],no_ModelG_draws=True,scope='Conditional synthetic development and already-seen saved trajectories; no rescore'))
def run():
    p=read(OUT/'plan.json');assert p['source']==binding() and not(OUT/'attempt.json').exists();RAW.mkdir(exist_ok=False);start=time.perf_counter();atomic_json(OUT/'attempt.json',dict(status='RUNNING'));checks=0
    try:
        with np.load(OLD/'prepared_inputs.npz') as z:ref=z['reference'];initial=z['initial'];omega=float(z['omega'])
        saved={}
        for job in read(PRIOR/'plan.json')['jobs']:
            path=OLD/'production'/job['id']/'observations.npz';assert hashlib.sha256(path.read_bytes()).hexdigest()==p['inputs'][str(path.relative_to(ROOT))]
            with np.load(path) as z:saved[job['id']]=z['counts'][:,:,:3]/omega
        rng=np.random.default_rng(p['seed']);models=[];calibrations=[];summary=[];regime_inputs=[]
        def deadline():
            if time.perf_counter()-start>p['cap']:raise TimeoutError('MG165 active600s cap')
        def verify(model,field,f):
            nonlocal checks
            m=model.x;y=field-field.mean(axis=1)[:,None];g=np.roll(model.bank[int(np.flatnonzero(m.scales==f['scale'])[0])],f['center'],axis=1);w=np.roll(m.window,f['center']);A=np.zeros((3,500,4));A[:,:,0]=g
            for i in range(3):A[i,:,i+1]=1
            normal=np.einsum('ina,ij,jnb,n->ab',A,model.inv,A,w);rhs=np.einsum('ina,ij,jn,n->a',A,model.inv,y,w);coef=np.linalg.solve(normal,rhs);err=y-np.einsum('ina,a->in',A,coef);res=np.einsum('in,ij,jn,n->',err,model.inv,err,w)
            np.testing.assert_allclose([coef[0],res],[f['amplitude'],f['residual']],rtol=1e-8,atol=1e-8);np.testing.assert_allclose((err*err)@w,f['component_residuals'],rtol=1e-8,atol=1e-8);checks+=1
        for regime,t in enumerate(p['regimes']):
            samples=np.array([saved[f'homogeneous_r{r}'][round(t*1000),0] for r in(0,1)]);samples-=samples.mean(axis=2)[:,:,None];ft=np.fft.rfft(samples,axis=2);spectra=np.einsum('rsf,rtf->fst',ft,np.conj(ft))/(2*500)
            spectra=np.stack([np.mean(spectra[max(0,i-2):min(251,i+3)],axis=0) for i in range(251)]);diag=np.zeros_like(spectra)
            for sp in range(3):diag[:,sp,sp]=spectra[:,sp,sp].real
            spectra=.9*spectra+.1*diag;eig,V=np.linalg.eigh(spectra);eig=np.maximum(eig,1e-8*eig.max(axis=1)[:,None]);L=V*np.sqrt(eig)[:,None,:];L[0]=0
            from mg161_local_shape import LocalShape
            smooth=LocalShape(ref[1],initial[1],omega).smoother;sf=np.fft.irfft(ft*smooth,n=500,axis=2);cov=np.einsum('rsn,rtn->st',sf,sf)/(2*500);cov=.95*cov+.05*np.diag(np.diag(cov));model=Joint(ref,initial,omega,cov);models.append(model);regime_inputs.append(dict(time=t,covariance=cov.tolist()))
            def sample(kind,index):
                white=np.fft.rfft(rng.normal(size=(3,500)),axis=1);noise=np.fft.irfft(np.einsum('fij,jf->if',L,white)*smooth,n=500,axis=1);amp=[.8,1.,1.2][index%3];scale=[.85,1.1,1.6][(index//3)%3];shift=(37*index+71*regime)%500
                if kind.startswith('intact_'):amp,scale=[(.8,.85),(1.,1.1),(1.2,1.6)][int(kind[-1])]
                signal=amp*model.signal(scale,shift)
                if kind=='null':signal[:]=0
                if kind.startswith('missing'):signal[int(kind[-1])]=0
                if kind.startswith('shift'):signal[int(kind[-1])]=np.roll(signal[int(kind[-1])],16)
                if kind=='split':signal+=amp*model.signal(scale,(shift+250)%500)
                return signal+noise
            def batch(kind,n,cal=None):
                fields=[];fits=[];labels=[]
                for i in range(n):
                    deadline();field=sample(kind,i);f=model.fit(field,already_smooth=True);fields.append(field);fits.append(f)
                    if i in (0,n-1):verify(model,field,f)
                    if cal is not None:labels.append(model.label(f,cal))
                np.savez_compressed(RAW/f'r{regime}_{kind}.npz',fields=fields)
                atomic_json(OUT/f'r{regime}_{kind}_fits.json',dict(fits=fits,labels=labels));return fits,labels
            null,_=batch('cal_null',256);good,_=batch('cal_shape',256);cal=dict(score=max(f['score'] for f in null),residual=max(f['residual'] for f in good),component=max(f['component_stat'] for f in good));calibrations.append(cal);atomic_json(OUT/f'r{regime}_calibration.json',cal)
            for kind,n in [('null',1024)]+[(f'intact_{i}',128) for i in range(3)]+[(f'{k}{i}',128) for k in('missing','shift') for i in range(3)]+[('split',128)]:
                fits,labels=batch(kind,n,cal);count=labels.count('candidate_compatible');positive=kind.startswith('intact');limit=.9 if positive else(.01 if kind=='null' else .05);passed=count/n>=limit if positive else count/n<=limit;lo=0 if count==0 else float(beta.ppf(.025,count,n-count+1));hi=1 if count==n else float(beta.ppf(.975,count+1,n-count));summary.append(dict(regime=t,case=kind,n=n,compatible=count,fraction=count/n,screen_passed=passed,descriptive_individual_CP95=[lo,hi]))
        qualified=all(r['screen_passed'] for r in summary);atomic_json(OUT/'benchmark.json',dict(conditional_development_screen_passed=qualified,cases=summary,regime_inputs=regime_inputs,calibrations=calibrations,limits='Raw finite benchmark screen, individual intervals not familywise qualification; estimated conditional Gaussian noise only. No real-trajectory held-out validation.'))
        # All saved paths are development/exploratory, regardless of benchmark failure.
        application=[]
        for job in read(PRIOR/'plan.json')['jobs']:
            features=[];labels=[]
            for side in range(2):
                fs=[];ls=[]
                for i,field in enumerate(saved[job['id']][:,side]):
                    deadline();ri=int(np.argmin(abs(np.array(p['regimes'])-i*.001)));f=models[ri].fit(field);fs.append(f);ls.append(models[ri].label(f,calibrations[ri]))
                    if i in(0,205,500,1000):
                        m=models[ri].x;smoothfield=np.fft.irfft(np.fft.rfft(field-field.mean(axis=1)[:,None],axis=1)*m.smoother,n=500,axis=1);verify(models[ri],smoothfield,f)
                features.append(fs);labels.append(ls)
            atomic_json(OUT/(job['id']+'_candidate.json'),dict(id=job['id'],fits=features,labels=labels,exploratory=True,qualified_persistence_verdict=False))
            application.append(dict(id=job['id'],condition=job['condition'],compatible_frames=[l.count('candidate_compatible') for l in labels],frames=1001))
        atomic_json(OUT/'analysis.json',dict(application=application,independent_GLS_checks=checks,conditional_screen_passed=qualified,new_persistence_verdict=None,elapsed_seconds=time.perf_counter()-start));atomic_json(OUT/'attempt.json',dict(status='COMPLETE',elapsed_seconds=time.perf_counter()-start))
    except BaseException:atomic_json(OUT/'attempt.json',dict(status='FAILED',elapsed_seconds=time.perf_counter()-start,error=traceback.format_exc()));raise
if __name__=='__main__':
    if __import__('sys').argv[1]=='register':register()
    else:run()

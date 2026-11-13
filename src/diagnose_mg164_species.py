"""Saved-only species diagnostics; no new model draws or joint pass/fail score."""
import os
for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):os.environ[k]='1'
from pathlib import Path
import time,json,gzip,hashlib,traceback
import numpy as np
from mg161_local_shape import LocalShape
from run_mg163_long import ROOT,BASE as RAW,PUBLIC,read,atomic_json
OUT=ROOT/'results/MG12H3/species_profile_diagnostic'

def run():
    start=time.perf_counter();OUT.mkdir(exist_ok=True);attempt=OUT/'attempt.json';assert not attempt.exists();atomic_json(attempt,dict(status='RUNNING'))
    try:
        plan=read(PUBLIC/'plan.json');scores=json.loads(gzip.decompress((PUBLIC/'candidate_scores.json.gz').read_bytes()));inputs={};checks=0;rows=[];profiles={};arrays={}
        with np.load(RAW/'prepared_inputs.npz') as z:reference=z['reference'][:3];initial=z['initial'][:3];omega=float(z['omega']);m=LocalShape(z['reference'][1],z['initial'][1],omega)
        far=abs(m.coords)>=m.n*m.h/4;raw=np.roll(reference,-m.anchor,axis=1);raw-=np.median(raw[:,far],axis=1)[:,None]
        def bank(scale):
            return np.array([np.fft.irfft(np.fft.rfft(np.interp(m.coords/scale,m.coords[m.order],r[m.order],left=0,right=0))*m.smoother,n=m.n) for r in raw])
        templates={float(s):bank(float(s)) for s in m.scales}
        def fit(field,f):
            smooth=np.fft.irfft(np.fft.rfft(field-field.mean(axis=1)[:,None],axis=1)*m.smoother,n=m.n,axis=1);k=f['center_site'];w=np.roll(m.window,k);g=np.roll(templates[f['scale']],k,axis=1);fm=smooth@w/m.ws;gm=g@w/m.ws;fc=smooth-fm[:,None];gc=g-gm[:,None];norm=(gc*gc)@w;amplitude=(fc*gc)@w/norm;residual=fc-amplitude[:,None]*gc;cov=(residual*w)@residual.T/m.ws;energy=(fc*fc)@w/m.ws
            distance=np.minimum((np.arange(m.n)-k)%m.n,(k-np.arange(m.n))%m.n)*m.h;outside=distance>m.background_exclusion;bg=smooth[:,outside];bg=bg-bg.mean(axis=1)[:,None];bcov=bg@bg.T/bg.shape[1]
            return amplitude,cov,energy,bcov,smooth,g,w
        initialfit=m.fit(initial[1]/omega);baseamp=fit(initial/omega,initialfit)[0];assert np.all(abs(baseamp)>1e-10)
        for job in scores['jobs']:
            path=RAW/'production'/job['id']/'observations.npz';inputs[str(path.relative_to(ROOT))]=hashlib.sha256(path.read_bytes()).hexdigest()
            with np.load(path) as z:counts=z['counts'];times=z['times']
            hist=np.zeros((2,1001,3,4));covhist=np.zeros((2,1001,2,3,3))
            for side,outcome in enumerate(job['outcomes']):
                for t,f in enumerate(outcome['fits']):
                    amp,cov,energy,bcov,*_=fit(counts[t,side,:3]/omega,f);hist[side,t,:,0]=amp/baseamp;hist[side,t,:,1]=np.diag(cov);hist[side,t,:,2]=np.divide(np.diag(cov),energy,out=np.full(3,np.nan),where=energy>1e-20);hist[side,t,:,3]=np.divide(np.diag(cov),np.diag(bcov),out=np.full(3,np.nan),where=np.diag(bcov)>1e-20);covhist[side,t]=[cov,bcov]
                    np.testing.assert_allclose(cov[1,1]*m.ws,f['residual'],rtol=1e-8,atol=1e-8)
                    if time.perf_counter()-start>600:raise TimeoutError('MG164 active600s cap')
                if job['condition']=='prepared':idx=next(i for i,l in enumerate(outcome['labels']) if l!='core_compatible')
                else:idx=next((i for i,l in enumerate(outcome['labels']) if l=='core_compatible'),1000)
                rows.append(dict(id=job['id'],condition=job['condition'],side=side,event_index=idx,event_time=float(times[idx]),event_label=outcome['labels'][idx],event_boundary=idx in outcome['boundary_advisory_frames'],initial=hist[side,0].tolist(),before=hist[side,max(0,idx-1)].tolist(),event=hist[side,idx].tolist(),final=hist[side,-1].tolist(),minimum_relative_amplitude=hist[side,:,:,0].min(axis=0).tolist(),maximum_relative_amplitude=hist[side,:,:,0].max(axis=0).tolist()))
            arrays[job['id']+'_metrics']=hist;arrays[job['id']+'_covariance']=covhist
            if job['condition']=='prepared':side=0;idx=next(r['event_index'] for r in rows if r['id']==job['id'] and r['side']==0)
            else:
                indices=[next((i for i,l in enumerate(o['labels']) if l=='core_compatible'),1001) for o in job['outcomes']];side=int(indices[1]<indices[0]);idx=min(indices[side],1000)
            for tag,t in [('initial',0),('before',max(0,idx-1)),('event',idx),('final',1000)]:
                f=job['outcomes'][side]['fits'][t];amp,cov,energy,bcov,smooth,g,w=fit(counts[t,side,:3]/omega,f)
                # independent weighted normal equations, including free offset
                for sp in range(3):
                    design=np.column_stack([g[sp],np.ones(m.n)]);coef=np.linalg.solve(design.T@(w[:,None]*design),design.T@(w*smooth[sp]));err=smooth[sp]-design@coef
                    np.testing.assert_allclose(coef[0],amp[sp],rtol=1e-9,atol=1e-9);np.testing.assert_allclose((w*err*err).sum()/m.ws,cov[sp,sp],rtol=1e-9,atol=1e-9);checks+=1
                profiles[job['id']+'_'+tag]=np.roll(smooth,-f['center_site'],axis=1)[:,m.order];profiles[job['id']+'_'+tag+'_time']=np.array([times[t],side])
        arrays['times']=times;arrays['metric_names']=np.array(['relative_amplitude','local_residual_mean_square','local_residual_energy_fraction','local_residual_to_background_variance']);arrays['species']=np.array(['G','X','Y']);arrays['initial_fitted_amplitudes']=baseamp
        np.savez_compressed(OUT/'species_metrics.npz',**arrays)
        # Only cropped selected profiles are public plot exports, not full trajectory tensors.
        select=abs(m.coords[m.order])<=12.8
        export={k:(v[:,select] if v.ndim==2 else v) for k,v in profiles.items()};export['coordinates']=m.coords[m.order][select];np.savez_compressed(OUT/'selected_profile_exports.npz',**export)
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        for condition in ['prepared','homogeneous']:
            jobs=[j for j in scores['jobs'] if j['condition']==condition];fig,axes=plt.subplots(4,3,figsize=(12,10),sharex=True,layout='constrained')
            for row,j in enumerate(jobs):
                for sp,name in enumerate(['G','X','Y']):
                    ax=axes[row,sp]
                    for tag,color in [('initial','black'),('before','tab:blue'),('event','tab:red'),('final','tab:green')]:
                        stamp=profiles[j['id']+'_'+tag+'_time'];ax.plot(export['coordinates'],export[j['id']+'_'+tag][sp],color=color,lw=1,label=f'{tag} T={stamp[0]:.3g}')
                    ax.set_title(f"{j['id']} / {name}");ax.grid(alpha=.2)
                    if sp==0:ax.set_ylabel('Smoothed concentration - spatial mean');ax.legend(fontsize=6)
                    if row==3:ax.set_xlabel('Distance from shared X-fitted center')
            fig.suptitle(condition+' profiles: common X-based position, no independent species alignment');fig.savefig(OUT/(condition+'_profiles.png'),dpi=150);plt.close(fig)
        atomic_json(OUT/'diagnostics.json',dict(rows=rows,input_sha256=inputs,independent_normal_equation_checks=checks,x_residual_checks=16016,source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),elapsed_seconds=time.perf_counter()-start,model_draws=0,joint_score_defined=False))
        atomic_json(attempt,dict(status='COMPLETE',elapsed_seconds=time.perf_counter()-start))
    except BaseException:atomic_json(attempt,dict(status='FAILED',elapsed_seconds=time.perf_counter()-start,error=traceback.format_exc()));raise
if __name__=='__main__':run()

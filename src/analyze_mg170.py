"""Approximate beam model selection using known conditional Fourier noise."""
import os
for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):os.environ[k]='1'
from pathlib import Path
import time,hashlib,traceback,json
import numpy as np
from scipy.special import logsumexp
from scipy.stats import chi2
from run_mg163_long import ROOT,BASE as OLD,read,atomic_json
from mg167_xy import XYCore
OUT=ROOT/'results/MG12H3/hypothesis_tracking';INPUT=ROOT/'local_raw/MG169_tracking';PRIOR=ROOT/'results/MG12H3/tracking_feasibility';BENCH=ROOT/'results/MG12H3/XY_core_instrument/benchmark.json'
CLASSES=['absent','single','Xonly','Yonly','multiple']
def binding():return {p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in ['src/analyze_mg170.py','spec/MG170_uncertainty_hypotheses.md','src/mg167_xy.py']}
def register():
    OUT.mkdir(exist_ok=True);assert not(OUT/'plan.json').exists();paths=[OLD/'prepared_inputs.npz',BENCH]+[OLD/'production'/f'homogeneous_r{r}'/'observations.npz' for r in(0,1)]+[p for p in INPUT.glob('*rho*.npz') if p.name.startswith(('r0_','r2_'))];atomic_json(OUT/'plan.json',dict(source=binding(),inputs={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},cap=600,replicates=4,frames=24,regime_indices=[0,2],rho=[0.,.95],beam_per_nonabsent_class=2,no_new_draws=True))
class Beam:
    def __init__(self,model,L,rho,seed):
        self.model=model;self.rho=rho;self.seed=seed;self.widths=np.arange(.5,2.0001,.125);self.ft=np.array([np.fft.rfft(model.signal(s,0),axis=1) for s in self.widths]);self.phase=np.exp(-2j*np.pi*np.arange(500)[:,None]*np.arange(251)[None,:]/500);S=np.einsum('fij,fkj->fik',L,np.conj(L))[:,1:,1:]*model.x.smoother[:,None,None]**2;S[-1]=(L[-1].real@L[-1].real.T)[1:,1:]*model.x.smoother[-1]**2;self.P=np.zeros_like(S);self.P[1:]=np.linalg.inv(S[1:]);self.freqweights=np.full(251,2.);self.freqweights[[0,-1]]=1.;self.freqweights[0]=0.;self.qcut=float(chi2.ppf(.999,998));self.norms=np.array([self.inner(f,f) for f in self.ft]);self.beam=[dict(kind=k,center=seed,si=4,second=None,mean=np.zeros((2,251),complex) if k=='absent' else self.mean(k,seed,4,None),log=-np.log(5)) for k in CLASSES if k!='multiple'];self.beam.append(dict(kind='multiple',center=seed,si=4,second=None,mean=self.mean('single',seed,4,None),log=-np.log(5)));self.previous=None
        assert np.all(self.norms>0) and np.allclose(self.inner(self.ft[4],self.ft[4]),self.norms[4])
    def inner(self,a,b):return np.einsum('...sf,fst,...tf,f->...',np.conj(a),self.P,b,self.freqweights).real/500
    def mean(self,kind,k,si,second):
        f=self.ft[si]*self.phase[k]
        if kind=='Xonly':f=f.copy();f[1]=0
        if kind=='Yonly':f=f.copy();f[0]=0
        if kind=='multiple' and second is not None:f=f+self.ft[second[1]]*self.phase[second[0]]
        return f
    def step(self,field):
        y=np.fft.rfft(field-field.mean(axis=1)[:,None],axis=1);rho=0 if self.previous is None else self.rho;innovation=y if self.previous is None else y-rho*self.previous;variance=1 if self.previous is None else 1-rho*rho
        # Spectral global proposals are discovery candidates, never true secondary coordinates.
        kernels=np.einsum('fst,itf->isf',self.P,self.ft);corr=np.fft.irfft((y[None,:,:]*np.conj(kernels)).sum(axis=1),n=500,axis=1);scores=corr-.5*self.norms[:,None];global_order=np.argsort(scores.ravel())[::-1]
        candidates=[];means=[];origins=[];logprior=[];secondary_cache={}
        for pi,parent in enumerate(self.beam):
            candidates.append(dict(kind='absent',center=parent['center'],si=parent['si'],second=None));means.append(np.zeros((2,251),complex));origins.append(pi);logprior.append(parent['log']+np.log(.98 if parent['kind']=='absent' else .005))
            locations=[((parent['center']+d)%500,si) for d in range(-2,3) for si in range(max(0,parent['si']-1),min(13,parent['si']+2))]
            for k,si in locations:
                seconds=list(secondary_cache.get(k,[]))
                for idx in ([] if k in secondary_cache else global_order):
                    ss,kk=np.unravel_index(int(idx),scores.shape);dist=min((kk-k)%500,(k-kk)%500)*.2
                    if dist>self.model.x.background_exclusion and all(min((kk-q[0])%500,(q[0]-kk)%500)*.2>self.model.x.background_exclusion for q in seconds):seconds.append((int(kk),int(ss)))
                    if len(seconds)==2:break
                secondary_cache[k]=seconds
                for kind in('single','Xonly','Yonly','multiple'):
                    choices=seconds if kind=='multiple' else [None]
                    for second in choices:
                        candidates.append(dict(kind=kind,center=k,si=si,second=second));means.append(self.mean(kind,k,si,second));origins.append(pi);logprior.append(parent['log']+np.log(.98 if kind==parent['kind'] else .005)-np.log(len(locations))-np.log(len(choices)))
        mu=np.array(means);parents=np.array([self.beam[i]['mean'] for i in origins]);prediction=mu-rho*parents;err=innovation[None,:,:]-prediction;q=self.inner(err,err)/variance;log=np.array(logprior)-.5*q
        retained=[]
        for kind in CLASSES:
            indices=[i for i,c in enumerate(candidates) if c['kind']==kind];order=sorted(indices,key=lambda i:log[i],reverse=True);seen=set()
            for idx in order:
                c=candidates[idx];key=(c['center'],c['si'],str(c['second']))
                if key in seen:continue
                seen.add(key);retained.append(idx)
                if len(seen)==(1 if kind=='absent' else 2):break
        retained=np.array(retained);normal=float(logsumexp(log[retained]));weights=np.exp(log[retained]-normal);classes={k:float(sum(w for w,i in zip(weights,retained) if candidates[i]['kind']==k)) for k in CLASSES};winner=max(classes,key=classes.get);best=max([i for i in retained if candidates[i]['kind']==winner],key=lambda i:log[i]);state=candidates[best];label=winner if classes[winner]>=.95 and q[best]<=self.qcut else 'uncertain'
        single=[candidates[i] for i in retained if candidates[i]['kind']=='single'];spread=max(min((a['center']-b['center'])%500,(b['center']-a['center'])%500) for a in single for b in single);widthspread=max(abs(self.widths[a['si']]-self.widths[b['si']]) for a in single for b in single)
        if label=='single' and(spread>4 or widthspread>.25):label='uncertain'
        self.beam=[dict(**candidates[i],mean=mu[i],log=float(log[i]-normal)) for i in retained];self.previous=y
        assert abs(sum(classes.values())-1)<1e-10 and q.min()>-1e-7
        return dict(label=label,class_weights=classes,center=state['center'],width=float(self.widths[state['si']]),winner=winner,innovation_q=float(q[best]),qcut=self.qcut,single_center_spread=spread,single_width_spread=float(widthspread),beam=[dict(kind=c['kind'],center=c['center'],width=float(self.widths[c['si']]),second=c['second'],log_weight=c['log']) for c in self.beam])
def run():
    p=read(OUT/'plan.json');assert p['source']==binding() and not(OUT/'attempt.json').exists();start=time.perf_counter();atomic_json(OUT/'attempt.json',dict(status='RUNNING'))
    try:
        for name,h in p['inputs'].items():assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==h
        with np.load(OLD/'prepared_inputs.npz') as z:ref=z['reference'];initial=z['initial'];omega=float(z['omega'])
        inputs=[]
        for r in(0,1):
            with np.load(OLD/'production'/f'homogeneous_r{r}'/'observations.npz') as z:inputs.append(z['counts'][:,0,:3]/omega)
        checker=Beam.__new__(Beam);checker.P=np.tile(np.eye(2),(251,1,1)).astype(complex);checker.freqweights=np.full(251,2.);checker.freqweights[[0,-1]]=1;checker.freqweights[0]=0;grid=np.arange(500);fixture=np.array([np.sin(2*np.pi*grid/500),np.cos(6*np.pi*grid/500)]);ftest=np.fft.rfft(fixture,axis=1);np.testing.assert_allclose(checker.inner(ftest,ftest),np.sum(fixture*fixture),rtol=1e-12)
        covariances=read(BENCH)['regime_inputs'];oldsummary=read(PRIOR/'analysis.json');rows=[]
        for ri in(0,2):
            when=[.1,.5,1.][ri];model=XYCore(ref,initial,omega,np.array(covariances[ri]['covariance'])[1:,1:]);samples=np.array([x[round(when*1000)] for x in inputs]);samples-=samples.mean(axis=2)[:,:,None];ft=np.fft.rfft(samples,axis=2);S=np.einsum('rsf,rtf->fst',ft,np.conj(ft))/(2*500);S=np.stack([S[max(0,i-2):min(251,i+3)].mean(axis=0) for i in range(251)]);diag=np.zeros_like(S)
            for sp in range(3):diag[:,sp,sp]=S[:,sp,sp].real
            S=.9*S+.1*diag;ev,V=np.linalg.eigh(S);ev=np.maximum(ev,1e-8*ev.max(axis=1)[:,None]);L=V*np.sqrt(ev)[:,None,:];L[0]=0
            for rho in(0.,.95):
                for oldrow in [r for r in oldsummary['rows'] if r['regime']==when and r['rho']==rho]:
                    case=oldrow['case'];tag=f'r{ri}_rho{rho:g}_{case}';previous=read(PRIOR/(tag+'.json'));paths=[];positive=case in('static','moving','deforming');onset=0 if case=='absent' else 10
                    with np.load(INPUT/(tag+'.npz')) as z:fields=z['fields'][:4,:24]
                    for r,field in enumerate(fields):
                        seed=previous['paths'][r]['initial_seed'];truth=previous['paths'][r]['truth'][:24];beam=Beam(model,L,rho,seed);records=[]
                        if case=='absent':
                            nf=np.fft.rfft(field[0]-field[0].mean(axis=1)[:,None],axis=1);qnoise=float(beam.inner(nf,nf));assert .2*998<qnoise<3*998
                        for frame in field:
                            records.append(beam.step(frame))
                            if time.perf_counter()-start>600:raise TimeoutError('MG170 active600s cap')
                        single=sum(x['label']=='single' for x in records);accurate=sum(x['label']=='single' and min((x['center']-tr[0])%500,(tr[0]-x['center'])%500)<=2 and abs(x['width']-tr[1])<=.12500001 for x,tr in zip(records,truth));false=records[onset+5]['label']=='single' if not positive else None;success=single>=22 and accurate>=22 if positive else not false
                        oldrecords=previous['paths'][r]['records']['local'][:24];oldsingle=sum(x['accepted'] for x in oldrecords);oldaccurate=sum(x['accepted'] and x['center'] is not None and min((x['center']-tr[0])%500,(tr[0]-x['center'])%500)<=2 and abs(x['scale']-tr[1])<=.12500001 for x,tr in zip(oldrecords,truth));oldfalse=bool(oldrecords[onset+5]['accepted']) if not positive else None;oldsuccess=oldsingle>=22 and oldaccurate>=22 if positive else not oldfalse
                        paths.append(dict(replicate=r,seed=seed,records=records,truth=truth,single_frames=single,accurate_frames=accurate,false_single_deadline=false,success=bool(success),old_success=bool(oldsuccess),old_false_single_deadline=oldfalse))
                    count=sum(x['success'] for x in paths);oldcount=sum(x['old_success'] for x in paths);rows.append(dict(regime=when,rho=rho,case=case,n=4,successful_paths=count,old_successful_paths=oldcount,false_single_deadline=sum(bool(x['false_single_deadline']) for x in paths) if not positive else None,screen_passed=count==4));atomic_json(OUT/(tag+'.json'),dict(paths=paths,regime=when,rho=rho,case=case))
        atomic_json(OUT/'analysis.json',dict(rows=rows,all_screens_passed=all(r['screen_passed'] for r in rows),new_persistence_verdict=None,G_status='UNRESOLVED',elapsed_seconds=time.perf_counter()-start));atomic_json(OUT/'attempt.json',dict(status='COMPLETE',elapsed_seconds=time.perf_counter()-start))
    except BaseException:atomic_json(OUT/'attempt.json',dict(status='FAILED',elapsed_seconds=time.perf_counter()-start,error=traceback.format_exc()));raise
if __name__=='__main__':
    if __import__('sys').argv[1]=='register':register()
    else:run()

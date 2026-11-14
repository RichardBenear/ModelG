"""Frozen MG170 replay with observational candidate-mass instrumentation."""
import os
for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):os.environ[k]='1'
from pathlib import Path
import inspect,textwrap,time,json,gzip,hashlib,traceback
import numpy as np
from scipy.special import logsumexp
import analyze_mg170 as original
from analyze_mg170 import Beam,ROOT,OLD,INPUT,PRIOR,BENCH,CLASSES,read,atomic_json
from mg167_xy import XYCore
OUT=ROOT/'results/MG12H3/weight_audit';RAW=ROOT/'local_raw/MG171_weight_audit';TRACK=ROOT/'results/MG12H3/hypothesis_tracking'
FILES=['src/audit_mg171.py','src/analyze_mg170.py','src/mg167_xy.py','spec/MG171_weight_audit.md']
def binding():return {p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in FILES}
STEP=textwrap.dedent(inspect.getsource(Beam.step));old='retained=np.array(retained);normal=';new='retained=np.array(retained);self.capture(log,candidates,origins,mu,parents,innovation,variance,y,q,logprior,retained);normal=';assert STEP.count(old)==1;INSTRUMENTED=STEP.replace(old,new);namespace=dict(vars(original));exec(compile(INSTRUMENTED,'MG171_observational_step','exec'),namespace)
class AuditedBeam(Beam):
    step=namespace['step']
    def begin(self,truth,case):self.truth=truth;self.case=case;self.frame=0;self.massrows=[];self.probes=[];self.logs=[];self.priors=[];self.kinds=[];self.groups=[];self.retained=[];self.offsets=[0];self.grids=[];self.residuals=[]
    def capture(self,log,candidates,origins,mu,parents,innovation,variance,y,q,logprior,retained):
        t=self.frame;logprior=np.array(logprior);assert abs(logsumexp(logprior))<1e-9;kind=np.array([CLASSES.index(c['kind']) for c in candidates],np.uint8);groups=[];lookup={}
        for c in candidates:
            key=(c['kind'],c['center'],c['si'],str(c['second']))
            if key not in lookup:lookup[key]=len(lookup)
            groups.append(lookup[key])
        groups=np.array(groups,np.int32);A=[];B=[];C=[]
        for ci in range(5):
            indices=np.flatnonzero(kind==ci);A.append(float(logsumexp(log[indices])));maxima={}
            for i in indices:maxima[int(groups[i])]=max(maxima.get(int(groups[i]),-np.inf),float(log[i]))
            B.append(float(logsumexp(list(maxima.values()))));C.append(float(logsumexp(log[retained[kind[retained]==ci]])))
        A=np.array(A);B=np.array(B);C=np.array(C);wa=np.exp(A-logsumexp(A));wb=np.exp(B-logsumexp(B));wc=np.exp(C-logsumexp(C));self.massrows.append(dict(frame=t,all_history_logmass=A.tolist(),max_history_logmass=B.tolist(),retained_logmass=C.tolist(),all_class_weights=wa.tolist(),collapsed_class_weights=wb.tolist(),retained_class_weights=wc.tolist(),retained_fraction=np.exp(C-A).tolist(),max_collapse_fraction=np.exp(B-A).tolist(),topK_fraction=np.exp(C-B).tolist(),multiple_single_odds_shift=float((C[4]-C[1])-(A[4]-A[1])),absent_single_odds_shift=float((C[0]-C[1])-(A[0]-A[1])),all_winner=CLASSES[int(wa.argmax())],retained_winner=CLASSES[int(wc.argmax())]))
        self.logs.append(log.copy());self.priors.append(logprior);self.kinds.append(kind);self.groups.append(groups);self.retained.append(retained);self.offsets.append(self.offsets[-1]+len(log))
        if t in(0,10,15,23):
            primary=max(np.flatnonzero(kind==1),key=lambda i:log[i]);c=candidates[primary];parent=origins[primary];rho=0 if self.previous is None else self.rho;residual=innovation-(mu[primary]-rho*parents[primary]);kernels=np.einsum('fst,itf->isf',self.P,self.ft);correlation=np.fft.irfft((residual[None,:,:]*np.conj(kernels)).sum(axis=1),n=500,axis=1);grid=(correlation-.5*self.norms[:,None])/variance;distance=np.minimum((np.arange(500)-c['center'])%500,(c['center']-np.arange(500))%500)*.2;eligible=np.broadcast_to(distance[None,:]>self.model.x.background_exclusion,grid.shape);N=int(eligible.sum());others=[i for i,d in enumerate(candidates) if d['kind']=='multiple' and origins[i]==parent and d['center']==c['center'] and d['si']==c['si']];assert len(others)==2;seconds=[candidates[i]['second'] for i in others];deltas=[float(grid[s,k]) for k,s in seconds];np.testing.assert_allclose(deltas,[(q[primary]-q[i])/2 for i in others],rtol=1e-7,atol=2e-4);assert all(eligible[s,k] for k,s in seconds)
            selected=float(logsumexp(deltas)-np.log(2));full=float(logsumexp(grid[eligible])-np.log(N));represented=float(logsumexp(deltas)-np.log(N));assert full+1e-10>=represented;tr=self.truth[t];accurate=min((c['center']-tr[0])%500,(tr[0]-c['center'])%500)<=2 and abs(self.widths[c['si']]-tr[1])<=.12500001;self.probes.append(dict(frame=t,primary_center=c['center'],primary_width=float(self.widths[c['si']]),primary_accurate=bool(accurate),parent_kind=self.beam[parent]['kind'],secondary_indices=seconds,eligible=N,variance=float(variance),selected_logBF=selected,full_uniform_logBF=full,represented_uniform_logBF=represented,selected_minus_full=selected-full,all_scene_single=self.case in('static','moving','deforming') or(self.case not in('absent',) and t<10),scene_split=self.case=='split' and t>=10));self.grids.append(grid);self.residuals.append(residual)
        self.frame+=1

def register():
    OUT.mkdir(exist_ok=True);assert not(OUT/'plan.json').exists();paths=[OLD/'prepared_inputs.npz',BENCH,PRIOR/'analysis.json',TRACK/'plan.json',TRACK/'analysis.json']+[OLD/'production'/f'homogeneous_r{r}'/'observations.npz' for r in(0,1)]
    for p in sorted(TRACK.glob('r*_rho*.json.gz')):paths.extend([p,PRIOR/p.name,INPUT/(p.name.replace('.json.gz','.npz'))])
    atomic_json(OUT/'plan.json',dict(source=binding(),instrumented_step_sha256=hashlib.sha256(INSTRUMENTED.encode()).hexdigest(),inputs={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},cap=600,workers=1,threads=1,new_draws=False,classes=CLASSES))
def readgz(p):return json.loads(gzip.decompress(p.read_bytes()))
def run():
    p=read(OUT/'plan.json');assert p['source']==binding() and not(OUT/'attempt.json').exists();RAW.mkdir(exist_ok=False);start=time.perf_counter();atomic_json(OUT/'attempt.json',dict(status='RUNNING'));rows=[];probes=[];frames=0
    try:
        for name,h in p['inputs'].items():assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==h
        with np.load(OLD/'prepared_inputs.npz') as z:reference=z['reference'];initial=z['initial'];omega=float(z['omega'])
        controls=[]
        for r in(0,1):
            with np.load(OLD/'production'/f'homogeneous_r{r}'/'observations.npz') as z:controls.append(z['counts'][:,0,:3]/omega)
        cov=read(BENCH)['regime_inputs'];summary=read(TRACK/'analysis.json')
        for ri in(0,2):
            when=[.1,.5,1.][ri];model=XYCore(reference,initial,omega,np.array(cov[ri]['covariance'])[1:,1:]);data=np.array([a[round(when*1000)] for a in controls]);data-=data.mean(axis=2)[:,:,None];ft=np.fft.rfft(data,axis=2);S=np.einsum('rsf,rtf->fst',ft,np.conj(ft))/(2*500);S=np.stack([S[max(0,i-2):min(251,i+3)].mean(axis=0) for i in range(251)]);diag=np.zeros_like(S)
            for sp in range(3):diag[:,sp,sp]=S[:,sp,sp].real
            S=.9*S+.1*diag;ev,V=np.linalg.eigh(S);ev=np.maximum(ev,1e-8*ev.max(axis=1)[:,None]);L=V*np.sqrt(ev)[:,None,:];L[0]=0
            for item in [r for r in summary['rows'] if r['regime']==when]:
                rho=item['rho'];case=item['case'];tag=f'r{ri}_rho{rho:g}_{case}';saved=readgz(TRACK/(tag+'.json.gz'));case_rows=[];case_probes=[]
                with np.load(INPUT/(tag+'.npz')) as z:fields=z['fields'][:4,:24]
                for r,field in enumerate(fields):
                    old=saved['paths'][r];beam=AuditedBeam(model,L,rho,old['seed']);beam.begin(old['truth'],case)
                    for t,frame in enumerate(field):
                        result=beam.step(frame);stored=old['records'][t];assert result['label']==stored['label'];np.testing.assert_allclose(list(result['class_weights'].values()),list(stored['class_weights'].values()),rtol=1e-10,atol=1e-10);assert [(s['kind'],s['center'],s['width'],s['second']) for s in result['beam']]==[(s['kind'],s['center'],s['width'],tuple(s['second']) if s['second'] is not None else None) for s in stored['beam']];frames+=1
                        if time.perf_counter()-start>600:raise TimeoutError('MG171 audit cap600s')
                    cache=RAW/f'{tag}_p{r}.npz';np.savez_compressed(cache,logs=np.concatenate(beam.logs),priors=np.concatenate(beam.priors),kinds=np.concatenate(beam.kinds),groups=np.concatenate(beam.groups),offsets=beam.offsets,retained=np.array(beam.retained),grids=beam.grids,residuals=beam.residuals,probe_variances=[q['variance'] for q in beam.probes],precision=beam.P,template_FT=beam.ft,freqweights=beam.freqweights)
                    for row in beam.massrows:row.update(regime=when,rho=rho,case=case,replicate=r);rows.append(row);case_rows.append(row)
                    for probe in beam.probes:probe.update(regime=when,rho=rho,case=case,replicate=r);probes.append(probe);case_probes.append(probe)
                atomic_json(OUT/(tag+'.json'),dict(masses=case_rows,probes=case_probes))
                if sum(p.stat().st_size for p in RAW.iterdir())>250000000:raise RuntimeError('MG171 audit cache reserve')
        atomic_json(OUT/'analysis.json',dict(frames_replayed=frames,probes=len(probes),source_labels_unchanged=True,no_new_draws=True,elapsed_seconds=time.perf_counter()-start));atomic_json(OUT/'attempt.json',dict(status='COMPLETE',elapsed_seconds=time.perf_counter()-start))
    except BaseException:atomic_json(OUT/'attempt.json',dict(status='FAILED',elapsed_seconds=time.perf_counter()-start,error=traceback.format_exc()));raise
if __name__=='__main__':
    if __import__('sys').argv[1]=='register':register()
    else:run()

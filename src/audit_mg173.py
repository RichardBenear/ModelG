"""Exact static reference and full two-class geometry-filter prior/compression audit."""
import os
for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):os.environ[k]='1'
from pathlib import Path
import time,json,gzip,hashlib,traceback
import numpy as np
from scipy.special import logsumexp,expit
from scipy.stats import norm
from run_mg172 import ROOT,OLD,OUT as PRIOR,RAW as INPUT,read,atomic_json
OUT=ROOT/'results/MG12H3/static_core_audit';RAW=ROOT/'local_raw/MG173_static_audit'
CASES=['static','moving','absent'];MODES={'uniform':0.,'sticky90':.9,'fixed':1.}
def binding():return {p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in ['src/audit_mg173.py','spec/MG173_static_core_audit.md']}
def register():
    OUT.mkdir(exist_ok=True);assert not(OUT/'plan.json').exists();paths=[PRIOR/'plan.json',PRIOR/'analysis.json']
    for ri in(0,2):
        paths.append(INPUT/f'r{ri}_rho0_static_snap0.npz')
        for rho in(0.,.95):
            for case in CASES:
                tag=f'r{ri}_rho{rho:g}_{case}';paths.extend([INPUT/(tag+'.npz'),PRIOR/(tag+'.json.gz')])
    atomic_json(OUT/'plan.json',dict(source=binding(),inputs={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},no_new_draws=True,workers=1,threads=1,cap=600,replicates=4,frames=24,case_list=CASES,mode_parameters=MODES))
def ls(a):return float(logsumexp(a))
class Exact:
    def __init__(self,ft,P,w,rho,seed,alpha,compressed):
        self.ft=ft;self.P=P;self.w=w;self.rho=rho;self.seed=seed;self.alpha=alpha;self.compressed=compressed;self.logs=np.full((2,13,500),-np.inf);self.logs[:,4,seed]=np.log(.5);self.previous=None;self.kernels=np.einsum('fst,itf->isf',P,ft);self.norms=np.array([self.inner(f,f) for f in ft]);self.cross=np.array([np.fft.irfft((f[None,:,:]*np.conj(self.kernels)).sum(axis=1),n=500,axis=1) for f in ft]);self.stay=.98/.985;self.switch=.005/.985
    def inner(self,a,b):return float(np.einsum('sf,fst,tf,f->',np.conj(a),self.P,b,self.w).real/500)
    def step(self,field):
        y=np.fft.rfft(field-field.mean(axis=1)[:,None],axis=1);rho=0 if self.previous is None else self.rho;var=1 if self.previous is None else 1-rho*rho;innovation=y if self.previous is None else y-rho*self.previous;Q=self.inner(innovation,innovation);C=np.fft.irfft((innovation[None,:,:]*np.conj(self.kernels)).sum(axis=1),n=500,axis=1);old=self.logs.copy();new=np.full_like(old,-np.inf);q00=Q/var;q10=(Q+2*rho*C+rho*rho*self.norms[:,None])/var;new[0]=np.logaddexp(old[0]+np.log(self.stay)-.5*q00,old[1]+np.log(self.switch)-.5*q10)
        for current in range(13):
            q01=(Q-2*C[current]+self.norms[current])/var
            for previous in range(max(0,current-1),min(13,current+2)):
                nloc=5*(min(13,previous+2)-max(0,previous-1))
                for displacement in range(-2,3):
                    probability=(1-self.alpha)/nloc+(self.alpha if previous==current and displacement==0 else 0.)
                    if probability==0:continue
                    base=np.log(probability);a0=np.roll(old[0,previous],displacement);a1=np.roll(old[1,previous],displacement);prevC=np.roll(C[previous],displacement);q11=(Q-2*C[current]+2*rho*prevC+self.norms[current]+rho*rho*self.norms[previous]-2*rho*self.cross[previous,current,displacement%500])/var
                    new[1,current]=np.logaddexp(new[1,current],a0+np.log(self.switch)+base-.5*q01);new[1,current]=np.logaddexp(new[1,current],a1+np.log(self.stay)+base-.5*q11)
        masses=np.array([ls(v) for v in new]);Z=ls(masses);weights=np.exp(masses-Z);post=new-Z;fractions=[]
        if self.compressed:
            post[:]=-np.inf
            for kind,K in((0,1),(1,2)):
                indices=np.argpartition(new[kind].ravel(),-K)[-K:];vals=new[kind].ravel()[indices];selected=ls(vals);fractions.append(float(np.exp(selected-masses[kind])));post[kind].ravel()[indices]=vals+masses[kind]-selected-Z
        else:fractions=[1.,1.]
        for kind in(0,1):np.testing.assert_allclose(np.exp(ls(post[kind])),weights[kind],rtol=1e-9,atol=1e-12)
        best=np.unravel_index(int(np.argmax(new[1])),new[1].shape);result=dict(present=float(weights[1]),absent=float(weights[0]),logodds=float(masses[1]-masses[0]),center=int(best[1]),width=float(.5+.125*best[0]),retained_geometry_fraction=fractions);self.logs=post;self.previous=y;self.snapshot=dict(old=old,new=new,post=post,innovation=innovation,rho=rho,variance=var,C=C,Q=Q,alpha=self.alpha,compressed=self.compressed)
        return result

def fixed_hmm(fields,ft,P,w,rho,seed):
    mu=ft[4]*np.exp(-2j*np.pi*np.arange(251)*seed/500);means=[np.zeros_like(mu),mu];logs=np.log([.5,.5]);out=[];last=None
    for field in fields:
        y=np.fft.rfft(field-field.mean(axis=1)[:,None],axis=1);r=0 if last is None else rho;var=1 if last is None else 1-r*r;innovation=y if last is None else y-r*last;new=[]
        for current in(0,1):
            vals=[]
            for previous in(0,1):
                e=innovation-(means[current]-r*means[previous]);q=float(np.einsum('sf,fst,tf,f->',np.conj(e),P,e,w).real/500/var);vals.append(logs[previous]+np.log((.98 if previous==current else .005)/.985)-.5*q)
            new.append(ls(vals))
        logs=np.array(new)-ls(new);out.append(float(np.exp(logs[1])));last=y
    return out

def oracle(fields,ft,P,w,rho,seed):
    mu=ft[4]*np.exp(-2j*np.pi*np.arange(251)*seed/500);inner=lambda a,b:float(np.einsum('sf,fst,tf,f->',np.conj(a),P,b,w).real/500);normmu=inner(mu,mu);BF=0.;last=None;curve=[];increments=[]
    for field in fields:
        y=np.fft.rfft(field-field.mean(axis=1)[:,None],axis=1);r=0 if last is None else rho;var=1 if last is None else 1-r*r;innov=y if last is None else y-r*last;contrast=(1-r)*mu;inc=(inner(innov,contrast)-.5*inner(contrast,contrast))/var;BF+=inc;curve.append(dict(logBF=BF,present=float(expit(BF))));increments.append(inc);last=y
    effective=1+(len(fields)-1)*(1-rho)/(1+rho);D=normmu*effective;return dict(curve=curve,increments=increments,norm_mu=normmu,effective_independent_observations=effective,expected_logBF_static=D/2,sd_logBF_static=float(np.sqrt(D)),conditional_static_probability_final_ge95=float(norm.cdf((D/2-np.log(19))/np.sqrt(D))),conditional_absent_false_present_probability=float(norm.sf((np.log(19)+D/2)/np.sqrt(D))))
def preflight():
    for alpha in MODES.values():
        for previous in range(13):
            nloc=5*(min(13,previous+2)-max(0,previous-1));total=sum((1-alpha)/nloc+(alpha if current==previous and d==0 else 0) for current in range(max(0,previous-1),min(13,previous+2)) for d in range(-2,3));assert abs(total-1)<1e-12

def run():
    p=read(OUT/'plan.json');assert p['source']==binding() and not(OUT/'attempt.json').exists();start=time.perf_counter();RAW.mkdir(exist_ok=False);atomic_json(OUT/'attempt.json',dict(status='RUNNING'));checks=0;summary=[];preflight()
    try:
        for name,h in p['inputs'].items():assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==h
        for ri in(0,2):
            when=[.1,.5,1.][ri]
            with np.load(INPUT/f'r{ri}_rho0_static_snap0.npz') as z:ft=z['FT'];P=z['precision'];w=z['weights']
            for rho in(0.,.95):
                for case in CASES:
                    tag=f'r{ri}_rho{rho:g}_{case}';old=json.loads(gzip.decompress((PRIOR/(tag+'.json.gz')).read_bytes()));paths=[]
                    with np.load(INPUT/(tag+'.npz')) as z:fields=z['fields']
                    for rep in range(4):
                        seed=old['paths'][rep]['seed'];o=oracle(fields[rep],ft,P,w,rho,seed);naive=oracle(fields[rep],ft,P,w,0.,seed);hmm=fixed_hmm(fields[rep],ft,P,w,rho,seed);filters={f'{mode}_{compression}':Exact(ft,P,w,rho,seed,alpha,compression=='compressed') for mode,alpha in MODES.items() for compression in('exact','compressed')};curves={k:[] for k in filters}
                        for t,field in enumerate(fields[rep]):
                            for name,f in filters.items():
                                result=f.step(field);curves[name].append(result)
                                if rep==0 and t in(0,10,23):np.savez_compressed(RAW/f'{tag}_{name}_t{t}.npz',**f.snapshot,FT=ft,precision=P,weights=w,seed=seed)
                            np.testing.assert_allclose([curves['fixed_exact'][-1]['present'],curves['fixed_compressed'][-1]['present']],[hmm[t],hmm[t]],rtol=1e-8,atol=1e-10);checks+=1
                            if time.perf_counter()-start>600:raise TimeoutError('MG173 active600s cap')
                        paths.append(dict(replicate=rep,seed=seed,oracle=o,independence_assumed_oracle=naive,fixed_switching_hmm=hmm,curves=curves,MG172_present_weights=[r['class_weights']['single'] for r in old['paths'][rep]['records']['corrected']]))
                    atomic_json(OUT/(tag+'.json'),dict(regime=when,rho=rho,case=case,paths=paths))
                    for method in list(filters)+['oracle','independence_assumed_oracle','fixed_switching_hmm','MG172']:
                        final=[];n95=[];logodds=[]
                        for path in paths:
                            if method in filters:vals=[r['present'] for r in path['curves'][method]];lo=path['curves'][method][-1]['logodds']
                            elif method in('oracle','independence_assumed_oracle'):vals=[r['present'] for r in path[method]['curve']];lo=path[method]['curve'][-1]['logBF']
                            elif method=='fixed_switching_hmm':vals=path[method];lo=None
                            else:vals=path['MG172_present_weights'];lo=None
                            final.append(vals[-1]);n95.append(sum(v>=.95 for v in vals));logodds.append(lo)
                        summary.append(dict(regime=when,rho=rho,case=case,method=method,final_present_weights=final,frames_present_ge95=n95,final_logodds=logodds))
                    if sum(p.stat().st_size for p in RAW.iterdir())>100000000:raise RuntimeError('MG173 auditcache0.1GB')
        atomic_json(OUT/'analysis.json',dict(rows=summary,fixed_grid_vs_scalar_hmm_checks=checks,geometry_probability_rows_checked=39,no_new_draws=True,new_tracker_qualification=False,elapsed_seconds=time.perf_counter()-start));atomic_json(OUT/'attempt.json',dict(status='COMPLETE',elapsed_seconds=time.perf_counter()-start))
    except BaseException:atomic_json(OUT/'attempt.json',dict(status='FAILED',elapsed_seconds=time.perf_counter()-start,error=traceback.format_exc()));raise
if __name__=='__main__':
    if __import__('sys').argv[1]=='register':register()
    else:run()

"""Oracle geometry diagnostic on fixed MG167 fields; no new draws."""
import os
for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):os.environ[k]='1'
from pathlib import Path
import time,hashlib,traceback,gzip
import numpy as np
from run_mg167_xy import ROOT,OLD,RAW,OUT as PRIOR,read,atomic_json
from mg167_xy import XYCore
OUT=ROOT/'results/MG12H3/known_geometry_XY'
FILES=['src/analyze_mg168.py','src/mg167_xy.py','spec/MG168_known_geometry_XY.md']
def binding():return {p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in FILES}
def register():
    OUT.mkdir(exist_ok=True);assert not(OUT/'plan.json').exists();paths=[OLD/'prepared_inputs.npz',PRIOR/'benchmark.json']+sorted(RAW.glob('*.npz'))+sorted(PRIOR.glob('*_fits.json.gz'));atomic_json(OUT/'plan.json',dict(source=binding(),inputs={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},active_cap=600,no_new_draws=True))
def geometry(index,regime,case,model,mode):
    scale=[.85,1.1,1.6][(index//3)%3]
    if case.startswith('intact_'):scale=[.85,1.1,1.6][int(case[-1])]
    if mode=='nearest':scale=float(model.x.scales[np.argmin(abs(model.x.scales-scale))])
    return (37*index+71*regime)%500,scale

def oracle(model,field,center,scale):
    m=model.x;f=field-field.mean(axis=1)[:,None];template=np.roll(model.signal(scale,0),center,axis=1);w=np.roll(m.window,center);g=template-(template@w/m.ws)[:,None];fc=f-(f@w/m.ws)[:,None];kernel=(model.inv@g)*w;norm=float((g*kernel).sum());dot=float((fc*kernel).sum());amp=dot/norm;err=fc-amp*g;res=float(np.einsum('in,ij,jn,n->',err,model.inv,err,w));own=(g*g)@w;individual=(fc*g*w).sum(axis=1)/own;deviation=abs(individual-amp)*np.sqrt(own/np.diag(model.cov));remainder=f-amp*template;dots=np.fft.irfft((np.fft.rfft(remainder,axis=1)[None,:,:]*model.kfft).sum(axis=1),n=500,axis=1);gains=np.maximum(dots,0)**2/model.norm[:,None];dist=np.minimum((np.arange(500)-center)%500,(center-np.arange(500))%500)*m.h
    return dict(amplitude=float(amp),scale=float(scale),center=center,score=max(dot,0)**2/norm,residual=res,individual_amplitudes=individual.tolist(),component_stat=float(deviation.max()),secondary_score=float(gains[:,dist>m.background_exclusion].max()),bank_edge=bool(scale in(.25,3.)))
def guard_flags(model,f,c):return dict(detection=f['score']<=c['score'],common_amplitude=f['amplitude']/model.base<.5,scale=f['bank_edge'] or not .5<=f['scale']<=2,component_amplitude=bool(np.any(np.array(f['individual_amplitudes'])/model.initial_components<.5)),component_consistency=f['component_stat']>c['component'],secondary=f['secondary_score']>c['secondary'],residual=f['residual']>c['residual'])
def run():
    start=time.perf_counter();p=read(OUT/'plan.json');assert p['source']==binding() and not(OUT/'attempt.json').exists();atomic_json(OUT/'attempt.json',dict(status='RUNNING'));checks=0;rows=[];allcal=[]
    try:
        for name,h in p['inputs'].items():assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==h
        with np.load(OLD/'prepared_inputs.npz') as z:ref=z['reference'];initial=z['initial'];omega=float(z['omega'])
        bench=read(PRIOR/'benchmark.json')
        for ri,regime in enumerate(bench['regime_inputs']):
            model=XYCore(ref,initial,omega,np.array(regime['covariance'])[1:,1:])
            for mode in('nearest','exact'):
                def fitcase(case):
                    nonlocal checks
                    with np.load(RAW/f'r{ri}_{case}.npz') as z:fields=z['fields']
                    fits=[]
                    for index,field in enumerate(fields):
                        center,scale=geometry(index,ri,case,model,mode);f=oracle(model,field,center,scale);fits.append(f)
                        if index in(0,len(fields)-1):
                            m=model.x;w=np.roll(m.window,center);template=model.signal(scale,center);A=np.zeros((2,500,3));A[:,:,0]=template
                            for sp in range(2):A[sp,:,sp+1]=1
                            normal=np.einsum('ina,ij,jnb,n->ab',A,model.inv,A,w);rhs=np.einsum('ina,ij,jn,n->a',A,model.inv,field,w);coef=np.linalg.solve(normal,rhs);err=field-np.einsum('ina,a->in',A,coef);res=np.einsum('in,ij,jn,n->',err,model.inv,err,w);np.testing.assert_allclose([coef[0],res],[f['amplitude'],f['residual']],rtol=1e-8,atol=1e-8);checks+=1
                        if time.perf_counter()-start>600:raise TimeoutError('MG168 cap600s')
                    return fields,fits
                _,null=fitcase('cal_null');_,good=fitcase('cal_shape');cal=dict(score=max(f['score'] for f in null),residual=max(f['residual'] for f in good),component=max(f['component_stat'] for f in good),secondary=max(f['secondary_score'] for f in good));allcal.append(dict(regime=regime['time'],mode=mode,thresholds=cal));atomic_json(OUT/f'r{ri}_{mode}_calibration.json',cal)
                for old in [r for r in bench['cases'] if r['regime']==regime['time']]:
                    case=old['case'];fields,fits=fitcase(case);flags=[guard_flags(model,f,cal) for f in fits];accepted=[not any(g.values()) for g in flags];count=sum(accepted);positive=case.startswith('intact');passed=count/len(fits)>=.9 if positive else count/len(fits)<=(.01 if case=='null' else .05)
                    baseline=read_gzip(PRIOR/f'r{ri}_{case}_fits.json.gz')['labels'];base=np.array([l=='candidate_compatible' for l in baseline]);now=np.array(accepted);rows.append(dict(regime=regime['time'],mode=mode,case=case,n=len(fits),accepted=count,fraction=count/len(fits),baseline_accepted=int(base.sum()),gain_only_passed=sum(not g['detection'] for g in flags),screen_passed=passed,guard_failures={k:sum(g[k] for g in flags) for k in flags[0]},baseline_rejected_oracle_accepted=int((~base&now).sum()),baseline_accepted_oracle_rejected=int((base&~now).sum())))
                    payload=dict(fits=fits,accepted=accepted,flags=flags);atomic_json(OUT/f'r{ri}_{mode}_{case}.json',payload)
        atomic_json(OUT/'analysis.json',dict(rows=rows,calibrations=allcal,weighted_normal_equation_checks=checks,mode_screens={mode:all(r['screen_passed'] for r in rows if r['mode']==mode) for mode in('nearest','exact')},no_new_draws=True,G_status='UNRESOLVED',new_persistence_verdict=None,elapsed_seconds=time.perf_counter()-start));atomic_json(OUT/'attempt.json',dict(status='COMPLETE',elapsed_seconds=time.perf_counter()-start))
    except BaseException:atomic_json(OUT/'attempt.json',dict(status='FAILED',elapsed_seconds=time.perf_counter()-start,error=traceback.format_exc()));raise

def read_gzip(p):return __import__('json').loads(gzip.decompress(p.read_bytes()))
if __name__=='__main__':
    if __import__('sys').argv[1]=='register':register()
    else:run()

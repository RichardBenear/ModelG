"""MG161 bounded instrument benchmark; fresh synthetic noise, no Model G runs."""
import os
for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):os.environ[k]='1'
from pathlib import Path
import time,json,hashlib,subprocess,traceback,math,sys,scipy
import numpy as np
from scipy.stats import beta
from mg161_local_shape import LocalShape
from mg160_local_archive import local_archive
from runner import seed_for
from mg10e2c1_pool import atomic_json
ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/'results/MG12H3/noise_aware_shape';RAW=ROOT/'local_raw/MG161_noise_shape';ARCHIVE=local_archive(ROOT/'results/MG12H3/pattern_persistence')
FILES=['src/mg161_local_shape.py','src/run_mg161_shape.py','spec/MG161_noise_aware_shape.md']
def read(p):return json.loads(Path(p).read_text(encoding='utf-8'))
def input_paths():return [ARCHIVE/'prepared_inputs.npz']+[ARCHIVE/'P1_production'/f'homogeneous_r{r}'/'observations.npz' for r in (0,1)]
def hashes():return {p:hashlib.sha256((ROOT/p).read_text(encoding='utf-8').replace('\r\n','\n').encode()).hexdigest() for p in FILES}
def setup():
    with np.load(ARCHIVE/'prepared_inputs.npz') as z:initial=z['initial'];reference=z['reference'];omega=float(z['omega'])
    model=LocalShape(reference[1],initial[1],omega);noise=[]
    for p in input_paths()[1:]:
        with np.load(p) as z:x=z['counts'][-1,0,1]/omega;noise.append(x-x.mean())
    psd=np.mean(abs(np.fft.rfft(noise,axis=1))**2/500,axis=0);psd=np.convolve(np.pad(psd,(2,2),mode='reflect'),np.ones(5)/5,mode='valid');floor=.01*np.median(psd[1:]);psd=np.maximum(psd,floor);psd[0]=0
    return model,psd,float(reference[1].mean())
def plan():
    assert not (BASE/'plan.json').exists();prior=set()
    def collect(x):
        if isinstance(x,dict):
            for v in x.values():collect(v)
        elif isinstance(x,list):
            for v in x:collect(v)
        elif type(x) is int and 0<=x<2**32:prior.add(x)
    for p in (ROOT/'results').rglob('*.json'):
        if 'plan' in p.name or 'registration' in p.name:collect(read(p))
    cases=[dict(id='cal_null',kind='null',N=1024),dict(id='cal_shape',kind='signal',N=1024,amplitude=1.,scale=1.)]
    cases.append(dict(id='null',kind='null',N=4096,primary=True,truth='reject'))
    for a in (.8,1.,1.2):
        for s in (.75,1.,1.5):cases.append(dict(id=f'retain_a{a:g}_s{s:g}',kind='signal',N=512,amplitude=a,scale=s,primary=True,truth='compatible'))
    for name,a,s,kind in [('weak',.2,1.,'signal'),('narrow',1.,.25,'signal'),('broad',1.,3.,'signal'),('inverted',-1.,1.,'signal'),('split',1.,1.,'split')]:cases.append(dict(id=name,kind=kind,N=512,amplitude=a,scale=s,primary=True,truth='reject'))
    cases.append(dict(id='double_noise_sensitivity',kind='signal',N=512,amplitude=1.,scale=1.,noise_scale=2.,primary=False))
    used=set(prior)
    for i,c in enumerate(cases):
        nonce=0;v=seed_for(12,i,0,16100,0,nonce)
        while v in used:nonce+=1;v=seed_for(12,i,0,16100,0,nonce)
        used.add(v);c.update(seed=v,seed_nonce=nonce,index=i)
    model,psd,_=setup();np.savez_compressed(BASE/'noise_model.npz',psd=psd,window=model.window,scales=model.scales,template_raw=model.raw,radius95=model.radius95,initial_amplitude=model.initial_amplitude)
    atomic_json(BASE/'plan.json',dict(stage='MG161',authorization='go ahead with the noise-aware shape measurement.',scope='bounded new-instrument development/known-truth synthetic benchmark; not a persistence rescore or new Model G simulation',source_hashes=hashes(),inputs={str(p.relative_to(ARCHIVE)):hashlib.sha256(p.read_bytes()).hexdigest() for p in input_paths()},cases=cases,cap_seconds=300,workers=1,threads=1,environment=dict(python=sys.version,numpy=np.__version__,scipy=scipy.__version__,rng='PCG64'),primary_bounds=15,bound_alpha=.05/15,raw_location='local_raw/MG161_noise_shape',threshold_calibration='maximum of1024 separate null-search scores and maximum of1024 correct-template residuals; fixed before held-out cases'))
def evaluate_case(c,model,psd,background,deadline):
    rng=np.random.default_rng(c['seed']);fields=[];features=[];names=['amplitude','scale','center_site','shift_sites','score','residual','secondary_score','secondary_ratio','bank_edge'];start=time.perf_counter()
    try:
        for r in range(c['N']):
            if time.perf_counter()>deadline:raise TimeoutError('MG161 benchmark cap')
            sigma=c.get('noise_scale',1.);noise=np.fft.irfft(np.fft.rfft(rng.normal(size=500))*np.sqrt(psd),n=500)*sigma;shift=(37*r+71*c['index'])%500;signal=np.zeros(500)
            if c['kind'] in ('signal','split'):signal=c['amplitude']*model.signal(c['scale'],shift)
            if c['kind']=='split':signal+=c['amplitude']*model.signal(c['scale'],(shift+250)%500)
            x=background+signal+noise;out=model.fit(x,sigma);fields.append(x);features.append([out[k] for k in names])
    finally:
        np.savez_compressed(RAW/(c['id']+'.npz'),fields=fields,features=features,feature_names=names,seed=c['seed'],case=json.dumps(c));atomic_json(RAW/(c['id']+'_metadata.json'),dict(case=c,completed=len(fields),seconds=time.perf_counter()-start))
    features=np.asarray(features);np.savez_compressed(BASE/(c['id']+'_features.npz'),features=features,feature_names=names,seed=c['seed']);return features
def run():
    p=read(BASE/'plan.json');assert p['source_hashes']==hashes() and not (BASE/'attempt.json').exists();RAW.mkdir(parents=True,exist_ok=False);atomic_json(BASE/'attempt.json',dict(status='RUNNING'));start=time.perf_counter();deadline=start+300
    try:
        model,psd,background=setup();cal=[]
        for c in p['cases'][:2]:cal.append(evaluate_case(c,model,psd,background,deadline))
        detection=float(cal[0][:,4].max());residual=float(cal[1][:,5].max());atomic_json(BASE/'calibration.json',dict(detection_threshold=detection,residual_threshold=residual,source='separate calibration samples; maxima frozen before any validation samples',null_rank_bound_unconditional=1/1025,qualification_limit='Rank bound holds only for exchangeable draws of the specified conditional synthetic null; shape residual calibration is at amplitude1/scale1.'))
        rows=[]
        for c in p['cases'][2:]:
            features=evaluate_case(c,model,psd,background,deadline);labels=[]
            for v in features:
                result=dict(zip(['amplitude','scale','center_site','shift_sites','score','residual','secondary_score','secondary_ratio','bank_edge'],v));labels.append(model.classify(result,detection,residual))
            k=labels.count('core_compatible');n=len(labels);row=dict(case=c,compatible=k,N=n,label_counts={label:labels.count(label) for label in sorted(set(labels))})
            if c.get('primary'):
                alpha=p['bound_alpha'];upper=float(beta.ppf(1-alpha,k+1,n-k)) if k<n else 1.;lower=float(beta.ppf(alpha,k,n-k+1)) if k>0 else 0.;row.update(lower=lower,upper=upper)
                row['passed']=lower>=.9 if c['truth']=='compatible' else upper<=(.01 if c['id']=='null' else .05)
            rows.append(row);atomic_json(BASE/'progress.json',dict(cases_complete=len(rows),required=len(p['cases'])-2))
        assert time.perf_counter()<deadline and hashes()==p['source_hashes'];atomic_json(BASE/'benchmark_results.json',dict(qualified_within_synthetic_scope=all(r['passed'] for r in rows if r['case'].get('primary')),cases=rows,elapsed_seconds=time.perf_counter()-start,source_hashes=hashes(),scope='Conditional stationary colored-Gaussian instrument benchmark; not real-trajectory survival qualification or MG160-P1 rescore.'));atomic_json(BASE/'attempt.json',dict(status='COMPLETE',elapsed_seconds=time.perf_counter()-start));print(json.dumps([dict(id=r['case']['id'],compatible=r['compatible'],N=r['N'],passed=r.get('passed')) for r in rows]))
    except BaseException:
        atomic_json(BASE/'failure.json',dict(error=traceback.format_exc(),elapsed_seconds=time.perf_counter()-start));atomic_json(BASE/'attempt.json',dict(status='FAILED',elapsed_seconds=time.perf_counter()-start));raise
if __name__=='__main__':
    if sys.argv[1]=='plan':BASE.mkdir(parents=True,exist_ok=True);plan()
    elif sys.argv[1]=='run':run()
    else:raise ValueError('Choose plan/run')

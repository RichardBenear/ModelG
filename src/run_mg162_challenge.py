"""MG162 exhaustive real-noise replay/off-grid challenge; no kinetic trajectories."""
import os
for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):os.environ[k]='1'
from pathlib import Path
import json,time,hashlib,traceback,sys
import numpy as np
from run_mg161_shape import setup,ARCHIVE
from mg10e2c1_pool import atomic_json
ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/'results/MG12H3/shape_real_noise_challenge';RAW=ROOT/'local_raw/MG162_shape_challenge'
FILES=['src/run_mg162_challenge.py','src/mg161_local_shape.py','src/run_mg161_shape.py','spec/MG162_shape_validation.md']
def read(p):return json.loads(Path(p).read_text(encoding='utf-8'))
def hashes():return {p:hashlib.sha256((ROOT/p).read_text(encoding='utf-8').replace('\r\n','\n').encode()).hexdigest() for p in FILES}
def inputs():return [ARCHIVE/'P1_production'/f'homogeneous_r{r}'/'observations.npz' for r in (2,3)]
def plan():
    assert not (BASE/'plan.json').exists();cases=[dict(id='null',kind='null',primary=True,truth='reject',amplitude=0.,scale=1.)]
    for a in (.8,1.,1.2):
        for s in (.85,1.1,1.6):cases.append(dict(id=f'interior_a{a:g}_s{s:g}',kind='signal',primary=True,truth='compatible',amplitude=a,scale=s))
    for name,a,s,kind in [('weak',.2,1.1,'signal'),('narrow',1.,.3,'signal'),('broad',1.,2.6,'signal'),('inverted',-1.,1.1,'signal'),('split',1.,1.1,'split')]:cases.append(dict(id=name,amplitude=a,scale=s,kind=kind,primary=True,truth='reject'))
    for a in (.4,.45,.49,.5,.51,.55,.6):cases.append(dict(id=f'boundary_amplitude_{a:g}',amplitude=a,scale=1.,kind='signal',primary=False))
    for s in (.4,.45,.49,.5,.51,.55,.6,1.9,1.95,1.99,2.,2.01,2.05,2.1):cases.append(dict(id=f'boundary_scale_{s:g}',amplitude=1.,scale=s,kind='signal',primary=False))
    atomic_json(BASE/'plan.json',dict(stage='MG162A',authorization='ok, go ahead with the recommended next step',scope='exhaustive finite empirical-noise replay challenge; no new Model G trajectories or stochastic draws',source_hashes=hashes(),cases=cases,parents=[2,3],reflections=[False,True],rotations=500,fields_per_case=2000,cap_seconds=300,workers=1,threads=1,inputs={str(p.relative_to(ARCHIVE)):hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs()},calibration_sha256=hashlib.sha256((ROOT/'results/MG12H3/noise_aware_shape/calibration.json').read_bytes()).hexdigest(),primary_positive_min=.9,primary_negative_max=.05,null_max=.01,raw_location='local_raw/MG162_shape_challenge',limits='Two calibration-excluded P1 parent fields were previously seen; rotations are not independent new real trajectories. Exact finite-orbit fractions, not population confidence intervals.'))
def run():
    p=read(BASE/'plan.json');assert p['source_hashes']==hashes() and not (BASE/'attempt.json').exists();RAW.mkdir(parents=True,exist_ok=False);start=time.perf_counter();deadline=start+300;atomic_json(BASE/'attempt.json',dict(status='RUNNING'))
    try:
        model,_,background=setup();calpath=ROOT/'results/MG12H3/noise_aware_shape/calibration.json';assert hashlib.sha256(calpath.read_bytes()).hexdigest()==p['calibration_sha256'];cal=read(calpath);noise=[]
        for path in inputs():
            assert hashlib.sha256(path.read_bytes()).hexdigest()==p['inputs'][str(path.relative_to(ARCHIVE))]
            with np.load(path) as z:x=z['counts'][-1,0,1]/model.omega;noise.append(x-x.mean())
        rows=[]
        for c in p['cases']:
            fields=[];features=[];labels=[];groups=[]
            try:
                signal=np.zeros(500)
                if c['kind']!='null':signal=c['amplitude']*model.signal(c['scale'],model.anchor)
                if c['kind']=='split':signal+=c['amplitude']*model.signal(c['scale'],(model.anchor+250)%500)
                for parent in (0,1):
                    for reflected in (False,True):
                        base=noise[parent][::-1] if reflected else noise[parent];group=[]
                        for rotation in range(500):
                            if time.perf_counter()>deadline:raise TimeoutError('MG162A cap')
                            field=background+signal+np.roll(base,rotation);r=model.fit(field);label=model.classify(r,cal['detection_threshold'],cal['residual_threshold']);fields.append(field);features.append([r[k] for k in ('amplitude','scale','center_site','shift_sites','score','residual','secondary_score','secondary_ratio','bank_edge')]);labels.append(label);group.append(label)
                        compatible=group.count('core_compatible');groups.append(dict(parent=parent+2,reflected=reflected,compatible=compatible,N=500,fraction=compatible/500))
            finally:
                np.savez_compressed(RAW/(c['id']+'.npz'),fields=fields,features=features,labels=labels,case=json.dumps(c))
            np.savez_compressed(BASE/(c['id']+'_features.npz'),features=features,labels=labels);row=dict(case=c,groups=groups,compatible=labels.count('core_compatible'),N=len(labels),label_counts={l:labels.count(l) for l in sorted(set(labels))})
            if c['primary']:row['passed']=all(g['fraction']>=.9 for g in groups) if c['truth']=='compatible' else all(g['fraction']<=(.01 if c['id']=='null' else .05) for g in groups)
            rows.append(row);atomic_json(BASE/'progress.json',dict(cases_complete=len(rows),required=len(p['cases'])))
        assert hashes()==p['source_hashes'];atomic_json(BASE/'challenge_results.json',dict(supported_within_finite_challenge=all(r['passed'] for r in rows if r['case']['primary']),cases=rows,source_hashes=hashes(),elapsed_seconds=time.perf_counter()-start,scope=p['scope'],limits=p['limits']));atomic_json(BASE/'attempt.json',dict(status='COMPLETE',elapsed_seconds=time.perf_counter()-start));print(json.dumps([dict(id=r['case']['id'],compatible=r['compatible'],N=r['N'],passed=r.get('passed')) for r in rows]))
    except BaseException:atomic_json(BASE/'attempt.json',dict(status='FAILED',elapsed_seconds=time.perf_counter()-start,error=traceback.format_exc()));raise
if __name__=='__main__':
    if sys.argv[1]=='plan':BASE.mkdir(parents=True,exist_ok=True);plan()
    elif sys.argv[1]=='run':run()
    else:raise ValueError('Choose plan/run')

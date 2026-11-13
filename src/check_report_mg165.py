"""Saved-only independent score search/labels and compact MG165 development report."""
import time,json,gzip,hashlib
import numpy as np
from run_mg165_joint_v2 import OUT,RAW,OLD,ROOT,read,atomic_json
from mg165_joint import Joint

def report():
    start=time.perf_counter();bench=read(OUT/'benchmark.json');analysis=read(OUT/'analysis.json');checked=0;searchchecks=0;cases=[]
    with np.load(OLD/'prepared_inputs.npz') as z:ref=z['reference'];initial=z['initial'];omega=float(z['omega'])
    for ri,regime in enumerate(bench['regime_inputs']):
        model=Joint(ref,initial,omega,np.array(regime['covariance']));cal=bench['calibrations'][ri]
        for name,field in [('cal_null','score'),('cal_shape','residual'),('cal_shape','component_stat')]:
            fs=read(OUT/f'r{ri}_{name}_fits.json')['fits'];value=max(f[field] for f in fs);key={'score':'score','residual':'residual','component_stat':'component'}[field];assert value==cal[key]
        for row in [r for r in bench['cases'] if r['regime']==regime['time']]:
            record=read(OUT/f"r{ri}_{row['case']}_fits.json");labels=[]
            for f in record['fits']:
                if f['score']<=cal['score'] or f['amplitude']/model.base<.5:label='weak'
                elif f['bank_edge'] or not .5<=f['scale']<=2:label='scale'
                elif f['component_stat']>cal['component']:label='component_mismatch'
                elif f['residual']>cal['residual']:label='shape_mismatch'
                else:label='candidate_compatible'
                labels.append(label);checked+=1
            assert labels==record['labels'];count=labels.count('candidate_compatible');assert count==row['compatible'];positive=row['case'].startswith('intact');passed=count/row['n']>=.9 if positive else count/row['n']<=(.01 if row['case']=='null' else .05);assert passed==row['screen_passed'];cases.append(passed)
            if row['case'] in('intact_1','missing0','null'):
                with np.load(RAW/f"r{ri}_{row['case']}.npz") as z:f=z['fields'][0];f=f-f.mean(axis=1)[:,None]
                rolls=np.array([np.roll(f,-k,axis=1).ravel() for k in range(500)]);dots=model.k.reshape(23,-1)@rolls.T;gain=np.maximum(dots,0)**2/model.norm[:,None];np.testing.assert_allclose(gain.max(),record['fits'][0]['score'],rtol=1e-9,atol=1e-9);searchchecks+=1
    assert all(cases)==bench['conditional_development_screen_passed']
    atomic_json(OUT/'separate_check.json',dict(passed=True,labels_recomputed=checked,full_search_direct_spatial_checks=searchchecks,benchmark_verdict_recomputed=True,scope='Direct spatial full-bank score maxima and separately branched labels/fractions; shared templates and fitted features. Same author, not external replication.'))
    # Preserve uncompressed diagnostics locally, publish lossless compact JSON.gz.
    for p in list(OUT.glob('*_fits.json'))+list(OUT.glob('*_candidate.json')):
        q=p.with_suffix(p.suffix+'.gz');q.write_bytes(gzip.compress(p.read_bytes(),mtime=0));assert gzip.decompress(q.read_bytes())==p.read_bytes()
    manifest={p.name:dict(bytes=p.stat().st_size,sha256=hashlib.sha256(p.read_bytes()).hexdigest()) for p in RAW.iterdir() if p.is_file()};atomic_json(OUT/'local_raw_manifest.json',dict(directory=str(RAW),files=manifest,total_bytes=sum(x['bytes'] for x in manifest.values()),full_raw_published=False))
    lines=['# MG165 - multivariate instrument development and saved-data analysis','',
    'The joint G/X/Y candidate has been implemented, but it is not qualified to decide whether the complete particle persists. Shared geometry and covariance weighting improve acceptance of many prepared profiles, yet the benchmark frequently accepts deliberately damaged structures and misses some weaker intact cores in later-time noise. Removing or displacing G is particularly difficult to detect because its reference spatial feature is weak. This is useful evidence about what the available data can resolve: adding species and noise weights does not by itself create a reliable particle test. The saved MG163 application is exploratory, not a revised persistence verdict.','',
    '## Instrument','',
    'All species share one position, one dilation and one positive amplitude, with three free local background offsets. The search covers23dilations/500positions. A generalized weighted least-squares residual uses full3x3species covariance, with species-specific component-consistency checks so averaging cannot silently erase a component mismatch. Three fixed noise regimes use already-seen homogeneous-control0/1coarse snapshots atT=0.1,0.5,1; saved observations use the nearest regime, with earlier ties. Spatial/cross-species Gaussian synthetic noise uses regularized cross-spectral matrices. Shrinkage, cases, threshold rules and seeds were fixed before their evaluation. This is a conditional development model, not demonstrated nonlinear-noise physics or a validated state-dependent estimator.','',
    '## Implementation correction retained','',
    'The first attempt incorrectly included a signal in its null-calibration samples. That invalidates every v1threshold/classification; v1source, plan, timing and local raw data are retained. Code inspection identified the bug before any v1outcomes were inspected or used for tuning. Version2corrected only null construction, froze new source hashes/registration, retained the same seed/cases/criteria and wrote separate outputs. Its null generator explicitly checks zero template contribution. All results below are v2; v1is not counted as a scientific failure or discarded valid result.','',
    '## Frozen conditional development screen: FAILED','',
    'Each regime used256null/256correct-shape calibration fields, then independent synthetic validation fields. Thresholds are calibration maxima; no thresholds were relaxed after validation. Required raw acceptance fractions: intact>=90%, damaged<=5%,null<=1%. These are development screening criteria, not confidence-qualified population error rates. Individual95%Clopper-Pearson intervals in benchmark.json are descriptive and not familywise guarantees.','',
    '| Case | T0.1 accepted | T0.5 accepted | T1 accepted |','|---|---:|---:|---:|']
    for name in ['null','intact_0','intact_1','intact_2','missing0','missing1','missing2','shift0','shift1','shift2','split']:
        rows=[next(r for r in bench['cases'] if r['case']==name and r['regime']==t) for t in(.1,.5,1.)];label={'missing0':'G removed','missing1':'X removed','missing2':'Y removed','shift0':'G displaced','shift1':'X displaced','shift2':'Y displaced','split':'Two separated cores','intact_0':'Intact a0.8/s0.85','intact_1':'Intact a1/s1.1','intact_2':'Intact a1.2/s1.6','null':'No core'}[name];lines.append('| '+label+' | '+' | '.join(f"{r['compatible']}/{r['n']}" for r in rows)+' |')
    lines+=['',
    'G-removed and G-displaced cases are often accepted. Component checks cannot recover information that is buried in noise. X removal is rejected reliably on this benchmark, but displacement can be accepted when another species supports the common fit. Some Y damage is also missed. Two far-separated cores often pass because this candidate tests one local core and lacks a successful second-core veto; it is not qualified for single-object identity. Later weak-intact cases also fail the recall requirement: broader or more permissive weighting is not automatically more accurate. These failures are preserved, not tuned away.','',
    '## Exploratory application to saved MG163','',
    '| Path pair | Coarse candidate-compatible frames | Fine candidate-compatible frames |','|---|---:|---:|']
    for r in analysis['application']:lines.append(f"| {r['id']} | {r['compatible_frames'][0]}/1001 | {r['compatible_frames'][1]}/1001 |")
    lines+=['',
    'More accepted prepared frames than the old X-only score do not validate this candidate. It has already demonstrated acceptance of damaged examples. Control crossings also remain. All these trajectories were seen during development; none supplies independent validation. There is no new particle-persistence classification, and MG163 remains NOT ESTABLISHED under its original score.','',
    '## Checks, compute and retention','',
    f"Separate weighted-normal-equation checks passed {analysis['independent_GLS_checks']} selected synthetic/saved fits, including per-species residuals. A further {searchchecks}full-bank direct spatial checks reproduced maximum scores, and {checked}validation labels and all33case fractions/verdicts were recomputed. Shared reference/templates/inputs remain dependencies. Same-author checks are not external replication.",'',
    f"Corrected benchmark plus exploratory application took {analysis['elapsed_seconds']:.6g}seconds on one worker/library-thread1; invalid v1took10.5741seconds. No new Model G simulations. All syntheticfields remain locally in local_raw/MG165_multivariate(v1) and local_raw/MG165_multivariate_v2; v2payload {sum(x['bytes'] for x in manifest.values())/1e6:.6g}MB. Local raw existence/hashes are recorded; no immutable raw Git archive is claimed. Public lossless gzip files contain fit/label decisions, not full fields/events. Planning, attempts and correction provenance are retained.",'',
    '## Recommended next step','',
    'Resolve identifiability before another long run: quantify the minimum recoverable G feature and how much observation averaging could improve it without erasing dynamics. Consider a separately registered X/Y core instrument with G as an explicitly unresolved consistency channel, rather than claiming a fully observed G/X/Y particle. Validate any such change against missing/displaced species, multiple cores, time-varying backgrounds and independent noise realizations. Add a separately validated competing-core veto. More synthetic training or a more permissive threshold cannot substitute for these tests. No further simulations or instrument revisions are authorized by this report.','']
    text='\n'.join(lines)
    for a,b in [('cover23','cover 23'),('full3x3species','full 3x3 species'),('used256null/256correct','used 256 null/256 correct'),('v1threshold','v1 threshold'),('v1source','v1 source'),('v1outcomes','v1 outcomes'),('Version2corrected','Version 2 corrected'),('v1is','v1 is'),('Individual95%Clopper','Individual 95% Clopper'),('all33case','all 33 case'),('}full-bank','} full-bank')]:text=text.replace(a,b)
    (OUT/'MG165_multivariate_results.md').write_text(text,encoding='utf-8');atomic_json(OUT/'report_attempt.json',dict(status='COMPLETE',elapsed_seconds=time.perf_counter()-start))
if __name__=='__main__':report()

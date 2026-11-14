"""Recompute X/Y decisions, full spatial maxima and residual secondary search."""
import time,json,gzip,hashlib,collections
import numpy as np
from run_mg167_xy import OUT,RAW,OLD,ROOT,read,atomic_json,binding
from mg167_xy import XYCore

def report():
    start=time.perf_counter();b=read(OUT/'benchmark.json');a=read(OUT/'analysis.json');checks=0;spatial=0;passed=[];reasons={}
    with np.load(OLD/'prepared_inputs.npz') as z:reference=z['reference'];initial=z['initial'];omega=float(z['omega'])
    for ri,regime in enumerate(b['regime_inputs']):
        model=XYCore(reference,initial,omega,np.array(regime['covariance'])[1:,1:]);cal=b['calibrations'][ri];null=read(OUT/f'r{ri}_cal_null_fits.json')['fits'];good=read(OUT/f'r{ri}_cal_shape_fits.json')['fits']
        assert cal==dict(score=max(f['score'] for f in null),residual=max(f['residual'] for f in good),component=max(f['component_stat'] for f in good),secondary=max(f['secondary_score'] for f in good))
        for r in [r for r in b['cases'] if r['regime']==regime['time']]:
            case=r['case'];d=read(OUT/f'r{ri}_{case}_fits.json');labels=[]
            for f in d['fits']:
                if f['score']<=cal['score'] or f['amplitude']/model.base<.5:l='weak'
                elif f['bank_edge'] or not .5<=f['scale']<=2:l='scale'
                elif min(np.array(f['individual_amplitudes'])/model.initial_components)<.5:l='component_weak'
                elif f['component_stat']>cal['component']:l='component_mismatch'
                elif f['secondary_score']>cal['secondary']:l='competing_core'
                elif f['residual']>cal['residual']:l='shape_mismatch'
                else:l='candidate_compatible'
                labels.append(l);checks+=1
            assert labels==d['labels'];n=labels.count('candidate_compatible');assert n==r['compatible'];positive=case.startswith('intact');goodcase=n/r['n']>=.9 if positive else n/r['n']<=(.01 if case=='null' else .05);assert goodcase==r['screen_passed'];passed.append(goodcase);reasons[f'r{ri}_{case}']=dict(collections.Counter(labels))
            with np.load(RAW/f'r{ri}_{case}.npz') as z:f=z['fields'][0];f=f-f.mean(axis=1)[:,None]
            rolls=np.array([np.roll(f,-k,axis=1).ravel() for k in range(500)]);dots=model.k.reshape(23,-1)@rolls.T;gains=np.maximum(dots,0)**2/model.norm[:,None];si,k=np.unravel_index(gains.argmax(),gains.shape);saved=d['fits'][0];np.testing.assert_allclose(gains.max(),saved['score'],rtol=1e-9,atol=1e-9);assert k==saved['center'] and model.x.scales[si]==saved['scale']
            remaining=f-saved['amplitude']*np.roll(model.bank[si],k,axis=1);rolls=np.array([np.roll(remaining,-j,axis=1).ravel() for j in range(500)]);dots=model.k.reshape(23,-1)@rolls.T;gains=np.maximum(dots,0)**2/model.norm[:,None];dist=np.minimum((np.arange(500)-k)%500,(k-np.arange(500))%500)*model.x.h;np.testing.assert_allclose(gains[:,dist>model.x.background_exclusion].max(),saved['secondary_score'],rtol=1e-9,atol=1e-9);spatial+=1
    assert all(passed)==b['conditional_development_screen_passed'];assert a['saved_application_blocked'] and not a['application'];assert binding()==read(OUT/'plan.json')['source']
    atomic_json(OUT/'separate_check.json',dict(passed=True,labels_checked=checks,full_spatial_primary_and_secondary_checks=spatial,weighted_normal_equation_checks=a['independent_GLS_checks'],reason_counts=reasons,source_binding_unchanged=True))
    for p in OUT.glob('*_fits.json'):
        q=p.with_suffix('.json.gz');q.write_bytes(gzip.compress(p.read_bytes(),mtime=0));assert gzip.decompress(q.read_bytes())==p.read_bytes()
    manifest={p.name:dict(bytes=p.stat().st_size,sha256=hashlib.sha256(p.read_bytes()).hexdigest()) for p in RAW.iterdir() if p.is_file()};atomic_json(OUT/'local_raw_manifest.json',dict(directory=str(RAW),files=manifest,total_bytes=sum(x['bytes'] for x in manifest.values()),full_fields_published=False))
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(1,2,figsize=(11,4.5),layout='constrained');x=np.arange(3)
    for i,name in enumerate(['intact_0','intact_1','intact_2']):axes[0].bar(x+(i-1)*.24,[next(r['fraction'] for r in b['cases'] if r['case']==name and r['regime']==t) for t in(.1,.5,1.)],width=.24,label=['a0.8/s0.85','a1/s1.1','a1.2/s1.6'][i])
    axes[0].axhline(.9,color='black',ls=':');axes[0].set_ylabel('Intact acceptance fraction');axes[0].legend(fontsize=8)
    for i,name in enumerate(['missing0','missing1','split']):axes[1].bar(x+(i-1)*.24,[next(r['fraction'] for r in b['cases'] if r['case']==name and r['regime']==t) for t in(.1,.5,1.)],width=.24,label=['X removed','Y removed','Double core'][i])
    axes[1].axhline(.05,color='black',ls=':');axes[1].set_ylabel('Damaged-structure acceptance fraction');axes[1].legend(fontsize=8)
    for ax in axes:ax.set_xticks(x,['0.1','0.5','1']);ax.set_xlabel('Conditional noise regime model time');ax.grid(axis='y',alpha=.2)
    fig.savefig(OUT/'MG167_validation.png',dpi=160);plt.close(fig)
    lines=['# MG167 - X/Y core detector development','',
    'The narrower X/Y detector rejects missing components more reliably, but it still cannot serve as a validated persistence measurement across the tested noise regimes. It misses many weaker intact cores under later noise and accepts too many double cores. The key problem therefore extends beyond the weak G feature: detecting and identifying a single coupled X/Y core against an evolving background remains difficult. The frozen validation failure blocked application to MG163. G remains explicitly unresolved; no complete-particle persistence claim or revised earlier score is issued.','',
    '![Validation fractions](MG167_validation.png)','',
    '## Candidate and independent test','',
    'The X/Y instrument uses a common position, dilation and positive amplitude, two background offsets, and the2x2X/Ycovariance. Each species must also retain at least half its own initial fitted amplitude, so one cannot disappear behind the other. A new competing-core check subtracts the fitted primary template and searches the remainder beyond2R95; its threshold is frozen at the largest residual second-core gain in256correct-single-core calibration examples. Existing residual/component-consistency guards remain. These are declared operational choices, not physical laws.','',
    'Source, specification, seed16700and input hashes were published prospectively atbe38180before evaluation. Synthetic draws are new and independent of the MG165realizations; the estimated noise regimes still come from previously seen MG163control0/1and therefore are not independently validated real dynamics. Three regimes T=.1,.5,1 use256null/256intact calibration fields each; validation uses1024null and128perintact/damaged case perregime. Full fields stay local. No new Model G trajectories were run.','',
    '## Frozen development screen: FAILED','',
    'Acceptance requirements: intact>=90%, damaged<=5%,null<=1%. These are raw finite-screen criteria. Individual95%binomial intervals are retained descriptively in benchmark.json; no familywise or population reliability claim. Thresholds were not altered after outcomes.','',
    '| Case | T0.1 accepted | T0.5 accepted | T1 accepted |','|---|---:|---:|---:|']
    names={'null':'No core','intact_0':'Intact a0.8/s0.85','intact_1':'Intact a1/s1.1','intact_2':'Intact a1.2/s1.6','missing0':'X removed','missing1':'Y removed','shift0':'X displaced by3.2','shift1':'Y displaced by3.2','split':'Two separated cores'}
    for name,label in names.items():
        rr=[next(r for r in b['cases'] if r['case']==name and r['regime']==t) for t in(.1,.5,1.)];lines.append('| '+label+' | '+' | '.join(f"{r['compatible']}/{r['n']}" for r in rr)+' |')
    lines+=['',
    'Missing X was rejected in every validation example; missing Y was accepted in only1-4of128examples perregime, within the finite5%screen. This is a useful improvement in component protection, not sufficient overall qualification. The weakest intact case failed even the early screen (114/128=89.0625%,below90%);13ofthose early failures came from the individual-component floor. AtT=.5and1,109/128and120/128weak-intact examples were rejected because total primary evidence did not exceed the frozen detection criterion (gain/amplitude guard), not because G was included. Widening acceptance would need an independently validated tradeoff, not an after-the-fact relaxation.','',
    'The residual second-core check rejects123/128early split examples directly; all128early split examples are rejected by some guard. Later,24/128and23/128double-core examples remain accepted, exceeding5%. Mid-regime displaced-X/Y examples also fail the damage screen at12/128each. Guard order means a rejection label records the first failed test, not every reason it might fail.','',
    '## Saved trajectories and interpretation','',
    'The prospectively declared dependency gate blocked saved-MG163classification because at least one validation screen failed. No new compatible-frame counts or persistence verdict were generated. MG163remains NOT ESTABLISHED under its original X score; MG165remains unqualified; G is UNRESOLVED. This failure concerns the measurement candidate and does not demonstrate disappearance of the physical or simulated particle.','',
    '## Verification and resources','',
    f"Separate weighted normal equations reproduced {a['independent_GLS_checks']}selected primary fits/per-species residuals. Independent direct spatial sums reproduced primary maxima and the residual secondary search for {spatial}validation fields across all23scales/500positions. All {checks}validation labels, calibration maxima and27case verdicts were recomputed. Frozen source bindings were verified. Same-author checks share inputs/templates; not external replication.",'',
    f"Evaluation took {a['elapsed_seconds']:.6g}seconds on one worker/library-thread1, below600seconds. Full synthetic payload {sum(v['bytes'] for v in manifest.values())/1e6:.6g}MB remains in local_raw/MG167_XYwith a verified SHA-256manifest; no immutable raw Git commit is claimed. Public fit/label exports are losslessly compressed; no full synthetic fields or event trajectories are published.",'',
    '## Recommended next step','',
    'Before modifying another classifier, quantify the best-case known-position/width X/Y discrimination against damaged and multiple-core alternatives. This separates lack of signal from the penalty of searching the full ring. If informative evidence exists at known geometry, prepare a separately registered temporal tracker with explicit motion/deformation assumptions and control false-track tests. If it does not, temporal tracking cannot simply be presumed to rescue the measurement. Do not change thresholds or launch another long persistence run on the strength of this failed candidate. No further test is authorized by this report.','']
    text='\n'.join(lines)
    for old,new in [('the2x2X/Ycovariance','the 2x2 X/Y covariance'),('beyond2R95','beyond 2R95'),('in256correct','in 256 correct'),('seed16700and','seed 16700 and'),('atbe38180before','at be38180 before'),('MG165realizations','MG165 realizations'),('MG163control0/1and','MG163 control0/1 and'),('use256null/256intact','use 256 null/256 intact'),('uses1024null and128perintact/damaged','uses 1024 null and 128 per intact/damaged'),('Individual95%binomial','Individual 95% binomial'),('only1-4of128examples perregime','only 1-4 of 128 examples per regime'),('finite5%screen','finite 5% screen'),('below90%','below 90%'),('13ofthose','13 of those'),('AtT=.5and1,109/128and120/128weak-intact','At T=0.5 and 1, 109/128 and 120/128 weak-intact'),('rejects123/128early','rejects 123/128 early'),('all128early','all 128 early'),('Later,24/128and23/128double-core','Later, 24/128 and 23/128 double-core'),('exceeding5%','exceeding 5%'),('at12/128each','at 12/128 each'),('MG163classification','MG163 classification'),('MG163remains','MG163 remains'),('MG165remains','MG165 remains'),('}selected','} selected'),('}validation','} validation'),('all23scales/500positions','all 23 scales/500 positions'),('and27case','and 27 case'),('}seconds','} seconds'),('below600seconds','below 600 seconds'),('}MB','} MB'),('MG167_XYwith','MG167_XY with'),('SHA-256manifest','SHA-256 manifest')]:text=text.replace(old,new)
    (OUT/'MG167_XY_core_results.md').write_text(text,encoding='utf-8');atomic_json(OUT/'report_attempt.json',dict(status='COMPLETE',elapsed_seconds=time.perf_counter()-start))
if __name__=='__main__':report()

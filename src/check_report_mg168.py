"""Saved-only oracle decision and independent secondary-search verification."""
import time,json,gzip,hashlib
import numpy as np
from analyze_mg168 import OUT,PRIOR,RAW,OLD,ROOT,read,atomic_json,read_gzip,binding
from mg167_xy import XYCore

def report():
    start=time.perf_counter();s=read(OUT/'analysis.json');b=read(PRIOR/'benchmark.json');ncheck=0;spatial=0;verdicts=[]
    with np.load(OLD/'prepared_inputs.npz') as z:ref=z['reference'];initial=z['initial'];omega=float(z['omega'])
    for ri,regime in enumerate(b['regime_inputs']):
        model=XYCore(ref,initial,omega,np.array(regime['covariance'])[1:,1:])
        for mode in('nearest','exact'):
            cal=read(OUT/f'r{ri}_{mode}_calibration.json')
            for row in [r for r in s['rows'] if r['regime']==regime['time'] and r['mode']==mode]:
                case=row['case'];d=read(OUT/f'r{ri}_{mode}_{case}.json');now=[];flags=[]
                for f in d['fits']:
                    g=dict(detection=f['score']<=cal['score'],common_amplitude=f['amplitude']/model.base<.5,scale=f['bank_edge'] or not .5<=f['scale']<=2,component_amplitude=min(np.array(f['individual_amplitudes'])/model.initial_components)<.5,component_consistency=f['component_stat']>cal['component'],secondary=f['secondary_score']>cal['secondary'],residual=f['residual']>cal['residual']);flags.append(g);now.append(not any(g.values()));ncheck+=1
                assert flags==d['flags'] and now==d['accepted'];assert sum(now)==row['accepted'];positive=case.startswith('intact');passed=np.mean(now)>=.9 if positive else np.mean(now)<=(.01 if case=='null' else .05);assert bool(passed)==row['screen_passed'];verdicts.append(dict(mode=mode,passed=bool(passed)))
                base=np.array([v=='candidate_compatible' for v in read_gzip(PRIOR/f'r{ri}_{case}_fits.json.gz')['labels']]);accepted=np.array(now);assert int((~base&accepted).sum())==row['baseline_rejected_oracle_accepted'];assert int((base&~accepted).sum())==row['baseline_accepted_oracle_rejected']
                with np.load(RAW/f'r{ri}_{case}.npz') as z:field=z['fields'][0];field=field-field.mean(axis=1)[:,None]
                f=d['fits'][0];m=model.x;k=f['center'];w=np.roll(m.window,k);template=model.signal(f['scale'],k);A=np.zeros((2,500,3));A[:,:,0]=template
                for sp in range(2):A[sp,:,sp+1]=1
                normal=np.einsum('ina,ij,jnb,n->ab',A,model.inv,A,w);rhs=np.einsum('ina,ij,jn,n->a',A,model.inv,field,w);coef=np.linalg.solve(normal,rhs);err=field-np.einsum('ina,a->in',A,coef);np.testing.assert_allclose(coef[0],f['amplitude'],rtol=1e-9,atol=1e-9);np.testing.assert_allclose(np.einsum('in,ij,jn,n->',err,model.inv,err,w),f['residual'],rtol=1e-9,atol=1e-9)
                remaining=field-coef[0]*template;rolls=np.array([np.roll(remaining,-j,axis=1).ravel() for j in range(500)]);dot=model.k.reshape(23,-1)@rolls.T;gain=np.maximum(dot,0)**2/model.norm[:,None];dist=np.minimum((np.arange(500)-k)%500,(k-np.arange(500))%500)*m.h;np.testing.assert_allclose(gain[:,dist>m.background_exclusion].max(),f['secondary_score'],rtol=1e-9,atol=1e-9);spatial+=1
    for mode in('nearest','exact'):assert all(v['passed'] for v in verdicts if v['mode']==mode)==s['mode_screens'][mode]
    assert binding()==read(OUT/'plan.json')['source'];atomic_json(OUT/'separate_check.json',dict(passed=True,guard_label_checks=ncheck,independent_primary_and_spatial_secondary_checks=spatial,weighted_normal_equation_checks=s['weighted_normal_equation_checks'],paired_decisions_and_verdicts_recomputed=True,limits='Same author/shared input fields/template generation; no external replication.'))
    for p in OUT.glob('r*.json'):
        if p.name.endswith('_calibration.json'):continue
        q=p.with_suffix('.json.gz');q.write_bytes(gzip.compress(p.read_bytes(),mtime=0));assert gzip.decompress(q.read_bytes())==p.read_bytes()
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(1,2,figsize=(11,4.5),layout='constrained')
    for ax,case,title in zip(axes,['intact_0','intact_1'],['Weaker intact core (a0.8 / s0.85)','Ordinary intact core (a1 / s1.1)']):
        x=np.arange(3);rows=[next(r for r in s['rows'] if r['mode']=='exact' and r['case']==case and r['regime']==t) for t in(.1,.5,1.)];near=[next(r for r in s['rows'] if r['mode']=='nearest' and r['case']==case and r['regime']==t) for t in(.1,.5,1.)]
        ax.bar(x-.24,[r['baseline_accepted']/r['n'] for r in rows],width=.24,label='Free search (MG167)');ax.bar(x,[r['fraction'] for r in near],width=.24,label='Known center, bank width');ax.bar(x+.24,[r['fraction'] for r in rows],width=.24,label='Known center, exact width');ax.axhline(.9,color='black',ls=':');ax.set_ylim(0,1.05);ax.set_xticks(x,['0.1','0.5','1']);ax.set_xlabel('Conditional noise regime model time');ax.set_ylabel('Intact acceptance fraction');ax.set_title(title);ax.legend(fontsize=7);ax.grid(axis='y',alpha=.2)
    fig.savefig(OUT/'MG168_oracle_comparison.png',dpi=160);plt.close(fig)
    lines=['# MG168 - known-position/width X/Y discrimination','',
    'Supplying the true core position and width substantially improves X/Y detection, showing that searching the ring is a major part of the measurement difficulty. Ordinary-strength intact cores pass the finite recall screen in every tested noise regime under exact geometry. Weaker intact cores improve substantially but still fail under middle and late noise, and displaced-Y/double-core exclusions remain inadequate there. Search uncertainty is therefore important but not the whole problem. These oracle results motivate a carefully bounded tracking feasibility study, not a persistence claim: a real tracker must estimate geometry and cannot be handed the synthetic truth. G remains unresolved and all prior verdicts are unchanged.','',
    '![Oracle comparison](MG168_oracle_comparison.png)','',
    '## What was held fixed and changed','',
    'The same saved MG167synthetic fields are used in all three procedures. Baseline is the original free23scale/500position search. The two new procedures are supplied the original injected center and either the nearest width on the frozen scale bank or the exact injected width. They retain the same X/Ycovariance, window, common amplitude, two background offsets, component floors/consistency checks, residual check and residual second-core search. Each known-geometry procedure calibrates its own thresholds from the existing calibration fields before validation; comparison therefore includes the reduced search penalty in calibration, not solely changed fitted coordinates. No threshold was tuned to improve validation counts.','',
    'This is a descriptive procedure comparison on already-seen conditional Gaussian fields, not fresh real-noise validation or a rescore of MG167. Null fields get prescribed coordinates independent of their noise realization; those coordinates would not be available for free in a real detector. Damaged cases retain the intended original coupled geometry, rather than separately realigning species. No new random draws or Model G trajectories.','',
    '## Intact-core acceptance','',
    '| Regime T | Intact case | Free search | Known center / bank width | Known center / exact width |','|---|---|---:|---:|---:|']
    for t in(.1,.5,1.):
        for case,label in [('intact_0','a0.8/s0.85'),('intact_1','a1/s1.1'),('intact_2','a1.2/s1.6')]:
            a=next(r for r in s['rows'] if r['regime']==t and r['mode']=='nearest' and r['case']==case);e=next(r for r in s['rows'] if r['regime']==t and r['mode']=='exact' and r['case']==case);lines.append(f"| {t:g} | {label} | {a['baseline_accepted']}/128 | {a['accepted']}/128 | {e['accepted']}/128 |")
    lines+=['',
    'The raw intact-recall requirement is>=90%(at least116of128). Exact geometry lifts the ordinary case atT=1from59/128to120/128and the weaker case from8/128to64/128. Most improvement is already present when the nearest bank width is supplied; correcting the small width-grid mismatch adds relatively little. This comparison does not separately isolate unknown-center versus unknown-width search penalties.','',
    '## Damage rejection and remaining limits','',
    '| Regime T | No core accepted | X removed accepted | Y removed accepted | X displaced accepted | Y displaced accepted | Double cores accepted |','|---|---:|---:|---:|---:|---:|---:|']
    for t in(.1,.5,1.):
        rr=[next(r for r in s['rows'] if r['mode']=='exact' and r['regime']==t and r['case']==c) for c in('null','missing0','missing1','shift0','shift1','split')];lines.append('| '+f'{t:g}'+' | '+' | '.join(f"{r['accepted']}/{r['n']}" for r in rr)+' |')
    lines+=['',
    'Exact geometry passes all nine finite case screens under early noise. Both known-geometry modes fail the combined screen across all regimes. Under later noise, the exact procedure still accepts21/128and37/128double cores, and12/128and13/128displaced-Y cases, above the damaged acceptance limit5%. Knowing the primary core position does not make a weak second core or a damaged component reliably identifiable by these checks.','',
    'For the weaker intact exact-geometry case atT=.5,20/128fail the evidence-gain threshold and21/128fail the component amplitude floor. AtT=1,62/128fail evidence gain and27/128fail the component floor. These failures overlap; full nonexclusive guard counts are retained in analysis.json and compressed case exports. This is not a proof of universal information impossibility: it identifies limitations of this detector, thresholds and conditional noise model. The individual-component floor remains an engineering safeguard, not a physical law.','',
    '## Verification, provenance and resources','',
    f"Primary weighted normal equations independently agreed on {s['weighted_normal_equation_checks']}selected fits during analysis. Separate primary normal equations plus full direct spatial residual secondary searches agreed for {spatial}validation examples. All {ncheck}guard vectors/accepted labels,54case summaries and paired baseline changes were recomputed. Source/input bindings are frozen. Shared template generation/fields and same-author work remain dependencies; no external replication.",'',
    'The source/specification and inputs were committed before evaluation and the registration upload was initiated before calculation. Upload completion was confirmed afterward rather than awaited before calculation, a publication-verification sequencing deviation. No source, criteria or inputs changed after outcomes. Prospective local registration and the completed public registration are retained in Git.','',
    f"Analysis took {s['elapsed_seconds']:.6g}seconds on one worker/library-thread1, below600active seconds. Full input fields remain local_raw/MG167_XY, verified against the registered SHA-256input manifest. Public exports contain compressed fit/guard decisions and summaries, not full trajectory or synthetic-field arrays. The supplied truth geometry is available only for these planted tests; MG163was not classified with it and no new persistence verdict was generated.",'',
    '## Recommended next step','',
    'Prepare a small tracking-feasibility study for ordinary-strength X/Y cores. Specify allowable displacement/width changes and predictive uncertainty before testing; compare it with the free detector on intact moving/deforming cores, missing/displaced species, competing cores, and no-particle controls using independent validation noise. Retain weak-core states as explicitly unresolved rather than forcing a pass/fail particle story. The oracle gap is a reason to test tracking, not evidence that tracking will work. Do not launch another long persistence simulation until that feasibility question is resolved. No tracker or fresh kinetic production is authorized by this report.','']
    text='\n'.join(lines)
    for a,c in [('MG167synthetic','MG167 synthetic'),('free23scale/500position','free 23-scale/500-position'),('X/Ycovariance','X/Y covariance'),('is>=90%(at least116of128)','is >=90% (at least 116 of 128)'),('atT=1from59/128to120/128and','at T=1 from 59/128 to 120/128 and'),('from8/128to64/128','from 8/128 to 64/128'),('accepts21/128and37/128double','accepts 21/128 and 37/128 double'),('and12/128and13/128displaced-Y','and 12/128 and 13/128 displaced-Y'),('limit5%','limit 5%'),('atT=.5,20/128fail','at T=0.5, 20/128 fail'),('and21/128fail','and 21/128 fail'),('AtT=1,62/128fail','At T=1, 62/128 fail'),('and27/128fail','and 27/128 fail'),('}selected','} selected'),('}validation','} validation'),('}guard','} guard'),('labels,54case','labels, 54 case'),('}seconds','} seconds'),('below600active','below 600 active'),('SHA-256input','SHA-256 input'),('MG163was','MG163 was')]:text=text.replace(a,c)
    (OUT/'MG168_known_geometry_results.md').write_text(text,encoding='utf-8');atomic_json(OUT/'report_attempt.json',dict(status='COMPLETE',elapsed_seconds=time.perf_counter()-start))
if __name__=='__main__':report()

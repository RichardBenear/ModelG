"""Independent tracker state-machine/guard/metric replay and feasibility report."""
import time,json,gzip,hashlib,collections
import numpy as np
from run_mg169_tracking import OUT,RAW,OLD,PRIOR,ROOT,read,atomic_json,binding
from mg167_xy import XYCore

def report():
    start=time.perf_counter();summary=read(OUT/'analysis.json');old=read(PRIOR/'benchmark.json');checks=0;primary=0;pathschecked=0;guardchecks=0;additional=[]
    with np.load(OLD/'prepared_inputs.npz') as z:ref=z['reference'];initial=z['initial'];omega=float(z['omega'])
    for ri,when in enumerate((.1,.5,1.)):
        model=XYCore(ref,initial,omega,np.array(old['regime_inputs'][ri]['covariance'])[1:,1:]);calrecord=read(OUT/f'r{ri}_calibration.json');good=calrecord['shape_fits'];null=calrecord['null_fits'];expected=dict(score=max(f['score'] for f in null),residual=max(f['residual'] for f in good),component=max(f['component_stat'] for f in good),secondary=max(f['secondary_score'] for f in good));assert expected==calrecord['thresholds']
        for row in [r for r in summary['rows'] if r['regime']==when]:
            tag=f"r{ri}_rho{row['rho']:g}_{row['case']}";d=read(OUT/(tag+'.json'));case=row['case'];positive=case in('static','moving','deforming');onset=0 if case=='absent' else 10
            with np.load(RAW/(tag+'.npz')) as z:firstfield=z['fields'][0,0]
            extra=dict(regime=when,rho=row['rho'],case=case)
            for algorithm in('local','free'):
                counts=[];false=[];early=0;aliveonset=0;failure_reasons=collections.Counter()
                for path in d['paths']:
                    alive=True;gap=0;center=path['initial_seed'];scale=1.;stop=None;acceptedcount=0;accurate=0
                    for t,(rec,tr) in enumerate(zip(path['records'][algorithm],path['truth'])):
                        original=(path['initial_seed']+(0 if case=='static' else t))%500;expectedscale=1.
                        if case=='deforming':original=(path['initial_seed']+t//2)%500;expectedscale=[1.,1.125,1.25,1.125][t//10]
                        if case=='jump' and t>=10:original=(original+20)%500
                        assert tr==[original,expectedscale]
                        if not alive:assert not rec['active'] and rec['center'] is None and not rec['accepted'];checks+=1;continue
                        assert rec['active'];f=rec['fit'];c=calrecord['thresholds'] if algorithm=='local' else old['calibrations'][ri]
                        flags=dict(detection=f['score']<=c['score'],common_amplitude=f['amplitude']/model.base<.5,scale=f['bank_edge'] or not .5<=f['scale']<=2,component_amplitude=min(np.array(f['individual_amplitudes'])/model.initial_components)<.5,component_consistency=f['component_stat']>c['component'],secondary=f['secondary_score']>c['secondary'],residual=f['residual']>c['residual']);accepted=not any(flags.values());assert accepted==rec['accepted'];guardchecks+=1
                        if algorithm=='local':
                            assert flags==rec['flags'];assert min((f['center']-center)%500,(center-f['center'])%500)<=2;assert abs(f['scale']-scale)<=.12500001 and .5<=f['scale']<=2
                        if accepted:center=f['center'];scale=f['scale'];gap=0;acceptedcount+=1
                        else:gap+=1
                        if gap==3:alive=False;stop=t
                        assert rec['center']==(center if alive else None) and rec['scale']==(scale if alive else None)
                        if alive and min((center-tr[0])%500,(tr[0]-center)%500)<=2 and abs(scale-tr[1])<=.12500001:accurate+=1
                        checks+=1
                        if t==0 and path['replicate']==0:
                            m=model.x;w=np.roll(m.window,f['center']);template=model.signal(f['scale'],f['center']);A=np.zeros((2,500,3));A[:,:,0]=template
                            for sp in range(2):A[sp,:,sp+1]=1
                            normal=np.einsum('ina,ij,jnb,n->ab',A,model.inv,A,w);rhs=np.einsum('ina,ij,jn,n->a',A,model.inv,firstfield,w);coef=np.linalg.solve(normal,rhs);np.testing.assert_allclose(coef[0],f['amplitude'],rtol=1e-9,atol=1e-9);primary+=1
                    o=path['outcomes'][algorithm];assert acceptedcount==o['compatible_frames'] and accurate==o['accurate_frames'] and stop==o['stop'] and alive==o['alive_final'];continued=stop is None or stop>onset+5;success=alive and acceptedcount>=36 and accurate>=36 if positive else not continued;assert bool(success)==o['success'];counts.append(success)
                    if not positive:
                        assert continued==o['false_continuation'];false.append(continued);early+=stop is not None and stop<onset;aliveonset+=stop is None or stop>=onset
                    elif not success:
                        failure_reasons['stopped' if not alive else 'too_few_compatible' if acceptedcount<36 else 'localization_width_error']+=1
                    pathschecked+=1
                assert sum(counts)==row[algorithm]['successful_paths'];assert row[algorithm]['screen_passed']==(sum(counts)>=15 if positive else sum(false)==0)
                if not positive:assert sum(false)==row[algorithm]['false_continuations']
                extra[algorithm]=dict(stopped_before_damage=int(early) if not positive else None,alive_at_damage_onset=int(aliveonset) if not positive else None,positive_failure_reasons=dict(failure_reasons))
            additional.append(extra)
    assert all(r['local']['screen_passed'] for r in summary['rows'])==summary['local_feasibility_screen'];assert binding()==read(OUT/'plan.json')['source'];atomic_json(OUT/'separate_check.json',dict(passed=True,state_records_checked=checks,guard_decisions_checked=guardchecks,path_outcomes_checked=pathschecked,independent_primary_fits=primary,additional_outcomes=additional,limits='Same author/shared fields and templates; independent state replay/branch aggregation, not external replication.'))
    for p in OUT.glob('r*.json'):
        if 'rho' not in p.name and not p.name.endswith('_calibration.json'):continue
        q=p.with_suffix('.json.gz');q.write_bytes(gzip.compress(p.read_bytes(),mtime=0));assert gzip.decompress(q.read_bytes())==p.read_bytes()
    manifest={p.name:dict(bytes=p.stat().st_size,sha256=hashlib.sha256(p.read_bytes()).hexdigest()) for p in RAW.iterdir() if p.is_file()};atomic_json(OUT/'local_raw_manifest.json',dict(directory=str(RAW),files=manifest,total_bytes=sum(v['bytes'] for v in manifest.values()),full_fields_published=False))
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(1,2,figsize=(11,4.5),layout='constrained')
    for ax,rho in zip(axes,(0.,.95)):
        for name,color in [('local','tab:blue'),('free','tab:orange')]:
            values=[sum(r[name]['successful_paths'] for r in summary['rows'] if r['regime']==t and r['rho']==rho and r['case'] in('static','moving','deforming'))/48 for t in(.1,.5,1.)];ax.plot([.1,.5,1.],values,'o-',color=color,label=name)
        ax.axhline(.9,color='black',ls=':');ax.set_ylim(0,1.05);ax.set_title(f'Synthetic temporal correlation rho={rho:g}');ax.set_xlabel('Noise regime model time');ax.set_ylabel('Intact successful paths (pooled descriptive fraction)');ax.legend();ax.grid(alpha=.2)
    fig.savefig(OUT/'MG169_tracking_results.png',dpi=160);plt.close(fig)
    lines=['# MG169 - X/Y tracking feasibility','',
    'A bounded tracker can recover ordinary-strength X/Y cores much more effectively than the free-ring detector in several tested conditions, confirming that motion constraints can reduce the search problem. It succeeds on every early uncorrelated intact path and rejects all early uncorrelated damage/no-core cases. However, it fails the full feasibility screen: temporally correlated noise loses intact tracks, while later displaced-Y and competing-core cases often continue despite ambiguity. Tracking is therefore a promising direction with demonstrated partial gains, not a qualified persistence instrument. G remains unresolved, MG163 is not rescored, and no new Model G simulations were run.','',
    '![Tracking results](MG169_tracking_results.png)','',
    '## Frozen tracker and denominator','',
    'Every path receives the same kind of initial position/width seed, including no-core controls. No true geometry is supplied afterward. The local tracker searches within2sites and one0.125width increment of the last accepted state, width bounds0.5-2. A failed observation holds state; three consecutive failures stop permanently, without reacquisition or expanding uncertainty. These are hard motion assumptions, not statistical confidence bounds. The jump control intentionally violates the permitted motion.','',
    'The new local per-frame thresholds are calibration maxima from256null and256single-core examples in each old noise regime. The free comparator uses the old MG167free detector/thresholds with the same stop rule. Thus this compares whole procedures, not just the displacement bound. The full3x3cross-spectrum generator supplies only the X/Ychannels; G is not measured. Noise regimes0.1/0.5/1 label fixed estimated background conditions, not a simulated evolution from one regime to the next.','',
    'Tests use40observations/path, spacing0.001modeltime,16independent synthetic noise sequences percase in each of3regimes and2temporal conditions (rho0or0.95). Temporal AR1is prescribed, not measured from Model G; correlations at fixed sites interact with motion. Positives are amplitude1static, moving1site/observation, or moving/deforming with width1->1.125->1.25->1.125. Negative cases are absent from start, or vanish/removeX/removeY/displaceY/addsecondcore/jump20sites atframe10. Seed16900/source/input hashes were published and push completion confirmed before draws.','',
    '## Intact-path success','',
    'A successful intact path must stay alive through all40observations, have>=36compatible frames, and have>=36live estimates within2sites and0.125width of injected truth. Gaps are not compatible observations. Require at least15/16successful paths percase; not a population reliability bound.','',
    '| Noise regime | rho | Static local/free | Moving local/free | Deforming local/free |','|---|---:|---:|---:|---:|']
    for when in(.1,.5,1.):
        for rho in(0.,.95):
            cells=[]
            for case in('static','moving','deforming'):
                r=next(r for r in summary['rows'] if r['regime']==when and r['rho']==rho and r['case']==case);cells.append(f"{r['local']['successful_paths']}/16 / {r['free']['successful_paths']}/16")
            lines.append(f'| {when:g} | {rho:g} | '+' | '.join(cells)+' |')
    lines+=['',
    'At late uncorrelated noise, local tracking succeeds on15/16static,14/16moving and16/16deforming paths, compared with0/16for each free procedure. The moving case still misses the15/16criterion. Correlated-noise performance is weaker and inconsistent across regimes, so the success in independent-noise conditions cannot be generalized. Pooled fractions in the figure are descriptive only; the frozen gate applies separately to each case.','',
    '## False continuation after damage or absent cores','',
    'The tracker must stop by onset+5observations (by frame5for absent controls,15for other negatives). False continuation requires0/16percase to meet the raw<=5%screen;16paths cannot certify a5%populationfalse-trackrate. Stopping before damage also satisfies this endpoint but is not successful damage recognition. Pre-damage losses and the number alive at damage onset are therefore retained in separate_check.json.','',
    '| Regime | rho | Absent | Vanish | X removed | Y removed | Y displaced | Second core | Jump |','|---|---:|---:|---:|---:|---:|---:|---:|---:|']
    for when in(.1,.5,1.):
        for rho in(0.,.95):
            vals=[next(r['local']['false_continuations'] for r in summary['rows'] if r['regime']==when and r['rho']==rho and r['case']==c) for c in('absent','vanish','missingX','missingY','displacedY','split','jump')];lines.append(f'| {when:g} | {rho:g} | '+' | '.join(f'{v}/16' for v in vals)+' |')
    lines+=['',
    'Late split cases are the strongest warning: local tracking continues in12/16uncorrelated and11/16correlated sequences after the second core appears. It follows the initialized object but does not reliably report the scene as ambiguous. Displaced-Y continuation occurs in5/16and4/16respectively. Conversely, the local jump controls stop in every test; the free detector often follows the jump under early/middle noise, consistent with its lack of a motion-bound identity rule. A stopped or lost free track before damage is not evidence of superior damage recognition. All stop/guard/state histories remain in lossless compact exports.','',
    'Across all 144 intact paths per temporal condition, the local tracker has 16 failures under independent noise (5 stopped, 11 with too few compatible frames). Under rho=0.95 it has 53 failures: 45 stopped, 2 had too few compatible frames, and 6 failed localization/width accuracy. These are post-result diagnostic categories, not changed criteria. Late split controls also show why baseline survival matters: the local tracker was alive at damage onset in 15/16 paths in each noise condition; the free tracker was alive in only 3/16 and 2/16. Its smaller false-continuation count therefore cannot be read as superior recognition of the added core.','', '## Verdict and limits','',
    '**Overall feasibility screen: FAILED for both procedures.** The local tracker passes all ten finite cases only in the early uncorrelated setting. No tracking criterion was tuned after outcomes. New random paths are independent conditional synthetic realizations; the estimated spectral/covariance inputs are previously seen, and this is not real-trajectory validation. Synthetic motion, deformation, noise correlation and stop rules need qualification against actual model dynamics before any persistence interpretation.','',
    '## Checks, compute and retention','',
    f"The primary evaluation checked {summary['independent_primary_checks']}selected fit amplitudes through direct weighted normal equations. Separate replay verified {checks}state records, {guardchecks}guard decisions and {pathschecked}path outcomes, including motion bounds, gap stopping, pre-damage losses and source binding. It independently checked {primary}additional primary fits. Shared templates/fields and same-author work remain dependencies; not external replication.",'',
    f"Synthetic evaluation took {summary['elapsed_seconds']:.6g}seconds on one worker/library-thread1 within600active seconds. Full sequence payload {sum(v['bytes'] for v in manifest.values())/1e6:.6g}MB remains in local_raw/MG169_tracking, below the0.5GBreserve. A local SHA-256manifest verifies retained physical files; no immutable raw Git archive is claimed. Public data are compact fit/state/truth metadata, summaries, plots and hashes; full fields remain local. G is UNRESOLVED; no new particle-persistence verdict or kinetic budget.",'',
    '## Recommended next step','',
    'Use the saved failures to separate frame-evidence losses from localization drift and unresolved competing cores before changing the tracker. A next prototype should represent uncertainty and competing hypotheses explicitly, with a separately validated ambiguity state, and calibrate sequential evidence under temporally correlated noise rather than merely allow more gaps. Preserve missing-component and no-core controls and make pre-damage tracking success explicit in damage tests. Freeze any revised rules before independent validation; another long Model G run is not justified by this failed feasibility screen. No revised tracker or production is authorized by this report.','']
    text='\n'.join(lines)
    replacements={'within2sites':'within 2 sites','one0.125width':'one 0.125 width','bounds0.5-2':'bounds 0.5-2','from256null and256single-core':'from 256 null and 256 single-core','MG167free':'MG167 free','full3x3cross-spectrum':'full 3x3 cross-spectrum','X/Ychannels':'X/Y channels','regimes0.1/0.5/1':'regimes 0.1/0.5/1','use40observations':'use 40 observations','spacing0.001modeltime,16independent':'spacing 0.001 model time, 16 independent','percase in each of3regimes and2temporal':'per case in each of 3 regimes and 2 temporal','rho0or0.95':'rho 0 or 0.95','AR1is':'AR1 is','amplitude1static':'amplitude 1 static','moving1site/observation':'moving 1 site/observation','width1->':'width 1->','jump20sites atframe10':'jump 20 sites at frame 10','Seed16900/source':'Seed 16900/source','all40observations':'all 40 observations','have>=36compatible':'have >=36 compatible','have>=36live':'have >=36 live','and0.125width':'and 0.125 width','least15/16successful':'least 15/16 successful','on15/16static,14/16moving and16/16deforming':'on 15/16 static, 14/16 moving and 16/16 deforming','with0/16for':'with 0/16 for','the15/16criterion':'the 15/16 criterion','onset+5observations':'onset+5 observations','frame5for':'frame 5 for','controls,15for':'controls, 15 for','requires0/16percase':'requires 0/16 per case','raw<=5%screen;16paths':'raw <=5% screen; 16 paths','a5%populationfalse-trackrate':'a 5% population false-track rate','in12/16uncorrelated and11/16correlated':'in 12/16 uncorrelated and 11/16 correlated','in5/16and4/16respectively':'in 5/16 and 4/16 respectively','}selected':'} selected','}state':'} state','}guard':'} guard','}path':'} path','}additional':'} additional','}seconds':'} seconds','within600active':'within 600 active','}MB':'} MB','the0.5GBreserve':'the 0.5 GB reserve','SHA-256manifest':'SHA-256 manifest'}
    for old,new in replacements.items():text=text.replace(old,new)
    (OUT/'MG169_tracking_feasibility_results.md').write_text(text,encoding='utf-8');atomic_json(OUT/'report_attempt_v2.json',dict(status='COMPLETE',elapsed_seconds=time.perf_counter()-start))
if __name__=='__main__':report()

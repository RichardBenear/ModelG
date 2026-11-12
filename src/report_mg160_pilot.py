"""MG160-P1 report preserves negative per-path result and labels exploratory diagnostics."""
import os
for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):os.environ[k]='1'
import time,json,re
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from run_mg160_persistence import BASE,read,hashes
from measure import smooth_periodic
from mg10e2c1_pool import atomic_json
from mg160_local_archive import local_archive
BASE=local_archive(BASE)
def report():
    started=time.perf_counter();attempt=BASE/'P1_report_attempt.json';number=1
    while attempt.exists():number+=1;attempt=BASE/f'P1_report_attempt_{number:03d}.json'
    atomic_json(attempt,dict(status='RUNNING'))
    try:
        a=read(BASE/'P1_analysis_results.json');assert a['passed'] and a['source_hashes']==hashes();d=read(BASE/'P1_exploratory_diagnostics.json');ledger=read(BASE/'P1_production/ledger.json');run=read(BASE/'P1_attempt.json');plan=read(BASE/'P1_plan.json');assert run['status']=='COMPLETE';prepared=[j for j in a['jobs'] if j['condition']=='prepared'];nulls=[j for j in a['jobs'] if j['condition']=='homogeneous'];pt=[t for j in prepared for t in j['trajectories']];nt=[t for j in nulls for t in j['trajectories']];assert len(prepared)==len(nulls)==4
        retained=[sum(j['trajectories'][s]['pattern']['retained_at_all_frames'] for j in prepared) for s in range(2)];triggers=[sum(j['trajectories'][s]['pattern']['positive_frames']>0 for j in nulls) for s in range(2)];assert retained==[0,0] and triggers==[0,0]
        rows=[]
        for j in prepared:
            pair=[]
            for t in j['trajectories']:
                bracket=t['pattern']['first_loss_interval'];pair.append(f'({bracket[0]:.6g}, {bracket[1]:.6g}]' if bracket else 'none')
                assert bracket==t['shape']['first_loss_interval'];assert t['minimum_relative_amplitude']>=.5
            rows.append(f"| {j['id']} | {pair[0]} | {pair[1]} |")
        with np.load(BASE/'prepared_inputs.npz') as z:initial=z['initial'];omega=float(z['omega'])
        with np.load(BASE/'P1_exploratory_mean_fields.npz') as z:means=z['mean_counts']
        times=np.arange(101)*.001;fig=plt.figure(figsize=(14,8));grid=fig.add_gridspec(2,3,height_ratios=[1,1.1]);axes=[fig.add_subplot(grid[0,i]) for i in range(3)];profile=fig.add_subplot(grid[1,:]);colors=['#245a9c','#c66c22'];keys=['relative_amplitude','correlation','relative_width'];titles=['Fitted amplitude / initial','Template correlation','Effective extent / initial']
        for job in plan['jobs']:
            metrics=read(BASE/'P1_production'/job['id']/'metrics.json')['values'];is_prepared=job['condition']=='prepared'
            for side in (0,1):
                for ax,key in zip(axes,keys):ax.plot(times,[v[side][key] for v in metrics],color=colors[side] if is_prepared else '#888888',lw=1 if is_prepared else .7,alpha=.65 if is_prepared else .3,ls='-' if side==0 else '--')
        for side in (0,1):
            for ax,key in zip(axes,keys):ax.plot(times,[f[side][key] for f in d['conditions'][0]['metrics']],color='#111111',lw=2,ls='-.' if side==0 else ':')
        axes[0].axhline(.5,color='red',ls=':',lw=1);axes[0].set_ylim(0,1.3);axes[1].axhline(.8,color='red',ls=':',lw=1);axes[1].set_ylim(0,1.05);axes[2].set_yscale('log');axes[2].axhspan(.5,2,color='#dcefdc',alpha=.7);axes[2].axhline(2,color='red',ls=':',lw=1)
        for ax,title in zip(axes,titles):ax.set_title(title);ax.set_xlabel('Post-preparation model time');ax.grid(alpha=.18)
        original=smooth_periodic(initial[1]/omega-initial[1].mean()/omega,.2);anchor=int(np.argmax(abs(original)));shift=250-anchor;x=(np.arange(500)-250)*.2
        for job in [j for j in plan['jobs'] if j['condition']=='prepared']:
            with np.load(BASE/'P1_production'/job['id']/'observations.npz') as z:final=z['counts'][-1]
            for side in (0,1):f=final[side,1]/omega;f=smooth_periodic(f-f.mean(),.2);profile.plot(x,np.roll(f,shift),color=colors[side],alpha=.14,lw=.7)
        profile.plot(x,np.roll(original,shift),color='black',lw=2,label='Initial prepared profile')
        for side in (0,1):f=means[0,-1,side,1]/omega;f=smooth_periodic(f-f.mean(),.2);profile.plot(x,np.roll(f,shift),color=colors[side],lw=2,label=f'Final unaligned ensemble-mean field, tick {plan["deltas"][side]:g}')
        profile.set_xlim(-25,25);profile.set_xlabel('Model-lengths about fixed initial feature anchor');profile.set_ylabel('Centered, smoothed X');profile.set_title('Final individual fields (faint) and ensemble means — exploratory, not per-path rescoring');profile.grid(alpha=.15);profile.legend(fontsize=8)
        fig.legend(handles=[Line2D([],[],color=colors[0],label='Prepared, tick0.0001'),Line2D([],[],color=colors[1],ls='--',label='Prepared, tick0.00002'),Line2D([],[],color='#888888',label='Homogeneous controls'),Line2D([],[],color='black',ls='-.',label='Ensemble-mean field diagnostics')],loc='upper center',bbox_to_anchor=(.5,.96),ncol=4,fontsize=9);fig.suptitle('MG160-P1: amplitude remains, while individual shape/coherence criteria are lost',fontsize=14);fig.tight_layout(rect=(0,0,1,.90));fig.savefig(BASE/'MG160_P1_pattern_diagnostics.png',dpi=160);fig.savefig(BASE/'MG160_P1_pattern_diagnostics.svg');plt.close(fig)
        amp=[t['final']['relative_amplitude'] for t in pt];rho=[t['final']['correlation'] for t in pt];width=[t['final']['relative_width'] for t in pt];prepstocks=[s for s in d['stocks'] if s['id'].startswith('prepared')];maxa=max(s['A_depleted_fraction'] for s in prepstocks);maxb=max(s['B_depleted_fraction'] for s in prepstocks);active=read(BASE/'P1_analysis_attempt.json')['elapsed_seconds']+read(BASE/'P1_diagnostic_attempt.json')['elapsed_seconds'];events=sum(r['events'] for r in ledger['records']);upper0=1-.025**.25
        atomic_json(BASE/'P1_budget_summary.json',dict(worker_seconds=run['charged_worker_seconds'],worker_core_hours=run['charged_worker_seconds']/3600,worker_cap=100000,wall_seconds=run['wall_seconds'],active_audit_diagnostics_seconds=active,parent_cap=1800,report_seconds_so_far=time.perf_counter()-started,completed_pairs=8,marginal_chains=16,failures=0,retries=0,event_count=events));atomic_json(BASE/'P1_outcome_summary.json',dict(prepared_retained_pairs_per_tick=retained,homogeneous_ever_triggered_pairs_per_tick=triggers,independent_pairs_per_condition=4,first_failure_reason='Upper effective-extent bound in every prepared marginal trajectory; amplitude threshold never failed.',exploratory_binomial_two_sided95_upper_zero_of_four=upper0,score='No continuous sampled initial-pattern retention through0.1 in observed prepared trajectories; complete pilot, not population/physical particle-death proof.'))
        ensemble=d['conditions'][0]['final'];ensemble_rows=[f"| {plan['deltas'][s]:g} | {m['relative_amplitude']:.6g} | {m['correlation']:.6g} | {m['relative_width']:.6g} |" for s,m in enumerate(ensemble)]
        lines=['# MG160-P1 — first pattern-persistence pilot results','',
        'The pilot completed successfully as a run, but it did not demonstrate uninterrupted retention of the prepared pattern under the criteria approved in advance. Every prepared trajectory at both tick sizes crossed the shape-retention limit before model time0.1, while all homogeneous controls remained below the detection criteria. Importantly, the fitted pattern amplitude stayed near its starting value: the first failures came from the global effective-extent measurement, followed later by reduced template correlation. Averaging the four prepared fields still reveals a recognizable feature, suggesting that background noise is an important measurement issue rather than establishing physical disappearance. The next step should separate noise-sensitive measurement from genuine deformation before funding a much longer run; this pilot neither proves particle death nor demonstrates long-lived stochastic particles.',
        '','## Registered findings','',
        '**PILOT COMPLETE; NO CONTINUOUS SAMPLED INITIAL-PATTERN RETENTION THROUGH T=0.1 IN THE OBSERVED PREPARED TRAJECTORIES.** No failed/replaced/retried jobs. No scoring thresholds, observation times or denominators changed. PI closure is pending.',
        '','| Tick | Prepared pairs retaining the criterion at every frame | Homogeneous pairs ever triggering the criterion |','|---|---:|---:|','|0.0001|0 / 4|0 / 4|','|0.00002|0 / 4|0 / 4|','',
        'Each condition has four independent seed pairs; the two tick paths within a pair are correlated. They are not eight independent replications. Nor are the101 snapshots per path independent survival trials. Four nulls do not establish a low false-detection rate. For scale, a descriptive exact two-sided95% binomial interval after0 successes in4 independent trials has upper endpoint about'+f'{upper0:.6g}'+', assuming the corresponding Bernoulli model; the zero count alone is weak population-level evidence.',
        '','## First sampled losses','',
        'The initial-pattern rule requires amplitude>=0.5 of initial, correlation>=0.8, and effective extent within[0.5,2] times initial. The first crossing was the upper extent bound in every prepared marginal path. Each bracket ends at the first failed sampled frame and begins at the preceding observed pass. Unobserved earlier losses/recoveries remain possible; these are not measured continuous-time lifetimes.',
        '','| Prepared pair | First sampled loss bracket, tick0.0001 | First sampled loss bracket, tick0.00002 |','|---|---|---|',*rows,'',
        f'At the final frame the prepared-path amplitude ratios ranged from{min(amp):.6g} to{max(amp):.6g}, correlation from{min(rho):.6g} to{max(rho):.6g}, and extent ratios from{min(width):.6g} to{max(width):.6g}. Amplitude did not cross its failure threshold anywhere. Some paths subsequently crossed back into the combined criterion temporarily; their first loss is retained and is not erased by later recovery. All final prepared frames failed both the extent and correlation thresholds. Both ticks gave the same overall retention count, but their crossing times sometimes differed; no formal T=0.1 tick-equivalence qualification was registered here.',
        '','![Pattern metrics and exploratory field averages](MG160_P1_pattern_diagnostics.png)','',
        '## What the exploratory field average suggests','',
        'After preserving the negative per-path outcome, we averaged raw fields across the four prepared pairs separately at each tick, without phase-aligning individual paths. These are post-result diagnostics with a different estimand, not substitute survival outcomes. At T=0.1 the ensemble-mean fields give:',
        '','| Tick | Relative amplitude | Template correlation | Relative effective extent |','|---|---:|---:|---:|',*ensemble_rows,'',
        'Thus the prepared feature remains visible in the ensemble mean even though the individual-path criterion is lost. This is consistent with noise outside the core inflating a global participation-length statistic and reducing cosine similarity. It does not uniquely separate noise from genuine deformation, prove that every original feature survives, or justify rescoring the registered result. The homogeneous ensemble fields also remain below the detection rule.',
        '',f'Bulk prepared A/B depletion through this interval was small: at most{maxa*100:.6g}% of initial A and{maxb*100:.6g}% of initial B across the prepared paths. That makes wholesale stock exhaustion an unsupported explanation for these early crossings; local stock redistribution and stochastic dynamics are not ruled out. Complete resource accounting is retained.',
        '',"The earlier accuracy study tested ensemble-mean agreement between resolutions at 0.005/0.01. This pilot tests an individual-path rule through 0.1. Pairing improves the precision of a difference while leaving each individual path's stochastic fluctuations intact.",
        '','## Scope and reproducibility','',
        'The frozen pilot has four prepared and four mass-matched homogeneous pairs,16 marginal paths, Omega=10,500sites, fixed ring/environmental links and finite A/B with retained Z/W. The initial state is a rounded prepared donor from historical absolute time10. Post-preparation time0.1 corresponds to a short continuation, not a newly formed particle or a20-unit lifetime experiment. Every0.001 segment continues the RNG state and ends at common tick boundaries. Approval and source/seed registration were published at`d73f1ef` before sampling; no later stage is authorized.',
        '',f"All{a['checkpoint_chain_endpoints']} chain checkpoint endpoints passed histogram reconstruction, nonnegativity and conservation checks. Both full numeric sentinels replayed, covering{a['sentinel_events']/1e6:.6g}million joint events and every retained reservation/reset boundary. RNG continuity, seeds, source hashes and full observation records passed. Separate Fourier smoothing/cross-correlation, amplitude/extent calculation and independently implemented sampled-loss classifications matched every stored metric/outcome. A/B depletion and Z/W production matched reaction counts exactly.",
        '',
        'The compiled replay validates integer reservations, event ordering, enabled participants, boundary state and histograms; it does not resample random clocks. Positivity of channel rates follows here from the checked equality of reservation and falling-factorial requirements and positive coefficients. Shared frozen channel/initial definitions and same-author analysis are disclosed; this is not external replication. The registered per-path result and exploratory ensemble-field analysis are retained separately. Earlier closed accuracy results remain unchanged.',
        '','## Compute and retention','',
        f"All eight jobs completed in{run['wall_seconds']/60:.6g}minutes ({run['wall_seconds']/3600:.6g}hours) wall time, using{run['charged_worker_seconds']:.6g}charged worker-seconds ({run['charged_worker_seconds']/3600:.6g}core-hours) against the100000-second ceiling. At most eight workers ran with numerical-library threads1. Audit and exploratory diagnostic analysis used{active:.6g}measured seconds; the active parent/control/audit allowance is1800seconds. Every analysis attempt records elapsed time, including failures. Production processed{events/1e6:.6g}million joint events (a shared event counts once). Full raw data and sentinel histories are retained locally under the PI storage policy; see P1_local_archive.json. Public artifacts contain the report, figures, code and compact decisions only. Raw-field/causal reproduction requires the local archive.",
        '','## Recommended next step','',
        'Qualify a noise-aware shape measurement before extending this pilot to much longer horizons. Use preserved fields and held-out synthetic/noisy controls to distinguish background noise, local core deformation, phase movement and actual disappearance. Keep the present criteria and negative result intact; any revised localization/filter/threshold is a separate prospectively validated version with its own control requirements. Only then decide whether T=1 or the historical T=15 window warrants its cost, whether finer-tick checks need expansion, and whether exact-law engine acceleration is necessary. No new simulations, particle-identity claim, or physical-vacuum interpretation is introduced by this report.',
        '',
        'Specification: `spec/MG12H3_pattern_persistence_plan.md`. Causal/measurement audit: `src/check_mg160_pilot.py`. Exploratory diagnostics: `src/diagnose_mg160_pilot.py`. Raw decisions: `P1_analysis_results.json`, `P1_exploratory_diagnostics.json`, `P1_outcome_summary.json`, `P1_budget_summary.json`. Do not rerun production to replace paths.','']
        # Presentation-only spacing; raw files retain full precision.
        text='\n'.join(lines)
        for old,new in [('time0.1','time 0.1'),('the101','the 101'),('two-sided95%','two-sided 95%'),('after0','after 0'),('in4','in 4'),('from0.','from 0.'),(' to0.',' to 0.'),('from1.','from 1.'),(' to1.',' to 1.'),('from2.','from 2.'),(' to2.',' to 2.'),('at`','at `'),('before sampling;','before sampling;'),('paths,16','paths, 16'),('pairs,16','pairs, 16'),('=10,500sites','=10, 500 sites'),('time10','time 10'),('a20-unit','a 20-unit'),('Every0.001','Every 0.001'),('All1600','All 1600'),('covering138','covering 138'),('in77','in 77'),('against the100000','against the 100000'),('threads1','threads 1'),('is1800','is 1800'),('about1.26GB','about 1.26 GB'),('the2GB','the 2 GB')]:text=text.replace(old,new)
        text=re.sub(r'\b([a-z]{2,})(?=\d)',r'\1 ',text)
        text=re.sub(r'(?<=\d)(?=(?:minutes|hours|core-hours|charged|million|measured|seconds)\b)',' ',text)
        text=text.replace('## What the exploratory field average suggests','The bottom figure panel displays the central50 of the100model-length ring. Correlation and effective extent are calculated on the entire ring.\n\n## What the exploratory field average suggests')
        text=text.replace('central50','central 50').replace('the100model-length','the 100 model-length')
        (BASE/'MG12H3_pattern_persistence_P1_results.md').write_text(text,encoding='utf-8');atomic_json(attempt,dict(status='COMPLETE',elapsed_seconds=time.perf_counter()-started));print('Report complete: registered per-path retention0/4 at both ticks; no control detections; exploratory ensemble result separate.')
    except BaseException:
        import traceback;atomic_json(attempt,dict(status='FAILED',elapsed_seconds=time.perf_counter()-started,error=traceback.format_exc()));raise
if __name__=='__main__':report()

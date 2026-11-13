"""MG162A finite empirical replay results and MG162B concrete readiness packet."""
import os
for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):os.environ[k]='1'
import json,time
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from run_mg162_challenge import BASE,read,hashes
from mg10e2c1_pool import atomic_json
def report():
    start=time.perf_counter();a=read(BASE/'challenge_results.json');v=read(BASE/'separate_recomputation.json');p=read(BASE/'plan.json');assert v['passed'] and a['source_hashes']==hashes();cases=a['cases'];primary=[r for r in cases if r['case']['primary']];assert a['supported_within_finite_challenge'];positive=[r for r in primary if r['case']['truth']=='compatible'];negative=[r for r in primary if r['case']['truth']=='reject'];boundaries=[r for r in cases if not r['case']['primary']];rows=[];edge_rows=[]
    for r in primary:rows.append(f"| {r['case']['id']} | {r['compatible']} / {r['N']} | {min(g['fraction'] for g in r['groups']):.6g}–{max(g['fraction'] for g in r['groups']):.6g} | {'PASS' if r['passed'] else 'NOT SUPPORTED'} |")
    for r in boundaries:edge_rows.append(f"| {r['case']['id']} | {r['compatible']} / {r['N']} | {100*r['compatible']/r['N']:.6g}% |")
    fig,axes=plt.subplots(1,3,figsize=(12,4),sharey=True)
    subsets=[([r for r in boundaries if 'amplitude' in r['case']['id']],'amplitude',.5,'Amplitude boundary'),([r for r in boundaries if 'scale' in r['case']['id'] and r['case']['scale']<1],'scale',.5,'Lower dilation boundary'),([r for r in boundaries if 'scale' in r['case']['id'] and r['case']['scale']>1],'scale',2.,'Upper dilation boundary')]
    for ax,(rows_,key,cut,title) in zip(axes,subsets):
        x=[r['case'][key] for r in rows_];y=[r['compatible']/r['N'] for r in rows_]
        for group in range(4):ax.plot(x,[r['groups'][group]['fraction'] for r in rows_],color='#888888',alpha=.5,lw=.8)
        ax.plot(x,y,'o-',color='#245a9c',lw=2);ax.axvline(cut,color='#b74127',ls='--');ax.set_title(title);ax.set_xlabel('Known planted '+key);ax.grid(alpha=.2);ax.set_ylim(-.03,1.03)
    axes[0].set_ylabel('Core-compatible fraction in finite replay set');fig.suptitle('MG162A: near-cutoff labels are sensitive to noise and fitted-bank resolution');fig.text(.5,.01,'Blue: all finite replays. Gray: each parent/reflection group. Two original noise snapshots; no population confidence bands.',ha='center',fontsize=8);fig.tight_layout(rect=(0,.05,1,.92));fig.savefig(BASE/'MG162_boundary_sensitivity.png',dpi=150);fig.savefig(BASE/'MG162_boundary_sensitivity.svg');plt.close(fig)
    budget=36912.6878852*1.25;forecast=budget/8/3600;raw=read(BASE/'local_raw_manifest.json');null=next(r for r in primary if r['case']['id']=='null');low=min(r['compatible']/r['N'] for r in positive);a04=next(r for r in boundaries if r['case']['id']=='boundary_amplitude_0.4');s201=next(r for r in boundaries if r['case']['id']=='boundary_scale_2.01')
    atomic_json(BASE/'fresh_readiness.json',dict(stage='MG162B',status='READY FOR CONCRETE CRITERIA/BUDGET APPROVAL',pairs=8,prepared_pairs=4,homogeneous_pairs=4,duration=.1,ticks=[.0001,.00002],observation=.001,forecast_worker_seconds=budget,forecast_wall_hours_at8=forecast,worker_cap=100000,parent_cap=1800,job_cap=12000,workers=8,threads=1,local_raw_reserve_GB=2,approval_record_absent=not (BASE/'fresh_approval.json').exists(),fresh_plan_absent=not (BASE/'fresh_plan.json').exists(),new_ModelG_trajectories=0))
    text=f'''# MG162A — empirical-noise/off-grid challenge and fresh-pilot readiness

The frozen noise-aware measurement handled clear, off-grid cores well in this additional challenge using actual control-noise snapshots, while revealing important limits near its decision cutoffs. All15primary finite-replay requirements passed, but fields close to the amplitude or dilation boundaries frequently changed label, and the0.125-wide fitted scale grid cannot sharply distinguish widths just below and above2. This supports a carefully scoped new-data pilot on clear core-versus-control separation, not a claim that noisy boundary decisions are exact. No new Model G trajectories were run; MG161 and the original P1 outcome remain unchanged.

## What passed, and what did not become qualified

**The prescribed finite interior challenge is supported. Sharp near-boundary classification remains unqualified.** The method, window, bank, noise calibration and classification rules were frozen throughout. Source/design were published before evaluation at`93b28d8`. No cases were removed or thresholds refitted after seeing these results.

The challenge exhaustively combined500rotations and both reflections of two final homogeneous P1 noise fields (replicates2/3). They were excluded from MG161 calibration but had previously been inspected; they are not unseen real trajectories. These transformations preserve their empirical values/tails and circular spatial structure. Every case has2000finite replays; they are not2000independent physical replications. No binomial/population confidence interval is claimed for their fractions.

Nine known retained-core cases used amplitudes0.8/1/1.2 and off-grid scales0.85/1.1/1.6. Five rejected shapes tested weak, narrow, broad, inverted and split features. The null had{null['compatible']}/{null['N']}compatible labels. Retained-case aggregate acceptance ranged from{100*low:.6g}% to100%; every individual parent/reflection group met the registered90%criterion. The weak case had8/2000compatible labels; all other primary rejected shapes hadzero. All rejected/noise group limits were met.

| Primary case | Core-compatible / finite replays | Fraction range across four groups | Result |
|---|---:|---:|---|
{chr(10).join(rows)}

## Boundary sensitivity matters

A true amplitude0.4, below the nominal0.5cutoff, was labeled compatible in{100*a04['compatible']/a04['N']:.6g}% of these replays. A true dilation2.01, just above2, was labeled compatible in{100*s201['compatible']/s201['N']:.6g}% because the fitted width commonly lands on the accepted2.0bank point. These are expected manifestations of noisy estimation and discrete resolution, but they prevent treating the point label as an accurate sharp test of a true parameter. The descriptive boundary cases were fixed before evaluation and are all reported; they do not count as passing interior cases.

![Empirical boundary sensitivity](MG162_boundary_sensitivity.png)

| Boundary case | Compatible / finite replays | Conditional finite fraction |
|---|---:|---:|
{chr(10).join(edge_rows)}

No confidence interval on a real physical width follows from the fitted bank index. The later pilot therefore records boundary-advisory flags, as well as all unchanged MG161 labels. The advisory is engineering sensitivity information, not a calibrated confidence interval or proof of true shape.

## Verification and limits

Separate direct weighted sums verified{v['scalar_fit_checks']}selected scalar fits/residuals; classification branches and all finite-group fractions were independently recomputed, agreeing with all15primary verdicts. Source and input hashes matched. The challenge evaluated72000test fields in{a['elapsed_seconds']:.6g}seconds under a300-second, single-worker/thread-one cap. Full replay fields occupy{raw['total_bytes']/1e6:.6g}MB and remain local under`local_raw/MG162_shape_challenge`; only compact features, report/figure, source and manifests are public. Same-author/shared-field and template-bank dependencies are disclosed; this is not external replication.

Only two previously seen noise parents were replayed. Planting a known signal into a saved noise field does not reproduce dynamic signal–noise coupling, time dependence or all possible covariance changes. Passing this finite challenge is therefore not independent validation of real-trajectory survival, a real-noise false-positive population rate, or a physical particle. Earlier P1's negative global-measurement score and MG161's conditional synthetic outcome stay intact. Closure of earlier studies is not inferred.

## Concrete next run: MG162B, approval pending

The code and criteria are prepared for eight entirely fresh seed pairs: four prepared and four mass-matched homogeneous, ticks0.0001/0.00002 through time0.1, with101frames per marginal path. All earlier P1 data are excluded from the new denominator. The instrument and MG161 calibration are frozen before new data; no PSD fitting or threshold selection on the fresh controls is permitted. No fresh seeds have been registered and no new kinetic trajectories exist.

The strict finite-pilot separation criterion requires every prepared marginal path to be core-compatible at every frame with no boundary-advisory frames, and no homogeneous marginal path to be compatible at any frame. Otherwise report not established, preserving all outcomes. This tests clear-case independent-data separation; four pairs per condition cannot establish small population error rates or long lifetimes. Control crossings may be false detections or new structures to investigate, not automatic physical births. All raw/count/event/RNG records must pass audits and the new metrics must be separately recomputed before reporting a result.

Measured prior pilot cost plus25% allowance forecasts{budget:.6g}combined worker-seconds ({budget/3600:.6g}core-hours), or{forecast:.6g}hours on eight workers. Proposed hard ceiling100000worker-seconds (27.7778core-hours),12000/job, plus1800active parent/audit seconds; max8workers/numerical-library threads1. Reserve2GB of local raw storage. Full data stay local; compact publication only. No retries, sample replacements, extensions or altered thresholds. This is a new numerical allocation requiring explicit PI criteria/budget approval under RESEARCH_PROTOCOL.md workflow3.

Specification of the evaluated challenge: `spec/MG162_shape_validation.md`. Detailed prospective fresh pilot: `spec/MG162B_fresh_pilot.md`. Entry: `src/run_mg162_fresh.py register`, then `run`, `audit`, `score` only after real approval, source/seed publication and the required checkpoints. The registration gate requires fresh_approval.json; it does not yet exist. No longerT=1/T=15 study is authorized.
'''
    import re
    text=re.sub(r'\b([a-z]{2,})(?=\d)',r'\1 ',text);text=text.replace('All15','All 15').replace('all15','all 15').replace('two2','two 2').replace('have2000','have 2000').replace('has2000','has 2000').replace('not2000','not 2000').replace('The0.','The 0.').replace('the0.','the 0.').replace('than2','than 2').replace('above2','above 2').replace('accepted2.0','accepted 2.0').replace('source/design','source/design').replace('at`','at `').replace('under`','under `').replace('evaluated72000','evaluated 72000').replace('with101','with 101').replace('plus25%','plus 25%').replace('ceiling100000','ceiling 100000').replace(',12000',', 12000').replace('plus1800','plus 1800').replace('Reserve2GB','Reserve 2 GB').replace('with all15','with all 15').replace('hadzero','had zero')
    (BASE/'MG162_validation_readiness.md').write_text(text,encoding='utf-8');print('MG162A reported; MG162B prepared, no fresh model registration or draws.')
if __name__=='__main__':report()

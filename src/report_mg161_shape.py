"""MG161 instrument benchmark report; P1 outcomes remain frozen and separate."""
import os
for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):os.environ[k]='1'
import time,json,hashlib
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from run_mg161_shape import BASE,RAW,read,hashes
from mg10e2c1_pool import atomic_json
def report():
    start=time.perf_counter();a=read(BASE/'benchmark_results.json');separate=read(BASE/'separate_recomputation.json');application=read(BASE/'P1_exploratory_application.json');assert separate['passed'] and a['source_hashes']==hashes();plan=read(BASE/'plan.json');assert a['qualified_within_synthetic_scope'];cal=read(BASE/'calibration.json');cases=a['cases'];null=next(r for r in cases if r['case']['id']=='null');positive=[r for r in cases if r['case'].get('truth')=='compatible'];negative=[r for r in cases if r['case'].get('primary') and r['case'].get('truth')=='reject' and r['case']['id']!='null'];sensitivity=cases[-1]
    rows=[]
    for r in cases:
        c=r['case'];bound=f"lower {r['lower']:.6g}" if c.get('truth')=='compatible' else (f"upper {r['upper']:.6g}" if c.get('primary') else 'descriptive only');rows.append(f"| {c['id']} | {r['compatible']} / {r['N']} | {bound} | {'PASS' if r.get('passed') else 'secondary'} |")
    real_rows=[]
    for j in application['jobs']:
        for o in j['outcomes']:real_rows.append(f"| {j['id']} | {(.0001,.00002)[o['side']]:g} | {o['compatible_frames']} / 101 | {o['final_amplitude']:.6g} | {o['final_scale']:.6g} | {o['final_label']} |")
    with np.load(BASE/'P1_exploratory_features.npz') as z:f=z['features'];times=z['times']
    fig,axes=plt.subplots(2,2,figsize=(12,7),sharex=True);colors=['#245a9c','#c66c22'];initial=float(np.load(BASE/'noise_model.npz')['initial_amplitude'])
    for j in range(8):
        prepared=j<4
        for s in (0,1):
            color=colors[s] if prepared else '#888888';style='-' if s==0 else '--'
            if prepared:
                axes[0,0].plot(times,f[j,:,s,0]/initial,color=color,ls=style,lw=1,alpha=.8);axes[0,1].plot(times,f[j,:,s,1],color=color,ls=style,lw=1,alpha=.8)
            axes[1,0].plot(times,f[j,:,s,4]/cal['detection_threshold'],color=color,ls=style,lw=1,alpha=.8 if prepared else .3);axes[1,1].plot(times,f[j,:,s,5]/cal['residual_threshold'],color=color,ls=style,lw=1,alpha=.8 if prepared else .3)
    axes[0,0].axhline(.5,color='red',ls=':');axes[0,0].set_ylim(.4,1.3);axes[0,0].set_title('Prepared cores: fitted amplitude / initial');axes[0,1].axhspan(.5,2,color='#e2f0e2');axes[0,1].set_ylim(.4,2.1);axes[0,1].set_title('Prepared cores: fitted local dilation');axes[1,0].axhline(1,color='red',ls=':');axes[1,0].set_title('Search gain / null threshold (higher supports detection)');axes[1,1].axhline(1,color='red',ls=':');axes[1,1].set_title('Local residual / calibrated limit (lower fits better)')
    for ax in axes.ravel():ax.grid(alpha=.18);ax.set_xlabel('Post-preparation model time')
    fig.suptitle('MG161: noise-aware core fits on previously seen P1 fields — exploratory',fontsize=14);fig.text(.5,.015,'Blue: tick0.0001. Orange dashed: tick0.00002. Gray: homogeneous controls. Original P1 retention scores are unchanged.',ha='center',fontsize=9);fig.tight_layout(rect=(0,.04,1,.94));fig.savefig(BASE/'MG161_saved_core_diagnostics.png',dpi=150);fig.savefig(BASE/'MG161_saved_core_diagnostics.svg');plt.close(fig)
    raw=read(BASE/'local_raw_manifest.json');total=sum(c['N'] for c in plan['cases']);minimum=min(r['compatible']/r['N'] for r in positive);source_binding=hashes();active=a['elapsed_seconds']+separate['elapsed_seconds']+application['elapsed_seconds'];atomic_json(BASE/'budget_summary.json',dict(benchmark_seconds=a['elapsed_seconds'],benchmark_cap=300,workers=1,numerical_threads=1,recorded_benchmark_check_application_seconds=active,synthetic_fields=total,local_raw_bytes=raw['total_bytes'],new_ModelG_trajectories=0,raw_published=False))
    text=f'''# MG161 — noise-aware local shape measurement results

The new measurement passed all15requirements of its prespecified synthetic-noise benchmark. It fits the central pattern's amplitude and local dilation while allowing the amount of background noise observed in homogeneous controls, rather than judging shape from whole-ring spread. Used exploratorily on the saved P1 trajectories, it finds compatible central cores in every prepared frame at both tick sizes and no compatible detections in the homogeneous controls. This strengthens the interpretation that noise-sensitive global measurements contributed to P1's apparent loss of shape. It does not prove physical-particle survival: the benchmark assumes a stationary Gaussian noise model estimated from two paths, the real P1 fields were already seen, and the original P1 negative retention result remains unchanged.

## Registered benchmark outcome

**QUALIFIED WITHIN THE PRESCRIBED CONDITIONAL SYNTHETIC BENCHMARK.** All15primary case requirements passed. This is an instrument-development result, not a new stochastic Model G persistence qualification. Source, specification and seeds were published before synthetic-noise draws at`3209a9b`. No benchmark cases, thresholds, sample counts or outcomes were changed after validation. A missing CLI import was corrected before plan generation/noise draws; no scientific retry occurred.

There were{total}independent synthetic test fields across separate calibration, validation and sensitivity cases. The null validation produced{null['compatible']}/{null['N']}core-compatible false detections ({100*null['compatible']/null['N']:.6g}%). Its family-adjusted exact upper bound is{null['upper']:.6g}, below the registered0.01limit. Across nine retained-core cases, point acceptance ranged from{minimum*100:.6g}% to100%; every adjusted lower bound exceeded0.90. Weak signals produced2compatible labels out of512; inverted, narrow, broad and split controls producedzero. The double-noise case was descriptive and did not affect the primary verdict.

| Known-truth case | Core-compatible classifications | Adjusted bound | Result |
|---|---:|---|---|
{chr(10).join(rows)}

Each primary bound uses alpha0.05/15 and exact Clopper–Pearson construction. The15one-sided requirements have nominal simultaneous95% conditional-binomial coverage. This covers independent draws from the fixed synthetic model, not uncertainty in the estimated spectrum, nonstationary/non-Gaussian real fluctuations or arbitrary deformation families. Max-score calibration covers the complete position/scale search under that same null; its rank interpretation is not a real-vacuum probability statement.

## How the measurement differs

The reference is the original unrounded X profile. A fixed initial core-energy radiusR95=3.2model-lengths defines Gaussian spatial weighting with sigma6.4. The method tries23dilations from0.25to3 and all500cyclic locations, fitting positive amplitude plus a local constant background. It records the best location, amplitude, dilation, search gain, residual and a distant second-core score. A bank-edge fit or multiple compatible cores is treated as deformation/ambiguity rather than quietly accepted as one persistent object. The amplitude and dilation limits are operational shape-family criteria, not physical particle identity.

Noise power was estimated only from the final coarse-tick homogeneous P1 fields of replicates0/1. Both had previously been seen; they are calibration inputs, not unseen test data. Fresh independent Gaussian fields were generated from the frozen smoothed periodogram with a fixed spectral floor. Independent calibration draws set maximum null-gain and correct-shape-residual thresholds. Separate validation draws then tested known retained/lost/split/inverted shapes. No kinetic equations, reservoirs, ticks, Omega or token states were changed and no new Model G trajectories were run.

## Saved P1 application — exploratory only

The new instrument labels all101frames compatible in each prepared marginal path: four independent seed pairs at each tick. Homogeneous paths havezero compatible frames. This is a post-result application of a newly tested instrument to previously inspected data, not independent real-trajectory validation and not a replacement survival score.

| Existing P1 path | Tick | Compatible frames | Final relative amplitude | Final fitted dilation | Final label |
|---|---:|---:|---:|---:|---|
{chr(10).join(real_rows)}

![Saved-core diagnostics, explicitly exploratory](MG161_saved_core_diagnostics.png)

The original P1 rule still has0/4prepared pairs retaining its global extent/coherence criterion continuously through0.1at each tick. The new local model estimates the central structure separately from diffuse background fluctuations; its different answer concerns a different measurement. The fitted dilation is a point estimate on a0.125grid, not a confidence interval on physical width. P1's four pairs cannot establish population survival, low real false-positive rates, long lifetimes or unique identity. A template-family test also leaves blind spots for other kinds of deformation and interference. The primary known-shape dilations lie on the fitted bank grid; off-grid shapes and cases near the decision boundaries still require validation.

## Verification, cost and storage

Noiseless translation/scale recovery, constant-offset and flat-field controls passed;92direct dot-product checks agreed with the FFT implementation before noise draws. Afterward,54direct weighted fits/residual sums agreed with stored features, all classifications were separately recomputed, and all15bounds were verified by inverting binomial tails rather than reusing beta quantiles. Calibration maxima and every input hash were rechecked. These are separate computations by the same author, sharing fields, the template bank and statistical libraries; no external replication is claimed.

The benchmark used{a['elapsed_seconds']:.6g}seconds against its300-second, one-worker/thread-one cap. Benchmark, verification and saved-field application used{active:.6g}recorded seconds in total. Full synthetic fields occupy{raw['total_bytes']/1e6:.6g}MB under`local_raw/MG161_noise_shape`; they remain local and ignored by Git. Public artifacts contain only compact features, thresholds, seeds, spectrum/template configuration, source/input hashes, report/figure and exact-bound summaries. The full P1 archive remains local and unchanged. The publication-size guard remains enabled.

## Recommended next step

Freeze this instrument as a candidate for a separately planned test on fresh real trajectories, including independent homogeneous and deformed/removed-core controls where their construction is scientifically appropriate. Quantify sensitivity to the estimated spectrum, Gaussian/stationary assumptions and the discrete fitted scale before treating its labels as survival evidence. Real-noise validation and any longer T=1/T=15 run require their own prospective scope and measured budget. Do not alter or reopen P1's frozen negative result to obtain a favorable score. MG160-P1 closure remains pending; no new model runs were launched.

Specification: `spec/MG161_noise_aware_shape.md`. Model/benchmark: `src/mg161_local_shape.py`, `src/run_mg161_shape.py`. Separate check: `src/check_mg161_shape.py`. Exploratory application: `src/apply_mg161_saved.py`. Compact features and local raw manifest support recomputation; full scalar-field checks use the local synthetic archive.
'''
    import re
    text=re.sub(r'\b([a-z]{2,})(?=\d)',r'\1 ',text);text=re.sub(r'(?<=\d)(?=(?:independent|seconds|model-lengths|core-compatible|recorded|MB|limit|at\b))',' ',text);text=text.replace('at`','at `').replace('under`','under `').replace('R95=3.2model-lengths','R95=3.2 model-lengths').replace('to100%','to 100%').replace('through0.1','through 0.1')
    text=text.replace('All15','All 15').replace('15requirements','15 requirements').replace('15primary','15 primary').replace('15one-sided','15 one-sided').replace('2compatible','2 compatible').replace('producedzero','produced zero').replace('sigma6.4','sigma 6.4').replace('0.25to3','0.25 to 3').replace('all500','all 500').replace('all101','all 101').replace('havezero','have zero').replace('a0.125grid','a 0.125 grid').replace('source14','source 14')
    (BASE/'MG161_noise_aware_shape_results.md').write_text(text,encoding='utf-8');print('Report complete; all15conditional benchmark requirements supported, P1 exploratory application separate.')
if __name__=='__main__':report()

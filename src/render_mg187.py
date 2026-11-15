"""Render fixed MG187 audit exports; no additional raw-data pass or tuning."""
import json,csv
from collections import Counter
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from audit_mg187 import OUT
from mg10e2c1_pool import atomic_json

def run():
 a=json.loads((OUT/'audit.json').read_text());rows=[v for v in a['loss_runs'] if v['X_detector_loss'] is not None];high=a['high_population_runs'];assert len(rows)==29
 # Separate arithmetic path recomputes timing categories from continuous differences.
 for key,stored in [('G_minus_X_half','counts_half_amplitude'),('G_minus_detector_loss','counts_detector')]:
  count=Counter('unresolved/censored' if v[key] is None else 'before' if v[key]<-2 else 'after' if v[key]>2 else 'together' for v in rows);assert dict(count)==a[stored]
 diffs=np.array([v['G_minus_X_half'] for v in rows]);windows=np.array([[[w['G_amplitude'],w['X_amplitude'],w['G_residual']] for w in r['windows']] for r in rows]);stats=dict(median_G_minus_X_half=float(np.median(diffs)),range_G_minus_X_half=[float(diffs.min()),float(diffs.max())],mean_windows=windows.mean(axis=0).tolist(),preloss_G_range=[float(windows[:,1,0].min()),float(windows[:,1,0].max())],preloss_G_above_half=int((windows[:,1,0]>.5).sum()),mean_detector_offset=float(np.mean([v['G_minus_detector_loss'] for v in rows])))
 atomic_json(OUT/'arithmetic_check.json',dict(passed=True,separate_category_recount=True,statistics=stats,external_replication=False))
 def f(v):return '' if v is None else f'{v:.6g}'
 with (OUT/'Omega1000_loss_timing.csv').open('w',newline='',encoding='utf-8') as file:
  writer=csv.writer(file);writer.writerow(['run','X_detector_loss','G_half_onset','X_half_onset','G_minus_X_half','G_vs_X_half','G_vs_detector','G_pre10','X_pre10','G_post10','X_post10'])
  for r in rows:writer.writerow([r['id'],f(r['X_detector_loss']),f(r['G_half_amplitude']['onset']),f(r['X_half_amplitude']['onset']),f(r['G_minus_X_half']),r['relative_to_X_half'],r['relative_to_detector'],*[f(x) for x in (r['windows'][1]['G_amplitude'],r['windows'][1]['X_amplitude'],r['windows'][2]['G_amplitude'],r['windows'][2]['X_amplitude'])]])
 for omega in (3000,10000):
  with (OUT/f'Omega{omega}_per_run_G.csv').open('w',newline='',encoding='utf-8') as file:
   writer=csv.writer(file);writer.writerow(['run','historical_X','strict_X','block','valid_frames','total_frames','G_amplitude_reference_ratio','G_amplitude_SD','mean_instantaneous_residual','averaged_profile_residual'])
   for r in high:
    if r['omega']!=omega:continue
    for b in r['blocks']:writer.writerow([r['id'],r['X_persistent'],r['strict'],b['name'],b['valid'],b['total'],*[f(b[k]) for k in ('amplitude','amplitude_SD','instantaneous_residual','averaged_profile_residual')]])
 fig,axes=plt.subplots(1,3,figsize=(15,4.2),layout='constrained');curve=a['loss_curve'];v=np.array(curve['mean']);axes[0].plot(curve['relative_time'],v[:,0],label='G');axes[0].plot(curve['relative_time'],v[:,1],label='X');axes[0].axvline(0,color='gray',ls='--');axes[0].axhline(.5,color='gray',lw=.6);axes[0].set(xlabel='Time relative to recorded X loss',ylabel='Mean amplitude / species reference',title='29 losses at Ω = 1000');axes[0].legend();axes[0].text(.03,.04,'After loss: anchor held at last X position.\nContributors vary near the initial boundary.',transform=axes[0].transAxes,fontsize=8)
 axes[1].axhspan(-2,2,color='gray',alpha=.15);axes[1].scatter(np.arange(1,30),diffs,s=20);axes[1].axhline(0,color='gray',lw=.6);axes[1].set(xlabel='Loss run (replicate order)',ylabel='G half-amplitude onset − X onset',title='22 together; 6 G earlier; 1 later')
 for omega,color in [(3000,'tab:green'),(10000,'tab:orange')]:
  rr=[r for r in high if r['omega']==omega and r['blocks'][1]['amplitude'] is not None];axes[2].scatter([r['blocks'][1]['amplitude'] for r in rr],[r['blocks'][1]['averaged_profile_residual'] for r in rr],s=22,color=color,label=f'Ω = {omega}, n = {len(rr)}')
 axes[2].set(xlabel='Final-block G amplitude / reference',ylabel='Averaged-profile orthogonal residual',title='Per-run G, time 900–1000');axes[2].axvline(1,color='gray',lw=.6);axes[2].legend(fontsize=8);fig.savefig(OUT/'MG187_G_loss_summary.png',dpi=160);plt.close(fig)
 lines=['# MG187 bounded G loss audit','',
 'G does not emerge as a consistently early weak link in these particle losses. Across the 29 lost Ω=1,000 trajectories, G and X usually faded together: 22 were within two model-time units of each other, six had G fade earlier, and one had G fade later. The median G lead was only 0.5 time units. G remained substantial immediately before the recorded X loss, then both species lost their template-shaped signal. At Ω=3,000 and 10,000, the per-run averaged G amplitude and shape remained similar from the beginning to the end, with appreciably less noise at 10,000. This points toward a coupled particle collapse rather than a general pattern of G disappearing first; temporal ordering alone cannot identify which reaction or fluctuation causes collapse.','',
 '## Scope and fixed rule','',
 'This was one completed, bounded pass through existing data: 30 prepared Ω=1,000 runs (29 losses plus one survivor), 30 at Ω=3,000 and 30 at Ω=10,000. No new simulations, detector changes, parameter fitting, threshold sweeps or instrument refinement. The choices were committed and published before computing these loss-aligned results, but the MG185 outcomes were already known: this is explicitly post hoc exploratory diagnosis, not a new persistence gate.','',
 'Each species amplitude is its centered reference projection divided by its own unrounded reference norm. Fading means a trailing ten-time-unit mean below half amplitude, remaining below for ten units. Timing labels allow ±2 units as together. Both species use the same rule. The first sustained crossing is a diagnostic of template amplitude, not physical extinction.','',
 '## What happens in the 29 Ω=1,000 losses','',
 '| Comparison | G before | Together | G after | Unresolved |','|---|---:|---:|---:|---:|',
 '| Same half-amplitude rule for G and X | 6 | 22 | 1 | 0 |',
 '| G half-amplitude versus recorded X detector loss | 0 | 2 | 27 | 0 |','',
 f"For the matched amplitude comparison, G minus X onset has median {f(stats['median_G_minus_X_half'])}, range [{f(diffs.min())}, {f(diffs.max())}] model-time units. The six earlier G crossings lead by only 2.5–3.5 units; the one later crossing lags by 3.5. These small leads are not a universal pre-collapse G disappearance. The surviving run r16 crossed neither sustained half-amplitude threshold through time 1,000 and is retained as censored, not one of the 29 losses.",
 '', 'The two rows compare different clocks. The X loss timestamp is the last matched detection, confirmed retrospectively after twenty unmatched units. The amplitude rule uses a trailing mean, which delays its reported crossing. Thus the 27 “after detector loss” cases do not establish that G outlives X physically. The matched G-versus-X amplitude row is the fairer ordering comparison.','',
 '| Window relative to X detector loss | Mean G amplitude | Mean X amplitude | Mean instantaneous G residual |','|---|---:|---:|---:|']
 for label,v in zip(['−20 to −10','−10 to 0','0 to +10','+10 to +20'],stats['mean_windows']):lines.append('| '+label+' | '+' | '.join(f(x) for x in v)+' |')
 lines+=['',f"All 29 runs had average G amplitude above 0.5 in the ten units before X loss; their values ranged from {f(stats['preloss_G_range'][0])} to {f(stats['preloss_G_range'][1])}. Both species declined across the loss transition. The full per-run timing table is [Omega1000_loss_timing.csv](Omega1000_loss_timing.csv). No timing threshold was adjusted after these results.",
 '', 'Until loss, alignment follows only the existing unique-original X track. Missing positions carry the last known position; after loss that anchor is frozen. Post-loss values therefore measure G/X overlap at the last known location. A moving remnant could evade that measurement, and this audit deliberately did not add a new tracker to resolve it.','',
 '## Per-run G amplitude and shape at higher populations','',
 'Amplitude is normalized to the reference. Shape residual is the norm of the component orthogonal to that reference, also divided by the reference norm; zero would mean a perfect scaled shape. “Averaged profile” first averages aligned fields within each run and block, then measures residual. It suppresses temporal noise. “Instantaneous residual” averages the residual of each separate snapshot and therefore retains noise. Neither is an individual particle pass/fail threshold.','',
 '| Ω | Block | Contributing runs | G amplitude median [range] | Averaged-profile residual median [range] | Instantaneous residual median |','|---:|---|---:|---|---|---:|']
 for row in a['cohort']:
  m=row['metrics'];amp=m['amplitude'];res=m['averaged_profile_residual'];lines.append(f"| {row['omega']} | {row['block']} | {amp['n']}/30 | {f(amp['median'])} [{f(amp['min'])}, {f(amp['max'])}] | {f(res['median'])} [{f(res['min'])}, {f(res['max'])}] | {f(m['instantaneous_residual']['median'])} |")
 lines+=['','Blocks are time 0–100 (excluding 100) and 900–1000 (including 1000). At Ω=3,000, r10 has no unique final-block track and is explicitly missing; r22 contributes 195/201 final-block frames but loses terminal uniqueness. Final totals are 5,823/6,030 valid frames at 3,000 and all 6,030 at 10,000. Full per-run amplitude, amplitude variability, both residuals and missing counts: [Ω=3,000 table](Omega3000_per_run_G.csv), [Ω=10,000 table](Omega10000_per_run_G.csv).','',
 'The averaged-profile residual is stable between initial and final blocks: about 0.300 at Ω=3,000 and 0.163 at 10,000. This gives no sign of widespread progressive G shape erosion among the trackable profiles. Instantaneous residuals remain much larger than the reference shape norm—about 3.82 and 2.09 in the final blocks. G is therefore still noisy in individual full-domain snapshots; near-unit mean amplitude must not be read as a noise-free instantaneous shape. No noise subtraction or residual recalibration was introduced.','',
 '![G loss timing and high-population per-run measurements](MG187_G_loss_summary.png)','',
 '## Checks and conclusion','',
 'All 9,000 input segments matched the previous local hash manifest. Raw-field projections matched saved valid-frame amplitudes to within 5.11e-15. The orthogonal-energy identity and all higher-population block means passed. Synthetic before/together/after/censored cases passed, and a separate arithmetic path recounted timing labels from continuous offsets. Shared reference, X track/loss events and FFT implementation remain dependencies; this is same-author checking, not external replication. The raw-data pass took about 22.4 seconds on one worker. Scientific histories remain unchanged and local.','',
 '**Answer:** G is noisy, but this audit does not show it systematically fading well ahead of X. The dominant observed behavior is G and X losing their coherent pattern together. The six small G leads are real under this fixed diagnostic but are insufficient to identify G as the causal weak link. Higher population preserves a stable averaged G component in the trackable particles, with Ω=10,000 clearly cleaner than 3,000.','',
 'The requested bounded audit ends here. No further instrument development or simulations are launched. MG185 persistence verdicts remain unchanged; PI closure remains pending. Inputs and exact choices: [MG187 specification](../../spec/MG187_bounded_G_loss_audit.md), [machine-readable audit](audit.json), [verification](verification.json), [arithmetic check](arithmetic_check.json).']
 (OUT/'MG187_results.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')
if __name__=='__main__':run()

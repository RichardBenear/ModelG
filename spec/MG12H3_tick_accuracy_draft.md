# MG158 — paired tick accuracy qualification, ready for budget approval

This is a draft, not production registration. PI asks to proceed to the next step. Preparation and saved-only validation are authorized; this packet presents the newly measured full cost for approval under RESEARCH_PROTOCOL.md workflow 3. No new simulation draws. MG156/MG157 closure remains separate.

## Scientific question and fixed scope

Do ticks 0.0001 and 0.00005 each agree sufficiently with 0.00002 on the canonical finite LOCAL/LOCAL Omega=10, 500-site state, through T=0.01? Test both primary times 0.005/0.01 and all G/X/Y projected contrasts plus G/X/Y high-band change powers: 24 mean differences overall. Finest tick is a numerical reference, not continuum truth or a calibrated physical time.

Retain the MG145v2 unrounded spatial templates, actual integer initial counts and contrast margins 0.05*abs(initial projected contrast): approximately G=0.0000723355, X=0.00332205, Y=0.00295473. Keep the previously approved absolute power margin 0.05 in its fixed reference scale. Retain the same Fourier band |k|>=125 on 500 sites, including Nyquist. This does not test full-field, particle-shape, physical-vacuum or long-term persistence equivalence. Retain full observation fields to expose projection blind spots.

## Fresh sampling and analysis amendment

Use 256 fresh independent coupled pairs at (0.0001,0.00002) and 64 at (0.00005,0.00002): 320 pairs/640 marginal chains. MG157 pilot samples are excluded. Fixed direct differences only, no fitted blend or compensator. Each pair follows the frozen MG156 strict-tick marginal-preserving generator; use the MG157 compact writer unchanged, with a separate study directory/source binding. Keep counts and checkpoints every 0.001; predeclare replicate zero per comparison as full numeric sentinels. No independent reference ensembles are needed to estimate the paired mean differences.

Prospective amendment to MG145v2: replace independent-ensemble Welch intervals with paired Student-t intervals on independent replicate differences. Each scalar uses mean +/- t(1-.05/48,N-1)*SD/sqrt(N). Bonferroni provides nominal 95% simultaneous coverage for 24 intervals if the individual t intervals have their stated coverage; independence between the 24 observables is not required. Coverage is approximate for non-Gaussian differences. Store all replicate observables and disclose skew/tail/leverage diagnostics without changing method, removing outliers or selecting estimators.

Per tick, equivalent requires all 12 intervals fully inside their fixed +/- margins; incompatible means any interval wholly outside its margin; otherwise unresolved. Any missing/failed job, nonfinite observable, unverified count/event archive, or changed source blocks a complete qualification verdict. Zero observed variance is flagged unresolved, rather than treating repeated values as proof of zero population uncertainty. No sample-size extension after viewing outcomes. No optional stopping or intermediate scientific verdicts. A detectable nonzero difference can still be acceptably small.

## Precision and uncertainty in planning

Saved MG157 paired differences and powers were reconstructed at both primary times. The preliminary 192/64 allocation appeared adequate using point SDs, but the coarser comparison's projected interval edge exceeded its margin when SD was inflated to its one-sided 95% chi-square upper estimate from eight pilot observations. Use 256/64 instead. This variance sensitivity assumes normal pilot differences and holds pilot means fixed; it does not bound the unknown true means. The new nominal worst interval edges are about 0.841466/0.588122 of margin; upper-SD sensitivity about 0.981564/0.782166. Neither conditional power estimates nor this sensitivity guarantee a pass. A true tick bias beyond the margin should yield incompatibility, not a relaxed tolerance.

No prior variance-benefit result is rescored as an accuracy verdict. The saved-only readiness script tests equivalent/incompatible/boundary-overlap/zero-variance and invalid or incomplete synthetic samples; high-band power is independently checked using Parseval normalization.

## Concrete proposed compute and retention budget

MG157 measured mean full-horizon pair costs about 582.484/577.836 worker-seconds under eight-worker contention. Use these measured values and 15% forecasting allowance. Total forecast is about 7.43 hours wall time at eight workers, not a subscription-usage estimate. Split the fixed denominator across three sessions: coarser 86/85/85 and finer 22/21/21 pairs. Each session is estimated below three hours at eight workers. This split is scheduling only; do not condition the second or third session on intermediate scientific outcomes. Integrity failures and budget stops still block continuation.

Proposed ceilings: 100000 charged worker-seconds per session, 300000 total (83.3333 core-hours); up to 1200 seconds per job; at most eight local workers and one numerical-library thread per worker. Forecast admission rule and no-retry policy remain. Allow 1800 separately measured active parent/control/audit seconds in total; time every analysis attempt including interruptions. Controller waiting is wall time, not an extra active core allocation. Report the actual time and any stop or overrun transparently. A short denominator is incomplete, not silently redefined.

Retain all decision checkpoints, histograms, RNG states/seeds, full observation fields, source/environment hashes, logs and budgets; full numeric event/boundary records for two predeclared sentinels. Use the corrected eager-loading audit rather than repeated NPZ access inside event loops. Archive raw data and separate recomputation before reporting. Expected storage is a few hundred MB; reserve 1 GB conservatively. Preserve all historical archives.

Freeze the final source/seed manifest and record PI approval of this concrete packet before launching. External review is not requested or claimed. This packet authorizes no subsequent altered study or physical claim.

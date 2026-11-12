# MG161 — noise-aware local shape measurement, bounded development benchmark

PI instruction: “go ahead with the noise-aware shape measurement.” Implement and test a separate instrument, preserving MG160-P1's negative registered per-trajectory result and its pending closure. No new Model G trajectories, reservoir modifications, threshold relaxation or persistence rescore. This is a bounded one-worker development/known-truth benchmark, not a new physical-survival qualification. Freeze/publish this specification, source and independent seeds before synthetic noise draws. Limit execution to300elapsed seconds, numerical threads1; preserve failures with no retries, sample substitution or benchmark tuning.

## Why a new instrument

P1's global effective extent can increase when background noise is added even if the central feature remains. Fit a local core while accounting for a separately estimated noise distribution. Location, amplitude and dilation are observable estimates, not physical identity or a proof of particle survival. A template-family fit can miss deformation outside that family.

## Fixed construction

Use the unrounded original X reference and actual rounded initial counts, not a template fitted to the later prepared trajectories. Center the reference on its largest absolute deviation and subtract its far-field median outside25model-lengths. Linearly interpolate the reference at dilation scales0.25..3 in steps0.125, zero outside its original coordinate range, then apply the existing sigma0.4 periodic Gaussian smoothing. Observation fields receive the same smoothing and mean centering. RadiusR95 encloses95% of the original smoothed core energy (3.2model-lengths on this input). The fixed Gaussian spatial weight has sigma2*R95. This replaces whole-ring extent with fitted local dilation; it is a new measurement, not a reinterpretation of P1's extent variable.

At each of23scales and all500cyclic positions, fit positive template amplitude and a free local constant background by weighted least squares. Select maximum explained weighted energy (gain). Record amplitude, scale, position, gain, residual weighted energy and a secondary-core score outside2*R95. The full search is used during null calibration, so the location/scale maximization is not ignored in false-positive testing. Known noiseless scale/shift recovery, constant-offset invariance and zero-field checks precede random draws. Computation uses cached FFT cross-correlations, with an independent direct weighted-sum verification.

## Conditional noise model and calibration

Estimate the X concentration noise spectrum from the final coarse-tick homogeneous fields of P1 replicates0/1 only. All P1 paths have already been inspected, so none is described as unseen or statistically held-out. Average their raw centered periodograms, smooth across5frequency bins by reflected padding, floor positive modes at1% of the positive-mode median, and set DC tozero. Generate fresh independent stationary colored-Gaussian fields by filtering seeded independent normal vectors. These are test fields, not token states or new kinetic simulations. The covariance estimate from two real paths is uncertain, and Gaussian/stationary assumptions are restrictive.

Calibration:1024null fields, plus1024known amplitude1/scale1 fields with independent seeds. Freeze the detection threshold at the maximum null search-gain, and the residual threshold at the maximum correct-shape residual. For an exchangeable draw from this exact synthetic null, the max-score rank gives unconditional rejection probability at most1/1025, ignoring ties; it is not a conditional real-vacuum false-positive guarantee. Residual calibration is only at its stated shape; varied-shape validation is separate. Save thresholds before generating any validation cases.

## New operational classification

Absent/weak if gain does not exceed its noise-calibrated threshold or relative amplitude<0.5. Ambiguous multiple cores if a distant secondary core also exceeds that detection threshold and has at least half the best gain. Deformed scale if the best fit is at a bank edge or outside[0.5,2]. Shape mismatch if the residual exceeds its calibrated threshold. Otherwise “core compatible.” This label refers to compatibility with this fitted core family under the specified noise benchmark; it is not the old per-path retention label.

## Fresh validation and outcome rules

Separate validation draws, excluded from calibration:

-4096null fields.
-Nine retained-core cases: amplitude0.8/1/1.2 crossed with scale0.75/1/1.5,512fields each.
-Five rejected cases: weak amplitude0.2; scale0.25; scale3; inverted template; two equal cores separated by half the ring.512fields each.
-One double-noise amplitude1/scale1 case,512fields, descriptive sensitivity only.

Locations follow a fixed cyclic schedule, independent of the noise and outcomes. The nine retained and five rejected shape cases plus the null form15primary cases. For each primary positive case require its one-sided Clopper–Pearson lower bound on “core compatible” frequency>=0.90. For each negative case require the corresponding upper bound<=0.05; the null upper bound must be<=0.01. Each bound uses alpha0.05/15, providing simultaneous nominal95% conditional-binomial coverage across15primary requirements. Gaussian validation fields are independent within each case conditional on the frozen spectrum; these bounds do not quantify the uncertainty of that spectrum or real model coverage.

Declare qualified within the synthetic benchmark only if all15requirements pass. Any failure means not qualified for the full prescribed benchmark; preserve failed cases, do not drop them or refit thresholds. No sample increase or scientific retries. Double-noise sensitivity does not change the verdict. Only after all primary requirements pass may the frozen instrument be applied descriptively to saved P1 fields; those labels remain exploratory and do not overwrite the P1 result or establish new-trajectory persistence. A benchmark failure blocks that dependent application and informs a separate future design.

## Retention and reporting

Full synthetic fields and per-case metadata stay under ignored local_raw/MG161_noise_shape. Commit only compact features, thresholds, seeds/configuration, source/input hashes, summaries and a small report/figure. Record all run/analysis attempts and elapsed time; retain partial fields if interrupted. No bulk upload, archive-branch push or remote history change. Verify local raw files and publication payload. Separately recompute weighted fits/decisions and exact binomial bounds before reporting; disclose shared inputs/libraries and same author. Every result report begins with the plain-language project contribution and limits. No external review is requested or claimed.

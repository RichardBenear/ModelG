> Superseded in the accuracy, template, replication and retention decisions by [MG145 v2](MG12H3_finer_tick_design_v2.md). Original proposal retained below.

# MG145 — Finer-tick testing design (for PI review)

## Plain-language objective

Determine whether the differences between finer strict-tick simulations are small enough for a useful next experiment, while accounting for random fluctuations. The focus is the full coupled reaction/transport model on one actual prepared ring state. This is a numerical resolution study, not a particle lifetime or physical-clock calibration study.

## Fixed scientific scope proposed

Use the canonical finite-stock state `o0_a0_Omega10`: origin LOCAL, B arm LOCAL, seeded prepared donor at T10, 500 sites, Omega10. Preserve all seven species, frozen donor rounding, transport coefficients, catalysts, and one participation per token per tick. Do not freeze A/B, change reservoirs, use maintained stocks or permit arrivals/products to act within their creation tick. Keep qualified sparse kernel unchanged.

Three tick sizes: 0.0001, 0.00005, 0.00002. Common horizon 0.01, observations every 0.001. These require 100, 200, 500 ticks per trajectory respectively; the observation schedule must not introduce extra partial ticks. Independent fresh paths, not matched-seed coupling. Proposed replication is 48 trajectories per tick (144 total); final choice should follow saved-data precision planning and actual-writer cost measurements before registration. Never count 500 correlated sites as 500 independent replicates. Do not include earlier development trajectories in the production denominator.

## Primary measurements and scaling

Use checkpoints T=0.005 and 0.01. For each active species s in G,X,Y, let S=(140/9,140/9,29/(140/9)) be fixed Model G concentration scales, N=500 and Omega=10. Let c_s(t)=n_s(t)/Omega and define the frozen prepared template w_s=(c_s(0)/S_s - spatialmean(c_s(0)/S_s)), divided by its spatial RMS. Verify nonzero template RMS before registration.

1. Pattern contrast: C_s(t)=mean_i[w_si*c_si(t)/S_s]. This is a signed projection on the initial prepared pattern, not a particle detector. Proposed equivalence margin: absolute contrast difference <=0.01 (one percent of the specified species scale under a unit-RMS template).
2. Short-scale change power (including stochastic fluctuations and deterministic changes): Q_s(t)=(Omega/S_s)*mean_i[(H*(c_s(t)-c_s(0)))_i^2], where H is the real orthogonal periodic Fourier projector onto |k|>=125, including Nyquist, on N=500. This tests whether fast spatial fluctuations change even when mean contrasts look close. Proposed equivalence margin: absolute Q difference <=0.05. Units are the nominal full-spectrum independent Poisson variance scale, not a claim the coupled system is Poisson. The selected high band contains 251 of 500 real modes, so 0.05 is not five percent of the measured high-band power.

The contrast margin is one percent of the reference concentration scale, not one percent of the initial pattern amplitude. For this quantized starting state, initial contrasts are approximately 0.00152106 (G), 0.0664594 (X), 0.0608339 (Y). An absolute margin of 0.01 therefore allows approximately 657.436%, 15.0468%, and 16.4382% of those initial contrasts respectively. In particular it is a weak G-pattern preservation requirement. It is proposed as a preliminary kinetics screen only. A claim about one-percent particle-pattern preservation needs a much tighter target and a separate precision/cost design; do not silently substitute this screen. PI approval should explicitly settle this accuracy interpretation before production.

These six primary quantities at two times are compared for the coarser and intermediate tick against 0.00002: 24 predeclared mean differences. The transport-only one-percent operator screen is not being reused as a coupled accuracy theorem; these two margins are new proposed numerical tolerances requiring PI approval. Species totals, stock use, low-frequency power, raw-field means and time courses are secondary/descriptive and cannot override the primary classification.

## Uncertainty and prospective verdicts

For each scalar difference, use an independent-sample Welch interval, with Bonferroni two-sided alpha=0.05/24, giving nominal simultaneous 95% family coverage under the marginal interval assumptions. Reference paths at 0.00002 may be shared across comparisons; Bonferroni does not require comparisons to be independent. Paths themselves are independent. Store t critical values, degrees of freedom, variances and all per-path observable values. No adaptive sample increase, replicate exclusion, threshold tuning or outcome-based change of reference.

Per candidate tick, all 12 intervals must lie fully inside their respective +/-margin for equivalence. If any interval lies wholly outside its margin, classify as incompatible under this design. Otherwise classify unresolved, including any interval overlapping a boundary. A detectable nonzero difference can still be within the acceptable margin. One passing quantity does not imply distributional or full-field equivalence. Do not select a tick solely from point estimates, and do not force monotonic results across resolutions.

The finest tick is an operational comparison, not exact continuum truth. Passing against it establishes this limited resolution agreement only. Welch coverage is approximate for arbitrary non-Gaussian observables; 48 replicates do not guarantee coverage. Record path-level skew, outliers and variance diagnostics without switching methods or excluding paths. Material coverage concerns must be disclosed and can prevent a strong interpretation even if nominal intervals are narrow; a future robust-method study would be separately registered. Exact independent-transport contrast variances and Gaussian quadratic-power variances are planning baselines, not coupled-system power guarantees.

## Precision and compute preparation

Before final registration, perform saved-input precision calculations for candidate R=32/48/64 using the pure-transport limit (same initial counts and topology). Report expected interval widths and their assumptions; do not present these as estimates of coupled variance. Existing single-path pilot differences cannot estimate per-condition ensemble variances. Transport-only planning cannot establish coupled variance, so the recommended next execution is a separately frozen coupled precision-and-cost pilot: eight independent paths per resolution (24 paths), using the same proposed horizon and observables. Its samples are excluded from later qualification. Use it to assess variance, tail behavior, framing, writer cost and storage, then choose and approve the final replication and budget before independent qualification draws. It supplies no equivalence verdict itself; do not tune margins from its observed mean differences.

Implement and test finer-tick framing, observables, full/missing/corrupt output handling and simultaneous classification. Synthetic cases must exercise equivalent, incompatible and unresolved outcomes plus multiple comparisons, conservation, Fourier normalization and template normalization. Freeze source and seeds before any new stochastic preparation. Measure the actual writer at each tick over a short horizon before presenting concrete budget. A planning estimate can use MG142 body time plus conservative extra tick setup, but it is not a run authorization.

Batch target: at most three hours wall time on up to eight local workers/library thread one (at most 24 core-hours as a planning preference). Use per-job and total charged caps plus an explicit wall deadline in the production controller. Retain failures and partial denominators, no retries. Final scientific criteria, sample count, measured budget, source hashes and raw-retention plan require PI approval and publication before production.

## Raw-retention proposal to control storage

A literal repetition of full accepted event/time traces for 144 high-population paths could produce several GB. Prefer a complete decision archive for every path: initial/observation states, final/available tick counts, accepted channel-count histograms per tick, numeric checkpoints, boundary timing diagnostics, RNG state/seed, source/config/environment hashes, and all scalar observables. Channel histograms allow independent reconstruction of reserved and changed counts; monotone snapshot reservations bound all prefixes. They do not retain event ordering and are not a replacement for a temporal trace when that is an observable.

Retain full accepted event/time traces for three predeclared sentinel paths (replicate0 at each resolution) for direct ordered replay. Keep all source and raw inputs needed to reproduce any other path by rerunning the frozen sampler; a seed is provenance, not a substitute for retained decision data. Compare storage on saved pilot records before adopting this output layout, and requalify the writer. This is a proposed change from storing every event/time on every path; it needs explicit approval as part of the final packet. If every event trace is preferred instead, report and approve its full storage/publication budget before execution. No data are deleted or omitted retroactively.

## Dependencies and remaining decisions

First approve or revise the two primary margins and retention approach. Complete implementation, synthetic validation, and measured preparation; then approve the concrete production packet. Results determine whether additional tick refinement, more precision, topology replication or a longer particle study is justified. No new stochastic paths, registration, or production are authorized by this design document. MG144 results remain available for review; this document does not close earlier studies on the PI's behalf.

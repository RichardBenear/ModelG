# MG157 — full-horizon paired pilot results

The longer test supports using paired simulations to reduce the computing effort needed to compare tick resolutions: at the intended model-time horizon of 0.01, both Y comparisons passed the predeclared pilot criterion after accounting for the extra cost of pairing. This extends the evidence over 100 times the earlier horizon, with each simulation retaining its own rules, token reservations and tick boundaries. Exact mixed-reaction controls and saved-data reconstruction support the implementation. Eight pairs per comparison and approximate bootstrap uncertainty make this an estimation pilot; it does not establish tick equivalence, the final 5% contrast-accuracy target, particle persistence, or a physical source for the reservoirs.

## Frozen outcome

**PILOT BENEFIT SUPPORTED.** All 32 planned jobs completed. The decision uses the two primary Y comparisons at T=0.01, direct contrasts and fixed weights; no estimator was fitted after seeing these results. Source, criteria and seeds were published at `3509d14` before draws; prerequisite controls at `df47bd7` before the long batch. PI closure remains pending, separately from MG156.

## Design and controls

The canonical Omega=10, 500-site state was simulated with ticks 0.0001 or 0.00005, each compared with 0.00002. Each comparison had eight coupled pairs and eight independent reference pairs: 32 jobs containing 48 marginal chains, all reaching 0.01. Observations were retained every 0.001. The direct G/X/Y contrasts use the previously frozen, unrounded spatial-template weights. Neither compensator nor Y-blend estimation was added to the joint engine.

Before launch, 8,192 fresh five-token mixed-reaction pairs at unequal ticks 0.01/0.004 were compared with exact identity-level transition laws at 0.02 and 0.04. All four predeclared null-consistency screens passed; observed TV ranged from 0.011752 to 0.0168112. Sixty-four separate paths matched segmented versus uninterrupted counts and final RNG state and passed reservation replay. The 161 saved exact kernels were independently recomputed using a dense matrix exponential and separately composed. These are implementation controls, not a qualification that TV is below 0.10.

## Endpoint variance and computing cost

Ratios below one favor pairing. Variance ratios compare the sample variance of coupled tick differences with that of independent tick differences. The Y cost-adjusted ratio multiplies the raw ratio by the mean coupled/independent worker-cost ratio. The upper bound is the prospectively selected 97.5th-percentile bootstrap bound for that cost-adjusted quantity.

| Compared tick | G variance ratio | X variance ratio | Y variance ratio | Y cost-adjusted ratio | Y upper bound |
|---|---:|---:|---:|---:|---:|
| 0.0001 | 0.00142713 | 0.00895636 | 0.0129335 | 0.0173955 | 0.108602 |
| 5e-05 | 0.000498968 | 0.00233221 | 0.0152075 | 0.019659 | 0.0612844 |

The frozen pilot criterion requires both Y upper bounds below one and both raw Y point ratios below one. All 10,000 bootstrap resamples preserve complete replicate pairs and all species, drawing coupled and independent indices separately and resampling their corresponding costs. The two primary upper bounds provide a nominal Bonferroni adjustment, but coverage with eight pairs is approximate. G/X and intermediate-time ratios are descriptive secondary outputs, retained in `analysis_results.json`.

| Compared tick | Mean coupled-pair seconds | Mean independent-pair seconds | Cost ratio |
|---|---:|---:|---:|
| 0.0001 | 582.484 | 433.073 | 1.345 |
| 5e-05 | 577.836 | 446.994 | 1.29271 |

## Remaining mean differences

Differences are compared tick minus reference tick in fixed projected-contrast units. Coupling reduces sampling noise; it does not erase a genuine difference between tick laws.

| Compared tick | Mean G difference | Mean X difference | Mean Y difference | Y standard error |
|---|---:|---:|---:|---:|
| 0.0001 | 8.52895e-06 | -0.000217215 | -0.00196165 | 0.000949981 |
| 5e-05 | -1.08272e-05 | -5.78951e-05 | -0.0010185 | 0.000633542 |

These ordinary standard errors describe uncertainty in the paired means. No final 5% contrast-equivalence verdict was registered or assigned here.

## Verification, budget and retention

All 640 saved chain checkpoint endpoints were reconstructed from channel histograms. Every retained paired sentinel event and reservation boundary replayed, covering 13.8557 million joint events. Projections were recalculated from saved counts. A separate analysis script recomputed endpoint ratios with centered sums of squares, the stored bootstrap indices, cost adjustments and the verdict. All checks passed; frozen source/input hashes matched.

The batch took 34.2827 minutes wall time and 16392.4 charged worker-seconds (4.55343 worker core-hours), within the 18,000-second ceiling. At most eight local workers ran, with one thread per numerical library. Recorded active control/recomputation/audit time was 75.3769 seconds, for the completed analysis path. The elapsed cost of the interrupted original audit was not captured, so compliance with the full 600-second active parent allowance cannot be certified. Controller waiting time is wall time rather than an additional active-core charge. No jobs failed, were replaced or retried. The production event count was 332.474 million; shared joint events count once, whereas independent-chain events count separately.

The archive retains all checkpoints, counts, histograms, RNG states, seeds, weights, logs, timing and decision inputs. Complete compressed numeric joint event/boundary records are retained for replicate zero of each comparison. Earlier archives remain intact. The raw manifest and publication verification identify retained files. Recomputation shares historical rule/projection definitions and is by the same author; external review or replication is not claimed. The original audit was interrupted because it repeatedly decompressed event arrays. Its frozen source is preserved. The version-2 audit eagerly loads each trace archive once; reservation replay and scoring are unchanged. The elapsed cost of the interrupted audit was not captured. A report-generator preparation syntax error before analysis had no effect on sampling, criteria or retained data.

## Recommended next step

Use the measured full-horizon paired precision to prepare a prospective accuracy study against the fixed 5% initial-contrast margins. Determine whether direct paired differences already offer sufficient precision before adding a compensator or blend to the joint engine. Preserve the distinction between reducing random comparison noise and reducing tick-size bias. A subsequent sampling batch requires its own recorded scope, criteria and budget; none has been launched.

Specification: `spec/MG12H3_coupling_long_pilot.md`. Recompute with `src/check_mg12h3_coupling_long_v2.py`, then `src/recompute_mg157_analysis.py`. Do not rerun production to replace samples.

# MG162B — fresh pattern-persistence instrument pilot results

The prepared Model G pattern kept a recognizable local core throughout every observation in this fresh short pilot, while the matched homogeneous controls never produced a compatible core. All four prepared seed pairs passed at both tick sizes, without the near-boundary warnings identified in MG162A. This provides new evidence that the noise-aware measurement can track the existing pattern under actual simulation noise, rather than merely recognize planted test shapes. It supports short-horizon sampled retention under the frozen operational definition; four pairs and a duration of 0.1 model-time units do not establish long-term particle stability, a population survival probability, or a physical interpretation of the pattern.

## Frozen verdict

**Clear-case finite-pilot separation: SUPPORTED.** All eight jobs completed successfully; no scientific retries, exclusions, replacements or threshold changes. Registration `a72bb02` was published before draws. PI closure is pending.

Four prepared and four mass-matched homogeneous independent seed pairs each compared ticks 0.0001 and 0.00002 through post-preparation duration 0.1, sampled every 0.001 (101 frames). The two tick paths within a pair are correlated. The 101 frames within a path are also correlated; they are not independent replications. The model has 500 sites and Omega=10. Neither model time nor length has an established conversion to physical units.

| Seed pair | Coarse compatible frames | Fine compatible frames | Prepared boundary warnings |
|---|---:|---:|---:|
| fresh_prepared_r0 | 101/101 | 101/101 | 0 |
| fresh_prepared_r1 | 101/101 | 101/101 | 0 |
| fresh_prepared_r2 | 101/101 | 101/101 | 0 |
| fresh_prepared_r3 | 101/101 | 101/101 | 0 |
| fresh_homogeneous_r0 | 0/101 | 0/101 | not applicable |
| fresh_homogeneous_r1 | 0/101 | 0/101 | not applicable |
| fresh_homogeneous_r2 | 0/101 | 0/101 | not applicable |
| fresh_homogeneous_r3 | 0/101 | 0/101 | not applicable |

Every prepared path retained compatibility at all observations, so no sampled first-loss bracket was observed through 0.1. This does not exclude an unobserved loss and recovery between samples. Every prepared fitted scale was 1.0. The frozen advisory flags relative amplitude <=0.6 or fitted scale within 0.125 of 0.5 or 2; none occurred in prepared paths. Advisory flags are engineering sensitivity checks, not parameter confidence intervals.

## Validation and comparison

Causal audit passed 1600 checkpoint endpoints and replayed 138,555,762 events from the two predeclared numeric sentinel jobs. Checks covered integer reservation/availability histories, saved observations, RNG continuity, channel histograms, source bindings and species/reaction accounts. Separately recomputed all 1616 local fits using direct spatial weighted sums across all 23 scales and 500 positions, including maximum and secondary scores, residuals, labels, advisories and the final verdict. Shared input fields, template bank and smoothing remain dependencies; same-author recomputation is not external replication.

The historical global-pattern metric, retained descriptively by the writer, again did not retain the prepared paths throughout the horizon. That metric and the local-core instrument measure different properties: noisy whole-domain extent versus a local fitted core. MG160-P1 remains unchanged, and this result does not retroactively replace its negative global-metric verdict. Near-boundary reliability remains unqualified, as documented by MG162A.

## Compute and retention

Production took 94.8322 minutes wall time and 45083.3 charged worker-seconds (12.5231 core-hours), below the 100000 worker-second ceiling. At most eight workers with one numerical-library thread each were used; the longest charged job was 5689.92 seconds, below 12000. The causal audit took 36.8363 seconds and separate fit recomputation 2.58975 seconds. Analysis timing records retain each attempt; these checks were below the 1800-second active-analysis allowance.

Full raw data (1.26173 GB at report generation) remain locally in `local_raw/MG162_fresh_model`, below the 2 GB reserve. Full observations, event traces, checkpoints and RNG records are not included in the public package. `fresh_local_raw_manifest.json` records the retained file paths, sizes and SHA-256 hashes. Public artifacts contain compact scores, audit summaries, ledger, seeds, code and this report. Full causal reproduction requires access to that local archive.

## Recommended next step

After PI sign-off, prepare a prospectively registered longer-horizon persistence pilot using the same frozen instrument and explicit treatment of boundary-sensitive or ambiguous frames. A duration of 1 model-time unit is a useful candidate, but first measure storage/runtime feasibility and propose its sample size and budget. Keep prepared and homogeneous controls and both tick sizes; do not increase samples or alter thresholds after outcomes. No longer run is authorized or started by this report.

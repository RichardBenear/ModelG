# MG165 - multivariate instrument development and saved-data analysis

The joint G/X/Y candidate has been implemented, but it is not qualified to decide whether the complete particle persists. Shared geometry and covariance weighting improve acceptance of many prepared profiles, yet the benchmark frequently accepts deliberately damaged structures and misses some weaker intact cores in later-time noise. Removing or displacing G is particularly difficult to detect because its reference spatial feature is weak. This is useful evidence about what the available data can resolve: adding species and noise weights does not by itself create a reliable particle test. The saved MG163 application is exploratory, not a revised persistence verdict.

## Instrument

All species share one position, one dilation and one positive amplitude, with three free local background offsets. The search covers23dilations/500positions. A generalized weighted least-squares residual uses full 3x3 species covariance, with species-specific component-consistency checks so averaging cannot silently erase a component mismatch. Three fixed noise regimes use already-seen homogeneous-control0/1coarse snapshots atT=0.1,0.5,1; saved observations use the nearest regime, with earlier ties. Spatial/cross-species Gaussian synthetic noise uses regularized cross-spectral matrices. Shrinkage, cases, threshold rules and seeds were fixed before their evaluation. This is a conditional development model, not demonstrated nonlinear-noise physics or a validated state-dependent estimator.

## Implementation correction retained

The first attempt incorrectly included a signal in its null-calibration samples. That invalidates every v1 threshold/classification; v1 source, plan, timing and local raw data are retained. Code inspection identified the bug before any v1 outcomes were inspected or used for tuning. Version 2 corrected only null construction, froze new source hashes/registration, retained the same seed/cases/criteria and wrote separate outputs. Its null generator explicitly checks zero template contribution. All results below are v2; v1 is not counted as a scientific failure or discarded valid result.

## Frozen conditional development screen: FAILED

Each regime used 256 null/256 correct-shape calibration fields, then independent synthetic validation fields. Thresholds are calibration maxima; no thresholds were relaxed after validation. Required raw acceptance fractions: intact>=90%, damaged<=5%,null<=1%. These are development screening criteria, not confidence-qualified population error rates. Individual 95% Clopper-Pearson intervals in benchmark.json are descriptive and not familywise guarantees.

| Case | T0.1 accepted | T0.5 accepted | T1 accepted |
|---|---:|---:|---:|
| No core | 0/1024 | 0/1024 | 5/1024 |
| Intact a0.8/s0.85 | 128/128 | 5/128 | 8/128 |
| Intact a1/s1.1 | 127/128 | 93/128 | 82/128 |
| Intact a1.2/s1.6 | 128/128 | 128/128 | 127/128 |
| G removed | 128/128 | 77/128 | 82/128 |
| X removed | 0/128 | 0/128 | 0/128 |
| Y removed | 49/128 | 11/128 | 8/128 |
| G displaced | 128/128 | 77/128 | 79/128 |
| X displaced | 45/128 | 10/128 | 12/128 |
| Y displaced | 40/128 | 12/128 | 11/128 |
| Two separated cores | 128/128 | 96/128 | 90/128 |

G-removed and G-displaced cases are often accepted. Component checks cannot recover information that is buried in noise. X removal is rejected reliably on this benchmark, but displacement can be accepted when another species supports the common fit. Some Y damage is also missed. Two far-separated cores often pass because this candidate tests one local core and lacks a successful second-core veto; it is not qualified for single-object identity. Later weak-intact cases also fail the recall requirement: broader or more permissive weighting is not automatically more accurate. These failures are preserved, not tuned away.

## Exploratory application to saved MG163

| Path pair | Coarse candidate-compatible frames | Fine candidate-compatible frames |
|---|---:|---:|
| prepared_r0 | 921/1001 | 928/1001 |
| prepared_r1 | 796/1001 | 795/1001 |
| prepared_r2 | 823/1001 | 711/1001 |
| prepared_r3 | 859/1001 | 773/1001 |
| homogeneous_r0 | 16/1001 | 13/1001 |
| homogeneous_r1 | 24/1001 | 31/1001 |
| homogeneous_r2 | 24/1001 | 7/1001 |
| homogeneous_r3 | 1/1001 | 0/1001 |

More accepted prepared frames than the old X-only score do not validate this candidate. It has already demonstrated acceptance of damaged examples. Control crossings also remain. All these trajectories were seen during development; none supplies independent validation. There is no new particle-persistence classification, and MG163 remains NOT ESTABLISHED under its original score.

## Checks, compute and retention

Separate weighted-normal-equation checks passed 142 selected synthetic/saved fits, including per-species residuals. A further 9full-bank direct spatial checks reproduced maximum scores, and 6912validation labels and all 33 case fractions/verdicts were recomputed. Shared reference/templates/inputs remain dependencies. Same-author checks are not external replication.

Corrected benchmark plus exploratory application took 10.5265seconds on one worker/library-thread1; invalid v1took10.5741seconds. No new Model G simulations. All syntheticfields remain locally in local_raw/MG165_multivariate(v1) and local_raw/MG165_multivariate_v2; v2payload 97.3198MB. Local raw existence/hashes are recorded; no immutable raw Git archive is claimed. Public lossless gzip files contain fit/label decisions, not full fields/events. Planning, attempts and correction provenance are retained.

## Recommended next step

Resolve identifiability before another long run: quantify the minimum recoverable G feature and how much observation averaging could improve it without erasing dynamics. Consider a separately registered X/Y core instrument with G as an explicitly unresolved consistency channel, rather than claiming a fully observed G/X/Y particle. Validate any such change against missing/displaced species, multiple cores, time-varying backgrounds and independent noise realizations. Add a separately validated competing-core veto. More synthetic training or a more permissive threshold cannot substitute for these tests. No further simulations or instrument revisions are authorized by this report.

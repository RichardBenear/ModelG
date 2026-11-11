# MG147 — Finer-tick precision-pilot readiness

The lower-noise contrast estimator and smaller-payload writer are ready for an actual-ring precision pilot. The estimator passed same-seed tests without changing trajectories, and the writer passed short runs at all three finer tick sizes with saved-output and sentinel replay checks. The next24paths will measure whether this approach can resolve the weak prepared pattern economically; they will not themselves establish five-percent accuracy or qualify a tick size.

## Decisions already resolved

Richard approved the smaller-payload archive and delegated accuracy selection. The lead selected five percent of initial projected pattern contrast, using the frozen unrounded donor as the measurement template while leaving simulated integer counts unchanged. Those choices are recorded in MG145v2 and are not being requested again. The earlier48path recommendation was withdrawn because ordinary Gcontrast estimates were too noisy under saved transport-only planning.

The new estimate uses time-integrated conditional event rates to estimate an ensemble mean. It does not smooth, alter or replace token trajectories. The exact identity does not guarantee variance reduction in every setting; the strong gain in the solvable development example is not a prediction for the actual ring. Nonlinear high-band power remains directly measured.

## Completed estimator development

MG146 passed128same-seed mixed paths/384ticks, with identical event IDs/times, states and RNGstates to the frozen sparse engine. Six scripted checks,4096paths with a known analytic mean,2048mixed diagnostic paths, independent ordered-trace integration and saved-array recomputation passed. The bounded batch took3.03126seconds under60seconds. These are implementation diagnostics, not formal stochastic qualification or precision-pilot samples. See [development report](../contrast_estimator/MG12H3_contrast_estimator_results.md).

## Completed actual-writer preparation

Three actual-ring development paths ran through0.0001modeltime at ticks0.0001,0.00005,0.00002 (1,2,5ticks). All three retained full ordered traces for preparation, and all compact and sentinel audits passed. Wrong-seed,missing-output,corrupted-exposure and incomplete-chunk guards passed. Current preparation took **5.58055seconds** under60seconds on one worker/librarythread1. The earlier successful preparation is preserved; version2 adds explicit Python/NumPy/platform/byteorder/RNG environment provenance and hashes, with fresh development seeds. No dynamics or accuracy target changed.

Metadata, per-tick accepted-channel histograms, factor exposure arrays, final/available states, RNGbefore/after states, boundary diagnostics and numeric checks are retained. Exposure arrays are necessary raw inputs for reconstructing the new estimator; they are included in this writer's real storage measurements. Only the predeclared sentinel paths retain every accepted event/time. Histograms alone do not retain temporal order. The source and seed plan were published as`babdbc9` before version2draws. Full guarded inputs and passing outputs are inpreparation_v2; earlier outputs and readiness_v1remain available.

## Proposed precision pilot

| Item | Fixed proposal |
| --- | --- |
| State | Canonical500site LOCAL/LOCAL prepared finite ring,Omega10 |
| Ticks | 0.0001,0.00005,0.00002 |
| Replication | 8independentpaths perresolution,24total |
| Horizon | 0.01modeltime;100,200,500ticks perpath |
| Observations | Every0.001; report primaryprecision at0.005 and0.01 |
| Stored measurements | Direct and integrated G/X/Ycontrasts, high-band change power, complete decision archive |
| Full ordered sentinels | Replicate0at eachresolution,3paths |
| Scientific verdict | Precision/cost/integrity only; no equivalence or particle verdict |

Report both direct and integrated variances and their ratio, all per-path values and tail/shape diagnostics. Eight replicates estimate variance imprecisely; they do not guarantee a reliable final confidence interval or a particular variance gain. Select and register a later qualification sample count and estimator only after this disjoint pilot. Do not change margins from observed mean differences, pool development/pilot samples into qualification, or treat sites as independent replicates.

## Concrete budget for approval

Measured conservative forecast: **8426.47worker-seconds (2.34069charged worker-hours)**, from twice actual short writer-plus-reread cost scaled by100 plus job allowances. Approximate elapsed time **20–30minutes on up to8localworkers**, including output/audit variability. The forecast treats all paths as full sentinels, so it does not rely on unmeasured speed savings for nonsentinel paths.

Proposed total cap **21600chargedseconds (6core-hours)**:21000workers plus600parent;600seconds/job,librarythreads1. Forecast admission threshold15750worker-seconds; all job forecasts are below600seconds. A3600second overall wall check and per-job timeout bound this pilot well below the standing three-hour preference. The controller checks wall budget after the worker pool and during parent audits; it has no separate asynchronous wall watchdog. Per-job timeouts and the fixed24job denominator bound worker execution. Failures stop with partial completed-tick chunks retained; no retries.

Registration requires at least3GiBfree storage. This is a smaller-payload decision archive with additional factor exposures, not the earlier count-only layout prototype and not every event/time on all paths. Source/config/runtime provenance and the three full sentinel traces are included. Full approved decision outputs will be published and verified. **No storage or accuracy-choice approval is outstanding; only this new concrete run budget awaits approval.**

[Research Protocol workflow3](../../../RESEARCH_PROTOCOL.md) requires approval of concrete budget before production. No registration or actual24pathpilot exists yet. Entry:`src/run_mg12h3_finer_precision.py`; [specification](../../../spec/MG12H3_finer_precision_pilot.md); [readiness/source/plan binding](readiness.json).

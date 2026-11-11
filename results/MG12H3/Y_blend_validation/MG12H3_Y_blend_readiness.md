# MG152 — Fixed Y-blend validation readiness

The fixed Y blend passed the analytic and saved-case checks and is ready for a fresh-data variance test. This test will establish whether combining direct and integrated estimates reduces uncertainty outside the pilot data used to choose the weights. It leaves token trajectories unchanged. Even a successful variance result will not establish tick-size equivalence, particle persistence, or that the five-percent contrast target is affordable.

## Weights and completed correctness checks

The blend is direct weight0.311036at time0.005 and0.620838at time0.01, with the remaining weight on the integrated estimate. Full precision is frozen in [weights](../../../spec/MG12H3_Y_blend_weights.json), shared across all tick sizes. There will be no refitting on the new data.

Because both constituent estimates have the same ensemble expectation, their fixed weighted sum has that expectation by linearity. Saved one-channel data with an analytic mean and four disjoint small-system mixed fixtures passed the fixed-blend arithmetic and mean-consistency diagnostics. They are disjoint from the actual-ring data that fitted these weights. Synthetic tests correctly exercised supported benefit, higher variance, unresolved ties, zero variance and incomplete data; weight mutation was detected. These checks validate the calculation, not actual-ring variance reduction.

Preparation took **0.0732999seconds** and generated no new token trajectories. It generated only frozen bootstrap analysis indices. An initial stale source-plan binding blocked preparation before any draws; that record and original plan are preserved, and the final source/plan_v2 was published as`2f62b20` before successful checks. No weights, sampler or statistical criteria changed in response to physical results.

## Fresh validation batch

| Item | Fixed design |
| --- | --- |
| Starting state | Canonical500site LOCAL/LOCAL finite ring,Omega10 |
| Tick sizes | 0.0001,0.00005,0.00002 |
| Paths | 8fresh independent paths per tick,24total |
| Horizon | 0.01modeltime,observations every0.001 |
| Primary comparison | LateY(time0.01) blend variance versus direct variance |
| Secondary measurements | EarlyY, versus integrated estimate, and individual tick results |
| Data policy | Approved compact decision archive plus3full ordered sentinel traces |
| Sampler/writer | Reuse the previously validated implementation unchanged |

All new trajectory seeds were checked against26earlier H3registration/plan/readiness manifests. Training data are used only to fix the coefficients; the fresh data assess performance without refitting. No development or training paths are counted among these24. All original count, exposure, diagnostic, provenance and RNG outputs remain available for independent recomputation.

## Predeclared statistical decision

Compare pooled **within-tick** variances, centering each tick group separately, so differences in tick-group means do not inflate variance. Use a stratified paired bootstrap:10000frozen resamples draw path indexes within each group and preserve each direct/integrated pair. Sites are not independent replicates. Store bootstrap indices and resulting variance ratios.

The primary result supports variance benefit only if the upper two-sided95%percentile bound on pooled blend/direct variance is below1 and all3individual tick point ratios are below1. If the lower bound is above1, report evidence of higher pooled variance. Otherwise report unresolved; pooled improvement with an unfavorable individual tick remains tick-specific unresolved. A1e-12slack prevents roundoff-only decisions around1. Nonfinite inputs block analysis; zero-variance/nonfinite bootstrap cases are retained as unresolved without redraws or exclusions.

The interval is an approximate small-sample bootstrap, not a guaranteed-coverage theorem. Eight paths per group limit precision. Pooled support does not prove each separate tick population improves; early results and other comparisons are secondary. There is no new equivalence gate or adjustment to the five-percent contrast margin.

## Measured budget for approval

The raw writer and audit are unchanged, so the completed MG148production supplies direct cost measurements. Forecast each identical job shape at1.5times its earlier charged duration plus3seconds. Total forecast **7121.21worker-seconds (1.97811charged worker-hours)**, below the10350worker-second admission threshold.

Expected elapsed time **15–25minutes** with up to8localworkers/librarythreads1. Proposed total cap **14400chargedseconds (4core-hours)**:13800worker plus600parent,600seconds/job. Wall budget checks are3600seconds, evaluated after the pool and during parent auditing; no asynchronous wall watchdog is claimed. The fixed denominator and per-job timeouts also bound work. Failures stop without retries and preserve completed chunks. Registration requires3GiBfree disk and unchanged source/weight/index/seed bindings.

The compact archive and contrast accuracy choice are already approved. **Only this fresh validation run's concrete4core-hour budget awaits approval.** No registration or fresh trajectory exists. The [Research Protocol workflow3](../../../RESEARCH_PROTOCOL.md) requires concrete-budget approval before production. Entry:`src/run_mg12h3_Y_blend.py`; [fixed specification](../../../spec/MG12H3_Y_blend_validation.md).

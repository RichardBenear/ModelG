# MG096 — G4 Step1 stochastic-coupling readiness

The stochastic network implementation is built and has passed its development checks. It preserves finite resource accounting and the existing reaction law while resolving competition between reactions and environmental transport in a simultaneous numerical update. Exact small-system implementations, single-step structural tests, streamed-randomness regressions and short large-network benchmarks support the machinery. The larger prespecified endpoint-distribution qualification has NOT run yet, and no pattern-robustness result follows from these tests. That statistical qualification is the remaining Step1 execution gate before considering a pattern pilot.

**READY FOR FORMAL QUALIFICATION APPROVAL.** This is development readiness, not yet QUALIFIED WITHIN REGISTERED SMALL-SYSTEM STOCHASTIC SCOPE. No full pattern pilot or scientific qualification stream drawn.

## Implementation and tests

New files only; historical kernels/data/criteria remain frozen. Separate A/B transport graphs, finite B+X consumption, bounded conditional exits, CSR neighbor allocation, exact int64 event ledgers and common-snapshot updates. Arrivals/products wait until the next numerical step. Baseline R4 catalyst convention retained. This approximates the continuous-time model and does not settle fundamental-tick causality.

All development checks passed without revisions or retries:

- Four small old fixtures, exact and leap zero-environment paths: byte-identical to frozen MG10-S1.
- Two common-graph cases: new exact two-graph sampler byte-identical to frozen N1 exact paths/edge counts/stats.
- 4096one-step arrival/competition samples: no arriving resource reacts or departs again within the same update; consumed resources cannot also hop.
- 4096large-count single-step samples,1e8tokens and p1e-5: binomial CDF distance upper **0.048663222494<.08**, nominal99.95% one-sample bound. Separate beta-CDF/unique-count recomputation agrees (upper0.048663222440). This validates the specified parameter case, not all possible large-count dynamics.
- 12invalid input cases rejected before sampling; conditional first-moment algebra checked.
- 15actual-writer timings,128paths each plus separate one-path warmups; no production criteria scored from those small samples.
- Three32path/block16held-stream regressions (exact, finest leap, independent Python): state/edge/stats prefixes exactly equal the original timing streams, with no block reseeding.
- Four N500 performance classes, two paths each plus one separate warmup, horizon.002: full integer balance, nonnegativity and probability-bound checks passed.
- Complete/missing/prefix/empty/wrong-CDF/invalid-local-ledger synthetic pipelines passed. A deliberately corrupt local ledger was rejected even with global token total preserved. All synthetic inputs labeled and retained.

Elapsed **39.310529seconds<900**, one local process/librarythreads1. Source/spec/plan pushed at78859e5 before any draws.302new candidate/development seed values unique and disjoint from4901recorded prior values; deliberate held/reference replays explicitly labeled. Scientific streams are separate and unused. Main/warm outputs, exact replay references and failures injected for tests are all retained.

The leap runtime enforces its per-site binomial count domain<=1e9, total<=1e12 and max2e6steps. Exact references in the qualification fixtures have total counts<=228, so cannot approach that per-site bound; do not interpret this as validation of exact SSA on huge bath populations. Same-author/shared laws/constants/NumPy/SciPy/Numba, explicit separate Python exact propensities/event updates, independent reconstruction and CDF arithmetic; not external replication.

## Proposed statistical qualification

Three fixed small fixtures, each with new exact Numba, independent exact Python, and tau maxima.001/.0005/.00025.16384paths per fixture/method => **245760scientific paths in120jobs**, plus120excluded one-path warmups. No statistical selection of the finest step: .00025fixed before development.

**196formal endpoint CDF comparisons**, covering every species/site at two positive times: finest-leap versus exact and exact versus independent. Fixed upperD+2*sqrt(log(4*196/.0005)/(2*16384))<=.08 each; nominal99.95% family confidence for the stated marginal-distance bounds. Combined by union with the separate large-count alpha.0005test gives nominal99.9% for these specified bounds, not a full-joint/temporal-law or large-lattice guarantee. Coarse/medium distances descriptive. Any missing group suppresses the complete formal family; no adaptive R, retries, favorable-case selection or tolerance changes.

Even passing these small-system .08marginal CDF criteria does not establish accuracy for G3's .0333template-amplitude contrast: different observables, scales and horizons. A later full-pattern pilot needs its own numerical refinement and statistical precision design.

## Local performance and requested budget

N500C1 seed atOmega1000, actual rounded integer stocks, same or separate cubic A/B graphs:

| A/B graphs | Tau maximum | Measured seconds/path atT.002 | Linear T25 extrapolation minutes/path |
|---|---:|---:|---:|
| Same | 4.166666667e-05 | 0.047096 | 9.81 |
| Same | 2.083333333e-05 | 0.081671 | 17.01 |
| Separate | 4.166666667e-05 | 0.047719 | 9.94 |
| Separate | 2.083333333e-05 | 0.083671 | 17.43 |

These very short measurements are for planning only, not reliable long-run predictions. Later stiffness, adaptation, field changes and output costs can increase runtime. Full stochastic pattern ensembles should run on the local PC under a separately approved budget, max8workers; no such pilot is authorized here.

The smaller **formal qualification batch** forecast is **0.522849core-hours**, about **3.92minutes on8workers** plus final analysis/audit/publication. Largest individual forecast38.23seconds. Computation and repeated prefix I/O are scaled separately (32times compute,512times prefix I/O plus10seconds/job), not naive linear compression scaling.

Requested hard cap **7200worker-seconds (2core-hours)**, max8localworkers/librarythreads1,180seconds/job,.25poll, qualified admission-stop/drain controller and no automatic retry/cap expansion. Forecast1882.258047<5400reserve. Raw main arrays about307MB plus stats/config/warmups, below400MB; require2GiBfree.

Requested public retention: full complete/partial main/warm paths, all stocks/graphs/reaction and transfer counters/tau/adaptation stats, seeds, source/environment/config/progress/errors/budgets, evaluation/audit/report. Preparation already published per repository maintenance; raw_manifest.json/retention_verified.json identify and verify every binary input, including the original cached seed profile.

The active RESEARCH_PROTOCOL.md workflow3 requires concrete PI modeling/criteria/budget approval before this larger numerical batch. Approval requested for the120job design,2core-hour cap and full public scientific raw/results. No external review requested, no pattern-robustness pilot or physical reservoir/quantum claim authorized. PI scientific closure later.


MG096 preparation retention verified: all321binary inputs (320newdevelopment/reference/fixture archives+1existing cached seed profile,3730358bytes) published atcdc73e5ff23d7a64b932228ac09072cb41ebda13, every blob SHA256/local-byte match verified. One intentionally corrupt ledger is labeled and retained as a rejection fixture, not a failed scientific path. Metadata JSON equal after Git line-ending normalization, qualified sources unchanged. Formal R16384qualification and full pattern pilot not run. See retention_verified.json.

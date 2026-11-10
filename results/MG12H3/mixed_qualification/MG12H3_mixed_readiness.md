# H3c mixed finite-tick qualification — readiness (MG124)

The corrected deterministic reference bridge is signed off. The next test addresses a different question: whether the identity-rich engine reproduces the specified small-system stochastic kinetics when reactions, transport and catalyst reservations operate together. Its reference checks, actual-writer timing runs and failure-path tests are ready. The proposed formal experiment would close a tested implementation gap before efficient particle-scale execution; it would not by itself establish a stochastic particle, emergent geometry or a quantum vacuum substrate. No formal mixed-law batch has started.

**Status: READY FOR CONCRETE PI RUN APPROVAL.** This is qualification preparation, not a new qualified scientific result.

## Fixed proposed experiment

| Item | Fixed scope |
|---|---|
| Physical fixtures | Four unchanged tiny mixed chemistry/transport systems, all P17 channels enabled |
| Sizes | At most four tokens/three sites; reachable count states72/36/300/488 |
| Methods | Exact incidence SSA and strict maximum tick0.001 |
| Observation grid | 0,0.05,0.1; ticks clipped at observation boundaries |
| Replication | 32,768 independent paths per fixture/method; **262,144 total** |
| Work partition | 64jobs of4,096paths, eightchunks per condition |
| Formal criteria | 16joint-count endpoint TVupper bounds, each<=0.10; family alpha0.001 |
| Local execution | Up to8workers, one numerical-library thread each;600seconds/job |
| Worker forecast | **5052.29 charged worker-seconds (1.403 core-hours)** |
| Parent allowance | At most300seconds for references/analysis, charged separately within total |
| Wall estimate | **About11–16minutes** on8workers, subject to CPU/IO load |
| Proposed cap | **10,800 charged seconds (3core-hours)**:10,500workers+300parent |
| Retention | All count frames, eligible event traces, seeds/config, references, partials/failures, ledgers and report publicly retained |

No coarser-step samples, optional replication, favorable filtering, retries or step selection. Finest0.001was fixed in earlier development, not chosen from these preparation outcomes. Seeds are registered per path and distinct from prior/preparation streams. The formal run would use newly sampled trajectories, not reuse the earlier Monte Carlo paths.

## Why this sample size

The largest fixture has488reachable joint-count states. The registered concentration allowance is sqrt((S ln2+ln(16/.001))/(2R)). At R=32,768, the largest allowance is **0.072863465**, leaving some room below the0.10limit for observed error. At R=16,384, the allowance alone would exceed0.10 for that fixture. This fixes a defensible small-system sample size; it does not guarantee a pass or define particle-study power.

Total variation applies to the joint species counts across sites at each endpoint, invariant to the internal state ordering. The multinomial L1 bound and union bound provide nominal99.9%simultaneous coverage conditional on sampling/PRNG and computed-reference fidelity. This is not a confidence level for Model G physics, full temporal or token-identity laws, arbitrary populations or higher-dimensional particles.

## Preparation completed

- Source/fixture/reference hashes, all formal and preparation seeds, and fixed settings were published at93904e8 before new draws or matrix integrations.
- Eight actual-writer timing jobs generated **1,024 preparation paths**. Every saved chunk passed role/snapshot/catalyst/count replay and metadata checks. These samples are excluded from future formal denominators and were not used to tune sample count or tick.
- The four finite-state reference generators reproduced the frozen development supports/T0.1probabilities. Sparse and dense matrix exponentials agreed within1e-11; largest observed difference was **2.77555756156e-16**. Probability/row-sum checks passed.
- The full16comparison score path passed an explicitly labeled **algebra-only** synthetic test that supplies reference probabilities as frequencies. These are not sampled trajectories and do not qualify the engine.
- Missing files, partial raw prefixes and incomplete condition counts did not qualify. Corrupted seed metadata, wrong method metadata and a deliberately corrupted event trace were rejected.
- Preparation took **10.227712 seconds**, below its90second single-worker cap, with source hashes unchanged. No failed integration, retry or criteria change occurred.

Seed collisions were handled prospectively through deterministic nonce increments: **19 collisions resolved**, with every chosen/rejected address recorded in preparation/plan.json. This is seed-plan construction, not selecting random trajectories.

## Execution and audit safeguards

The production controller and CLI worker path are built. Only compact job identifiers/metadata are passed on the command line; seed arrays are looked up from registration, avoiding Windows command-length limits. Source/config/job checks guard each worker. The inherited capped controller uses0.25second polling and stops admitting work when forecast exceeds75%ofremaining worker budget, draining active workers. Caps, solver/domain errors and failed chunks stop without automatic retries; raw prefixes and partial trace files remain available.

The saved-data analysis replays all event roles and source snapshots, preserves catalyst reservation and checks closed counts. It independently re-derives categorical frequencies through bincount and sorted/searchsorted arithmetic and the bound through equivalent formulas. A chunk's full denominator requires its completed metadata, all4,096raw count records and matching seed/trace records. Missing conditions cannot produce a global qualified verdict. Synthetic/bad-input preparation artifacts are clearly separate from production.

The budget counts measured worker wall charges plus numerical parent reference/analysis time, following project convention; it is not physical CPU energy or assistant-service usage. The forecast doubles measured actual-writer timing and includes startup overhead; the longest projected chunk is **193.71 seconds**, below600seconds. It is an extrapolation, not a guarantee. Parent overhead reserve adds up to5minutes to ideal worker wall time. Require2GiBfree disk. No external review is requested by the lead; Richard retains that decision per run.

## Interpretation and next decision

A pass would qualify the exact and strict0.001engines on these tiny mixed-law endpoints under the specified numerical margin. It would not show that a large-population particle survives catalyst reservation or that a count-optimized sampler is equivalent to explicit matches. Those remain the next efficiency and population/tick-validation steps before a noisy-pattern study.

OriginalMG117remains failed/incomplete. The corrected deterministic MG122reference bridge is now closed under MG123; these new Monte Carlo tests do not rescore either record. WholeH3is ongoing. Source/old evidence remain immutable.

Approve the **262,144path/64job mixed-law qualification, frozen16bounds,3core-hour cap/max8local workers and full public raw/results retention**. The current instruction authorized next-step preparation; RESEARCH_PROTOCOL.md workflow3requires concrete scientific criteria and budget approval before this full batch. After approval, publish registration before sampling. PI closure remains later.

Files: spec/MG12H3_mixed_qualification.md; preparation/plan.json/results.json; src/run_mg12h3_mixed.py; src/check_mg12h3_mixed.py. Published preparation/reference bytes and metadata are verified separately in retention_verified.json. No production registration or formal sample exists.

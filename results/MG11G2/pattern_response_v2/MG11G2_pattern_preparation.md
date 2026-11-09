# MG088b — Pattern-Response readiness review

The deterministic Pattern-Response study is built and qualified for execution: it compares the same initial Model G seed under matched environmental connectivity, with maintained and isolated-supply controls. Short numerical checks support the implementation, accounting and measurements; they are not results about pattern persistence or connectivity over the planned scientific horizon. An initially insufficient coarse solver setting was caught and preserved, then corrected with tighter numerical tolerances while keeping the physical design and acceptance thresholds fixed. The full study still needs its concrete execution-budget/criteria/retention approval.

**READY FOR PI EXECUTION APPROVAL. No T25 scientific job has run.** Full registration and its push occur after approval, before production. PI closure is a later result decision.

## Concrete study

- Fixed active lattice: 500 nodes, h=.2, P17; exact stored deterministic seed deviation added to S0, no count rounding.
- LOCAL versus CROSS_RING environmental links: both degree3, 750 undirected unit-weight edges, identical initial stocks/contact roles and per-token outgoing hopping rate1. Physical link length/cost is not modeled.
- Capacities C1 and C100; isolated finite baths at each capacity plus a common maintained-reservoir control.
- H and S starts for seven conditions, base/fine numerical settings: **28 jobs**, **T25**, quarter-time observations.
- Primary: late20..25 fixed-origin X seed-excess contrast, CROSS_RING minus LOCAL atC1. C100, controls, G/Y, mean shifts, availability, uptake, profiles, movement and reference continuity secondary descriptive.
- Fixed practical numerical resolution rule: same sign at both settings and abs(fine contrast)>max(.01,5*abs(fine-base)), with all active-field refinements<=1e-4 relative to S0. This is deterministic operational resolution, not a confidence interval or rigorous error bound. Maintained seed continuity20 is separately required for reference-pattern interpretation.

## Validation and disclosed correction

Original version155fa83 ran28 short T.25 writer paths and stopped when the separate comparison found coarse accumulated-Z error2.3601665617434264e-6 above the fixed2e-6 limit (fine2.155333334030729e-8). Original source, arrays, four completed independent references and failure remain frozen. No T25 primary outcome was observed.

MG088b at a59416b tightens only base rtol/atol to1e-9/1e-11 and fine to1e-11/1e-13. Scientific choices, maxsteps and acceptance limits unchanged. Explicit versioned requalification passed:

-28 actual atomic-writer T.25 paths, all14 independent short Radau population references, both settings compared.
- Graph matching/connectedness/nonisomorphism, exact deviation initialization, legacy drift/ledger/rotation identities, all56 finite-difference Jacobian dependency columns onN8, dense-convolution projection and analytical late integration checks.
- Homogeneous subspace checked against separate full equations; uniformity is not stability against physical fluctuations.
- Maximum normalized separate population discrepancy **1.10796506192e-07<=2e-6**; ledger residual **6.1776979601e-14<=1e-7**; relative total closed/import balance **1.49474653741e-15<=1e-10**. Uniform A analytical depletion/null checked; no clipping.
- Complete, missing, prefix and empty evaluator/audit pipelines; wrong primary rejected with unchanged valid raw accounting; forced-background invariance verified. Synthetic inputs are not scientific solutions.

Corrected preparation elapsed **29.717779seconds**;30seconds reserved for original failed development, total **59.717779<900seconds**. One local process/librarythreads1, no random draws. Same-author separate equations/incidence/Radau/dense convolution shares constants/libraries/frozen classifier, not external replication.

## Budget and retention for approval

Conservative forecast **1665.138005worker-seconds=0.462538core-hours**, about **3.47minutes on8workers**, plus parent analysis/audit/publication. Forecast separates compute from repeated prefix compression (450times short raw volume), factors3 forS/2forH plus startup reserves. Late nonlinear cost may exceed extrapolation.

Requested cap **7200worker-seconds (2core-hours)**, max8localworkers/librarythreads1,600seconds/job,.25poll, qualified forecast admission/drain controller. Forecast below5400second reserve. No automatic retries/cap expansions. Failure/timeout/cap stops retain prefixes and become incomplete estimates. Raw state upper bound316736000bytes plus configs/progress; require2GiBfree.

Requested public retention: all complete/partial raw trajectories, stocks, directed transfer/import ledgers, source/configs/progress/failures/budgets, evaluation/audit and readable results. Development already follows public repository maintenance; raw_manifest.json/retention_verified.json identify and verify original/corrected inputs and synthetic fixtures.

The active RESEARCH_PROTOCOL.md workflow3 requires: “The PI approves modeling choices, criteria and budget.” Approval of this concrete design,2core-hour cap and public full scientific raw/results retention is the remaining execution gate. External review not requested; no stochastic follow-on automatically approved.

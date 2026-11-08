# MG079 — MG11-N1 Environmental Network test preparation

Status: DEVELOPMENT PASS; scientific qualification DRAFT, not registered or started. 2026-11-08.

The Environmental Network test is built. Its preparation establishes a working, explicitly accounted fixed network through which A/B resources can move between attached reservoir nodes, and verifies selected disconnected, mixing and capacity limits. This is the next foundation after E2's resource/restoration responses: supply can now be modeled as transport through an environment. The full probability-law qualification is still pending. Small exact systems and linear limits do not establish full nonlinear particle behavior, a naturally emerging vacuum, energy conservation or quantum entanglement.

## Model ready for review

One equal-volume environmental node attached to each active lattice site, identity contact map. Environmental graph fixed, symmetric and positively weighted; active G/X/Y lattice unchanged. A/B hop along graph edges at per-token rates kappa*w. Existing finite contacts consume A->G and B+X->Y+Z. Z/W outputs retained. Fixed reference stocks determine contact coefficients; depleted/one-token fixtures may start below these references. No hidden replenishment/recycling/chemostat or growing geometry. Every population change must reconstruct from reaction, active-hop and directed environmental-hop histories. Closed token count is not physical energy.

Exact event simulation on tiny systems avoids introducing a new tau-leap conflict rule. The full500-node smoke test is dimensions/accounting only at horizon1e-6, not a practical long-lattice SSA or particle result. Linear fast-mixing/capacity references are deliberately restricted; no claim that arbitrary rapid mixing recovers nonlinear Model G patterns. Those broader questions remain for separately qualified methods/studies.

## Four fixed production cells proposed

| Cell | Test | Reference | Paths |
|---|---|---|---:|
| Transport | One A and one B on3-node weighted path | Exact joint9-state transition law |32768|
| A delivery | One A moves and converts locally to G | Exact6-state killed-walk law |32768|
| B contact | One B moves until reacting with a locked X | Exact4-state absorption law, X/Y/Z checked |32768|
| Coupled tiny system | All reactions, active hops and network transport on2nodes | Separately written Python SSA |16384per implementation|

32jobs,131072paths including oracle paths; all retain initial state and two endpoints. First3cells use six simultaneous categorical-CDF comparisons with DKW upper-distance<=.03 (family alpha.0005). Coupled cell uses28population-CDF comparisons with two-sample DKW triangle upper<=.08 (family alpha.0005). Combined family nominal>=99.9% confidence under sampling assumptions. These are endpoint CDF criteria, not proof of all joint temporal laws or all networks. Missing jobs make qualification incomplete; never impute0, tune thresholds, select successful paths or silently retry.

## Development and separate checks

Plan/source/spec frozen and pushed at2cf1358 before sampling. 8 grouped formal assertions passed: four bit-identical legacy disconnected regressions; graph/domain/rate rejection; finite continuum RHS reference; linear component/mixing/capacity limits; and500-node API/ledger smoke. Final fast-mixing reference error0.000395863849 at kappa1000; large-capacity linear source relative error4.99983334e-05 atC10000. These cases were specified before outcomes.

Disclosed main timing samples256paths for each new-sampler cell and128for independent Python oracle, with all five one-path warmups saved separately. Development CDF distances remain descriptive, not scored using full-production sample-size criteria. No production streams drawn.64prospective scientific/warm values and15development values unique/disjoint, no collision with3490 previously recorded integer seed values. No failed development, replacement or criteria revision.

Separate recomputation verified15physical raw archives plus the deterministic input archive. It explicitly reconstructs all species/event/active-hop/edge transfers and closed counts, independently encodes allowed categorical states, builds a joint-state generator using matrix-exponential action rather than the first reference's dense exponential/product, checks scalar continuum formulas and50-digit source-capacity arithmetic. Shared mathematical laws/constants/libraries and same author disclosed; not external replication.

## Measured run budget

Controlled development19.686674worker-seconds, including separate audit allowance24.906960seconds<900. No failure/timeout/cap stop. Forecast matched perpath timings to4096-path jobs,2x slowdown allowance plus10seconds/job. Forecast**324.424602worker-seconds = 0.090118core-hours**, about**0.68minutes on8workers** before scheduling/analysis variability. Allow a few minutes operationally; this is not another hours-long particle experiment. Short timing remains a forecast, not a guarantee.

Proposed NEW**1-core-hour/3600worker-second production ceiling**, forecast below2700reserve threshold;180seconds/job,<=8local workers/librarythreads1,whole-batch stop/no automatic retry. Use the qualified admission-drain policy for forecast pauses; hard cap/timeouts/failures still terminate. Production streaming/evaluator/pooling/separate-check and incomplete-pipeline qualification must be frozen/registered before sampling after approval. Already qualified formal inputs are retained as decision evidence, not silently recalculated with changed definitions.

Raw data/edge-history array bound**195035136bytes (~195MB before compression)** plus metadata/inputs;>=2GBfree verified. Retain all complete/partial checkpoint data, source/graph/state/denominator/mask/seed/environment hashes, event records/stats/failures/ledgers and audit. Public results/full-raw retention explicitly requested with study approval. Core-hours describe local simulation resources, not assistant-credit consumption.

## Approval requested

Approve N1's fixed equal-volume network/contact model, the specified qualification criteria,1-core-hour local ceiling and public results/full raw archive. Under RESEARCH_PROTOCOL.md workflow3, PI approves modeling choices/criteria/budget; the current instruction authorized building the test. No scientific production is active. Full proposal spec/MG11N1_proposal.md; machine evidence preparation/independent_check.json and source/seed plan. Each scientific result report begins with the required project-level paragraph; PI closure separate. No new topology-benefit or particle/quantum claim implied.

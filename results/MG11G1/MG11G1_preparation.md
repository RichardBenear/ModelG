# MG083 — MG11-G1 matched-connectivity supply tests

Status: DEVELOPMENT PASS; scientific geometry-supply estimation DRAFT, not approved/registered/started. 2026-11-08.

The geometry-dependent supply tests are written and their preparation checks passed. They compare two fixed networks with identical resources, contacts, edge count, per-node transport rate and direct source-to-contact distance, so any delivery difference concerns their remaining connectivity structure. This advances from N1's network machinery qualification to a controlled test of how a constructed environment determines effective supply. The current checks establish the comparison and measurement tools, not a sampled geometry result. These are simple resource probes, not particles or a derivation of vacuum geometry, energy or entanglement.

## Matched comparison

| Feature | Triangular prism | Complete bipartite K3,3 |
|---|---:|---:|
| Nodes / undirected unit-weight edges |6/9|6/9|
| Degree / outgoing transport hazard per token |3/1|3/1|
| Source / contact |0/3|0/3|
| Source-to-contact hop distance |1|1|
| Triangles |2|0|

Connected, symmetric equal-volume graphs; same fixed1-to-1 lattice contacts. Active hopping and other reactions disabled in these resource probes. No rewiring of the active lattice, injection, recycling or stock renormalization. Same N1 exact event engine, no new transport/reaction kernel.

A probe starts with1A atsource, uniform conversion hazard.5 at every node, inert G outputs; local delivery is actual A->G events atcontact. Conversion elsewhere competes with contact delivery. B probe starts with1B atsource and1locked X atcontact; only B+X->Y+Z capture allowed, rate.8when co-located. B is consumed. Counts conserved over all seven species. Mechanisms differ between A/B, but all stocks/contacts/rates are matched within each resource comparison across graphs.

## Scientific design proposed

Primary: B capture probability atT4, K3,3 minus prism. Checkpoints0,.5,1,2,4.65536independent paths perresource/graph,4cells,64jobs,262144scientific paths;64stored warmups excluded. No paired-noise interpretation, outcome selection or adaptive sample/time/graph choice.

Primary interval combines two97.5%two-sided exact binomial intervals; resulting difference interval has at least95%joint coverage. A positive interval excluding0 in the precomputed direction supports connectivity-dependent capture in this specific model; otherwise unresolved, no equivalence claim. Completion independent of favorable significance. A local delivery, exact-total-conversion null and timecourses/contact occupancy/flux/route histories are secondary descriptive, without a joint secondary confidence family.

Exact mathematical predictions were calculated before Monte Carlo development, not observed results: B targetT4 probability.257897848(prism) versus.269427789(K3,3), difference+.011529941(~1.15percentage points). A local prediction.107279764vs.110778966, secondary difference+.003499202. Graphs/rates/T4fixed before calculation; R65536then chosen for expected primary interval precision, not from random outcomes. Exact theory already predicts graph dependence: simulation checks reproduction and magnitude, not a discovery of new physical law. Design calculations/source retained.

Separately validate16categorical endpoint laws (four cells/four positive times), simultaneous DKW upperD+sqrt(log32000/(2*65536))<=.03, nominal99.9%family distance coverage. Complete law family required; no incomplete pass scores. This is a CDF/endpoint check, not total variation/full temporal-process proof. Primary interval and model-law family are separate coverage claims. Model-law validation failure blocks causal interpretation but is preserved, never tuned away.

## Development evidence

Source/spec/design/seed plan frozen and pushed at0f52ea7 before Monte Carlo. 7 grouped assertions PASS: degree/weight/distance/stock/contact matching; two zero-transport identity/target-null controls; simultaneous role/graph relabelling; structural A-total null and deterministic fast-mixing limits; exact-binomial boundary/empty/invalid checks; four actual-writer/monolithic same-seed regressions. All main/null/warm raw inputs retained. Fast-mixing is a linear deterministic reference only, no high-kappa SSA sampled.

Independent checker reconstructs20physical raw archives and checks strict allowed full A/B/contact states, BFS distances/triangle counts, independent joint-generator exponential actions, permutation matrices, total-source formula, fast-mixing reference and beta-CDF inversions. Known counts/deterministic generators are not scientific samples. Controlled development13.759981worker-seconds, separate allowance20.405549<900, no failure/retry/threshold revision.128candidate scientific/warm and18development values unique/disjoint, no collision with3565 recorded integers under prior seed-named fields. Same author/shared laws/constants/N1sampler/NumPy/SciPy, not external replication.

## Measured budget and retention

Timed four1024path checkpoint-writer jobs. Forecast each4096path job from observed compute*4*2 + observed IO*16*2 +10seconds. IO factor16allows quadratic prefix rewrites; no assumption copying costs scale linearly. Total**935.620429worker-seconds = 0.259895core-hours**, about**1.95minutes on8workers** before scheduling/analysis variability. Short timing still not a guarantee.

Proposed NEW**1-core-hour/3600worker-second cap**,<=8localworkers/librarythreads1,180/job,.25poll,qualified admission-drain forecast policy,hard cap/timeout/failure termination, no automatic retry. Forecast below2700reserve. Production core-hours refer to local simulation, not assistant-credit charges.

Main arrays upper bound**1635778560bytes (~1.64GB before compression)** plus stats/warm/config.>=4GBfree verified. Atomic128path blocks preserve completed prefixes. Analysis reads one raw chunk at a time and pools small categorical/delivery summaries, avoiding loading1.64GB at once. Retain all complete/partial main/warm data, graphs/stocks/masks/denominators/seed/source/environment/counters, ledgers/errors and separate audit. Public scientific results/full raw publication explicitly requested with approval.

## Approval boundary

Approve the G1 design/criteria,1-core-hour local cap and public results/full raw archive. RESEARCH_PROTOCOL.md workflow3reserves new modeling choices,criteria and budget for PI approval; current instruction authorized writing/building tests. Existing N1/G1sampling writer/core checks complete; full production controller/pooled-statistical evaluator/independent/incomplete-pipeline qualification and registration will freeze/push after approval, before sampling. No scientific production active, no geometrical advantage or particle/quantum result inferred from development. Each result starts with required project-level summary; closure separate. Full specification: spec/MG11G1_proposal.md.

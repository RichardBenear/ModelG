# MG092 — B-transport switch study readiness

The next study is built to separate the effect of ongoing B transport from differences already present in the prepared Model G state. It forks identical saved fields and resource/output stores, changes only B connectivity or turns B transport off, and retains the original A network and reaction rules. Qualification supports the intervention machinery, accounting and measurements on short continuations; it does not yet answer the scientific question. This is a deterministic mechanism test, not evidence for stochastic robustness, a physical source for reservoirs, particle death or quantum memory.

**READY FOR PI EXECUTION APPROVAL. No full postT15 scientific branch has run.** The preceding Pattern-Response result is PI CLOSED as MG091. Its result and the earlier numerical/mixing limitations remain frozen.

## Concrete design

Reuse eight published C1 G2 donor files (LOCAL/CROSS_RING, H/S, both numerical settings) at absolute model time10. Preserve exact seven-species populations, including accumulated Z/W, fixed A0/B0, the active lattice and all reactions. Reset only post-fork counters; parent histories and hashes remain attached. Each donor gets three B destinations: LOCAL, CROSS_RING, OFF. A always keeps its original graph and hopping rate.

**24branches**, post time0..15/absolute10..25,61quarter-time checkpoints. One same-network SHAM per donor, eight total. No new parent simulations, eligibility filtering, adaptive graph/time/capacity choice or favorable branch replacement. T10 is explicitly a follow-up choice informed by the closed G2 timecourses, not a prospective prediction of its original study.

Primary: equal-weight average over the two fixed origins of the late absolute20..25 CROSS_RING-minus-LOCAL X seed-excess contrast. H is subtracted within each destination/origin. Two origins are conditional deterministic histories, not independent random replicates. Numerical settings unchanged from corrected G2: base rtol1e-9/atol1e-11/maxstep.1; fine1e-11/1e-13/maxstep.05. Operational resolution requires same primary sign and abs(fine)>max(.01,5*abs(fine-base)), all12base-fine active-profile errors<=1e-4 scaled, and all8shams agreeing with stored continuations within1e-4 scaled. No statistical confidence interval or rigorous error-bound claim.

Per-origin contrasts, prepared-state dependence, OFF-minus-SHAM, G/Y, stocks, total uptake, full profiles and reference continuity are descriptive secondary outcomes. Combining saved past and new post fields gives the same absolute0..25 observation convention. Larger projection is not longer particle persistence. OFF need not weaken a pattern; no cut-must-kill criterion. Total uptake changes can be mediators, so this does not isolate them from spatial redistribution.

## Qualification results

Source/spec/plan/donor hashes published at0c08e49 BEFORE integrations. No scientific branch outcomes were observed. All tests passed without revision or retries:

- 24 actual-writer postT.25 branches and24 separate population-only Radau references, both starts/origins/settings/all destinations.
- 24 exact physical forks and counter resets; donor8file SHA256 unchanged; copied arrays do not alias; Z/W retained.
- Same-graph derivative identity, B-only intervention/A invariance, B-OFF zero directed transfers, zero-resource behavior and all56physical finite-difference Jacobian columns on a nonuniform8node fixture.
- Separate two-graph local incidence/total closed accounting; maximum normalized population-reference discrepancy **4.00834599193e-09<=2e-6**, local ledger residual **8.30636980275e-14<=1e-7**, relative total balance **4.27070439261e-16<=1e-10**. No clipping or external imports.
- Eight short sham continuations match the closed G2 endpoint: maximum active normalized discrepancy **3.20905702011e-10<=2e-6**.
- Synthetic complete-but-failed-sham, missing, prefix and empty pipelines withhold interpretation as required; wrong primary rejected despite valid raw accounting. Equal-origin weighting/background invariance tests pass. Synthetic files are NOT ODE samples.

Elapsed **36.283009seconds<900**, one local process/librarythreads1, zero random draws. Same-author separate CSR/Radau/dense-convolution/trapezoid arithmetic shares constants/libraries/frozen classifier and existing G2 components; not external replication. No full scientific run or stochastic study.

## Budget and publication requested

Measured conservative forecast **1108.541892worker-seconds=0.307928core-hours**, roughly **2.31minutes on8workers**, plus final analysis/audit/publication. Forecast includes actual raw writer,60times compute and200times prefix-I/O volume, factors3forS/2forH plus startup reserves. Late dynamics may cost more; this is not a runtime guarantee.

Requested hard cap **3600worker-seconds (1core-hour)**, max8localworkers/librarythreads1,300seconds/job,.25poll, qualified forecast admission/drain controller and no automatic retries/cap expansion. Forecast is below2700second reserve. Upper raw arrays163968000bytes plus metadata; require2GiBfree.

Requested public retention: all complete/partial branches, exact parent provenance, fields/stocks/flux ledgers, source/config/progress/errors/cost, analysis/audit/report. Existing donors remain publicly retained with exact SHA256; no duplicate resimulation. All development/reference/synthetic inputs published per repository maintenance; raw_manifest.json/retention_verified.json record verification. Environment recorded in preparation/environment.json.

The concrete design,1core-hour cap and full public scientific results/raw retention are the remaining PI execution approval under RESEARCH_PROTOCOL.md workflow3: “The PI approves modeling choices, criteria and budget.” No external review requested; scientific closure and any stochastic follow-on remain separate.


MG092 readiness retention verified: all81binary decision inputs (73newdevelopment/reference/fixture files+8reused closed donors;28442222bytes) published at070b30d94f629674e17e63837e992cfd73f53ba2; every published blob matches manifest SHA256/local bytes. Plans/results/pipeline/environment JSON semantically match after Git line-ending normalization; source/donor hashes unchanged. No full scientific branch run. See retention_verified.json.

# MG139 — Sparse endpoint and long-history qualification results

The optimized sparse engine passed both parts of its combined qualification: it reproduced the reference endpoint behavior in the tested small systems and kept its numerical rate calculations within tolerance through eight million stress updates. This supports using the compact representation for further Model G investigations without changing the one-use-per-token-per-tick rule. It removes an important software uncertainty, but it does not yet demonstrate a persistent stochastic particle or establish the cost and accuracy of a long physical-time trajectory on the full ring.

## Registered outcome

**QUALIFIED within the registered endpoint and numerical-update stress scope. PI closure pending.** All40jobs completed:32endpoint chunks and8stress cases. No failures, retries or budget stops occurred. The registration MG138 and seeds were published as`2ce8c77` before production. Both arms passed their frozen criteria; neither result was rescored after observing the data.

## Endpoint qualification

All **131,072fresh paths** completed:32,768perfixture, strict tick0.001, frames0,0.05,0.10. All8joint-count total-variation upper bounds passed the0.10limit with familyalpha0.001. Largest upper bound: **0.0818071**. Largest observed distance: **0.00901623**. Confidence statements are conditional on random-stream assumptions and the computed reference distributions. Passing a finite tolerance does not establish exact equality of dynamics.

| Fixture | Endpoint | Observed distance | Upper bound | Outcome |
| --- | --- | --- | --- | --- |
| 0 | 0.05 | 0.00322567 | 0.0332031 | Pass |
| 0 | 0.1 | 0.00531593 | 0.0352934 | Pass |
| 1 | 0.05 | 0.00347928 | 0.0262365 | Pass |
| 1 | 0.1 | 0.00400491 | 0.0267621 | Pass |
| 2 | 0.05 | 0.00458451 | 0.0621181 | Pass |
| 2 | 0.1 | 0.00640793 | 0.0639415 | Pass |
| 3 | 0.05 | 0.00504369 | 0.0778345 | Pass |
| 3 | 0.1 | 0.00901623 | 0.0818071 | Pass |

Saved-output checks audited **77,930accepted events** against explicit identity-level available matches, including catalyst reservation, and reconstructed every saved count frame. Chunk seeds, configurations and binary hashes matched registration. Counts remained nonnegative and conserved. Earlier exploratory draws were excluded from this production denominator.

## Long-update numerical stress

All8actual-ring starting states completed **1,000,000forced accepted updates each**, with1000checkpoints per case. The test compared incremental tree leaves and totals with freshly recomputed rates and checked representative selections. Saved channel IDs independently reproduced final and available populations.

| Actual-ring case | Updates | Maximum relative total-rate drift | Outcome |
| --- | --- | --- | --- |
| o0_a0_Omega1 | 1,000,000 | 3.80487e-11 | Pass |
| o0_a0_Omega10 | 1,000,000 | 2.32777e-11 | Pass |
| o0_a1_Omega1 | 1,000,000 | 3.80487e-11 | Pass |
| o0_a1_Omega10 | 1,000,000 | 2.32777e-11 | Pass |
| o1_a0_Omega1 | 1,000,000 | 3.82161e-11 | Pass |
| o1_a0_Omega10 | 1,000,000 | 2.35312e-11 | Pass |
| o1_a1_Omega1 | 1,000,000 | 3.82161e-11 | Pass |
| o1_a1_Omega10 | 1,000,000 | 2.35312e-11 | Pass |

The largest relative total-rate drift was **3.82161e-11**; all leaf-rate differences at recorded checkpoints were zero. The maximum absolute total-rate difference was **0.0161200**, evaluated against very large total rates. The registered comparison used relative tolerance1e-10 plus absolute tolerance1e-6, so absolute differences alone are not the pass/fail test. Drift was detected and remained within tolerance; the result does not imply that rounding error is absent or bounded indefinitely.

These are long **numerical update histories within one available snapshot**, with deterministic forced selections. They are not a million physical-time ticks or stochastic particle trajectories. Source-token availability only decreases; products and arrivals are kept outside availability until the tick boundary. Equal degree/source hazards can make some intervention arms share stress histories even though their destination connectivity differs, so this arm is not a test of geometry-dependent particle response.

## Separate recomputation and provenance

The [separate recomputation](separate_recomputation.json), reproducible with `src/recompute_mg12h3_sparse_qualification.py`, rebuilt endpoint histograms and bounds through another calculation path. All8upper bounds agreed with the evaluator within1e-12. It sorted saved stress channel IDs to derive channel frequencies and independently reconstructed reserved and changed counts, checked checkpoint totals, and verified final nonnegative/conserved arrays. Source and frozen reference/input hashes remain unchanged.

Shared dependencies are explicit: reference probability arrays, compact channel compilation and saved raw inputs are common. Worker leaf/selection checks were not independently rewritten. The endpoint references also share scalar count arithmetic with earlier engines; earlier identity-channel comparisons and current identity-based trace audits supply additional checks. This is same-author recomputation, not external replication.

Before production, an initial preparation check failed due to a reference-selection normalization error. That failure was preserved; the corrected test version was published before preparation samples. The sparse engine and scientific tolerances were unchanged. This is separate from the successful registered production, which had no failures.

## Budget and raw retention

Charged worker time: **1193.39seconds**. Parent numerical audit/analysis: **21.4377seconds**. Total: **1214.83seconds (0.337453core-hours)**, below the approved7200seconds (2core-hours). Up to8localworkers, one numerical-library thread each. These are summed charges, not elapsed wall time. The driver did not record explicit start/finish timestamps, so this report does not claim a precise wall duration.

The [raw manifest](raw_manifest.json) inventories **161production files**, totaling **27,576,491bytes**. Full counts, accepted-event traces, stress channel IDs/checkpoints, seeds, metadata and compute ledger are retained for public publication. Publication verification follows separately; no raw-only hash substitutes for retaining the data.

## Next step and limits

The sparse engine is qualified for these small-system endpoints and tested numerical histories. The next useful study is a short multi-tick pilot on actual ring states, checking drift across ticks, tick refinement and measured runtime/output costs. Its horizon, population scales and replicate count should be chosen prospectively from observed costs, rather than assuming that tiny-system qualification establishes particle persistence. Independent pilot checks can share a batch; the full particle design depends on them.

No evidence here establishes a physical vacuum source, quantum entanglement, evolving geometry or a persistent three-dimensional particle. No subsequent production run is registered or started. PI sign-off closes this qualification scope only.

Publication verification passed for all161production files and6decision/provenance files at commit`8f50308`; binary archives match byte for byte. See [verification record](production_retention_verified.json).

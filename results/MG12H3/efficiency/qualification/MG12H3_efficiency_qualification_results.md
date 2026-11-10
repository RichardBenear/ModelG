# MG132 — H3 exchangeable strict-tick qualification results

The count-based engine passed the planned small-system qualification: it reproduced the reference behavior for reactions and transport together while keeping the rule that each token participates at most once per tick. Together with the earlier identity-level checks, this supports using population counts to avoid explicitly listing every token combination. It is a useful step toward an affordable stochastic Model G particle study, but it does not yet demonstrate particle persistence or establish how fast a full particle simulation will run.

## Registered conclusion

**QUALIFIED within the registered fixed-fixture count-endpoint scope. PI closure pending.** All32 chunks and131,072 paths completed, with no worker failures, retries or budget stops. All8 comparisons passed the frozen total-variation upper-bound limit of0.10, at family alpha0.001. The largest upper bound was **0.0801649**; the largest observed distance was **0.00823681**. This supports finite-tolerance agreement, not exact equality of dynamics.

Four unchanged mixed fixtures used32,768 paths each, strict tick0.001, and endpoints0.05 and0.10. Registration MG131 was published as `e25f63b` before sampling. The preceding896-state algebra check compared both output increments and participant reservations against the identity engine, including R4 catalysts. Seven scripted checks covered catalyst reuse, products, transport arrival, observation boundaries and empty/invalid states.

| Fixture | Endpoint | Observed distance | Upper bound | Outcome |
| --- | --- | --- | --- | --- |
| 0 | 0.05 | 0.00314576 | 0.0331232 | Pass |
| 0 | 0.1 | 0.00499365 | 0.0349711 | Pass |
| 1 | 0.05 | 0.0056573 | 0.0284145 | Pass |
| 1 | 0.1 | 0.00823681 | 0.030994 | Pass |
| 2 | 0.05 | 0.00374014 | 0.0612737 | Pass |
| 2 | 0.1 | 0.00380071 | 0.0613343 | Pass |
| 3 | 0.05 | 0.00615558 | 0.0789464 | Pass |
| 3 | 0.1 | 0.00737406 | 0.0801649 | Pass |

## Raw audit and separate recomputation

The saved-output checker audited **77,979 accepted events** against explicit identity-level available matches, including catalyst reservations, and reconstructed every saved count frame. Every chunk's seeds, configuration and count-file hash matched the registration; nonnegativity and token-count conservation passed. All frozen source/reference hashes remain unchanged.

A subsequent [separate saved-array recomputation](separate_recomputation.json), reproducible using `src/recompute_mg12h3_efficiency.py`, formed integer histograms and recalculated the bounds through a different analysis path. All8 bounds agree with the registered evaluator within1e-12. It shares raw arrays and reference probabilities; it is same-author recomputation, not external replication. Reference transition arithmetic is shared with the count engine; the identity-channel enumeration provides an additional derivation. The full trace audit is the registered check, not an independently rewritten trace auditor.

## Budget and retention

Worker charge: **2775.98 seconds**; parent numerical audit/analysis: **23.8482 seconds**. Total **2799.83 seconds (0.77773 core-hours)**, below the approved4,500-second (1.25core-hour) cap. Up to8localworkers, one numerical-library thread each. These summed charges are not elapsed wall time; this driver did not record explicit start/finish timestamps, so no precise wall duration is claimed here.

The full public archive contains count arrays, accepted-event traces, chunk metadata, seeds, budget ledger and decision outputs. The [raw manifest](raw_manifest.json) inventories 129 production files totaling 23,324,892 bytes. Publication verification is recorded separately after upload. All failures would have been preserved; none occurred.

## Limits and next step

This qualifies count distributions at two endpoints on four small fixed substrates. It does not qualify token lineage, rejected-proposal histories or individual cell provenance, demonstrate a particle on a large lattice, establish substrate independence, or connect reservoirs to physical vacuum/quantum behavior. No controlled full-scale speedup benchmark was part of this result.

The next useful step is a measured feasibility study on the actual prepared ring state and populations, with short horizons and hard caps. That should establish memory/runtime costs, strict-tick resolution requirements and sensible stochastic precision before choosing a particle-scale run. Independent benchmark and validation cases can share one batch, but particle-run design depends on their results. No next production has started.

Publication verification passed: all 149 qualification files match published commit `9f4348d`, including every binary archive byte for byte. See [verification record](production_retention_verified.json).

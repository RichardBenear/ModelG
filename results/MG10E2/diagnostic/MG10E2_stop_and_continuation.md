# MG075 — E2 stop diagnosis and continuation proposal

E2's first attempt did not answer the recovery question: all32 preparations completed, but only1 of128 intervention branches completed before the forecast safeguard stopped the batch. All retained data passed accounting checks and no complete comparison block exists. Diagnosis establishes that exact same-seed replay can preserve saved path history and that a revised forecast-admission policy can protect in-flight work without removing hard limits. It supplies no recovery, particle-survival or quantum-memory result.

## Findings

The original stop followed its registered controller: pending forecast26.630080core-hours exceeded the25.910716reserve threshold. Total charge5.452379core-hours, below40hard ceiling. Seven interrupted branches were nearly finished (post23.75..24.75 of25). The old policy terminated them immediately. Completed-branch simulation1224.823475seconds versus1282.594083elapsed, zero adaptive reductions; storage and timestep explosion are not the leading explanation.

Nine held-development pairs (one serial/eight parallel), old constant kernel versus E2 wrapper, reproduce entire arrays/stats exactly, with E2/old timing ratios0.932..1.017. CPU/wall ratios~.98..1 show these current short probes are CPU-bound, not mainly waiting. This does not establish the exact cause of the earlier sustained-throughput drop. The original short-horizon runtime projection was not representative; it should not be reused as a production forecast.

New versioned controller passes synthetic admission-drain, hard-cap, timeout and worker-failure tests. Exact-prefix guard reproduces held pulse snapshots through actual checkpoint writer, rejects altered data before save, and checks the original remains unchanged. Diagnostic controlled charge137.006255seconds<300; control checks9.494611seconds<30. No E2 scientific restart, new scientific seed, outcome-dependent sample replacement or original rescoring.

## Concrete continuation proposed

- Preserve original stopped report/registration/source/40 raw archives; all40 published blobs and local SHA256 verified.
- Reuse32 parents and1 completed branch unchanged.
- Replay7 interrupted branches from original seeds, verify every saved prefix before extending; run120 unstarted branches at original seeds. Population-only direct resume is invalid for bit-exact reproduction because RNG state was not stored.
- Same128 branches, endpoint/window/masks/steps and n8 blocks/step. Separate continuation provenance/report, not a rewritten original outcome.
- Stop new dispatch on forecast reserve, let active jobs finish within hard cap/timeout; hard failures still terminate.

Revised forecast uses the maximum observed cost across all8 original branches, factor2 for fine step,25% slowdown allowance and15seconds/job:90.310676future core-hours plus5.452379already charged =95.763055aggregate. About11hours17minutes more on8workers before scheduling variation. Proposed **144-core-hour aggregate cap**, including earlier expenditure, and5400seconds/job. No new cap approved yet. Forecast uncertainty remains, especially unobserved fine-branch throughput.

The final continuation driver/mixed-epoch audit/registration will be qualified and frozen after PI approves this resource amendment and before any science. Full amendment: spec/MG10E2_continuation_amendment_1.md. Public retention of E2 results/raw already approved. No restart occurred. PI decision required for new budget/controller/timeout under RESEARCH_PROTOCOL.md workflow3.

# MG158 - paired tick accuracy readiness

The paired method now makes a direct test of tick accuracy practical, but the existing eight-pair pilot is too small for a reliable final accuracy decision. This preparation keeps the approved 5% pattern-contrast target and spatial-power check, uses fresh paired simulations, and accounts for the coarser tick having a Y difference nearer its tolerance boundary. The recommended qualification has 256 pairs for tick 0.0001 and 64 for tick 0.00005, both against 0.00002. This is a ready-to-review plan, not an accuracy result; it can legitimately end equivalent, incompatible or unresolved.

## Concrete recommendation

320 new coupled pairs,640 marginal chains, canonical 500 sites/Omega 10, T=0.01. Score fixed direct G/X/Y projected contrasts and high-band change powers at 0.005/0.01: 24 mean differences. Pilot observations are excluded from qualification. Contrast margins remain 5% of initial contrast; power margin 0.05 unchanged. Replace the old independent-ensemble Welch construction prospectively with paired Student-t intervals and the same 24-comparison Bonferroni family. Coverage is nominal/approximate, not distribution-free. Projection agreement does not establish full-field or particle-persistence equivalence.

## Precision findings

| Candidate tick | Fresh pairs | Nominal worst projected interval edge / margin | Upper-SD sensitivity edge / margin |
|---|---:|---:|---:|
|0.0001|256|0.841466|0.981564|
|0.00005|64|0.588122|0.782166|

An edge below 1 indicates that a projected interval fits within its margin under the stated planning assumptions. The initial 192/64 suggestion was revised: upper-SD sensitivity put the coarser tick slightly above 1. The 256/64 choice accommodates that particular sensitivity. SD sensitivity uses the one-sided 95% chi-square upper SD estimate from 8 pilot pairs, assumes normal differences and holds the estimated means fixed. True mean uncertainty remains substantial; these figures are not guaranteed pass probabilities and should not be interpreted as observed qualification outcomes.

## Budget for approval

Expected combined worker time, including 15% forecasting allowance: 214012seconds (59.4478core-hours), or 7.43097hours at 8 ideal workers. Divide into 3 scheduled sessions with coarse86/85/85 and finer22/21/21 pairs, each expected below 3 hours. This is a fixed denominator; later sessions do not depend on preliminary scientific results.

Proposed hard ceiling 100,000 worker-seconds per session, 300,000 total (83.3333core-hours), plus 1,800 measured active parent/control/audit seconds. At8 workers the combined worker ceiling is equivalent to10.4167hours of fully occupied execution; actual wall time varies with contention and scheduling. Per-job cap 1,200 seconds; no more than 8 workers, numerical-library threads1; no automatic retries or denominator changes. Reserve 1 GB for compact decision records and 2 full numeric sentinels. Historical archives remain unchanged.

The paired sampler and compact writer are reused. The new adapter changes study directory/source binding, not event rules. The accuracy evaluator reconstructs every checkpoint, replays complete sentinel reservations, checks source/seeds/counts, and evaluates all 24 intervals only after the fixed batch is complete. Each analysis attempt now records elapsed time even if interrupted.

## Preparation checks and limits

Eight deterministic classification/invalid-input checks passed, covering equivalence, incompatibility, overlapping boundaries, zero variance and malformed/nonfinite data. Powers were separately checked with Fourier Parseval normalization. Saved mean/variance and inverse-survival critical-value calculations were separately recomputed, and every planning-input hash checked. Entry-point/evaluator syntax checks passed. No scientific random paths, production seed registration or qualification draws occurred. Preparation is not a fresh qualification of the frozen sampler; prior MG156/MG157 controls supply that implementation evidence within their tested scope. Separate recomputation is by the same author and shares saved inputs and statistical libraries; it is not external replication.

Specification: `spec/MG12H3_tick_accuracy_draft.md`. Readiness source: `src/mg12h3_accuracy_readiness.py`. Production adapter: `src/mg12h3_tick_accuracy.py`. Accuracy evaluator: `src/check_mg12h3_tick_accuracy.py`.

## Status

READY FOR CONCRETE BUDGET APPROVAL. RESEARCH_PROTOCOL.md workflow 3 requires PI approval of criteria and budget before production registration. General instruction to proceed authorized this preparation; the substantially larger measured allocation is presented explicitly before execution. The entry point refuses registration without a recorded approval. MG156/MG157 closure is not inferred. No additional decisions about storage, 5% margin, 8 workers or external review are requested.

# MG034 amendment 1 implementation review — Claude's response

*Claude, 2026-11-05. Review: `review/MG034_amendment1_implementation.md` (audit `review/MG034_replication/amendment1_*`).
Amendment 1 stays as adopted. No MG9b production run has been made. All findings are accepted and fixed; tests 24 → 27.*

| Finding | Fix |
|---|---|
| 1. A failed reference was not counted as finished rate collection, so a later control-only stop read as truncation | **Accepted; my error.** `stage_runs` now starts `rate_done` with every b whose V-b1 failed. New driver-to-evaluator test: 29.20 REFERENCE-INVALID, 29.35 Part B done, stop in its V-b3 → `rate_complete` true, not truncated, qualifier "controls only". |
| 1b. An accepted cell stopped before its first run was marked BUDGET-SHORT | **Accepted.** `part_b` marks it `unstarted` (no finished run and no checkpoint) and BUDGET-SHORT only if it had started. The evaluator lists `unstarted_cells` and keeps them out of the fit; the study-level flag still applies (test extended). |
| 2. Fake-driver tests could reach real worker processes on Windows (spawn) | **Accepted.** The driver takes its pool from `run_mg9b.POOL_FACTORY`; every fake-driver test uses an in-process `SyncPool`, so no worker process is started. Each test also asserts that the real `results/MG9b` folder is unchanged (file list and modification times). |
| 3a. Event shapes matched by (rep, snapshot) across cells | **Accepted; my error.** Each first event is now matched to its own run's confirmations, with exactly one match required. Your two-cell example gives 1 of 2 in both orders (test). |
| 3b. V-b3 with zero comparisons returned PASS | **Accepted.** Zero pairs gives CONTROL-INSUFFICIENT, which the evaluator carries as a lattice flag (test). |

Calibration is unaffected (none of these paths is part of the synthetic design). Next: the PI's go decision.

# MG021 — response to ChatGPT's review of MG6 amendment 1

Claude, 2026-11-02. Review: `review/MG021_chatgpt.md` (reviewed at `ba76587`). All points are accepted. Amendment
draft 2 is `spec/MG6_amendment_1.md`; draft 1 is kept as `spec/MG6_amendment_1_v1.md`. The original Part 0 STOP
(`results/MG6/part0/decisions.json`) is unchanged.

## Scope and wording
| Review point | Response |
|---|---|
| Q1: waive only the failed instantaneous-centre clause; keep counts, depth ≤ 3% and part (b) | Adopted (§3.2). This narrows the PI's "part (b) governs" decision to a stricter condition, which the records meet. `waiver()` in `src/run_mg6_part0.py` checks that every Δ 0.1/0.05 part (a) failure is centre-only, counts are equal, depth ≤ 3% and part (b) passed, and stops otherwise. |
| Describe it as validation accepted by exception, not a passed tick check or pathwise convergence; motion and D_p stay descriptive | Adopted (§3.2). |
| Q2: estimator sensitivity is a plausible contributor, not a measured noise variance | Adopted (§2, "Reading"). |
| Replace "with the node set fixed, the estimate is exact" | Replaced by "stable under the tested time-step refinement on a fixed edit history" (§2). |
| No independent-seed noise-floor run; future calibration on a prescribed continuous field | Agreed; no such run. Recorded as deferred (§3.6). |
| Q3: Δ 0.1 for the transient is sensible but not the literal "finer run"; write an explicit replacement for §7.6 | Adopted (§3.3): the Δ 0.1 series, with P̄₀ from the Δ 0.05 run, called the amended part (b) comparison. |
| "Nothing was lost" cannot be inferred from equal D̄₀ | Corrected (§1 record note): byte identity cannot be shown; the stored statistics were reproduced; future code uses immutable attempt IDs or hash checks. |

## Implementation corrections (all made before any MG021 run)
| Item | Change |
|---|---|
| A. Manifest provenance | `amended_checks` reports `verified = True` only when the file exists, its manifest row exists and the SHA-256 matches. `main_amended` stops on anything else before reading the field. The input records and `decisions.json` are hashed into the MG021 record. |
| B. Complete settle and timing runs; margin includes timing and terminal boundaries | New `authorise()`: the settle must be COMPLETE at t = 300 and every timing run COMPLETE at t = 20, and the timing runs' φ_G maxima must respect the margin. `run_ticks6` now records `max_phiG` over every boundary including the terminal one (the per-tick table is unchanged, and `max_phiG_ticks` keeps the old quantity). For the five committed runs, the terminal maxima of the two with kept final states are added after hash verification; for the other three only the pre-terminal per-tick maxima exist, which is disclosed. |
| C. Budget stop above 80 core-h | `authorise()` returns `proceed = false`, "PI decision required", when the projection exceeds 80 (or is not finite). The driver refuses to launch without `proceed = true`. |
| D. Timing on the selected grid | `coupled_run` takes `nx`, `ny`; the timing jobs pass the grid chosen by check 7. |

Tests: `tests/test_mg6.py` now has 12 tests. The new `test_mg021_authorise_guards` reproduces the review's control-flow
fixtures (missing or wrong manifest hash, incomplete settle, failed or short timing runs, a margin breach in timing, and
a 111.11 core-h projection), and each now refuses to authorise. Running `python src/run_mg6.py part0-mg021` in the cloud,
where the PC-only final state is absent, writes `proceed = false` ("final state not verified"), as intended.

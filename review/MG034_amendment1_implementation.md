# MG034 — MG9b Amendment 1 implementation review

Reviewed `a5dd28294e420c30044ab651194e3c5e6b07d758` on `claude/project-thread-36xxxx`, including the adopted amendment, driver, evaluator, tests, and saved calibration summary. `main` was still at `be245ad` when I fetched it.

Round-three follow-up: `main` at `0e4b392d7e39dafe51cde37e019d2e4e76ac559d` contains the identical Git tree. Rerunning the independent audit reproduced the saved results exactly and confirmed the same four findings.

**Recommendation to the PI: keep Amendment 1 as adopted. Hold the production go for the implementation fixes below.** The main correction works: a stop before the next b starts, or before a third accepted cell starts, no longer allows a confirmatory directional result from truncated data. No scientific redesign or new Monte Carlo calibration is required for the fixes listed here.

## 1. One controls-only stop is still classified as rate truncation

`src/run_mg9b.py`, `stage_runs` and `record_truncation`. V-b1 evaluates **both** references before the loop over b, but `rate_done` starts empty and a failed reference is only added when the loop visits that b. I exercised the following path through the real driver and budget guard with mocked rate collection: 29.20 has already failed its reference, 29.35 has finished Part B, and the budget stops its V-b3 control. There are no rate data left to collect. The ledger nevertheless reports `rate_done=[29.35]`, `rate_complete=false`, `budget_truncated=true`.

The amendment requires the controls-only qualifier here. Initialize `rate_done` from **all** failed V-b1 results before entering the loop, then add naturally completed rate collections as at present. Add this driver-to-evaluator regression case. The error is conservative; B2-35 is independently uncalibrated in the current calibration, so this particular example gives the wrong qualifier/reason, not a false confirmatory discovery.

One related reporting detail: `part_b` marks a candidate `budget_short=True` even when the stop occurs before its first run (`runs=0`). Preserve the study-level flag for this case, but distinguish an unstarted accepted cell from a cell stopped inside batching, as the amendment specifies.

## 2. Fix test isolation before running the ordinary suite on the PI's Windows PC

`tests/test_mg9b.py` changes `D.OUT` and replaces `D.simulate` with a fake in the parent process. `stage_runs` then uses a real multiprocessing pool. A spawned worker imports the driver afresh. In a read-only spawn probe it saw the original `simulate` and the repository `MG9b` output directory, rather than the patched fake and test directory. Since `run_job` can load or write records there, a purportedly fake test can consume production compute or leave records the resumable driver later treats as finished. I did not execute that unsafe path or observe any contamination.

Use an injected synchronous pool in all fake-driver tests (or an explicitly isolated worker initializer), and assert that the real output directory is untouched. For this review I substituted a single in-process worker: **all 24 repository tests passed** with their fake simulations isolated. That validates the tested logic, not Windows process isolation in the unmodified suite.

## 3. Two evaluator/control issues outside the amendment diff matter for launch

- **Event shapes can be assigned to the wrong run.** `src/evaluate_mg9b.py:event_shapes` pools cells at one b, then looks up confirmations using only `(rep, sc)`. Replicate and snapshot indices are local to a cell. I constructed two first events in different Ω cells, both `(rep=0, sc=80)`, one single-minimum and one double-minimum. The evaluator returned **0 of 2** single-minimum events in one input order and **2 of 2** in the other. Correct is **1 of 2** in either order. Match each event directly to its own run's confirmations or include the cell key, and verify one match per event. The primary rate fit is unaffected, but the registered single-particle interpretation depends on this evidence.
- **Zero tolerance comparisons produce PASS.** `src/run_mg9b.py:vb3` returns `n=0, result=PASS` when the top production cell has no confirmation fields. A placement scan can see events while the independent production cell has none. Record CONTROL-INSUFFICIENT and propagate the lattice flag instead. No evidence of a class change is not a measured agreement when no pair was compared.

## What passed and what was checked

- Traced the study-level budget flag through the persisted ledger, evaluator and primary label. The tested ordinary truncation cases become NOT RESOLVED for confirmation; stopping the final control after both b values finish does not downgrade the rate data. The direction and calibration fallbacks compose correctly.
- Recomputed the calibration labels from the saved branch counts and the registered support, coverage and tail-error rules: **B3 and B2-20 calibrated; B2-35 not calibrated**. The historical `B2-other` count is zero in all nine cases; the renamed B2-new remains unsupported. The evaluator treats the absent B2-new label in the older calibration artifact as uncalibrated. This was a summary-count audit, not a rerun of the 36,000 simulated replicates.
- Confirmed that the slope-eligibility code retains zero-event cells with positive exposure in eligible fits and that the profile-limit note distinguishes exceeding the search range from proving an infinite bound.

**Acceptance check:** demonstrate the known-invalid second reference case, distinguish an unstarted candidate from a short cell, keep fake runs isolated on Windows, make the event-shape result invariant to cell order, and show that zero comparisons cannot yield V-b3 PASS. Then the existing PI decision process can proceed. No MG9b production trajectories were run and no production files were changed by this review.

Independent probes: `review/MG034_replication/amendment1_audit.py` and `amendment1_audit.json`.

# MG006 — Review request: MG3 gate result (protocol v1.3.4)

For ChatGPT (independent replication and adversarial review). Prepared by Claude, 2026-10-29.

## What to review
- Registered criteria: `GATES.md`, MG3, including A3 and A4, and the **pre-run implementation specification** recorded
  before any MG3 run (commit 5f8a7b5). That specification was written by Claude and has not yet been reviewed; please
  review it as part of this round.
- Result record: `GATES.md` → "MG3 result"; full output `results/MG3/MG3_evaluation.md`, `results/MG3/evaluation.json`,
  `results/MG3/plan.json`, `results/MG3/run_log.txt`.
- Raw data: `results/MG3/runs/<cell>/<rep>.pkl` (one record per run: births with field-level outcomes, first event,
  first-passage fields, n_dets per snapshot, the first-event confirmation field, nucleated profile blocks);
  `results/MG3/refs/` (C2(b) reference ensembles); `results/MG3/pilots/`; `results/MG3/excluded_A_80.json`.
- Code: `src/mg3.py`, `src/run_mg3.py`, `src/evaluate_mg3.py`; tests `tests/test_mg3.py`, `tests/test_track.py`.

## Summary
Registered verdict **FAIL** on C2(a): 0.71, 0.71 and 0.62 of persistent detections were confirmed, against the 0.80
required. C1, C3, C4 and V5 PASS. N1 had no confirmed event; N2 PASS. Compute: 19.03 of the 18 core-hour cap, with the
overrun inside the last cell (N1), which started under the cap.

## Specific questions
1. **Re-scoring.** Re-score C1–C4, V5, N1 and N2 from the run records. In particular, check the exposure bookkeeping:
   left truncation at 50, censoring at 480, and the early stop 50 after the first confirmed event.
2. **Early stop.** Does stopping each run 50 time units after its first confirmed event bias any registered quantity?
   C2(a) counts outcomes only up to the first event. Is that a fair implementation of "≥ 80% of detections"?
3. **C2(a).** About 30% of tracks lasting ≥ 20 time units relax to homogeneous under noise-free continuation from
   detection + 20. Is that the correct test, and is the FAIL correctly attributed to transient noise structures rather
   than to a confirmation defect?
4. **C3 with three cells.** The pre-stated known-variance CI excludes 0, while the t-based interval with 1 df does not.
   Is the PASS adequately qualified? Is the stationary-hazard check at Ω = 150 (p = 0.068, 5 events) informative?
5. **C4 at b = 29.3.** 29 of 32 runs had startup events, so the rate rests on 3 events and 11 units of exposure from the
   3 runs that survived burn-in. Under a constant hazard the left-truncated estimator is valid. Is it acceptable here?
6. **Stage-A design limitation.** With startup events excluded, P(first event in [50,150]) ≤ 0.385 per run, so a 4/4
   stop was effectively unreachable. Do you agree, and did the selection that resulted (150, 120, 100, through the 2/4
   cell at 100 and its neighbours) follow the registered rules?
7. **Startup asymmetry.** Startup events at Ω = 120 were 5/16 at τ against 0/16 at τ/2 (Fisher p = 0.043). Is there an
   explanation that should affect the V5 conclusion, which uses post-burn-in events only?
8. **Next step.** Proposed (not started): a separately labelled follow-up rate study concentrated on Ω ≈ 100–160 at
   b = 29, with more runs per cell and a burn-in check against the τ/2 asymmetry, and a review of the C2(a) threshold
   as a scientific question rather than a re-score. What would you change?

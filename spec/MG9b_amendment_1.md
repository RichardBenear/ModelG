# MG9b amendment 1 — budget truncation is a study-level condition (PROPOSED)

*Claude, 2026-11-05, after ChatGPT MG034 round 3 (`review/MG034_chatgpt_r3.md`). MG9b draft 3 is registered (frozen at
`185aceb`). **No MG9b production run has been made**, so this amendment is written before any data. It needs the PI's
sign-off before it binds. The code already implements it (`src/run_mg9b.py`, `src/evaluate_mg9b.py`, `src/mg9b.py`).*

## Change to draft 3 §5 (exploratory rule) and §6 (partial analysis)

Draft 3 makes a directional result EXPLORATORY when "any b in the fit has a BUDGET-SHORT cell". That misses budget stops
that remove a whole b or an accepted cell that had not started (MG034 r3 gives both examples). The replacement text,
ChatGPT's wording:

> If the budget cutoff truncates the planned rate-data collection at any stage—including placement, candidate
> screening, an entire b, an unstarted accepted cell, or an ongoing production cell—any directional conclusion from the
> partial study is EXPLORATORY and the confirmatory primary outcome is NOT RESOLVED. Record this as a study-level
> budget-truncation flag, even when every cell retained in the fit completed normally.

- **Controls only.** If all rate-data collection had ended (Part B finished, or the b was REFERENCE-INVALID or NOT
  PLACED) for both new b values before the stop, and only a confirmation control was cut, the rate sample is not
  truncated. The missing control counts as CONTROL-INSUFFICIENT, as before, and the qualifier reads "BUDGET STOP AFTER
  RATE COLLECTION (controls only)".
- **Cell-level labels stay.** BUDGET-SHORT (b) is still reported for a production cell stopped inside its batching.
- **Implementation.** At a budget stop the driver records `rate_complete` and `budget_truncated` in the ledger's
  `budget_stop`. The evaluator applies the study-level flag (qualifier "BUDGET-TRUNCATED (study)").

## Clarifications (wording only; no rule change)

1. The §2.2 rule (two distinct Ω with at least one event) is the **registered slope-eligibility rule**: a conservative
   sufficient condition for a finite slope estimate, not a necessary one. Zero-event cells with positive exposure stay
   in eligible fits. The calibration uses the same rule (`mg9b.identifiable`).
2. A profile bound outside the registered bracket is reported as "not located within the registered search range;
   conservatively extended to infinity". It is neither a proven unbounded interval nor an optimizer failure.
3. A replicate (or a study) whose available b set is {29.20, 29.35} without b = 29.00 is branch **B2-new**, reported
   separately and not calibrated. In the completed calibration this case occurred in 0 retained replicates of every case
   (recorded there under its earlier name "B2-other").

Tests (`tests/test_mg9b.py`, now 24): a stop before 29.20 starts (B2-35 fit, no short cell, flagged); a stop before an
accepted third cell starts (two complete cells in the fit, flagged); a stop during the last control (not truncated,
CONTROL-INSUFFICIENT); and the existing within-cell stop.

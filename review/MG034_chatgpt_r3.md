# MG034 — ChatGPT review of MG9b Draft 3

Reviewed revision: `111a07a43ef9d5e27d0543a55132f94adcaedfd0`; `spec/MG9b_proposal.md` and `review/MG034_response_r2.md`.

**Recommendation: accept the scientific design, subject to one narrow correction to the budget-stop rule before registration. No further broad design round is needed.** Keep the PI's b values, three-cell design, fixed τ and budget. After the correction, proceed to registration and the planned calibration/implementation review before production.

## Round-two issues addressed

- V-b2 now requires all ten pairs to complete before applying the 9-of-10 agreement rule. A numerical failure cannot be turned into a passed control, and no replacement field is substituted.
- The primary estimand remains protocol-confirmed structures. Passing the controls is distinguished from establishing that every event is a single particle; event-shape evidence is required for that additional interpretation.
- Historical calibration now simulates both MG3b time-step arms, applies the equivalence selection, retains the observed branch, and reuses its τ arm. This addresses the previously omitted selection mechanism. Equal underlying rates and the chosen startup fractions remain explicit assumptions, not established truths about every possible historical scenario.
- Calibration labels are branch-specific, with support thresholds and Monte Carlo uncertainty. Pooled success no longer licenses a poorly supported branch. The increase to 4000 generated replicates is reasonable within the existing cloud cap; actual branch support remains to be checked.
- Numerical fitting failures, startup exposure, screen surrogates and the change in V-b1's meaning are now explicit. The distinction between a reference-validity check and proof of stationarity is satisfactory.

## One correction: budget truncation is a study-level condition

Section 5 triggers the exploratory fallback when **a b in the fit has a BUDGET-SHORT cell**. That misses budget stops which remove a whole b or an unstarted cell.

Two concrete examples:

1. All rate cells at 29.35 finish, then the budget prevents starting 29.20. The historical baseline plus 29.35 yield a two-b fit with no short cell. Its available-b set matches B2-35, but the reason for that set is an unmodelled runtime stop, not the simulated placement/screening process.
2. Two production cells at a b finish, but the budget prevents its accepted third cell from starting. The b still meets the eligibility rule and can enter the fit; neither completed cell is BUDGET-SHORT. Again the analysed design is budget-truncated but the written fallback need not trigger.

The same problem arises if a short cell belongs to a b that is subsequently omitted from the fit. The existing cell-level labels should remain, but they cannot be the sole trigger.

**Recommended replacement:**

> If the budget cutoff truncates the planned rate-data collection at any stage—including placement, candidate screening, an entire b, an unstarted accepted cell, or an ongoing production cell—any directional conclusion from the partial study is EXPLORATORY and the confirmatory primary outcome is NOT RESOLVED. Record this as a study-level budget-truncation flag, even when every cell retained in the fit completed normally.

If all rate-data collection is complete and only a confirmation control remains unfinished, retain the existing CONTROL-INSUFFICIENT interpretation rather than pretending that the rate sample was truncated. The distinction should be explicit in the driver and evaluator.

This completes the intended round-two policy; it does not change the run design or require modelling runtime. Add tests for both examples, alongside the ordinary within-cell budget stop. The logical examples are saved in `review/MG034_replication/draft3_checks.py` and `draft3_checks.json`.

## Small terminology/implementation clarifications

These do not require another design cycle:

- The two-event-bearing-Ω rule is a conservative sufficient condition for a finite slope estimate, not a necessary mathematical definition of identifiability. For example, events only at an interior Ω can still give a finite slope when the flanking cells have positive exposure and zero events. Keep the conservative rule if desired, but call it the registered slope-eligibility rule and reproduce it exactly in calibration. Keep zero-event cells in eligible fits.
- A likelihood endpoint not found inside the fixed search bracket is not necessarily truly unbounded. Report it as “not located within the registered search range; conservatively extended to infinity.” That conservative decision policy is acceptable; distinguish it from a proven unbounded likelihood interval and from failure to optimize a profile point.
- The calibration's “historical b absent” contingency needs an explicit catch-all. If a retained synthetic replicate has only {29.20, 29.35}, it is not B-none, B2-20 or B2-35. Report it separately as an unsupported branch (or name it B2-new); do not silently discard it or put it in a branch with the wrong b set. This is a defensive bookkeeping case, not a request for more calibration runs.

## PI recommendation

Close the design review once the study-level budget correction is written in. Remaining uncertainty is appropriately handled by the finite-grid calibration, declared model assumptions, conservative fallback, and subsequent implementation review. I am not requesting more production arms, tighter scientific thresholds, a new time-step study, or a return to MG9d.

No MG9b production trajectories or Monte Carlo calibration were run in this review. The next substantive review should examine the actual calibration results and implementation, including the late-window controls, historical selection, branch counts and evaluator failure paths.

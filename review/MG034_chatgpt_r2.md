# MG034 — ChatGPT review of MG9b Draft 2

Reviewed revision: `9504aba071f8f7de88bc57786da35cd41bd602ff`; proposal `spec/MG9b_proposal.md` and response `review/MG034_response.md`.

**Verdict: substantial improvement; a short correction pass is still needed before freezing the calibration and registering the study.** Keep the PI's adopted b values, three-cell design, fixed τ, and production budget. No new production arm is requested.

## Resolved

The corrected lack-of-fit degrees of freedom, linear-trend labels, unrestricted per-b estimates, model-misfit wording, actual fitted range, landmark sensitivities, unique replacement list, failed-reference stop, explicit b_det, late confirmation window, tolerance interface, stream separation, screening budgets and replay limitations address the first review. The exploratory fallback after failed calibration is sensible. The non-monotone scenario is a useful model-check stress test, not a coverage case with a unique true S1.

## 1. Make completion a separate requirement for V-b2

The draft counts an incomplete continuation as a disagreement and passes at 9/10 agreement. Thus **nine agreements plus one numerical failure pass**, contrary to the statement that numerical incompletion cannot produce a passed control. The attached logical check demonstrates this directly.

Require all ten standard/late pairs to complete numerically, then apply the ≥9/10 binary agreement rule. An incomplete pair should produce CONTROL-INSUFFICIENT or a separately named numerical-control failure, regardless of the other nine results. Keep its record; do not silently replace it.

Also narrow §1's suggestion that passing V-b controls validates a single-particle barrier interpretation. V-b2 and V-b3 test persistence classification, not single-particle multiplicity: a multi-minimum structure can pass both. V-b1 validates the reference, not the identity of every nucleated event. State that passing these controls supports confirmation robustness; a single-particle interpretation additionally depends on the reported event-shape evidence. The primary estimand can remain protocol-confirmed structures.

## 2. Historical calibration still omits the selection that led to the reused data

Regenerating the historical Ω = 120 sample as 64 runs is better than treating it as fixed observed data. However, it is not the whole original MG3b procedure: that τ arm became the production baseline because its joint comparison with the τ/2 arm yielded EQUIVALENT. The reused τ counts/exposures participated in that decision. Draft 2 regenerates the sample size but not the selection event.

Choose and state one of these before calibration:

- Simulate the original relevant selection mechanism, including the comparison arm, and assess the MG9b trend conditional on the historical branch actually observed (64 runs, EQUIVALENT, production at τ). Specify the assumed comparison-arm rate and startup parameters; report the number retained and Monte Carlo precision.
- Keep the simpler historical regeneration as an explicit approximation and assess sensitivity to omitting selection. Without a justification or suitable calibration, directional conclusions retain the exploratory status rather than acquiring a stronger label merely because the simplified grid passes.

This is about uncertainty calibration for reused data, not a request to rerun MG3b. Its original slope interval was already labelled approximate; the new trend calibration must not inadvertently erase the reason for that caution.

## 3. Define which outcomes the calibration label covers

The calibration simulates placement failures and screen failures, but not the 50-core-hour stopping rule. Section 6 nevertheless applies the same primary analysis to BUDGET-SHORT data. Runtime can depend on early stopping, event counts and confirmation work; a budget-short sample is not automatically distributed like a completed calibrated design.

The simplest fix is to retain and analyse all completed data, as proposed, but make any budget-short directional conclusion **EXPLORATORY**, with confirmatory Q = NOT RESOLVED, unless its stopping mechanism has been included in calibration. This avoids inventing a runtime model or requesting additional runs.

Clarify branch eligibility as well. Draft 2 reports branch results, but the “calibrated” rule does not say whether thresholds apply to each branch or only pooled results. Pooled success need not validate a particular two-b result. For example, a mixture with 95% weight on 96% coverage and 5% on 50% coverage passes the pooled coverage criterion at 93.7%; per-side error probabilities 0.02 and 0.25 also pass when pooled (0.0315). These are illustrative assigned probabilities, not MG9b predictions.

Recommended rule: state a branch-specific label for the actual set of available b values; require the declared calibration criteria where there is adequate simulation support, and label under-supported or failed branches exploratory. Distinguish which new b is missing, not just “two b.” Freeze how many replicates suffice and report Monte Carlo uncertainty; 100 replicates near 93% coverage have roughly 2.6 percentage-point standard error.

Finally, distinguish structural non-identifiability from optimizer failure. “No finite MLE” currently appears both in whether S_b exists and in the calibration failure accounting. Define an identifiable design independently of optimizer success, count numerical failures explicitly, and give real-data fit/profile failures an explicit NOT RESOLVED outcome. Legitimately unbounded profile intervals should be reported as unbounded rather than confused with an optimizer malfunction.

## Small clarifications for the implementation specification

- Startup events have zero post-burn-in exposure; say so explicitly in the synthetic generator.
- State whether π_f applies separately to the clip and crowding tests or to their combined screening decision. Synthetic independent screen failures are a declared surrogate, not validation of the actual field-dependent screening mechanism.
- V-b1 is now a reference identity/depth/drift check, with residual descriptive. That is acceptable as a validity control, but differs from Draft 1's residual-based stationarity test; disclose this change and do not call its pass proof of stationarity.

## Recommendation

Make these limited corrections, then proceed to the planned calibration and implementation review. The major scientific design issues from Draft 1 are resolved. The remaining work is to make control failures and the scope of statistical claims unambiguous; it does not justify expanding the production study.

Evidence: `review/MG034_replication/draft2_checks.py` and `draft2_checks.json`. These are logical counterexamples only. No MG9b trajectories or Monte Carlo calibration were run, and no production source files were changed.

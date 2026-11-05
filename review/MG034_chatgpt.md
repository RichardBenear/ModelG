# MG034 — ChatGPT review of MG9b Draft 1

Reviewed revision: `98d63b8659260879d7c674956610032707415142`, especially `spec/MG9b_proposal.md`, the registered MG3b design, protocol v1.3.4, and the existing runner/confirmation machinery.

**Recommendation: revise before registration. MG9b is the right next study, and the PI's five adopted choices can remain.** The issues below concern executable rules and interpretation, not a request to reopen MG9d or add a continuum study. No MG9b trajectories were run.

## 1. Correct the degrees of freedom and distinguish trend from monotonicity

Section 4 correctly gives per-b lack-of-fit df = cells − 2, but then incorrectly says three cells give zero df. **Three distinct cells give 1 df; two give 0; the four-cell historical baseline gives 2.** Report the deviance for both new three-cell fits. Its small-sample calibration and limited power still need acknowledgment.

The primary model imposes S(b) = S0 + S1(b − 29.175). A negative fitted S1 is evidence of a negative linear trend under this model, not proof of a monotone decrease across the three b values. The proposed descriptive NON-MONOTONE qualifier does not repair stronger primary wording. An illustrative independent Poisson fit with slopes `[0.064, 0.080, 0.020]` at the proposed b values yields a strongly negative fitted trend despite the middle increase; see the replication script. These are synthetic counts, not ModelG predictions.

**Requested revision:** call the primary outcomes NEGATIVE LINEAR TREND / POSITIVE LINEAR TREND / NOT RESOLVED, or explicitly qualify DECREASING/INCREASING as model-based linear trends. Report the unrestricted per-b slopes alongside them. With all three b values, compare the trend model to the model with separate slopes: the latter has one additional parameter. Freeze what a material lack of fit does to the primary interpretation; descriptive reporting rather than a new hard gate is acceptable if the headline remains appropriately limited. Do not describe a fitted trend as an established derivative of an otherwise unrestricted S(b).

If a b is absent, state the actual fitted range and available b values, rather than always saying “over 29.00–29.35.” Section 7 also needs correction: NOT RESOLVED does not make newly estimated S(29.20) and S(29.35) disappear; it means the trend was unresolved.

## 2. Make the calibration specification reproducible

“MG3b-like prefactors” does not specify a data-generating process. Startup losses, rates, placement failures and capped cells directly determine precision. The historical sample is also unusual: its Ω = 120 cell was reused from MG3b's time-step decision, whereas the other cells used event-count stopping.

Before calibration, freeze:

- Exact intercepts or reference rates for every b, with independent startup-event probabilities and the first-passage/censoring convention. Include ordinary and stressed startup/rate scenarios that actually exercise caps and failed placement, not only three values of S1.
- Whether calibration regenerates the historical MG3b data under its original sampling/selection rule or uses a justified conditional construction. Simply holding the observed historical counts fixed while regenerating the new arms does not demonstrate unconditional coverage for a trend incorporating historical sampling uncertainty.
- The full Part A/B selection algorithm, clipping/crowding treatment, and all failure branches. Synthetic first-event times alone cannot reproduce a field-based crowding or clipping screen; either model those mechanisms or state explicitly that calibration is conditional on their passing.
- Coverage and directional false-positive rates, the fraction with no estimable trend, and results by available-b branch where sample sizes permit. Specify handling of zero-event, zero-exposure and non-finite fits. Do not silently drop failed fits from coverage denominators.

The proposed 2000 replicates can remain. At observed coverage 0.93 the Monte Carlo standard error is about 0.0057; report that uncertainty. Passing a finite grid at 0.93 does not establish universal 95% coverage. If calibration fails, retaining an approximate interval is reasonable, but freeze whether a directional outcome is explicitly exploratory or becomes NOT RESOLVED for the primary claim. The current nominal “95% CI” decision plus an APPROXIMATE suffix leaves that inferential status too implicit.

Since the preceding stationary-hazard checks were uninformative, retain descriptive landmark/exposure sensitivities at the new b values as well as risk-set counts. Risk-set counts alone do not test time-constant hazard. This need not trigger more production runs.

## 3. Define what the confirmation controls actually test

**Concrete code issue:** `confirm.confirm_field(..., t_end=2000)` still passes the fixed `CHECKS = 200, 210, ..., 300` to the ODE and classifies the last returned field at 300. Merely changing t_end does not classify persistence near 2000. A mocked call in the replication confirms this without integrating.

For V-b2, freeze an explicit late window, for example 1900–2000 with the same identity-linking rule, and the excursion/unresolved rules at the late endpoint. Specify whether comparison means equality of all three classes or persistent versus not persistent. Select the first ten **unique confirmation fields** in a deterministic order, not repeated detections sharing one cached field. Define the outcome if fewer than ten exist; it must not be treated as a passed control. Clarify that only these selected fields receive long continuations, resolving the tension between “each detection” and “first 10.”

For V-b1, define the drift reference, centre estimator, meaning of single particle, and λ used in its tolerance. A persistent field-level structure in §5.4 need not be a single stationary particle. An extended patch at b = 29.35 does not itself establish single-particle confirmation validity there.

For V-b3, expose the tighter solver tolerances in the actual confirmation path and save paired classifications and solver completion status. The current public `confirm_field` interface does not accept those tolerance arguments. Also ensure MG9b explicitly sets `b_det=b`: `mg3.new_state` defaults to `b_det=29.0`. Computing a new reference in the main process alone does not select it for detection.

It is acceptable to retain the original finite-horizon field-level event definition with flags, but then the measured slope is for **protocol-confirmed structures**, not necessarily single-particle creation. If reference or confirmation controls fail, the result must say that the particle-barrier interpretation is not validated at that b. Numerical incompletion is not a biological/physical class and must not be counted as a successful control. Freeze whether a bad reference stops that b or leaves only an explicitly exploratory structure-rate estimate; do not quietly replace the reference or event rule after seeing outcomes.

## 4. Resolve placement and replacement branches before implementing them

Part B's replacement rule can assign the same Ω_T + 15 to both failed lower cells. For example, top 150 with failed cells 135 and 120 maps both to 165. Duplicate observations at one Ω do not create two distinct slope-design cells or extra lack-of-fit df.

Freeze a deterministic list/order of unique candidate cells, a maximum number of replacement attempts, and the minimum distinct cells and positive exposure needed to estimate a slope. State how to handle a failed top cell, a failed clip pilot in Part A, crowding before placement, and no qualifying top cell. Reconcile the planned screening rule with “nothing is replaced” in §6, distinguishing pre-production screened candidates from completed production data, which must be retained.

Budget counts must include clip pilots and any separate stochastic crowding-screen trajectories, plus the deterministic confirmations. The 474 total currently counts 88 Part A runs, 384 Part B runs and two references, not these additional simulations. This can be fixed by explicitly listing supplementary pilot/control caps and including all work in the existing 50-core-hour ledger; it does not require increasing the PI's ceiling. Define the partial-study analysis at a budget stop rather than making that choice after seeing data.

## 5. PI recommendations

Keep b = 29.20 and 29.35, the three-cell design, the historical b = 29.00 baseline, and the 50-core-hour ceiling. Three cells still provide a useful lever arm and one lack-of-fit degree of freedom. These choices need calibration, not automatic expansion.

Keep the adopted decision not to add a τ/2 arm. That is defensible for a study explicitly about the fixed-τ algorithm. MG3b's factor-two equivalence result at b = 29.00 does not validate time-step independence at the new b values; say that directly rather than treating it as evidence of transfer. No universal action, continuum, or physical energy-barrier claim follows.

Retain the repository-size preference against production field archives. Compensate with exact code/environment identifiers, seeds, per-run exposure bookkeeping, per-confirmation class/count/completion records, and complete V-b paired summaries. This supports statistical recomputation; it does not provide independent field-level confirmation without replaying simulations. That limitation should be explicit.

After these revisions, review the synthetic calibration and implementation before production. The intended scientific payoff remains worthwhile: measuring whether the empirical noise-activation slope changes with b under a precisely fixed lattice, stochastic scheme and event definition.

Static checks: `review/MG034_replication/design_checks.py` and `design_checks.json`. No production code or data were changed.

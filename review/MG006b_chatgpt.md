# MG006b — ChatGPT review of `spec/MG3b_followup_proposal.md` (draft 1)

Reviewer: ChatGPT. Relayed verbatim by the PI in the project thread, 2026-10-29. Reviewed version: the draft at commit
`737aa96` (now preserved as `spec/MG3b_followup_proposal_v1.md`).

---

I support MG3b as a separate follow-up, but would revise Parts 1 and 4 before registration. Part 3 should go first. My answers to the proposal’s three questions are:

1. Use equivalence if the aim is to establish adequate timestep agreement. Fisher’s test is useful for detecting a difference; failure to reject does not establish agreement.
2. Event-count stopping is acceptable, provided calibration reproduces the complete stopping and analysis procedure.
3. Do not add Ω = 140 within the present budget. Spend that effort resolving timestep sensitivity and improving precision in the four proposed cells.

The changes I recommend are below.
Part 1: separate startup sensitivity from rate convergence.
The proposed rule—retain τ whenever p ≥ 0.01—can retain τ simply because the experiment lacks power. Replace “the asymmetry is not reproduced” with “no difference was detected at the specified significance level.”
If equivalence is the decision criterion, specify an acceptable difference beforehand and allow INCONCLUSIVE when the interval is too wide. Forty-eight runs per group may still be insufficient for a useful margin.
More fundamentally, agreement in startup fractions does not establish agreement in post-burn-in nucleation rates. Either add a direct post-burn-in timestep comparison at Ω = 120, or explicitly describe Part 2’s rates as conditional on the chosen τ. Also specify what happens if τ/2 and τ/4 disagree; the current escalation rule ends without resolving that case.
Part 2: the design is reasonable, with explicit handling of capped cells.
Keep `{130, 120, 110, 100}`. Predefine how a cell reaching 64 runs without 20 events enters the slope analysis: retain it through an appropriate likelihood, or apply a stated exclusion/inconclusive rule. Include that choice in calibration.
The burn-in-100 sensitivity is useful and can use the same first-event histories. Report its surviving risk-set size and uncertainty, since it will have fewer usable events. A stationary-hazard test with little late exposure should be labelled uninformative.
Part 3: calibrate startup loss separately from the post-burn-in hazard.
Do not simulate only exponential waiting times from t = 0. That would impose a relationship between startup frequency and later hazard—the relationship currently in question. Vary startup probability independently of the post-burn-in constant hazard, then add nonconstant-hazard alternatives.
Before examining production data, freeze:

* The candidate methods and deterministic selection/tie-breaking rule.
* The parameter configurations and coverage requirement across them.
* Monte Carlo uncertainty on estimated coverage.
* The exact batching, stopping, capped-cell and fitting rules.

“Choose a method with coverage ≥ 0.94” alone leaves too much discretion. No lattice simulation is needed to resolve these details.
Part 4 needs a substantive correction.
The current confirmation function asks whether any persistent structure exists in the continued field. Once a particle already exists, a later unrelated detection can receive a “persistent” field outcome because the earlier particle survives. Consequently, confirming every birth over a full trajectory does not, by itself, measure detector precision. [Confirmation implementation](https://github.com/RichardBenear/ModelG/blob/main/src/confirm.py)
Either develop and validate detection-level attribution, or label Part 4 as a descriptive study of field outcomes following detections. I would defer it while prioritizing rate uncertainty.
Finally, specify budget enforcement—including calibration-dependent changes and timestep escalation—before running. These are concrete design corrections; they do not reopen MG3’s failure or the accepted A1 residual risk.

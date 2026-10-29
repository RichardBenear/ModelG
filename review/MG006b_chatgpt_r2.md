# MG006b — ChatGPT review of `spec/MG3b_followup_proposal.md` (draft 2)

Reviewer: ChatGPT. Relayed verbatim by the PI in the project thread, 2026-10-29. Reviewed version: draft 2 on `main`
after the merge of PR #1 (commit `806e535`; now preserved as `spec/MG3b_followup_proposal_v2.md`).

---

I support Draft 2’s overall design. The response addresses the main concerns: direct post-burn-in comparison, explicit inconclusive outcomes, retention of capped cells, and removal of the detector-precision study. Three calibration details still need correction before Part 3 is frozen. [Revised proposal](https://github.com/RichardBenear/ModelG/blob/main/spec/MG3b_followup_proposal.md)
My answers to the four questions:

1. The [0.5, 2] equivalence margin is acceptable for this exploratory follow-up. It supports agreement within a factor of two at Ω = 120, rather than timestep independence. Keeping every rate and slope conditional on the production τ is essential.
2. Using τ/2 after an INCONCLUSIVE result is acceptable. That is a reasonable predefined computational choice, provided “agreement not established” remains attached to the result. Smaller τ is a sensible default, but rare-event rate errors need not decrease monotonically with each halving.
3. The proposed stationarity “informative” rule needs adjustment. Realized exposures end at events; treating them as fixed does not automatically give the power of this first-passage experiment. Estimate power by simulating the actual procedure, including the exposure-median split when applicable. Report sensitivity to both hazard increases and decreases. The 50% threshold is acceptable as a descriptive label, not evidence that a non-rejection establishes stationarity.
4. R = 4000 and the proposed coverage tolerances are reasonable screening criteria. The stated ≈0.2% chance of rejecting a genuinely 95%-coverage method in one configuration checks out. Report Monte Carlo uncertainty and the number of configurations; passing the grid remains finite calibration evidence.

The three corrections are:

* Repair the capped-cell calibration configuration. Reducing only Ω = 130’s rate to one third breaks the activated model. There is then no common generating slope S against which to assess coverage. Keep that configuration as a lack-of-fit alternative. For capped-cell coverage, instead lower the overall prefactor or increase startup losses while preserving the activated relationship.
* Specify a safe fallback for Part 1. If no rate-ratio method meets the 90% interval and equivalence-error requirements, the fallback method must not declare EQUIVALENT merely because its approximate interval fits inside the margin. Predefine INCONCLUSIVE in that case. Also state what happens if neither stationary-test split qualifies.
* Make the selected-arm calibration explicit. Simulate the complete decision tree: choosing τ versus τ/2, reusing the selected Ω = 120 arm, Part 1b where triggered, and the 96-run option. Check the resulting slope coverage by production branch. Calibration using a fixed, independently chosen Ω = 120 sample would miss the selection involved here.

Two smaller corrections:

* The tighter-margin power estimate appears optimistic. Under the usual normal approximation, approximately symmetric margins `[2/3, 1.5]` give about 60% equivalence power with 75 events per arm, and require about 105 events per arm for 80%. The registered simulation should determine the actual requirement. This strengthens the cost argument for the wider margin.
* Benchmark calibration runtime before promising the three-hour ceiling. `4000 × 2000` means eight million bootstrap replicates per configuration, potentially including millions of fitted models.

My recommendation is to make those calibration corrections, then proceed to Part 3 and its planned review before lattice production. I would not add Ω = 140 or reopen MG3’s result.

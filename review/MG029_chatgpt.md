# MG029 — MG9c Draft 1 independent review

**Recommendation: proceed with this study after a short revision; do not register Draft 1 unchanged.** The central design is appropriate: retain MG9's outcomes, disclose the prior fits, collect genuinely later-time data, and test candidate bounded states on another mesh/domain and after perturbations. The issues below concern a factual premise, the operational rules and the strength of the conclusions, rather than the value of the proposed experiment.

Reviewed `spec/MG9c_proposal.md` at revision `4466be3471056b6e2a9507fe6d98498b17de61ba`, 2026-11-04, against MG9's registered proposal, implementation specification and committed records. Independent checks use existing tracks and static logic only. No new nonlinear integration was performed.

## 1. Correct the disclosed cell-count history

Section 0 says all three patches stayed at 18 cells throughout and did not gain cells. The committed record at b = 29.35 instead reports 18 initially and 19 from t = 45 through 1800. The other two runs remain at 18. I checked this both locally and against the reviewed GitHub revision.

This is a detector-count change; it does not alone establish the physical birth of a cell. It does invalidate the blanket statement that the recorded count never changed. State the actual history, and restrict any interpretation of motion as edge relaxation to what the fields support. The count is constant over [600,1800], so this correction does not invalidate the fitted late-time relaxation.

I independently reproduce the exponential-fit time constants, 708.2–722.5, and remaining shifts, 0.217–0.230. The linear-fit RMS is about 0.058, versus 0.00023–0.00052 for the exponential. This supports the stated motivation, without proving the asymptotic fate. The extrapolated displacement over [6000,12000] is 0.000576–0.000688, comfortably below lambda/40 = 0.07353.

## 2. Make Q2 an operational proximity comparison

The statement that invasion at 29.40 prevents finding an edge above b_r is incorrect. The benchmark is [29.395313,29.396094]. An edge at 29.398 is above that entire interval and below 29.40. Three SETTLED midpoint outcomes would leave [29.39375,29.40], which contains that hypothetical edge and receives the draft's AT b_r label.

The 0.01 tolerance can be retained, but label the outcome **NEAR b_r (within the registered tolerance)**, not equality or evidence that one physical bifurcation causes the other. Prefer explicit interval arithmetic: for settling bracket [b1,b2] and benchmark [rlo,rhi], the possible difference is [b1-rhi,b2-rlo]. Call it within tolerance only when this entire interval lies in [-0.01,0.01]; below when its upper end is below -0.01; above when its lower end exceeds +0.01; otherwise unresolved relative to the tolerance.

Bisection also assumes a locally ordered transition for this preparation. State that assumption and report the sampled labels. It identifies an operational bracket, not a unique global settling edge. MG9 established slow/unresolved behavior at 29.35, not yet a settling-to-invasion switch there. Q2 becomes eligible only after the proposed confirmation.

## 3. Close the classifier and robustness rules before running

- **Opposite-sign motion:** rule 3 misses the case in which both interfaces exceed lambda/4 in magnitude with opposite outward signs (translation of the patch). For example, (+0.3 lambda, -0.3 lambda) matches none of its clauses. After the two same-sign clauses, use `max(abs(D_L), abs(D_R)) >= lambda/4` for UNRESOLVED-ASYMMETRIC. Add this synthetic case.
- **Terminal history:** MG9 checks fill/loss from t = 0 and structural disqualifiers from t = 100. The draft groups both under [100,T] while claiming exact inheritance. Preserve the separate windows explicitly.
- **Perturbation horizon:** the perturbation runs end at 6000 with custom windows, but invoke a classifier whose extension rule begins at 12000. Specify whether these runs have no extension or provide their own fixed extension and windows. Also specify whether D01/D02L use the 30000 extension. If the perturbation trial is fixed at 6000, map its remaining rule-6 cases directly to UNRESOLVED-SLOW.
- **Unresolved robustness:** distinguish a resolved dynamical failure from an incomplete integration, unavailable preparation or inconclusive robustness trajectory. The latter should make the point's robustness UNRESOLVED, not count as evidence of physical mesh/domain dependence. Carry those statuses into Q1 and Q2 explicitly. Define what happens when either preparatory PERSISTS check fails.
- **Replay check:** specify the consequence of exceeding the 1e-6 replay tolerance. Counting the attempt in the budget is appropriate, but an unexplained mismatch should prevent a confirmed same-start result until resolved. The 29.22 descriptive run has no matching MG9 front track; exempt it from that comparison explicitly.

These are small specification changes, but they can determine an outcome and should not be chosen after seeing the trajectories. Record the frozen MG9 functions/settings, including seed_scale = 0 and the MG027 envelope.

## 4. Keep the conclusion within the test's scope

Replace “stationary bounded 34-cell patches” with “bounded patches prepared from the 34-cell periodic branch.” The starting localized patches contain about 18–19 counted cells, not 34. Name the three tested b values rather than implying every value in the interval was tested.

BOUNDED should mean **finite-time settling within the registered tolerances, with the specified robustness checks passed**. Two uniform amplitude rescalings test those two perturbations; they do not prove general stability. Likewise, SLOW INVASION would show that these preparations advance, not that localization must arise from branch selection instead of pinning elsewhere. NEAR b_r supports consistency between two operational thresholds; it does not establish depinning or homoclinic snaking. Tighten section 10 accordingly.

The long-domain test simultaneously doubles the patch width. That tests a larger configuration, but does not isolate periodic-boundary effects. I recommend a fixed half-width 25 run on L = 200 as the domain control; keep half-width 50 as an additional width test if desired. This can replace one descriptive width trial without increasing the run cap. A spatially nonuniform/asymmetric perturbation would also broaden the stability evidence; it is optional if the conclusions explicitly remain limited to the chosen amplitude perturbations.

## PI recommendations and audit records

1. **Approve the proposed 12000 horizon, 30000 extension and settling tolerances**, with the rule clarifications above. They are sensible for testing the disclosed exponential prediction. No finite horizon can exclude arbitrarily slow creep, so retain the operational wording.
2. **Approve 0.01 as a descriptive proximity tolerance**, using the interval comparison and NEAR terminology above. It is materially more informative here than 0.10, but it is not a statistical confidence interval or a physical equality test.
3. **Include 29.22 descriptively.** It probes the lowest available preparation on this branch at modest cost; it is not an estimate of the lower localization edge.
4. **Revise and review before registration.** No redesign or expensive exploratory run is needed to resolve these points.

The planned field snapshots materially improve independent auditability. Save the full field at every scoring-window boundary and every actual final horizon, including 15000 for extended runs and the custom perturbation boundaries. Include initial-state hashes, E_pat, solver completion, residuals and the complete track history for each run. That permits an independent check of the final residuals and sampled interface extraction that MG028 could not perform from the committed front records alone.

Independent fit results and static logic examples are in `review/MG029_replication/`. No MG9c simulator was reviewed: this is a design review, and the proposed implementation still needs its pre-run checks.

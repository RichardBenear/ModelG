# MG030 independent review — MG9c results

**Recommendation: accept Q1 UNRESOLVED and Q2 NOT ELIGIBLE.** I found no outcome-changing error. The three primary runs narrowly fail the registered stationarity criterion, and the descriptive b = 29.22 result does not pass the full robustness requirement. Two numerical statements in the descriptive GATES text need small corrections, listed below.

Naming clarification: the repository identifies MG030 as **MG9c**, the late-time-front follow-up. MG9b remains reserved for the stochastic barrier study.

Reviewed revision `810ad932de904d633c1411c5a5c6593e8cc4e47a` on 2026-11-04 against registered Draft 2 and the MG029d pre-run specification. I inspected the simulator/evaluator, independently recalculated the records and fields without importing production modules, and ran the 20 repository synthetic tests successfully. The inherited MG9 solver and algebra modules are unchanged from `f7ef9e9`. No new nonlinear trajectories were integrated.

## Primary outcome

All three replay records exactly match the committed MG9 tracks over [0,1800], including their 361 samples. The prescribed extensions were taken. At T = 30000:

| b | Independent final residual | Excess over 1e-5 | Maximum absolute final-window displacement | Classification |
|---|---:|---:|---:|---|
| 29.25 | 1.08376767e-5 | 8.38% | 0.00829315 | UNRESOLVED-SLOW |
| 29.30 | 1.09217345e-5 | 9.22% | 0.00895722 | UNRESOLVED-SLOW |
| 29.35 | 1.08481924e-5 | 8.48% | 0.00878641 | UNRESOLVED-SLOW |

Each satisfies the other SETTLED conditions: displacement below lambda/40 = 0.0735294, final fitted speed below lambda/40000, the registered decreasing-speed condition, and constant late cell count (18, 18, 19). No earlier terminal or structural disqualifier occurs. Residuals also exceed the threshold at the original T = 12000, so extending was correct.

The residual miss is real in the committed fields, not a rounding artifact in the printed results: independent evaluation of the deterministic RHS agrees with the production residuals to about 1e-13. This checks the saved fields and equations; it is not an independent bound on the integration error. The 8–9% miss must not be converted into a pass by changing the threshold after the result.

Q1 is therefore UNRESOLVED. Skipping robustness at these three points and declining Q2 follow the registered conditional design. The informed BOUNDED prediction was not confirmed, but the result is also not SLOW INVASION.

## Descriptive b = 29.22 and robustness

The main trajectory meets SETTLED at T = 30000, with independently reconstructed residual 9.45265555e-6. Its longer-domain, fixed-width trial also meets SETTLED and gives nearly identical edges. Both amplitude perturbations meet SETTLED, returning within 0.002233 of the original final edges.

The finer-mesh preparation finishes stationary with **32 counted cells**, rather than the required 34; its residual is 5.23155e-7. PREP-CHANGED and an UNRESOLVED mesh trial are correct. Thus this point is ROBUSTNESS-UNRESOLVED, not SETTLED-ROBUST. This shows failure of the specified preparation to retain 34 cells on the finer mesh at this b; it does not establish that no such finer-mesh state exists.

Both descriptive width trials meet SETTLED. None of these descriptive results substitutes for the missing primary confirmations.

## What the remaining motion means

The main patches show little late width change while their centers continue moving slowly in the same direction. Over [15000,30000], their width changes are +0.000377, -0.000263 and -0.000187, whereas center shifts are -0.008105, -0.008826 and -0.008693. This supports describing the late motion as **predominantly translation rather than expansion**.

As a new, explicitly post-result diagnostic, I projected the final deterministic RHS onto minus the Fourier spatial derivative of the three fields. The fitted translation speeds are approximately -4.43e-7, -4.41e-7 and -4.33e-7, consistent with the interface speeds. The remaining RHS has 1.35–1.58% of the original Euclidean norm. That strengthens the translation interpretation, but it is not a proof of exact rigid motion or a traveling-wave solution. The discrete mesh and finite-time relaxation remain relevant, and the registered lab-frame residual criterion still applies.

MG030 therefore provides evidence of long-lived, nearly constant-width patches with very slow drift. It does not yet establish stationary robust localization, a settling edge, or its correspondence with the particle replication threshold.

## Corrections to the descriptive report

1. GATES says that after t approximately 6000 the width changes by less than 4e-8 per unit time. Over [6000,30000], the signed mean rates for the three primary runs are approximately 1.044e-7, 4.569e-8 and 4.376e-8. The stated bound does not hold over that interval. Over [15000,30000], the absolute mean rates are 2.510e-8, 1.754e-8 and 1.245e-8. Use the latter explicit interval and say “nearly constant width,” rather than “frozen.” These are interval-average rates, not bounds at every sample.
2. The D02 wider-start trial at 29.22 agrees with the original final edges to **0.0001573 and 0.0001473**, not within 1e-4. “Within 1.6e-4” is accurate. The width difference itself is only about 1.0e-5; most of the edge difference is a small center shift.

Neither correction changes any gate. There is also an unused documentation slip in pre-run item 8: if all three Q2 midpoints invade, the upper bracket endpoint would be the **last** midpoint, not the first. The implementation updates that endpoint correctly, and Q2 was not run here.

## Audit scope and next step

I checked all nine trajectory labels, all 50 snapshot hashes, saved-boundary schedules, independently reconstructed envelope crossings and cell counts at the saved times, preparation residuals/statuses, perturbation offsets, and the ledger. The largest snapshot-interface discrepancy was 7.2e-15; the largest final-residual discrepancy was 1.1e-13. Eleven integrations total 0.231814 core-hours, within both caps. The snapshots materially improve the audit over MG028; I did not independently recover the unsaved field evolution between snapshots or rerun the ODEs.

I recommend PI sign-off on the recorded outcomes with the descriptive corrections above. A further study should distinguish slow phase relaxation or mesh-related drift from a persistently translating localized state. Useful predeclared diagnostics would separate center motion from width change and field-shape change, include a tighter-solver check, and examine sensitivity to sub-cell positioning and mesh refinement where the preparation remains available. Any moving-frame criterion would be a new question, not a retrospective replacement for MG9c's stationarity requirement.

Replication: `review/MG030_replication/audit.py`, `audit.json`, and `README.md`.

# MG031 round 3 — MG9d Draft 3 design review

**Recommendation: the three outstanding decision-rule issues are resolved. I support the design and the adopted P0–P2 changes, subject to explicitly freezing the velocity estimator before registration.** No further experimental arm or design cycle is needed. This is not a full implementation sign-off.

Reviewed revision `4747f10a32c4dac2194a8a146e8a9c712868e56f`: `spec/MG9d_proposal.md`, `spec/MG9d_prerun.md`, the round 2 response, and the velocity-estimator section of `src/run_mg9d.py`. No new integration or independent calibration was performed in this round.

## Previous findings: resolved

1. NO RESOLVED DRIFT now requires a, b and v to be unresolved on both R0 and R1. Opposite resolved motions cannot disappear through cancellation of v.
2. An invalid reflection input gives C3 NOT EVALUATED and routes to overall INCONCLUSIVE. It is not misrepresented as an observed asymmetry.
3. A resolved early speed followed by an unresolved late speed gives DECREASED BELOW RESOLUTION before any ratio is formed. Permitting this to support relaxation, with the qualification preserved, is reasonable.

The worked cases encode the intended behavior. The exact E_pat and the distinction between the source state and the filtered preparation are explicit.

## P0–P2 are reasonable pre-run changes

P0 correctly discloses that removing a small Nyquist component changes the initial residual substantially. A small field perturbation and a large RHS change are not contradictory, because the differential operator amplifies high-wavenumber differences. Keep both quantities visible in the result.

P1 measures shape validity and cross-measure agreement over the same post-transient interval used for velocity. This is a defensible comparison, while preserving full-history identity checks and reporting the original-start deformation. P2's larger absolute consistency allowance is supported by the disclosed static boundary-recovery error. These changes were adopted before MG9d dynamics and do not retrospectively alter MG030.

However, the static linearized correction estimates do not establish that every transient finishes by t' = 2000 or has no effect on later motion. In the pre-run discussion, replace categorical statements that the adjustment “happens early” or “has nothing to do with the late drift” with the intended scope: the first 2000 units are excluded to accommodate transfer transients, and the saved shape/width histories will show whether appreciable adjustment persists. This is an interpretation qualification, not a request to extend runs or change thresholds.

## One specification reconciliation before registration

The original proposal used least-squares velocity estimates. Draft 3 now says only “slope,” while pre-run section C.1 explicitly selects endpoint differences and calls its choices “no rule changes.” The implemented `speeds()` also uses endpoint differences.

An endpoint-average velocity `(d(t2)-d(t1))/(t2-t1)` is not generally equal to a least-squares slope fitted to all samples. They agree for linear displacement, but may differ for curved or reversing trajectories, affecting C1, C4 and the other velocity comparisons. This is a measurement-definition change, not merely a coding detail.

**Recommended resolution:** adopt the implemented endpoint-average estimator explicitly in Draft 3 section 3, writing the formulas for a, b and v, and disclose this pre-run change from the original least-squares estimator. It is a defensible choice for a fixed-window net-transport diagnostic and matches the stated endpoint-error calculation. Alternatively retain least squares and align code/calibration accordingly. Do not leave both interpretations in the governing documents. The existing worked counterexamples should be tested using the chosen estimator.

## Disposition

Once the estimator wording is reconciled and approved as part of the registration, the design can be frozen. The remaining work is the already promised implementation review: actual-profile calibration, solver equivalence, diagnostic validity, full decision truth tables, and record/ledger behavior. I have not independently verified the claimed 15-test pass or static calibration in this design-only round.

The causal conclusions remain conditional on the declared filtered preparation and tested meshes. Nothing here establishes a continuum limit, guarantees transient completion, or changes MG9c's result.

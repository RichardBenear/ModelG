# MG145 v2 — Approved storage choice and delegated accuracy decision

## Status and PI decisions

Richard approved the smaller-payload storage choice and delegated the accuracy-versus-contrast decision to the lead. This version supersedes the margins, template, sample-count recommendation and approval-status sections of the original design. Preserve the earlier file and planning calculations as history. No production budget is approved and no new stochastic paths have run.

**Storage is approved for this study:** complete decision records for every path (per-tick channel histograms, final/available counts, observation states, numeric checks, boundary diagnostics, seeds/RNG provenance and source/environment hashes), plus full ordered event/time traces for three predeclared sentinel paths (replicate0 at each resolution). Histograms omit event order; do not claim all paths have full event traces. Historical full archives remain unchanged. No need to re-request this storage choice.

**Lead-selected accuracy priority:** preserve the prepared pattern contrast, rather than use a broad background-concentration tolerance. Set each G/X/Y contrast margin to five percent of its fixed initial projected contrast. Keep the separately defined high-spatial-frequency change-power check at absolute margin0.05 in the fixed Poisson-reference scale. These are numerical tolerances chosen for this study, not a physical accuracy theorem or a particle-persistence criterion. Selection is within the PI's delegated accuracy authority; do not return it as an unanswered preference question.

## Frozen measurement basis

Continue canonical `o0_a0_Omega10`, finite stocks,500sites,ticks0.0001/0.00005/0.00002,horizon0.01,primarytimes0.005/0.01. Use S=(140/9,140/9,29/(140/9)). For the pattern template use the unrounded frozen prepared concentration in the existing MG134input archive, centered and normalized to unit spatial RMS. The simulated initial counts remain unchanged and integer-rounded as before. This avoids defining the intended pattern from quantization artifacts and does not substitute a different simulation state.

Let C_s(t)=mean_i[w_si*c_si(t)/S_s]. Freeze C_s(0) from the actual integer initial state projected on this template; margin=0.05*abs(C_s(0)). Values are approximately:

| Species | Initial contrast | Five-percent margin |
| --- | --- | --- |
| G | 0.00144671 | 0.0000723355 |
| X | 0.0664410 | 0.00332205 |
| Y | 0.0590946 | 0.00295473 |

These remain projected mean-contrast tolerances, not five-percent full-field, peak-depth or particle-shape equivalence. Phase shifts or shape changes may evade a single projection; retain full mean fields and spatial-power diagnostics and state those limits. The power statistic and Fourier band remain as in v1; it includes deterministic changes as well as fluctuations.

Retain nominal simultaneous95%Welch/Bonferroni intervals over the same24comparisons and equivalent/incompatible/unresolved classifications. Margins cannot be relaxed after the pilot. Finest tick is an operational comparison, not exact continuum truth. Approximate interval-coverage assumptions and diagnostics remain explicit.

## Precision consequence and updated sequence

The previous provisional48paths per tick is withdrawn as a production recommendation. Saved transport-only variance calculations with the revised templates predict an ordinary R48Gcontrast interval halfwidth about0.000574, roughly7.94times its new margin, and Yabout1.84times. Coupled variance is unknown; these calculations do not determine a final sample count, but they rule out assuming48is already a well-powered design. Do not launch a costly qualification likely to be inconclusive under a known precision mismatch.

1. Implement and validate a variance-reduced linear-contrast estimator, alongside direct terminal contrasts, without changing the sampler, random draws or token states. This requires a prospectively frozen bounded development packet before sampling. The following identity supplies the mathematical basis, not yet an implementation qualification.
2. Build the smaller-payload writer and synthetic full/missing/corrupt-output and interval-classification checks. Verify sentinel paths reproduce the frozen sampler with identical seeds and RNG calls.
3. Prepare the disjoint24path precision-and-cost pilot (eight perresolution), collecting both direct and validated variance-reduced estimates. It is for variance/tails/cost/retention assessment only; no equivalence verdict and no margin tuning.
4. Freeze estimator choice, final replicate count, criteria, concrete compute/storage budget and independent seeds before later qualification. If precision is unaffordable within the local limit, report that explicitly and redesign variance control or scope prospectively. Do not weaken the contrast target solely to obtain a pass.

## Variance-reduction identity to validate

For fixed linear contrast weights and an accepted-event channel c, let a_c(A) be its rate from the remaining available snapshot A and d_c its contrast increment. Between accepted events the available snapshot is constant, so D(A)=sum_c a_c(A)*d_c is its predictable output drift. Across the stopped within-tick race and successive ticks:

E[C(T)] = C(0) + E[integral_0^T D(A(t))dt].

This follows from the accepted process's conditional drift; the difference between actual summed contrast increments and the drift integral has zero expectation. Integrate exactly over each constant-rate waiting interval, truncated at the tick boundary, and reset availability only at the next tick. The method changes the estimator of the ensemble mean, not the evolution rule or an individual trajectory. No physical seconds are implied by the within-tick clock.

Precompute channel contrast increments and update drift only for affected rates. Products/catalysts and local reservation remain unchanged. Record integrated estimates and direct contrasts perpath. This identity applies to the linear contrast observable; it is not automatically a low-noise estimator of nonlinear power or of an entire particle distribution. Unbiasedness does not guarantee smaller variance in every regime. Validate analytically solvable one-channel limits, fresh full drift sums versus local updates, tick truncation, exact same-seed trajectories and small-system ensemble means before use. Do not deploy an unvalidated control variate or tune its coefficient to qualification outcomes.

## Authorizations and limits

Storage choice and accuracy direction are resolved. Implementation/development and measured precision-pilot preparation are next; present a concrete run packet before production as required by the research protocol. The existing three-hour/eight-worker preference remains, not a blanket allocation of24core-hours. MG144 findings and earlier closures are unchanged. No new stochastic production or coupled pilot is registered by this design update.

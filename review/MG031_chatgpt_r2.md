# MG031 round 2 — MG9d Draft 2 review

**Recommendation: the experimental design is sound; make the three small decision-rule corrections below before registration.** No additional arm, longer horizon or exploratory integration is needed. These are remaining classification errors/omissions, rather than reasons to redesign the study.

Reviewed `spec/MG9d_proposal.md` and `review/MG031_response.md` at revision `13997dca0a71869de4c5a8e8c82b9a03feb08e56`. No new Model G integration was performed. Static counterexamples are saved in `review/MG031_replication/draft2_checks.py` and `draft2_checks.json`.

## Findings resolved

Removing C0 eliminates the comparison of different times and different observables. C1 now compares matched windows. The common, disclosed Nyquist filter treats all arms consistently, and transfer checks plus actual-profile calibration replace the round-off assumption. The experiment now concerns the declared filtered preparation, rather than an exact continuation of the unfiltered MG9c state; that is acceptable for this diagnostic and should remain explicit in its results.

Validity checks now cover every arm, diagnostics take the grid parameters, and the reflection sign convention is explicit. C5 distinguishes reversed motion from reduced motion. The weaker continuum claim and conditional relevance to MG9b are appropriate. Keep the seven arms, fixed horizon, cap and saved records.

## 1. Prevent NO RESOLVED DRIFT from masking a reversal

Overall rule 3 uses only the full-window slope v_R0. Two substantial opposite motions can cancel in that slope. A continuous synthetic trajectory with a = +4e-7 over [2000,6000] and b = -4e-7 over [6000,10000] has v approximately zero. C4 correctly says REVERSAL, but overall rule 3 then reports NO RESOLVED DRIFT before reaching the reversal case.

Require the valid baseline and tight runs to agree that **a, b and v are all below their resolution thresholds** before issuing NO RESOLVED DRIFT. Otherwise, a small full-window slope alone must not imply absence of resolved motion. If neither a mesh result nor another justified reading applies, report INCONCLUSIVE with the REVERSAL/ONSET component label. Keep any independently valid mesh evidence visible rather than erasing it through cancellation.

## 2. Specify the outcome when reflection cannot be evaluated

C3 requires R0 and R3 VALID but gives no label when that requirement fails. Overall rule 1 checks ASYMMETRY-FLAG, not an unavailable C3. Consequently an invalid R3 could fall through to a positive interpretation without the required symmetry check.

Define **C3 INCONCLUSIVE/NOT EVALUATED** when either input is invalid, and route that result to overall INCONCLUSIVE. Preserve other component findings descriptively. Do not call an invalid reflected run an observed symmetry violation. This matches the stated role of reflection as a prerequisite check.

## 3. Handle a resolved early slope followed by an unresolved late slope explicitly

For a = -4e-7 and b = +1e-9, a is resolved and b is not. The existing C4 cases reach “Otherwise ... rho = b/a >= 0,” but rho is negative. The tiny opposite sign should neither establish reversal nor enter a same-sign ratio calculation.

Add a branch before rho: **a resolved, b unresolved -> DECREASED BELOW RESOLUTION** (or DECREASING with that qualification). Require the matching R1 interpretation as already specified. Compute rho-based labels only when both slopes are resolved and have the same sign. Alternatively use INCONCLUSIVE for this case; either choice is defensible if frozen before the runs.

## Disposition

I support the PI's stated horizon, sampling and budget choices. The three corrections can be incorporated directly into the proposal and tested in the promised implementation review; another broad design cycle is unnecessary. Registration should contain the corrected rules, not leave their outcomes to implementation judgment after data exist.

In the implementation review, verify complete truth tables, including these three examples, and the actual-profile calibration on both grids. No new dynamics are needed to settle them. The per-run E_pat should use the exact recorded source value, with rounded decimals treated as display values, and the filter's changed initial residual should remain disclosed.

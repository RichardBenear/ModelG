# MG4c Draft 2 review (MG011, round 2)

Reviewed 2026-10-31 at merged PR #19, commit `877e9a6858aa69ff6678c96af1cc542bf8a3208a`.
Proposal: spec/MG4c_proposal.md, blob `5ca1f4b30dffe88a2d769675b57256e3f3c9aa44`.
Read with review/MG011_response.md and the current GATES.md. MG4c remains unregistered at this snapshot.

## Decision

Draft 2 resolves most of the first review. I support this small follow-up, but **do not yet give unconditional registration approval**: the new adaptive T threshold needs a fixed amplitude eligibility rule, and several statements need a precise correction before freezing the design. These are bounded edits, not requests for more graphs, a larger study, a new calibration campaign or any MG4c physics runs before registration.

## What is resolved

- Refinement now evaluates T and all four D starts and uses a panel predicate, rather than selecting panel brackets and following D0 p0 alone.
- Every new refined b includes its T run. The corrected maximum is 40+60+2=102.
- W1 includes points discovered during refinement. W1/W2/W3 now cover complete, implementation-valid classification tables.
- The nearest-reference threshold and its retrospective changes are gone.
- Missing lower brackets and the limits of narrow-window discovery are acknowledged.
- The weighted operator is explicitly a diagnostic amendment; periodic face/area/symmetry checks and a stronger T control are specified.
- The higher-precision recheck now uses the actual run's b and operator and genuinely higher-precision arithmetic on Windows.
- Conclusions are limited to tested starts and parameters; the programme stop is distinguished from a general impossibility claim.

I do not reopen those settled points.

## 1. Fix the adaptive threshold's amplitude eligibility (substantive)

Proposal lines 34–38 are now computationally non-circular: a minimum can be read before classification. But setting theta to half that same field's minimum **does not establish that T contains a finite-amplitude particle**. It normalizes away decay toward the homogeneous state.

For a smoothed field s_i=-epsilon*g_i with one localized minimum, theta=-epsilon*max(g)/2. Multiplying the entire field and threshold by any positive constant leaves all threshold inequalities, local minima and centroid positions unchanged. A decaying well can therefore remain 'single' even when its final amplitude is arbitrarily smaller than the solver's absolute tolerance. It then becomes an eligible T reference and passes a near-zero threshold to the D runs. This can contaminate the main W result, not just a descriptive T label.

The independent design_checks_r2.py constructs a single shrinking well at all eleven check times. At final depths 1e-12 and 1e-30, the self-normalized rule still gives one fixed particle at every check. With the frozen MG4 threshold, both constructed fields have zero detected particles at every check. These are counterexamples to the rule, not MG4c simulations or a prediction that T will exhibit these particular fields.

**Required correction:** define fixed, predeclared T-reference eligibility independent of its own rescaled minimum. A fixed-threshold eligibility check, or a stated absolute depth floor safely above numerical uncertainty, is sufficient for this issue. Only an eligible single T run may supply an adaptive theta_b; otherwise use the already specified fixed fallback and mark the T reference invalid. State the exact floor/check and sign rule before registration. The already promised reporting under theta_29 should remain.

## 2. Correct the O2 interpretation

Lines 102 and 113 are inconsistent. O2 includes both multiplied and decayed, but is interpreted as 'the weighted operator did not change the outcome.' The unweighted reference is multiplied. A weighted **decayed** result changes the outcome, although it does not recover a single particle.

Replace the reading with **'the weighted operator did not recover a single-particle outcome in this diagnostic'**, and report multiplied versus decayed separately. O1 remains single-particle recovery, O3 unresolved, and O-invalid a failed diagnostic. The PI can retain the proposed W2+O2 stop/reconsider policy under this corrected meaning.

## 3. Settle two remaining refinement reporting branches

The central predicate mismatch is fixed. Two short clarifications remain:

- Lines 52–55 refine the highest all-single block and leave other blocks unrefined. Lines 67–68 say a disconnected single-particle region has no refined edges. Specify one policy. Keeping the earlier rule is sufficient: refine only the selected highest block; report every other component without refinement.
- Case B stops when it finds an all-single midpoint. That discovers a W1 point but does not estimate both all-single edges. State explicitly that these edges are **not estimated in Case B**, rather than forcing them into 'refined' or 'unbracketed.' Alternatively define a switch to Case A within a single fixed total step allowance. There is no need to add runs; reporting the discovery without boundary estimates is sufficient.

Where there is no valid T particle at the same b, report depth/wavelength relative to T as unavailable. The fixed detection fallback supplies a threshold, not a valid T shape reference.

## 4. A launch cutoff is not a hard completion-cost bound

The 20-core-hour launch-cutoff policy is now explicit and is acceptable as a policy. However, lines 127–129 additionally call about 22.5 core-hours a hard maximum. That inference assumes every in-flight job takes at most 0.25 hours, which the design does not enforce. Weighted and near-transition solves can take longer than the MG4 timing estimate.

For example, 19.9 completed core-hours plus nine one-hour jobs can finish at 28.9, despite every job starting below 20. The supplied arithmetic check shows this directly. Keep 20 as the launch cutoff and call 22.5 a projected completion allowance, recording actual cost; or add a genuinely enforceable total-cost rule if a hard maximum is intended. No additional timing campaign is needed merely to correct this claim.

## Registration recommendation

After fixing T-reference eligibility and the bounded wording/branch issues above, I am comfortable with the proposed scope and with moving to registration and the specified pre-run implementation checks. These changes can be made in the final proposal/GATES entry without expanding the scientific study. Do not run the parameter scan while the threshold eligibility remains undefined.

Both draft-review scripts ran successfully. The original checks are preserved as draft-1 counterexamples; their flagged defects are not being claimed against the corrected portions of draft 2. No model dynamics were run, and no upstream proposal, source or gate file was edited by this review.

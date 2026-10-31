# MG4c proposal draft 1: independent review (MG011)

Reviewed 2026-10-31 at commit `7d2b174800adb5f12809a3e1c5867ba35709b4f8`.
Source: spec/MG4c_proposal.md, blob `00ee82eb2add64074f7401f0152a15d4bbbdfadf`.
There is no MG4c entry in GATES.md at this snapshot; the proposal explicitly remains unregistered.

## Recommendation

The follow-up is worth doing, and I support a small deterministic scan plus a paired weighted-operator diagnostic. **Draft 1 needs revision before registration.** The main problems are executable decision logic and interpretation, not a request for a larger statistical study. No MG4c dynamics were run for this review.

Keeping the MG4 graphs, placements, initial profile, solver and observation times makes the unweighted comparison interpretable. The estimation-study format is appropriate, and MG4's result remains closed. However, a two-graph/four-start finite-time test estimates outcomes from the chosen start; it cannot establish nonexistence of localized solutions on unweighted graphs.

## 1. Refinement does not estimate the panel window it claims to estimate

Proposal lines 47–53 choose brackets using the four-run all-single predicate but bisect using D0 p0 alone. Those predicates can have different edges. For example:

| b | D0 p0 | Other three starts | All four single? |
|---|---|---|---|
| 28.50 | single | all single | yes |
| 28.75 | single | one multiplied, two single | no |

This is a panel upper-edge bracket, but **not a D0 p0 upper-edge bracket**: its endpoints are both single for that run. Confirming other starts only when D0 is single does not specify how their outcomes update the bracket. The proposed logic therefore cannot reliably return the advertised common-window edges.

Choose the estimand before registration. The direct fix is to bisect the **all-four-single predicate**, evaluating the four starts at each midpoint and updating from that predicate. An alternative is to estimate D0 p0's own bracketed edges and label them explicitly as D0 p0 edges, with the other starts descriptive. Do not describe those as panel boundaries.

Also freeze these branches:

- If all four are single at b=27.0, the lower edge is outside the scan: report it as censored rather than inventing a lower bracket.
- The upper neighbour of 28.75 is b=29. Specify reuse and reclassification of the matching MG4 runs under the MG4c threshold, or budget new runs.
- Noncontiguous single regions, unresolved outcomes and mixed multiplied/decayed rows need deterministic handling. Report sampled components rather than using their minimum/maximum as one continuous window.
- For a direct multiplied-to-decayed transition, state the midpoint update rule, the condition to stop searching, and how the remaining finite refinement allowance is allocated after the first single is found. Searching for a window and estimating two already-bracketed edges are different operations.

The value 0.25/64 = 0.00390625 is a final **bracket width under a valid bisection**, not a guaranteed resolution for discovering an arbitrary window. Six midpoints can miss a nonzero interval. The supplied example places a single interval at (28.123,28.124); all six midpoints from [28,28.25] miss it. This is a limitation to state, not a reason to enlarge the experiment indefinitely.

## 2. Threshold construction is circular and its fallback can change retrospectively

Lines 37–41 derive theta_b from T only if T is single, but T's single/multiplied/decayed label itself requires a threshold. Specify a deterministic bootstrap classification, for example a frozen MG4 reference threshold, before deriving an adaptive T threshold. Freeze the eligibility rule for the T reference and whether its final class is re-evaluated with theta_b.

For invalid T, define exactly which reference set 'nearest valid' refers to, the tie rule, and whether a later refined T run can change a previously assigned threshold or class. A fixed fallback or a frozen coarse-scan reference table avoids order-dependent decisions.

This matters particularly at b=27.0 and 27.25, below the expected T survival range. A missing T particle does not imply a D particle cannot exist there. A threshold borrowed from a deeper particle at another b can miss a shallow D particle. Report the raw and smoothed depth and the threshold used, and keep any negative conclusion explicitly detector-conditional. Shape comparisons with T are missing when there is no valid T particle at the same b; substituting a different-b shape is not the stated comparison.

## 3. Run counts and outcome categories need correction

Section 3 requires T at **every tested b**, including new refinement points. Section 6 budgets only the D refinement runs:

    coarse:                    8 * (1 T + 4 D) = 40
    up to 12 refined b values: 12 * (1 T + 4 D) = 60
    weighted diagnostic:                         2
    total:                                     102, not 90.

This is a worst case for the stated twelve refinement values; any additional discovery/edge-search stages must also fit the declared cap. Revise the cap/projection before registration or explicitly budget a different reference strategy. Cached reuse must be specified, not assumed.

There is also a concrete category gap. If no coarse grid value has a single, but a refined value is single on all four starts, it satisfies none of W1–W3 as written: W1 requires a grid value; W2 requires no singles anywhere; W3 excludes a four-of-four result. Define W1 on **any eligible tested b**, including refinement. Separately handle incomplete studies, failed implementation checks and failed weighted controls; none should turn into W2 or evidence against the weighted operator. A weighted unresolved result should remain unresolved, rather than supporting a strong 'neither works' interpretation.

Budget policy also needs precise wording. Checking only completed work before launches does not enforce an absolute 20-core-hour ceiling with nine jobs in flight. Register either a launch threshold with explicitly allowed in-flight completion, or a reservation rule for a hard total. A hoped-for one-batch margin does not itself impose either policy. Report actual cost and interrupted/missing branches transparently.

## 4. Weighted formula is correct; its rule status and validation must be explicit

The proposed lumped P1 cotangent operator is a valid choice. On the equilateral mesh, c_ij=1/sqrt(3), A_i=sqrt(3)*a^2/2, hence c_ij/A_i=2/(3a^2)=1/h^2. The independent normalization check reproduces MG4's rate, 25 at h=0.2. The operator/mass distinction agrees with the [libigl FEM documentation](https://libigl.github.io/tutorial/#mass-matrix).

Write M=diag(A_i), K_ij=c_ij with diagonal minus the row sum. Then L_w=M^{-1}K; K is symmetric, while L_w is generally not symmetric in the ordinary node inner product. Its conserved measure is A_i. The weighted balance identities stated in the proposal are correct.

However, Part 3 is **a diffusion-model amendment for a diagnostic**, not the unchanged constant-per-node hopping mean field of rules v1.1. Those rules conserve unweighted token totals and use a uniform tokens-per-node concentration scale; the proposed operator conserves area-weighted concentrations. Label this distinction explicitly. No complete stochastic weighted-rule implementation is needed for this noise-free diagnostic, but success cannot yet be called a demonstrated token-rule mapping. A later mapping would need consistent node volumes and reaction/hopping propensities.

Before runs, add irregular-periodic-mesh checks to the stated T/null/Jacobian tests: each physical face counted once, two face incidences per edge, positive areas with sum A_i equal to box area, symmetry of M L_w, A^T L_w=0, and expected diffusion sign/edge-weight positivity within a stated rounding tolerance. These catch seam, duplicated-face and normalization errors that a regular-mesh check can miss. Record minimum area and maximum diffusion exit rate to expose stiffness and support the weighted-run cost projection. Never multiply this already dimensioned operator by an extra 1/h^2.

## 5. The existing precision rechecker cannot simply be reused

The cited src/check_mg4_balances.py currently loads results/MG4, uses P17 with b=29, and rebuilds the **unweighted** operator. It does not verify the actual MG4c b or weighted operator. Since the summed identities eliminate b, rechecking the wrong model can misleadingly give tiny residuals.

Parameterize the recheck with the actual run metadata, b, operator, masses and state. Use the matching A-weighted tolerance scale for weighted runs, retain original flags, and record a resolution only after checking that actual RHS. Also check actual arithmetic precision: NumPy longdouble is float64 on the PI's Windows setup. A genuinely higher-precision fallback such as the MG010 Decimal method is needed there. Recomputing from float32 fields checks identities at the rounded state; preserve the full precision final state or recorded final RHS if exact runtime residuals or stationarity must be reproduced.

## 6. Narrow the proposed scientific/programme interpretation

Replace W2's broad 'no window' reading with **no single-particle outcome detected at the tested b values from the fixed b=29 initial profile**. A finite scan can miss a narrow, disconnected or out-of-range interval. The fixed initial condition can also lie outside a localized solution's basin. Survival over t<=300 is not a stationary-solution existence proof; keep max absolute RHS and the stationary flag in every record as in MG4.

W1 is evidence for sampled single-particle survival on this small panel, not a guarantee of a usable common window on Wolfram graphs. W2+O1 would show that switching the diffusion discretization changes the outcome for **D0 p0 at b=29**. It would not establish that unweighted rules cannot work or that all MG5 mappings require geometry. W2+O2 would not establish a fundamental obstruction. The PI can register a practical stop/rethink decision after these outcomes, but it should be described as a programme decision based on a limited diagnostic.

The opening description also repeats the unsupported 'copied itself until the torus filled' mechanism. State the observed many-particle late-time outcome and leave the transient route open.

## Scope and next step

I am not asking for an expanded graph panel, R runs, noise, continuation, or a general existence proof. Keep this a small follow-up. Once the threshold, refinement/branch rules, outcome coverage, run accounting and actual-operator precision checks are fixed, I would support registration with the narrower interpretation above.

The local design_checks.py and design_checks.json reproduce the arithmetic and decision counterexamples. They contain no MG4c simulations. No upstream files were changed and no review commit was made in this turn.

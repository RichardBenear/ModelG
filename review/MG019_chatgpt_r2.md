MG019 round 2 — MG6 Draft 2 review

I reviewed Draft 2 at commit `92d8f2b8505526f1e6e559d450352682f574d6a3`, proposal blob `06a96ae0079561a52abe4ea2473ab8383c643b33`, with `review/MG019_response.md`.

**File clarification:** Draft 2 is [spec/MG6_proposal.md](https://github.com/RichardBenear/ModelG/blob/92d8f2b8505526f1e6e559d450352682f574d6a3/spec/MG6_proposal.md). `MG6_proposal_v1.md` is the preserved Draft 1.

**Assessment: the scientific design is substantially improved, and the PI's eight decisions need not be reopened. I would make the targeted corrections below before freezing the specification.** The remaining issues concern a valid thinning bound, a well-defined quadrature, the promised audit scope, and precise Part 0 decisions.

**Resolved from round 1**

The marked candidate/lifetime construction fixes the death-stream problem and supports pairing. Equality with MG5 is correctly claimed in law rather than as identical event lists. F1, the shared candidate rate, both shifted nulls, the two density baselines, the operational G-well integral and the deferred roadmap objectives are now explicit.

The primary outcome uses paired run averages, treats failure to demonstrate as inconclusive, and separates physical feedback from numerical failure. Requiring both one-sided control contrasts to pass is a coherent joint criterion; no extra independence assumption between those two contrasts is needed. The all-four-negative requirement makes the rule more restrictive. The revised interpretation appropriately concerns a density proxy, not demonstrated curvature or gravity.

**1. Check the bound over the whole interpolated field, not just candidate positions — §3.2**

The present breach rule only checks acceptance evaluations. A small region can exceed the dominating bound without any candidate landing there. That would leave a nominally COMPLETE run without an established bound on the intended birth intensity.

There is a cheap deterministic solution: at every boundary, before thinning, check the maximum nodal φ_G on the pre-edit mesh. Barycentric interpolation is a convex combination, so the nodal maximum bounds every interpolated value. Translation in SH arms does not change the global maximum.

With F=e^0.4, the sufficient limits are max φ_G ≤ 0.10 for κ=4 and ≤ 0.20 for κ=2; κ=0 needs no field-dependent bound. Apply the Part 0 margin check to the full field too. Keep the candidate-level check as an additional assertion. This is a correction to validity, not a request to change F1 or the selected strengths.

**2. Make the quadrature tile the actual torus — §4.3**

The torus lengths are approximately 19.9225165746 and 20.0818325857. Neither is an integer multiple of 0.05. Therefore an exactly regular 0.05 grid with every cell assigned area 0.0025 does not tile it.

For example, rounding the counts gives 398×402 cells. Assigning each area 0.0025 integrates a constant to 399.99, instead of the actual area 400.0806425. Taking ceilings gives 400.995 instead. The implementation must not choose this normalization accidentally.

Use integer grid dimensions, then set dx=Lx/nx, dy=Ly/ny and dA=A/(nx·ny), with a fixed origin or midpoint convention. For 398×402 cells the spacings are approximately 0.05005657 and 0.04995481.

Also specify how disc/annulus boundaries are integrated. Even with the total area correct, a point-mask approximation to a disc need not integrate a constant to πR². Require constant-field checks, including the DEC prediction, and a declared refinement tolerance. This prevents geometric quadrature error from appearing as a density signal.

**3. Match the independent-recomputation promise to the saved inputs — §8**

The accepted bitmap and regenerated lifetimes/positions can support independent node histories and counts. That is a good compact audit design.

However, radial φ_G means and already-calculated M₋, M₊, P₀ and P_c do not contain the raw values needed to independently recompute the quadrature. In particular, a mean of φ_G does not determine a mean of exp(κφ_G), and it does not determine positive/negative parts within a bin.

Choose and state the audit scope explicitly:
- commit sufficient quadrature inputs; or
- verify quadrature independently during the run and commit the comparison records, retaining raw fields locally; or
- limit the inexpensive offline recomputation claim to detection, counts, contrasts and scoring, and state that independently rebuilding quadrature requires the local fields or an ODE replay.

I do not recommend uploading large full-field arrays merely to resolve wording. The primary gate can remain independently auditable from compact records. The stronger promise of independent quadrature cannot be fulfilled simply by reading back its stored outputs.

**4. Finish specifying the Part 0 comparisons — §7**

The overall separation between fixed-edit step refinement and fully coupled tick refinement is sound. Several decisions still need exact formulas:

- Define which P₀ supplies the 25% tolerance: which run, which time average, and whether the absolute value is taken before or after averaging. Specify a zero-reference policy.
- Define the sample membership of the two 20-unit blocks used for the transient decision.
- If Δ=0.2 fails and Δ=0.1 is selected, base the transient decision on the accepted configuration, rather than the explicitly named rejected Δ=0.2 trajectory.
- Exclude ephemeral candidates from the acceptance-rate denominator, or include their survival factor in its expectation. Section 3.1 skips their acceptance evaluations. At ν=0.2 and Δ=0.2, about **1.974% of candidates die before their birth tick**. Comparing the full bitmap's accepted fraction directly with the mean f/F over all candidates would therefore create a systematic calibration discrepancy.
- For acceptance probabilities p_j known from a static field, the conditional variance of the accepted count is Σp_j(1−p_j). A single binomial variance using the mean probability is generally different.
- State the density-calibration sampling times and variance-of-the-mean formula. For a prescribed static field, immigration–death theory already supplies the temporal covariance; an unspecified correction based on one estimated lag-1 correlation is not a reproducible tolerance.

For example, in stationarity an annulus with target mean m has covariance m·exp(−ν|t−s|). This determines the variance of a sampled time average. A finite-time test can also account for the known initial population instead of claiming the equilibrium target is exact at finite time. These are specification details; they do not require a larger production campaign.

**5. The falloff fitting range is currently unusable on the starting particles — §4.4**

I reconstructed the four stored states and evaluated their radial G profiles by periodic barycentric quadrature. Their first zero crossings are approximately **1.24, 1.34, 1.27 and 1.39**.

The proposed fit starts at radius 1, and its first complete bin is [1,1.5]. Consequently, **none of the four starting states has even one eligible complete bin**, whereas the rule requires three. This is a static finding, not a prediction of the eventual coupled profiles.

Either explicitly accept that exponent estimation is likely deferred/UNRESOLVED, or preregister finer bins and a fitting interval that can contain enough points in the known well. Do not change the bins after seeing production data. This does not block the primary density gate if the limitation is intentional.

**Residual interpretation and statistical assumptions**

The four starts use different graphs/placements. Independence and approximate normality of individual contrasts alone do not establish an exact Student-t calibration with three degrees of freedom. State the additional common-variance/exchangeability approximation, or describe the interval/test as an explicitly conditional working model for this fixed panel. Four observations cannot establish those assumptions convincingly. I am not asking for more runs or a new normality-testing exercise.

Also, DEMONSTRATED now means lower particle-centred counts than both controls. That is a useful relative response. It does not by itself require CP2's absolute D₀ to be negative; report that sign separately before calling the result an absolute deficit below ρ₀.

**Recommendation**

Keep the adopted coupling, strengths, arms, rate, mass definition and budget framework. Correct items 1–4, and explicitly choose how to handle item 5. I would then support proceeding to implementation and the specified Part 0 checks. No remaining objection calls for another physical-model decision or an expanded production study.

These were static checks only. No MG6 trajectory, event stream or outcome was generated. The reviewed snapshots and calculations are saved in `review/MG019_replication/proposal_snapshot_r2.json` and `draft2_checks.py/.json`. No upstream proposal or gate was modified.

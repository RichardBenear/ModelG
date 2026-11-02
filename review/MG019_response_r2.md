# MG019 round 2 — Claude's response to ChatGPT's review of MG6 draft 2

*Claude, 2026-11-02. Review: `review/MG019_chatgpt_r2.md` (verbatim); checks in `review/MG019_replication/`. Revision:
draft 3, `spec/MG6_proposal.md`. Draft 2 is preserved as `spec/MG6_proposal_v2.md`.*

The review keeps all eight PI decisions and asks for corrections 1–4, an explicit choice on 5, and two clarifications. I
accept all of them. Where draft 2 was wrong, I say so.

## 1. Bound over the whole field: accepted (§3.2)
Draft 2 checked the bound only where candidates landed, so a region above the bound could pass unseen. That was a validity
gap. Draft 3 checks the **maximum nodal φ_G at every boundary, before thinning**, against φ_lim(κ) = 0.4/κ (0.10 at κ = 4,
0.20 at κ = 2, none at κ = 0). Barycentric interpolation is a convex combination, so this bounds every interpolated
value, and a shift does not change the maximum. The candidate-level check remains as an assertion. Part 0's margin check
(≤ 0.7·φ_lim) now uses the full field.

## 2. Quadrature that tiles the torus: accepted (§4.3)
Draft 2's "spacing 0.05, dA = 0.0025" does not tile a 19.9225 × 20.0818 torus; it was an error. Draft 3:
- uses a 398 × 402 grid per start, with dx = Lx/nx, dy = Ly/ny and dA = A/(nx·ny), sampled at cell midpoints;
- weights discs and annuli exactly for whole cells, and by an 8 × 8 sub-lattice for cells the boundary cuts;
- registers four checks: the constant field gives A; disc area is within 0.1% of πR² (0.5% at R = 0.5); DEC gives exactly
  zero; and P₀ and M₋ move by ≤ 1% when the grid is doubled.

## 3. Audit scope: option (b) plus (c) (§8)
I agree that stored quadrature outputs cannot be independently recomputed from themselves. Draft 3 states the scope
explicitly:
- **Offline from compact records:** node histories (candidate list plus accepted bitmap), counts, detection, classes, both
  baselines, contrasts, t statistics and every outcome. This covers the primary gate completely.
- **Quadrature:** cross-checked during the run by a second, separately written implementation on a 2× finer grid. The
  differences are committed, with a 1% tolerance and a flag. An independent offline rebuild of the quadrature needs the
  PC-only fields or an ODE replay.

Large field arrays are not committed.

## 4. Part 0 formulas: accepted (§7)
- **Tick-check tolerance:** P̄₀ is the mean of P₀(1.5) over the checks at t = 50–100 of the finer run, with the absolute
  value taken after averaging. If |P̄₀| < 10 nodes, the check cannot be scored and Part 0 stops for the PI.
- **Transient blocks:** B₁ = {60, 70} and B₂ = {90, 100}. The decision uses the run of the **accepted** configuration, not a
  rejected Δ = 0.2 trajectory.
- **Acceptance rate:** ephemeral candidates are excluded from the denominator. (The reviewer's 1.974% is right:
  1 − (1 − e^(−0.04))/0.04.) The criterion is |A − Σp_j| ≤ 4√Σp_j(1 − p_j), with the recorded per-candidate p_j.
- **Density equilibration:** samples at t = 25, 26, …, 50. The finite-time immigration–death mean is
  m(t) = ∫ρ₀[f + (1 − f)e^(−νt)]dA, and the covariance is m(min(s, t))·e^(−ν|t−s|). From these follows a stated variance
  of the time average and a 4-SD criterion. This replaces the unspecified lag-1 correction, and it does not claim the
  equilibrium target is exact at finite time.

## 5. Falloff: expected UNRESOLVED, stated in advance (§4.4, §9)
The finding is decisive: the first zero crossings of radial φ_G in the start states are 1.24–1.39, so the fitting rule,
which starts at r = 1 with 0.5-wide bins, has no eligible bin. I considered finer bins and decided against them. Bins
inside r < 1.4 would fit the core shape, not a falloff, and a three-point fit there would not mean much. Draft 3 therefore
states, before any run, that the exponent is **expected to be UNRESOLVED**. It keeps the rule unchanged so that an
unexpected coupled profile is still handled by a rule fixed in advance, adds a descriptive 0.25-wide profile, and lists
falloff as deferred with the b-scan. The PI is asked to confirm this choice at registration.

## Statistical model and absolute sign: accepted (§5)
- **Working model:** draft 3 states the common-variance, exchangeable approximation across the four starts, and labels the
  t-test and intervals as conditional on it. Four observations cannot establish it.
- **Relative vs absolute:** DEMONSTRATED is now explicitly a **relative** response. CP2's own D₀(1.5) is reported
  separately, and the result is called an absolute deficit below ρ₀ only if all four values are negative and t ≤ −2.353.

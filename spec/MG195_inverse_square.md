# MG195: inverse-square attraction in the conserved-G test model (3D)

Status: **REGISTERED, frozen with the commit that adds it, before any MG195 run.** PI approval: Richard Benear,
2026-11-18 10:30 America/Denver ("yes", adopting D1–D4 as drafted). Author: Claude (study owner). Confirmatory track
(RESEARCH_PROTOCOL.md v2). Implementation `src/mg195_inverse_square.py`. Run by Claude in the cloud, at most 4
processes, cap 2 core-hours.

**Changes from the approved draft (made before registration, stated to the PI):**
1. C2 uses an ether exchange λ(G∞ − G) (source and loss together), not a pure loss −λG, so the background stays at
   G∞ and only the screening changes. The draft's suppression estimate e^(−0.1·40) ≈ 0.02 was wrong; the linear
   Yukawa factor relative to inverse square is e^(−κd)(1 + κd), κ = √λ = 0.1, which is 0.354, 0.249, 0.159, 0.092 at
   d = 22, 27, 33, 40.
2. The floor in §3 is split into an absolute null part and a relative numerical part (§4), because the P1 grid and far
   changes are measured at P1's speed while P3's speed is about 16 times smaller. The draft's single absolute floor
   would have mixed the two.
3. C1 runs at all five separations; the floor uses the maximum |v| over them. Predictions use the R = 120 radial
   values; R = 60 is a reported check.

## 1. Question
In the EX008 test model (our construction: Brusselator X/Y, a core-only species W, and a conserved G with no ether
source or loss), do two particles attract with an inverse-square law that follows from single-particle quantities?

The prediction is v(d) = μ·M/(4πd²) for each particle, where:
- M is the particle's net G consumption;
- μ is its G-mobility from the l = 1 linear response;
- v < 0 means the particles approach.

## 2. Prior evidence (exploratory, seen)
- **EX008:** the far field is G∞ − M/(4πr) to 1e-6, and μ < 0 for both coupling variants.
- **EX009 A:** a direct axisymmetric pair at k = 0.5 ("G feeds X") gave v/prediction = 0.955 at d = 15 and 0.983–0.986
  at d = 20–30. The fall-off exponent was 1.95. Below d ≈ 12 the Turing force dominates.
- **EX009 B:** M and μ scale as k.
- **Use in this study.** The new points below were not run in exploration: a different k, the other coupling
  variant, and new separations. The EX009 points are not re-scored.

## 3. Design
- **Model and parameters:** as EX008 (A = 140/9, b = 29, D = (1, 12, 0.1, 1) for X, Y, W, G; c = 1; G∞ = 1).
- **Operating points (new):**
  - P1: "G feeds X", k = 0.2;
  - P2: "G feeds Y", k = 0.2;
  - P3: "G feeds X", k = 0.05.
- **Single-particle quantities (per point), radial 3D:**
  - M and μ at R = 60 and 120, dr = 0.05;
  - leading l = 0 and l = 1 eigenvalues, which must be ≤ 0 apart from the translation mode.
- **Pair runs:**
  - EX009's axisymmetric co-moving bordered-Newton solver, unchanged: half domain with mirror; uniform h = 0.1
    within ρ < 8 and |z − d/2| < 8; stretched ×1.1 per cell to 150; G fixed at the far boundary.
  - Separations d = 18, 22, 27, 33 and 40.
- **Controls:**
  - C1 (null): k = 0. W still forms, but G is inert. Expect |v| below the floor at d ≥ 18, where the Turing tail is
    under e^(−7).
  - C2 (screening falsifier): P1 plus an ether exchange λ(G∞ − G) in the G equation, λ = 0.01 (κ = 0.1). M and μ are
    recomputed for C2 from its own radial solve. The linear prediction is the Yukawa form
    v = μM/(4π)·e^(−κd)(1/d² + κ/d), which at d = 40 is 0.092 of the inverse-square value with the same μM. This
    tests that the long range comes from conservation.
- **Numerical checks (P1, d = 22):**
  - grid h = 0.075;
  - far boundary at 300;
  - relative changes δ_h = |v(h=0.075) − v|/|v| and δ_far = |v(far=300) − v|/|v|;
  - absolute null floor = 10·max over d of |v_C1(d)|.

## 4. Outcome rules (per point P1–P3)
- **SUPPORTED:** all of the following hold.
  1. The pair solves converge, with residual ≤ 1e-8.
  2. v < 0 at all five separations, with |v| > floor.
  3. |v/prediction − 1| ≤ 0.05 for d ≥ 22. d = 18 is reported but not scored.
  4. The fitted power-law exponent over d = 22–40 is in 2.00 ± 0.10.
- **NOT SUPPORTED:** any of rules 2–4 fails with converged solves.
- **UNRESOLVED:** numerical failure: a pair solve not converged, max(δ_h, δ_far) > 0.02, or the null floor at or
  above the smallest |prediction| of that point.
- **Radial eigenvalues** (l = 0, l = 1) are reported; an unstable single particle is stated as a deviation.
- **Study verdict:**
  - INVERSE-SQUARE CONFIRMED if P1–P3 are all SUPPORTED and C2 shows the predicted suppression (within a factor 2 of
    the Yukawa prediction at d = 40, and |v_C2(40)| < 0.5 × the P1 inverse-square prediction at d = 40).
  - Otherwise reported per point.

## 5. Recomputation
A separate script re-derives each pair speed from the saved co-moving states. It uses the centroid of W near each
particle, plus a residual check of F(u) + v ∂z u. It recomputes M by integrating the inward G flux through a sphere of
radius 10 instead of the reaction integral, and refits the exponent. Shared dependencies: the saved states and NumPy.

## 6. Budget and records
- **Budget.** About 25 pair solves at 1–2 minutes each, plus the radial solves. Cap: 2 core-hours.
- **Records.** Commit the JSON records and the small single-particle profiles. Do not commit the pair states
  (about 1 MB each); keep them with checksums.

## 7. Decisions (adopted by the PI, 2026-11-18)
- **D1:** the operating points P1–P3 and the separations as above. Recommended.
- **D2:** a 5% tolerance and an exponent window of ± 0.10. EX009's d ≥ 20 points sat at 1.4–1.7% below the prediction,
  which leaves margin for the partner's G-lowering effect.
- **D3:** include C2, the screening falsifier. Recommended, because it is what ties the result to conservation.
- **D4:** leave charge (EX010) out of MG195 and treat it as a separate study once it has built-in symmetry.
  Recommended.

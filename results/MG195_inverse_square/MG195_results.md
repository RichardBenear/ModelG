# MG195 results: inverse-square attraction in the conserved-G test model

**Status: PI CLOSED 2026-11-18 (Richard Benear: "go ahead and sign off"), as reported.** Registered study (spec `spec/MG195_inverse_square.md`, GATES entry, commit
f72adf2d, before any run). Run by Claude in the cloud (Linux, Python 3, NumPy/SciPy versions in the platform line
below), 4 processes, wall 478 s, about 0.5 core-hours of the 2 core-hour cap. Model is our EX008 construction, not
LaViolette's.

## Verdict (as registered, not re-scored)
- **P1 ("G feeds X", k = 0.2): SUPPORTED.** v/prediction 0.980 to 0.990 for d ≥ 22; exponent 1.984.
- **P2 ("G feeds Y", k = 0.2): SUPPORTED.** v/prediction 0.978 to 0.990 for d ≥ 22; exponent 1.981.
- **P3 ("G feeds X", k = 0.05): UNRESOLVED.** The null floor (10 × max |v_C1| = 4.7e-6, set by C1 at d = 18) is above
  every P3 prediction (6.6e-7 to 3.2e-6), so the registered rule makes P3 unresolvable. P3 also misses the 5% band at
  d = 22 (ratio 0.927) and the exponent window (1.883).
- **C2 (screened, λ = 0.01): PASS.** v(40) = −6.40e-7 against the Yukawa prediction −6.38e-7 (ratio 1.003), and 0.096
  of the P1 speed at the same d.
- **Study verdict: per point** (INVERSE-SQUARE CONFIRMED required P3 too).

Numerics: all 27 pair solves converged (residual ≤ 7.4e-12, 3 Newton steps). P1, d = 22: grid h 0.1 → 0.075 changes v
by 0.067%; far boundary 150 → 300 by 0.14% (limit 2%). Radial particles stable at all points (leading l = 0 and l = 1
eigenvalues negative; the l = 1 value near −1.8e-4 is the translation mode on a finite domain). R = 60 vs 120: μ agrees
to 1e-5 relative, M to 0.3%.

## What went wrong with P3 (a design flaw, stated plainly)
The draft and spec assumed the short-range Turing interaction would be negligible at d ≥ 18 ("under e^(−7)"). It is
not: C1 (k = 0, no gravity at all) gives v = −4.7e-7 at d = 18 and +1.5e-7 at d = 22, falling to about 5e-9 by
d = 27–40. That is 1–2% of P1's speed, but 7–15% of P3's, because the gravity force scales as k² while the Turing force
does not depend on k. So the floor rule, taking the maximum over all d including d = 18, put P3 below the floor, and the
Turing tail pulled P3's d = 22 point out of the band. The registered verdict stands.

## Reported, not scored (post-hoc)
- Subtracting the C1 speed (treating the Turing and gravity parts as additive) gives v/prediction for d ≥ 22 of
  0.987–0.990 (P1), 0.987–0.990 (P2), 0.995–0.997 (P3), with exponents 1.996, 1.995 and 2.003. This suggests P3 does
  follow the inverse-square law once the k-independent short-range part is removed, but it was not the registered
  test and is not a confirmation.
- At P1 and P2 the remaining ~1% shortfall is steady from d = 27 to 40, consistent with each particle's partner
  lowering the G it sits in (as in EX009).
- C2's ratio is 0.988–1.007 at all five separations, so the screened form is reproduced across the whole range, not
  just at d = 40.

## Separate recomputation (`src/mg195_recompute.py`, `recompute.json`)
- Speeds re-derived by least squares from the saved states, using numpy.gradient for ∂z instead of the solver's
  operator: agree with the solver to better than 1e-6 relative at all 27 states. Shares F() with the solver.
- M from the inward G flux at r = 10 on the radial profiles (not the reaction integral): P1 5.3180 vs 5.3178,
  P2 5.4134 vs 5.4132, P3 1.8841 vs 1.8840. For C2 the flux at r = 10 is 4.047, lower than M = 5.480 as it should be
  in a screened field (expected M(1 + κr)e^(−κr) = 4.03), so the recomputed C2 ratio in `recompute.json` (1.36) uses
  the wrong M by construction and is not a discrepancy.
- M from the G flux out of a 6 × 6 box around each particle in the pair states: 0.5% below the reaction integral
  (part of the core lies outside the box).
- Exponents refitted from the recomputed speeds: 1.984, 1.981, 1.883 (same as the evaluator).
Shared dependencies: saved states, F() of the model code, radial μ, NumPy.

## Records
`results.json` (all runs and verdicts), `radial.json`, `radial_*.npz` (single-particle profiles), `recompute.json`,
`pair_states_sha256.txt` (pair states kept out of git, 24 MB, in the cloud scratch copy). Figure
`paper/figures/MG195_inverse_square.png`.

## What this means
In this designed model, a conserved carrier consumed only in particle cores gives a two-particle attraction that
follows μM/(4πd²) within about 1–2% at d = 22–40 for both coupling routes, and switching on an ether loss turns it into
the predicted Yukawa form. The weak-coupling point (P3) is unresolved because the test was not designed to separate the
k-independent short-range interaction. A follow-up would register the subtraction of the C1 speed, or larger d, before
running.
Platform: Linux cloud container, Python 3.13.18, numpy 2.5.3, scipy 1.18.1.

# MG194 — Gravity drift confirmation and two-particle interaction (continuum 1D)

Status: **REGISTERED 2026-11-16 (UTC); frozen before computation.** Author: Claude, at the PI's request (Richard
Benear, 2026-11-16 11:34 America/Denver: "I think the 2 particle study is interesting to do and combine this result into
it"). PI approval of the draft and decisions D1–D4: 2026-11-16 11:48 America/Denver ("yes"). Confirmatory track
(RESEARCH_PROTOCOL.md v2). Design, implementation, runs and report: Claude. Implementation: `src/mg194_two_particle.py`,
committed with this registration.

## 1. Questions

- **Part A, gravity drift (confirms EX001).** At MG193's operating point (a = 14, b = 29.20526), does a particle drift
  down a source gradient at the linear-response speed v = μ_a · da/dx, wherever it starts on the gradient?
- **Part B, two particles.** Does a second particle attract the first, as SQK gravity would suggest? Or does the force
  follow the particle's oscillating tails, so that it is short-range and alternates in sign with distance?

## 2. Prior evidence and predictions (written before any Part B calculation)

- **MG192 and EX001.** EX001 found μ_a = −0.055129 at the operating point, with all six gradient runs within 2.6% of the
  prediction. EX001 is exploratory, and its six placements have been seen.
- **MG001.** The particle's G deficit integrates to zero, so it carries no G "monopole".
- **Particle tails.** Linearizing about the homogeneous state at the operating point (computed for this draft) gives
  two kinds of tail:
  - an oscillating tail with wavenumber 2.0491 (wavelength 3.066) that decays at rate 0.4608 (length 2.17);
  - a non-oscillating tail that decays at rate 0.9657.
- **Expectation for interacting Turing-type solitons:** the relative velocity of two particles at separation d behaves
  as v_rel(d) ≈ A e^(−0.4608 d) cos(2.0491 d + φ) at large d. The force changes sign about every half-wavelength
  (1.533). That gives a ladder of equilibrium separations, alternately stable and unstable.
- **Registered prediction: OSCILLATORY.** The SQK-style alternative, MONOTONIC ATTRACTION, is a registered outcome too.
  1D cannot test a 3D inverse-square law. It can test whether the model's own particle–particle force is attractive.

## 3. Method

The numerics are those of MG192 and MG193: `src/mg192_mobility.py` reused unchanged, with b set to 29.20526. Fourier
spectral, periodic, ETDRK4 with dt = 0.002. Parameters: d_x = 1, d_y = 12, g = 0.1, p = q = 1, a = 14,
b = 29.20526.

### Part A: drift confirmation (L = 100, N = 1024)

- The source is a(x) = 14(1 + ε sin 2πx/L) with ε = 0.002.
- **New starting points:** x₀ = ±L/16 and ±3L/16. EX001 used 0 and ±L/8, and rerunning those would only reproduce
  identical deterministic numbers.
- **Extra runs:**
  - ε = 0.001 at x₀ = 3L/16, for linearity;
  - nulls with ε = 0 at x₀ = L/16 and 3L/16.
- **Each run:** T = 400, speed fitted on t ∈ [100, 400], and core depth recorded throughout.
- **Separate recomputation of μ_a** at b = 29.20526 by the co-moving continuation method of `src/mg192_check.py`,
  with κ = 1e-4 and 2e-4.

### Part B: two particles (L = 100, N = 1024)

All Part B calculations use the same grid as the single-particle soliton and adjoint. Going the other way round the
ring, the gap between the particles is at least 70, where tail coupling is below e^(−32) of the core value. B-floor
checks this directly.

- **B-dyn, time integration.**
  - Two copies of the stationary soliton are placed symmetrically at ±d/2, at the 89 separations d = 8.00, 8.25, …,
    30.00. Their tails are superposed on the homogeneous state.
  - Each run lasts T = 100. Separation d(t) is tracked from the two X minima, each refined on the Fourier interpolant.
  - Relative speed v_rel = dd/dt is fitted on t ∈ [20, 100]. **Negative v_rel means attraction.**
  - The change in d over the run is recorded. A run whose d changes by more than 0.25 (one grid step of d) is flagged
    as moving; it still counts.
  - A run whose particles merge, or where one dissolves (depth above −1, or fewer than two separate cores), is
    recorded as MERGED and is not resolved.
- **B-floor.**
  - The single-particle null speed (one particle, L = 100).
  - The difference in v_rel at d = 30 between L = 100 (N = 1024) and L = 200 (N = 2048).
  - The resolution floor is the larger of 10× the null speed and 2× that L difference. A separation is resolved when
    |v_rel| exceeds the floor.
- **B-lin, linear-response curve (reported).**
  - With the adjoint ψ† of the single particle (normalized ⟨ψ†, U′⟩ = 1) and the partner's tail
    δu_d(x) = U(x − d) − h₀, the velocity of the particle at the origin is v(d) = −⟨ψ†, (J(U) − J(h₀)) δu_d⟩. This is
    the same rule as MG192's mobility, with the partner's tail as the forcing. Positive v(d) means motion toward the
    partner, so the pair's v_rel = −2 v(d).
  - Evaluated on d ∈ [6, 40] in steps of 0.05.
- **B-eq, exact pair equilibria (reported).**
  - Newton in the even subspace, started from a symmetric pair at each zero of the B-lin curve in [8, 30]. The
    separation is free in that subspace, so Newton finds the nearby two-particle stationary state if one exists, and
    reports its separation d*.
  - Stability comes from the full spectrum. Only the exact translation mode (the odd eigenvalue nearest 0) is removed.
    The relative (breathing) mode is kept, because it decides whether the separation is stable.
  - A zero with no converged equilibrium is reported as such.

## 4. Gates and outcome rules

| Gate | Criterion | Role |
| --- | --- | --- |
| **A0 numerics** | Stationary residual ≤ 1e-10; continuation μ_a agrees with adjoint μ_a within 1e-3 relative; nulls below 1e-12 | STOP if fails |
| **A1 drift** | Every ε > 0 run has measured speed within 10% of μ_a · da/dx, and the core depth stays within 5% of its start | Gate: CONFIRMED / NOT CONFIRMED |
| **B0 numerics** | Pair Newton residuals ≤ 1e-10; single-particle null below 1e-12; B-floor determined | STOP if fails |
| **B1 sign pattern** | Uses only separations with \|v_rel\| above the floor. **OSCILLATORY**: at least 3 sign changes in d ∈ [8, 30], and every gap between consecutive zeros is 1.533 ± 25%. **MONOTONIC ATTRACTION**: all resolved v_rel < 0. **MONOTONIC REPULSION**: all > 0. **UNRESOLVED**: fewer than 10 resolved separations, or none of the above | Primary. Prediction: OSCILLATORY |
| **B2 decay law** | A least-squares fit of v_rel = e^(−κd)(A cos kd + B sin kd) over the resolved points, with each residual weighted by e^(0.4608 d) so all separations count about equally, gives κ within 10% of 0.4608 and k within 5% of 2.0491 | Gate: SUPPORTED / NOT SUPPORTED |
| B3 linear response | Zeros of B-lin compared with B-dyn and B-eq; median relative error of v_rel where resolved | Reported |
| B4 equilibria | List of d*, stable or unstable, and the alternation pattern | Reported |
| B5 net bias | Fit an added non-oscillating term C e^(−0.9657 d); report its sign and size relative to A | Reported |

**Overall labels.**
- Part A: CONFIRMED if A0 and A1 pass. If A0 fails, A1 is not scored.
- Part B: B1's label, plus B2's.

**Plain-language reading.**
- OSCILLATORY with B2 SUPPORTED means the model's particle–particle force is the ordinary short-range, sign-alternating
  interaction of Turing solitons, not a gravity-like attraction, at least in 1D.
- MONOTONIC ATTRACTION would support an SQK-style attraction in 1D, and would call for 2D and 3D follow-up.

No result is rescored after the data are seen. Deviations are recorded.

## 5. Limits

- 1D continuum, deterministic, one operating point.
- No noise, no lattice, no token engine.
- Separations are limited to [8, 30]. Below 8 the cores overlap and may merge, which is reported but not scored.
- A 1D force cannot test inverse-square falloff. The question is only attraction versus oscillating sign.

## 6. Budget, execution and records

- **Where:** Claude's cloud container, at most 4 processes, NumPy and SciPy threads set to 1.
- **Cap: 3 core-hours.** Estimate about 2:
  - Part A: about 10 minutes.
  - B-dyn: 89 runs at about 25 s each.
  - B-floor, B-lin and B-eq (about 14 dense eigenvalue solves): about 40 minutes.
- **Code:** `src/mg194_two_particle.py`, plus `src/mg194_check.py` as a separate recomputation of the B1 and B2
  classification and fit from the saved trajectories.
- **Records:** `results/MG194_two_particle/` (JSON results, trajectories, logs and report). Figure in `paper/figures/`.

## 7. PI decisions (all adopted as recommended, 2026-11-16)

- **D1.** Should Part A use new starting points (±L/16, ±3L/16) plus a separate μ_a recomputation, rather than repeating
  EX001's placements? *Recommended: yes.*
- **D2.** Should OSCILLATORY be registered as the prediction, with SQK-style monotonic attraction as the alternative?
  *Recommended: yes.*
- **D3.** Should the separation grid be d = 8 to 30 in steps of 0.25 (89 runs, about 7 wavelengths, about 6 points per
  half-wavelength)? Coarser steps make the ±25% zero-spacing test too blunt. *Recommended: yes.*
- **D4.** Is a cap of 3 core-hours, in the cloud, acceptable? *Recommended: yes.*

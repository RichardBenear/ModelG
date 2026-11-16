# MG192 results — SQK gravity mobility by linear response (continuum 1D)

**Registered outcome: UNRESOLVED.** The registered calculations were all done, but one gate failed for a reason that is
itself informative.

**What passed.**
- A 1D Model G soliton is pushed toward lower G by a source gradient, with mobility μ_a = −0.059847, as SQK predicts
  (P1 SUPPORTED).
- The value matches Pulver & LaViolette's published −0.0599 to 0.09% (P3 CONSISTENT).
- A separate calculation by a different method agrees to 4 parts in 10⁸.
- Time-dependent runs started at the point of steepest slope moved at the predicted speed, within 0.6%.

**What failed (P2).** Runs started an eighth of a period away, where the source is only 0.07% or 0.14% higher, did not
drift at all: the particle dissolved within about 50 time units. A post-hoc check shows why. At the published
parameters (a = 14, b = 29), the continuum soliton exists only up to about 0.011% above a = 14. The particle sits almost
exactly at the edge of its existence window.

Registered specification: `spec/MG192_gravity_mobility.md`, frozen and committed at `a2ba0fa` before any calculation. PI
approval: Richard Benear, 2026-11-15 21:27 America/Denver. The study was run by Claude at the PI's request.

## Registered results

| Gate | Result | Key numbers |
| --- | --- | --- |
| P0 numerical qualification | **PASS** | Residual 4.9e-12 at N = 1024 and 1.8e-11 at N = 2048; profiles agree to 9.9e-11. Translation eigenvalue −2.5e-12; least-stable other eigenvalue −0.0422. Core amplitudes within 0.08% of the published values. Mobilities agree to 7e-13 across N and 1e-10 across L = 100/200. |
| P1 sign (gate) | **SUPPORTED** | μ_a = −0.0598475 (v = μ_a · da/dx), numerical uncertainty 6.3e-12, so the sign is resolved about 10¹⁰ times over |
| P2 dynamic confirmation (gate) | **FAIL** | x₀ = 0: measured speed within 0.59% (ε = 0.001) and 0.60% (ε = 0.002) of the prediction; linearity ratio 1.9997. x₀ = L/8: particle dissolved, measured "speed" about 0. Null runs below 1e-17. Halving the time step changed the speed by 8e-8 relative. |
| P3 published benchmark (reported) | **CONSISTENT** | μ_a = −0.059847 against the published −0.0599 (difference 0.09%) |
| P4 specificity (reported) | See below | μ_w/μ_a = 1.0000000; μ_b/μ_a = −0.146 |

The rule registered in advance gives UNRESOLVED: P0 and P1 pass, but P2 fails. That label stands. The analysis below
explains the failure; it does not rescore it.

| ε | x₀ | Predicted v | Measured v | Relative error |
| ---: | ---: | ---: | ---: | ---: |
| 0.001 | 0 | −5.2645e-05 | −5.2336e-05 | 0.59% |
| 0.002 | 0 | −1.0529e-04 | −1.0466e-04 | 0.60% |
| 0.001 | 12.5 | −3.7225e-05 | −2.5e-08 | particle dissolved |
| 0.002 | 12.5 | −7.4451e-05 | −6.7e-10 | particle dissolved |
| 0 | 0 | 0 | 6.5e-19 | null |
| 0 | 12.5 | 0 | 7.2e-18 | null |

The small remaining gap at x₀ = 0 (about 0.6%) is of the order expected for a finite-wavelength correction from the
sinusoid's curvature. I estimated that order of magnitude; it was not computed.

## Separate recomputation

`src/mg192_check.py` recomputes the mobilities by a different route. It solves for a travelling soliton under a small,
windowed uniform-slope forcing, with the velocity as an unknown, using its own residual code and a finite-difference
Jacobian. The velocity derivative is taken by central differences at κ = 1e-4 and 2e-4.

| Mobility | Adjoint method | Continuation, κ = 1e-4 | Continuation, κ = 2e-4 |
| --- | ---: | ---: | ---: |
| μ_a (G source) | −0.0598474726 | −0.0598474703 | −0.0598474633 |
| μ_w (X source) | −0.0598474726 | −0.0598474703 | −0.0598474638 |
| μ_b (criticality) | +0.0087584477 | +0.0087584474 | +0.0087584473 |

These agree to better than 2e-7 relative, against the 1e-3 required. The check shares the model equations, the
reference profile file and NumPy/SciPy with the main code. It is same-author work, not external replication.

## Why the off-centre runs failed (post-hoc diagnostics, not registered tests)

1. **The particle dissolves.** `src/mg192_posthoc_decay.py` reran the four ε > 0 cases to T = 60 and recorded core depth
   (`posthoc_core_depth.json`):
   - x₀ = 0 (local a = 14.000): depth −8.358 at the start, still −8.358 at T = 60.
   - x₀ = 12.5 (local a = 14.0099 or 14.0198): depth falls from −8.358 to about 0 by T = 60. The particle is gone,
     not stuck.
2. **There is a fold just above a = 14.** `src/mg192_posthoc_fold.py` continued the uniform soliton in a, in steps of
   2e-5 relative (`posthoc_fold.json`):
   - Going up, Newton converges to +0.010% and fails at +0.012%. Meanwhile the least-stable even eigenvalue rises
     toward zero: −0.042, −0.038, −0.034, −0.029, −0.023, −0.014. That pattern is consistent with a saddle-node fold at
     about a = 14 × 1.00011.
   - Going down, the soliton continues smoothly to at least −0.4%, and deepens (core depth −9.43 at −0.4%).

So the published 1D soliton at b = 29 is about 0.011% from the edge of existence in a. Raising a moves X₀ up, which
pushes the ether further below the Turing threshold, so this is effectively the lower edge in b. The x₀ = L/8
placements sat 0.07% and 0.14% above a = 14, beyond the fold. **This is a flaw in my design:** I did not check the
existence window before choosing those placements. It is recorded here as such.

## What the specificity numbers mean (P4)

- **μ_w = μ_a exactly**, and this is an identity, not a coincidence. The G row of the adjoint equations is
  (∂ₓₓ − q)ψ_G + p ψ_X = 0. Taking its first moment gives ⟨ψ_X, x⟩ = (q/p)⟨ψ_G, x⟩, and q = p = 1 here. (This
  derivation was made after seeing the number, and the numerics confirm it to 4e-10.)
  - Consequence: a G-source gradient and an X-source gradient of the same slope push the particle equally hard. Per
    unit of background **G** slope, the X-source push is 10 times larger, because an X source raises G₀ only by
    g = 0.1 times as much.
  - So in linear response the drift is **not specific to G**. It follows source gradients that move G and X together,
    and it tracks the background X slope at least as closely as the G slope.
- **μ_b/μ_a = −0.146:** a criticality gradient, with G₀ and X₀ flat, moves the particle toward higher b (more
  supercritical ether), about 15% as strongly per unit of parameter slope.
- **Mobility per unit background G slope (a arm):** v = 0.0539 · (−dG_bg/dx). This agrees with MG191's raw,
  unqualified coarse-lattice median of about 0.050.

## Interpretation for the programme

1. **SQK's gravity direction holds in the continuum Model G equations.** A particle is pushed toward lower G by a source
   gradient. The mobility is exact to numerical precision and reproduces LaViolette's published coefficient. The
   pinning in MG191 was the obstacle, not the physics.
2. **The push is not uniquely a G effect.** Model G couples G and X so that any source gradient moves both, and an X
   gradient pushes just as hard. Calling the drift "gravity, a response to G" is an interpretation that this model
   does not single out.
3. **Fragility at the published parameters may matter more than gravity.** At a = 14, b = 29 the continuum soliton
   sits about 0.01% from the edge of its existence window. Moving uphill, or any small increase in the source, makes
   it dissolve; moving downhill makes it deeper and more robust. This plausibly ties together several earlier
   findings: MG191's fine-grid failures at a + 0.05%; MG1's slow 1D settling; and the small loss barrier in MG188
   (lattice b = 29 is only about 0.25% above the lattice's lower edge, 28.93). These links are hypotheses, not results.
4. **Not shown:** that the particle has a gravitational monopole of its own (MG001's zero G integral still stands);
   anything beyond 1D; or anything about the token engines.

## Suggested follow-ups (not authorized here)

- **Map the continuum existence window in (a, b)** by continuation. This takes minutes and tells every later study
  where a robust operating point is.
- **Repeat the P2 dynamic check** at a b inside that window, or with gradients that stay within it. Placements should
  be chosen with the window known.
- Then move on to study 2 (two-particle force against separation) on the same spectral setup.

## Resources and records

- **Compute:** one core in Claude's cloud container. Platform: Linux x86_64, Python 3.13.18, NumPy 2.5.3, SciPy 1.18.1.
  - Main run: 1413 s, including the eigen-decompositions and seven dynamic runs.
  - Separate check: 55 s.
  - Post-hoc diagnostics: about 15 minutes in total.
  - Total: about 0.7 core-hours, against the 2 core-hour cap.
- **Minor deviation:** the Newton loop's stopping tolerance (1e-12) sat below the round-off floor, so it ran its full
  40 iterations at a residual of about 5e-12 (2e-11 at N = 2048). The final residuals still meet the registered 1e-10
  limit, so no result is affected.
- **Files**, all committed:
  - code: `src/mg192_mobility.py`, `src/mg192_check.py`, `src/mg192_posthoc_decay.py`, `src/mg192_posthoc_fold.py`;
  - records: `results.json`, `run_log.txt`, `check_continuation.json`, `posthoc_core_depth.json`, `posthoc_fold.json`,
    `soliton_adjoint_L100_N1024.npz` (soliton and adjoint), `trajectories.npz` (position series).

PI closure pending.

## Addendum 2026-11-16 (after PI closure sign-off; no result changed)

The closeness of b = 29 to the 1D continuum lower edge was already on record before MG192. The Phase 2 protocol (MG003,
`spec/phase2_protocol_v1.3.4.md` §2) lists an operational continuum lower edge b_s = 28.993 in 1D, 0.024% below b = 29,
and notes that b = 29 is near the edge only in 1D (2D 27.343, 3D 26.972). The fold in a found post hoc here is the same
edge, approached in a different direction. I should have cited this, and checked it when designing P2.

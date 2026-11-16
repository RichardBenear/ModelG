# MG192 — SQK gravity mobility by linear response (continuum, 1D)

Status: **REGISTERED 2026-11-16 (UTC); frozen.** Author: Claude (reviewer, running this study at the PI's request).
PI approval: Richard Benear, 2026-11-15 21:27 America/Denver: "I approve." This file and the `GATES.md` entry were
committed to `main` before any scientific calculation. Later changes go in a dated deviation note, not edits here.

## 1. Question and motivation

MG191 was UNRESOLVED. At the strongest gradient, the coarse lattice drifted downhill in 8/8 phases, with raw mobility
about 0.05 at every ε. That was only 8 times the null drift, against the 10 times required, because lattice pinning
(MG9d) is as large as the signal at these speeds. The fine lattice could not hold the particle at the upper source
levels. The obstacle is the measurement method, not necessarily the physics.

MG192 asks the same physical question in the setting where SQK makes its claim, the continuum Model G equations. It uses
a method free of lattice pinning: **does a stationary 1D Model G soliton drift toward lower G under a weak imposed
source gradient, and with what mobility?**

## 2. Method

**2a. Discretization.** Fourier pseudo-spectral collocation on a periodic domain of length L = 100, with N = 1024 (primary)
and N = 2048 (convergence check) points. Published parameters (eq. 17): d_x = 1, d_y = 12, a = 14, b = 29, g = 0.1,
p = q = 1, s = u = w = 0, with A and B held constant. A resolved Fourier discretization has no lattice pinning beyond
spectrally small aliasing; the null runs in 2e measure what remains.

**2b. Stationary soliton.** Solve F(U) = 0 by Newton's method within the even (cosine) subspace, which removes the
translation zero mode. The initial guess is the cached h = 0.2 lattice reference (`results/MG003/ref_profiles.npz`,
key `h0.2_b29.0`), interpolated. Positivity is required.

**2c. Linearization and adjoint mode.** Build the full-space Jacobian J at U. Compute its spectrum, the translation mode
U′, and the adjoint (left) null vector ψ† (J^T ψ† = 0), normalized so that ⟨ψ†, U′⟩ = 1.

**2d. Linear-response mobility (solvability condition).** A perturbation of the right-hand side, P(x) = κ (x − x₀) e_S,
added to the equation of species S, moves the soliton at first order with speed v = −κ ⟨ψ†_S, x − x₀⟩. This is the
standard solvability result for localized patterns, valid when the gradient's length scale is much larger than the soliton
(here L = 100 against a core of order 1). Coefficients computed:

| Coefficient | Perturbation | Meaning |
| --- | --- | --- |
| μ_a | source a(x) in the G equation | Main SQK test; the same forcing as the Pulver & LaViolette gradient benchmark (source γ = m·x in the φ_G equation) |
| μ_b | b(x), entering X as −b·X and Y as +b·X | Criticality-gradient comparison (MG191 arm C) |
| μ_w | source w(x) in the X equation | X-source comparison: G₀ changes 10 times less than X₀ |

Velocities are reported as v = μ_a · da/dx. For slow modulation the background follows G₀ = a/0.9, so
dG_bg/dx = da/dx / 0.9. **SQK predicts motion toward lower G, i.e. μ_a < 0** (MG190 D3; Pulver & LaViolette eq. 28,
v = −0.0599 m).

**2e. Dynamic confirmation.** Time-dependent spectral integration (diffusion treated exactly in Fourier space, reactions
by fourth-order exponential time differencing, ETDRK4) with a(x) = 14 · (1 + ε sin(2π x/L)):
- the soliton is placed on the particle-free modulated background at x₀ = 0 (maximum slope) and x₀ = L/8 (slope
  reduced by cos 45°);
- ε = 0.001 and 0.002; plus null runs (ε = 0) at both positions; horizon T = 400.
- Position is the minimum of φ_X found by Fourier interpolation, sampled every 1 time unit. Velocity is the
  least-squares slope over t = 100–400; earlier times are a fixed transient.
- Step dt = 0.002, checked against dt = 0.001 on one run.

## 3. Registered predictions and gates

**P0 — numerical qualification (STOP if failed).**
- (a) Newton residual (max norm) ≤ 1e-10 at N = 1024 and N = 2048, and the two profiles agree within 1e-8 (max norm, at
  the shared points).
- (b) Translation eigenvalue |λ₀| ≤ 1e-7. Every other eigenvalue has Re λ ≤ −1e-4 (the soliton is linearly stable in
  the continuum).
- (c) Core amplitudes within 2% of the published 1D values (φ_X −8.36, φ_G −0.161, φ_Y 0.930). The published values
  used Dirichlet boundaries at ±50, so this is a sanity check, not a convergence test.
- (d) μ_a, μ_b and μ_w agree between N = 1024 and N = 2048 within 1e-4 relative, and between L = 100 and L = 200 within
  1% relative.

If P0(b) fails because the continuum soliton is unstable, the study stops and reports that. That would itself be a
finding.

**P1 — sign (gate).** SUPPORTED if μ_a < 0 and |μ_a| exceeds 100 times its P0(d) numerical uncertainty. CONTRADICTED if
μ_a > 0 by the same margin.

**P2 — dynamic confirmation (gate).** In each of the four ε > 0 runs, the measured velocity matches the solvability
prediction v = μ_a · da/dx(x₀) within 10%. Each null run's speed is below 1% of the smallest predicted speed. The
ε = 0.002 speeds are 2.0 ± 0.2 times the ε = 0.001 speeds (linearity).

**P3 — published benchmark (reported, not gated).** Report μ_a against Pulver & LaViolette's −0.0599. The label is
CONSISTENT if within 15%, otherwise DIFFERENT. Differences in protocol (Dirichlet boundaries, ramp source, unstated
velocity method) mean disagreement does not imply an error on either side.

**P4 — specificity (reported, not gated).** Report μ_b/μ_a and μ_w/μ_a, and express each as a drift direction relative
to the background G and X gradients that perturbation produces. This shows whether the drift follows G specifically or
follows criticality and X. No threshold is set, because there is no prior expectation.

**Overall outcome.**
- **SUPPORTED:** P0, P1 and P2 pass with μ_a < 0.
- **CONTRADICTED:** P0 and P2 pass with μ_a > 0.
- **UNRESOLVED:** P0 passes, but P1 or P2 fails.
- **STOPPED:** P0 fails.

## 4. Scope and limits (stated in advance)

- This is continuum 1D and deterministic. It tests the PDE behind SQK, not the token engines. MG2 C6 and MG9d suggest
  the lattice approaches the continuum as h → 0, but this study does not test the lattice.
- A supported result shows that a particle is pushed down an imposed source gradient. It does **not** show that the
  particle has a gravitational monopole of its own. MG001's zero-integral identity for φ_G still stands; that is
  study 3.
- A source gradient changes G and X together; μ_w and μ_b give partial separation only.
- Linear response holds for weak gradients. Large-gradient behaviour and arrival at a well minimum are not tested.

## 5. Execution, records and verification

- **Where:** Claude's cloud container (Linux, Python 3, NumPy and SciPy versions recorded), one core. Estimated well
  under one core-hour, with a cap of 2 core-hours.
- **Committed:** this specification; the code (`src/mg192_mobility.py`); a separate checker (`src/mg192_check.py`);
  the soliton profiles and adjoint vectors at N = 1024 (compact `.npz`, about 0.1 MB); trajectory position series; and
  the results report `results/MG192_gravity_mobility/MG192_results.md`.
- **Separate recomputation:** the checker recomputes μ_a by a different route, namely finite-difference continuation:
  solve for the stationary soliton under a small uniform-slope forcing in a co-moving frame (unknowns U and v), at two
  slope sizes, and take dv/dκ. It must agree with the adjoint value within 1e-3 relative. Same author, so this is not
  external replication; an independent recheck by ChatGPT is optional at the PI's request.
- No result is rescored after the data are seen. Any deviation is recorded in the report.

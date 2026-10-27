# Phase 2 Protocol v1.0: Stochastic Token Rules on a Fixed Lattice (MG003)

Status: PROPOSED 2026-10-27. Becomes frozen on PI sign-off. Gate definitions are in `GATES.md` (MG2, MG3).
Rules: `spec/rules_v1.1.md`, unchanged. Reference data: `results/MG003/`.

## 1. Substrate

- **Primary lattice:** 1D periodic ring, spacing h = 0.2, length L = 100 (N = 500 nodes). This is the
  simplest Wolfram-style spatial graph: a cycle graph with relations {v, v+1}.
- **Secondary lattice:** h = 0.1, L = 100 (N = 1000), used only for the MG2 lattice-independence check.
- **Why 1D and h = 0.2:** cost. Tau-leaping costs ≈ 0.5 s per unit time at h = 0.2 and ≈ 3 s at h = 0.1
  (2-core container). A 2D lattice at an adequate h costs ~10²× more and is deferred to Phase 3.
- **Noise parameter:** Ω = tokens per node at unit concentration. The h-independent noise strength is
  Ω_len = Ω/h (tokens per unit length). Both are reported.
- **Boundary:** periodic. It removes boundary effects on nucleation and has no reservoir artifacts.

## 2. Lattice model of Model G (deterministic references, computed in MG003)

On a lattice, Model G is a slightly different model from the continuum one. O(h²) discretization shifts
it, but its Ω → ∞ limit is exactly the lattice ODE with the same Laplacian. MG2 therefore compares
stochastic results to the **same-lattice ODE**. The lattice → continuum link is deterministic and
documented separately.

| Quantity (1D, b = 29) | continuum (MG002) | h = 0.1 | h = 0.2 |
|---|---|---|---|
| φ_X(0) | −8.3645 | −8.6165 (+3.0%) | −9.0225 (+7.9%) |
| λ₀ (registered estimator) | 3.0754 | 3.0751 | 3.0619 |
| deterministic seed threshold (scale) | 0.975 | 0.968 | 0.949 |
| Turing threshold b_c | 30.0441 | 30.0476 | 30.0458 |
| **single-soliton window [b_s, b_r]** | b_s = 28.9932 (1D); 27.3434 (2D); 26.9718 (3D); b_r not computed | [28.9782, 29.4119] | [28.9286, 29.3967] |
| LNA Ω·Var(c) per node, homogeneous (G, X, Y) | — | 15.57, 63.46, 2.29 | 15.59, 112.7, 2.73 |

Sources: `results/MG003/lattice_reference.txt`, `existence_window.txt`, `replication_edge.txt`.

**New finding (MG003):** at the published parameters, the stationary single soliton exists only in a
narrow band of b:
- Below **b_s**, a soliton continued from b = 29 decays to the homogeneous state.
- Above **b_r**, it self-replicates and the pattern invades the domain. At b = 29.5, 21 minima form by
  t = 400.
- In the 1D continuum, b_s = 28.993, so the published b = 29 sits **0.02%** above the lower edge. That fits
  MG1's slow 1D relaxation (critical slowing near a saddle-node) and the 0.975 seed threshold.
  "Exists" means φ_X(0) < −1 at t = 300 after seed-free continuation. Slow decay near the saddle-node
  can only bias b_s low, so 28.993 ≤ b_s ≤ 29.
- In the 2D and 3D continuum the window is much wider: b_s = 27.343 (2D) and 26.972 (3D), so b = 29 sits
  6.1% and 7.5% above the lower edge there. Being near the edge is a **1D-specific** feature, which
  matters for Phase 2 being 1D. In higher dimensions, the published point sits comfortably inside the
  window.

This window is the deterministic backbone of MG3. The homogeneous state and the soliton are both stable
only inside it, and that bistability is what noise-driven nucleation requires.

## 3. Stochastic simulation

- **Rules:** R1, R1r, R2, R2r, R3, R4, R5 and hopping D_S exactly as in rules v1.1 (propensities in
  dimensionless units).
- **Seed forcing (MG2-C5 only):** an external rule R_seed {v,X} → ∅ with zero-order propensity
  Ω·|χ(x_v, t)|, χ from eq. 18 (scaled). Its mean field is the published seed. It is a protocol forcing,
  not an ether rule.
- **Algorithm:** binomial/Poisson tau-leaping (`src/lattice.py::tau_leap_1d`), Lie splitting
  (reactions, then hops).
  - First-order losses are binomial, so counts can't go negative.
  - R4 firings are Poisson, capped at n_Y.
  - Hops: Binomial(n, 2p_S), split 50/50, with p_S = d_S τ/h².
- **Step:** τ = 0.4 h²/d_Y (p_Y = 0.4 per direction): τ = 1.333×10⁻³ at h = 0.2.
  - **Contingency (pre-declared):** if validation V2 fails, all gate runs use τ/2.
- **Exact reference for validation:** Gillespie direct-method SSA on a small ring (V1).
- **Random seeds:** run k uses seed 1000·(stage id) + k. They are recorded in the result files.

## 4. Pre-gate validations (must pass before any MG2/MG3 gate run)

- **V1: tau-leap vs exact SSA.** Ring N = 30, h = 0.2, Ω = 3 (strongly discrete), homogeneous start,
  T = 500, first 50 discarded.
  - Measure single-node mean and variance of G, X, Y, and the nearest-neighbour covariance of X.
  - Run tau-leaping at τ and τ/2 against the SSA.
  - **Pass:** every statistic agrees with the SSA within max(3%, 3 batch-means standard errors).
- **V2: step-size convergence.** MG2 configuration (b = 29, h = 0.2), Ω = 1000, 4 runs each at τ and τ/2.
  - **Pass:** aligned-mean φ_X profiles differ by < 0.5% (L2, |x| < 10), and far-field Ω·Var(φ_X) differs
    by < 5%.
- **V3: LNA check of the homogeneous state** (part of V1 data at Ω = 3, and MG2-C3 at large Ω). Reported
  only.

## 5. Measurement definitions

- **Smoothing:** Gaussian, σ = 0.4 (2 nodes), applied to φ_S = n_S/Ω − S₀ for detection and alignment only.
  Profiles are averaged unsmoothed.
- **Particle position:** sub-node location of the minimum of smoothed φ_X (parabolic fit through 3 nodes).
- **Aligned mean profile:** each recorded snapshot is shifted (linear interpolation, periodic) so the
  particle sits at x = 0, then averaged over snapshots and runs.
- **Profile error:** E(Ω) = max over S ∈ {G, X, Y} of ‖⟨φ_S⟩ − φ_S*‖₂ / ‖φ_S*‖₂ on |x| < 10, where φ*
  is the same-lattice ODE soliton. The statistical part σ_E is estimated by bootstrap over runs, and the
  bias part is E_b = √max(E² − σ_E², 0).
- **Particle (for detection):** a local minimum of smoothed φ_X below −0.5·|A*|, with A* the lattice-ODE
  φ_X(0) (−9.02 at h = 0.2), persisting ≥ 20 time units. Nucleation time = first time the persistent
  minimum crosses the threshold.
- **Far-field variance:** Ω·Var(φ_X) over nodes with |x − x_particle| > 20, time-averaged.
- **Recording:** snapshots every 0.5 time units.

## 6. Compute budget (2 cores)

| Block | Runs | Approx. core-hours |
|---|---|---|
| V1 + V2 | SSA + ~12 tau runs | 1.0 |
| MG2 (2 b-values × 7 Ω × 6 runs × T = 100) | 84 | 1.2 |
| MG2-C5 seeded nucleation (3 Ω × 20 runs × 2 scales, T = 60) | 120 | 1.0 |
| MG2-C6 (h = 0.1) | 8 | 0.7 |
| MG3 Stage A + B + C4 + null | ≤ 200 | ≤ 6 |

Total ≤ 10 core-hours, about 5 hours of wall time. MG3 is capped by its own stopping rules.

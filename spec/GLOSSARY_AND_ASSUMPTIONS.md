# Glossary, Definitions and Assumptions

Living reference document (not a frozen spec). Created MG004f, 2026-10-28, at the PI's request.
It records what every symbol means in this project, and every assumption made so far, with where it was
made and its status. When this file and a frozen spec disagree, the frozen spec governs, and the
disagreement is a bug to report.

---

## 1. Notation warning: two different "Ω"

| Symbol | Meaning | Where |
|---|---|---|
| **Ω (system size)** | Tokens per node at unit concentration (§2). **This is the Ω in all Phase 2 work, gates and plots.** | Rules v1.1 onward |
| **Ω (LaViolette sink species)** | The product of reaction X → Ω (rate k₅) in Model G. It is a sink: produced X is simply removed. Its reverse rate k₋₅ appears in the dimensionless constant w = √k₄·k₋₅Ω/K^{3/2}, which is 0 in the reference set. | `spec/model_g_published.md` §1–3 |

The sink species never appears as a token, because tokens that reach it are deleted (rule R5). Wherever
the context is the simulation, Ω means system size.

---

## 2. Tokens, etherons and Ω

**Etheron.** LaViolette's elementary ether entity: types A, B, G, X, Y, Z and the sink Ω. In SQK they
react and diffuse, and their concentration patterns are physical particles and fields.

**Token.** In the rewrite-rule simulation, **one token = one etheron** of type G, X or Y, sitting on a
graph node.
- Formally it is a unary relation {v, S} attached to node v (rules v1.1 §1).
- Reactions turn tokens of one type into another; diffusion moves a token to a neighbouring node.
- A and B are **not** tokens. They are uniform reservoirs whose constant concentrations are absorbed
  into rate constants (assumption M3).
- Z and the sink Ω are not tokens either; tokens sent there are deleted.

**Ω (system size).** The conversion factor between LaViolette's continuous concentrations and whole
token counts: n_S(v) = Ω · c_S(v).
- n_S(v) is the number of species-S tokens at node v; c_S is the dimensionless concentration in
  LaViolette's equations.
- Ω is **not** the total number of tokens, and it is **not** a quantity LaViolette defines. It is the
  knob this project adds when it makes the continuum ether discrete.
- At the uniform steady state (b = 29), each node holds on average:

  | Species | Tokens per node |
  |---|---|
  | G | G₀·Ω ≈ 15.6 Ω |
  | X | X₀·Ω ≈ 15.6 Ω |
  | Y | Y₀·Ω ≈ 1.86 Ω |
  | **Total** | **≈ 33 Ω** |

- On the 500-node ring that is ≈ 16,500 Ω tokens in all: about 1.65 million at Ω = 100, and
  ≈ 1.65 × 10¹² at Ω = 10⁸.

**Noise strength** [corrected MG004g]. In the weak-noise (linear-noise) regime, concentration
fluctuations scale as **Ω^(−1/2)**. That is the general statement here. "√n" fluctuations are only the
Poisson benchmark: this interacting system's homogeneous X variance is substantially *larger* than
Poisson (LNA Ω·Var(X) ≈ 113 vs Poisson 15.6 at h = 0.2).
- **Small Ω** means strong intrinsic noise.
- **The limits are separate** [corrected MG004g]:
  - at fixed h and τ, Ω → ∞ gives the **forward-Euler lattice map**;
  - recovering the continuous-time **lattice ODE** additionally needs τ → 0;
  - recovering LaViolette's **continuum PDE** additionally needs spatial refinement (h → 0).
  - Stationary states are shared exactly between the Euler map and the lattice ODE (when no clipping
    or capping is active), which is why stationary profiles are compared directly.
- MG2 tests the Ω → ∞ convergence; MG1 covers the continuum PDE.

**Ω/h (tokens per unit length).** Ω is per node, and the node spacing h is our choice. The
lattice-independent noise parameter is Ω/h (in 1D; Ω/h^d in d dimensions). Both are reported.

**Ω values used so far:**

| Where | Ω values | Why |
|---|---|---|
| V1 (exact SSA comparison) | 3, 30 | strongly and moderately discrete |
| V3 (time-step / clip pilot) | 1, 2, 3, 5, 7, 10, 15, 20, 30, 50, 100, 300, 1000, 3000, 10⁴ | every Ω used anywhere |
| V2 (time-step equivalence) | 300 (v1.3.1, aborted / FAIL record); **10⁸** (v1.3.2, PASS) | the decision point |
| V4 (alignment) | 100, 10³, 10⁴ (v1.3.1); **10⁴, 10⁶, 10⁸** (v1.3.2, PASS at 10⁸) | the decision point |
| V6 (assumption-A1 pilot) | 300, 3000 | real survival-conditioned runs |
| V7 (detector characterization) | 1, 2, 3, 5, 7, 10, 15, 20, 30, 50, 100, 300 | MG3's range |
| MG2 ladder (registered, v1.3.2) | 10, 30, 100, 300, 10³, 3×10³, 10⁴, 10⁵, 10⁶, **10⁸** (decision point for C1/C6) | convergence and bias scaling |
| MG3 grid (registered v1.3) | 300, 100, 50, 30, 20, 15, 10, 7, 5, 3, 2, 1 | nucleation scan |
| MG3 grid (Amendment 3, approved; v1.3.3) | 300, 250, 200, 150, 120, 100, 80, 65, 50 | V7's transitional range, with crowding exclusion |

**Rough physical calibration of Ω.** This is an order-of-magnitude estimate, not a project result, and
it has heavy caveats.
- LaViolette's nucleon paper gives an etheron density of over 10²⁵ per cubic fermi inside a nucleon.
  This is as transcribed from that paper during MG001's source review; re-verify against the paper
  before citing it.
- Suppose the Turing wavelength λ₀ = 3.08 (dimensionless) equals the nucleon Compton wavelength,
  1.32 fm. Then the length unit is L ≈ 0.43 fm, a lattice cell (h = 0.2 L) is ≈ 0.086 fm, and a 3D
  cell has volume ≈ 6 × 10⁻⁴ fm³.
- That gives ≳ 6 × 10²¹ etherons per cell, spread over a total dimensionless concentration of ≈ 33,
  so **Ω_physical ≳ 10²⁰**.
- **Caveats:** it mixes a 3D physical density with our 1D lattice; it assumes the Compton-wavelength
  identification; and it treats all etheron species alike.
- **Implication if it holds:** real-world intrinsic noise would be utterly negligible in this model.
  The intermediate-Ω regime where spontaneous nucleation occurs (≈ 100–300 at b = 29, per V7) would
  require etherons far sparser than LaViolette's estimate, or an effective noise source other than
  etheron counting. This is an interpretive constraint on MG3, not a test result.

---

## 3. Model and lattice quantities

| Term | Definition | Value or source |
|---|---|---|
| **Model G** | LaViolette's 3-variable reaction–diffusion ether model (G, X, Y), with the G ⇌ X back-reaction | `spec/model_g_published.md` |
| **Reference parameters (eq. 17)** | d_x = 1, d_y = 12, a = 14, b = 29, g = 1/10, p = q = 1, s = u = w = 0 | Pulver & LaViolette 2013 |
| **Dimensionless units** | Time T = 1/(k₋₂ + k₅); length L = √(D_G T); concentration C = 1/√(k₄ T) | §3 of the published spec |
| **φ_S** | Potential: φ_S = S − S₀ (deviation of species S from the uniform steady state) | eq. 12 |
| **S₀** | Uniform steady state: G₀ = X₀ = 140/9 ≈ 15.556, Y₀ ≈ 1.864 (b = 29) | |
| **b** | Control parameter (k₃B/K): how strongly B drives X → Y. **The main dial in Phase 2.** | 29 published |
| **b_c** | Turing threshold: above it the uniform state is linearly unstable | 30.044 (continuum); 30.046 (lattice h = 0.2) |
| **b_s** | Lower edge of the single-soliton window: below it a particle decays | 28.993 (1D), 27.343 (2D), 26.972 (3D) continuum; 28.929 (lattice h = 0.2). *Operational estimates.* |
| **b_r** | Upper edge: above it a particle self-replicates and fills space | 29.397 (lattice h = 0.2); 29.412 (h = 0.1) |
| **b_mid** | Centre of the h = 0.2 single-soliton window | 29.16 |
| **Soliton / particle** | Stationary localized dissipative structure: an X-well with Y bump and G-well, surrounded by decaying rings | MG1 |
| **λ₀ / λ*** | Turing wavelength, i.e. twice the mean zero spacing (estimator-dependent at about 0.3%) | 3.08 published; 3.06 lattice h = 0.2 |
| **A*** | Lattice-ODE φ_X(0), the unsmoothed particle depth | −9.02 (b = 29, h = 0.2) |
| **A*_s** | Depth of the same reference after detector smoothing (σ = 0.4) | −5.54 (b = 29); detection threshold −0.5·\|A*_s\| = −2.77 |
| **Seed** | Published Gaussian kick χ(r, t) added to the X equation to nucleate a particle (eq. 18) | Threshold scale: 0.975 (1D), 0.572 (2D), 0.486 (3D) continuum |
| **Ring / lattice** | 1D periodic cycle graph, the simplest Wolfram-style spatial graph | N = 500 nodes, spacing h = 0.2, length L = 100 |
| **h** | Node spacing in units of L | 0.2 primary; 0.1 secondary |
| **τ** | Simulation time step | 6.667 × 10⁻⁴ (h = 0.2); 1.667 × 10⁻⁴ (h = 0.1); τ/2 for Ω ≤ 5 (V3) |
| **Lattice ODE** | The deterministic Ω → ∞ limit of the token rules on the same ring. MG2's reference. | `src/lattice.py` |
| **Integral identity** | With ordinary diffusion, a stationary localized soliton has ∫φ_G = ∫φ_X = 0 in infinite space (if q − gp ≠ 0) | MG001a–b; confirmed numerically in MG002 |

---

## 4. Simulation and measurement terms

| Term | Definition |
|---|---|
| **Rules R1–R5, R1r, R2r** | Token rewrite rules for each Model G reaction step (rules v1.1). R4 is the autocatalytic step 2X + Y → 3X, with propensity n_X(n_X − 1)n_Y/Ω². |
| **Tau-leaping (unsplit multinomial)** | `tau_leap_euler_1d`: every token independently picks an exit (react, hop left, hop right, stay) each step, with probability rate × τ. As Ω → ∞ it is exactly forward Euler of the lattice ODE. |
| **Clip** | A node-step where exit probabilities summed to more than 1 and were renormalized. Allowed rate ≤ 10⁻⁴; τ is halved per Ω if exceeded (V3). |
| **SSA** | Gillespie's exact stochastic simulation algorithm, used only to validate tau-leaping (V1) |
| **Detection** | A local minimum of smoothed φ_X (Gaussian σ = 0.4) below −0.5·\|A*_s\| |
| **Track / events** | Detections linked over time. **Loss**: the original particle is gone ≥ 20 time units. **Birth**: a new track persisting ≥ 20. **Replication**: a birth within 2λ* of an existing particle. **Merge**: two tracks join. |
| **Eligible run** | A run with no loss, replication or merge before t = 100, and ≥ 20 single-particle snapshots in each of the 6 blocks. Profile statistics use eligible runs only (survival-conditioned). |
| **Survival** | Whether the original particle persists. Reported separately from profiles. |
| **Aligned profile** | Snapshot shifted (Fourier interpolation) so the particle sits at x = 0; window \|x\| < 10 |
| **Block profile** | Mean aligned profile over one of 6 consecutive 15-time-unit blocks of a run |
| **Far-field variance** | Ω·Var of φ_S at nodes more than 20 from any detected particle; compared with the LNA |
| **LNA** | Linear noise approximation: the analytic prediction of fluctuation size about the uniform state (continuous-time and the scheme's exact discrete-time version) |
| **Confirmation** | Field-level deterministic test (v1.3.3, `src/confirm.py`). The seed-free lattice ODE runs from the field at detection + 20; the outcome is "deterministically persistent particle-depth structure over continuation times 200–300" (identity-linked across 11 checks), excursion, or unresolved. Interpreted as basin membership, with a finite-horizon limitation. Validated by controls against long-run truth. |
| **Excursion** | A detection that relaxes back to the uniform state, i.e. noise that did not create a particle |
| **Stochastic Turing pattern** | A noise-sustained Turing-like pattern below the deterministic threshold. A **candidate interpretation** of V7's crowded detector regime (Ω ≤ 50), **pending E1**. Established in the Brusselator (Biancalani, Fanelli & Di Patti, arXiv:0910.4984): a precedent, not verification for Model G. |

## 5. Statistical terms

| Term | Definition |
|---|---|
| **E_total** | Relative L2 difference between the run-mean profile and the reference. Inflated by noise. |
| **U-statistic Ê²_b** | Unbiased estimate of the squared relative bias, from cross-products of different runs. Can be negative. |
| **Bias bounds** | [lower, upper] on the true bias via the triangle inequality and a t-radius. Valid under A1. |
| **Three-way decision** | PASS if upper < tolerance; FAIL if lower > tolerance; otherwise INCONCLUSIVE. Insufficient precision is never a pass. |
| **Welch two-sample bound** | Bounds on the difference between two ensembles, normalized by a deterministic reference norm |
| **α** | Exponent of the bias decay, bias ∝ Ω^(−α). Expected 1 asymptotically, if the leading coefficient is non-zero. |
| **A1 diagnostics / pooled audit** | Per-cell normality and outlier checks, and a gate-wide outlier check across all runs |

## 6. Process terms

| Term | Definition |
|---|---|
| **Stage (MGxxx)** | One unit of work, logged in `STAGE_LOG.md`. Lettered sub-stages (MG004a…) are corrections or continuations. |
| **Gate (MG1–MG7)** | A pre-registered pass/fail test with criteria committed before its run (`GATES.md`) |
| **Validation (V1–V7)** | Pre-gate checks of the tools. V7 is descriptive. |
| **Amendment** | A PI-approved change to a frozen protocol, applied as a new version (v1.3.1, v1.3.2, …). Earlier versions and records are preserved. |
| **PASS by exception** | A gate failing one registered criterion, accepted by the PI with the failure kept on record (MG1) |

---

## 7. Assumptions register

Status: **Adopted** (in force), **Checked** (tested and supported), **Open** (untested or under
review), **Retired**.

### Modelling assumptions

| # | Assumption | Introduced | Status / evidence |
|---|---|---|---|
| M1 | Model G as published by Pulver & LaViolette (2013), eq. 17 parameters, is the reference physics | MG001 | Adopted; reproduced in MG1 (14 of 15 criteria; 1D stationarity by PI exception) |
| M2 | Only the X → G reverse contribution is active; the other reverse terms vanish as combinations (k₋₃Z = k₋₅Ω_sink = 0, k₋₁ = k₋₄ = 0) | MG001a | Adopted (forced by the published parameters) |
| M3 | Reservoirs A and B are uniform and constant, absorbed into rate constants; Z and the sink are not tracked | Rules v1.0 | Adopted |
| M4 | One token = one etheron; a reaction = a local rewrite of tokens on one node; diffusion = a token hop to a neighbour | Rules v1.0 | Adopted (definition of the mapping) |
| M5 | Rate constants and the R4 falling-factorial propensity n_X(n_X − 1)n_Y/Ω² follow standard chemical-master-equation (mass-action) kinetics | Rules v1.1 | Adopted; it gives LaViolette's equations exactly as Ω → ∞ (verified MG003a) |
| M6 | Ordinary (conservative) diffusion, as in standard Model G, not LaViolette's δ(r)∇² core-broadening variant | MG001a | Adopted. Consequence: ∫φ_G = 0 for a stationary particle, so no net G "mass" |
| M7 | Space is a fixed 1D periodic ring (h = 0.2, L = 100) in Phase 2; irregular and evolving graphs come in Phases 3–4 | MG003 | Adopted (chosen for cost) |
| M8 | Continuous-time stochastic dynamics define the model; tau-leaping is a numerical approximation of them | MG003 | Adopted; V1 (exact SSA) and V2 (time-step check at Ω = 10⁸) PASS |
| M9 | **Ω is a free parameter** of the discretization, not taken from LaViolette | MG003 | Adopted. Physical scale only roughly estimated (§2, Ω ≳ 10²⁰), with caveats |
| M10 | The single-particle state at b = 29 is the target object; b_mid = 29.16 is used where robustness is needed | MG003 | Adopted. b = 29 is near the window's lower edge in 1D (not in 2D/3D) |
| M11 | **Node spacing h = 0.2 is a free parameter**, chosen for cost (MG003): about 15 nodes per wavelength; the deterministic particle is 8% off the continuum (3% at h = 0.1, 1% at 0.05); measured cost about 4× per unit time at h = 0.1. **Two readings:** (i) numerical discretization of continuum Model G (truth at h → 0), or (ii) the physical node spacing of the space graph (Wolfram reading; a model in its own right). **Effects:** the window edges move with h (b_s = 28.929 / 28.978 / 28.993 at h = 0.2 / 0.1 / continuum), so near-edge results at b = 29 (e.g. particle survival) are h-dependent; seed threshold 0.949 vs 0.975; noise comparisons across lattices should use Ω/h. MG2 is h-neutral by design (same-lattice reference), with C6 checking h = 0.1 | MG003; documented MG005 | Adopted. Suggested mitigation: report b relative to each lattice's window, and repeat key MG3 cells at h = 0.1 |

### Numerical assumptions

| # | Assumption | Introduced | Status / evidence |
|---|---|---|---|
| N1 | MG1 continuum solutions are converged (grid halving < 0.5%; tolerance check) | MG002 | Checked; independently replicated by ChatGPT (within 0.0012%) |
| N2 | τ chosen by bounded selection: exact fixed points, discrete LNA within 2%, exit-probability margin | MG003a | Checked (`tau_selection*.txt`) |
| N3 | Clipping is negligible (≤ 10⁻⁴ per node-step), with τ/2 for Ω ≤ 5 | MG004 | Checked (V3) |
| N4 | Persistence of a particle-depth structure over continuation times 200–300 indicates basin membership (confirmation) | MG003a; revised MG004g | **Checked on controls** (all agree with t = 2000 truth). Finite-horizon limitation stated; field-level, not per-detection |
| N5 | Window edges b_s and b_r found by finite-time continuation | MG003 | Adopted as operational estimates, not bifurcation points |

### Statistical assumptions

| # | Assumption | Introduced | Status / evidence |
|---|---|---|---|
| S1 (= A1) | Eligible runs' aligned profiles are i.i.d., approximately Gaussian, with no rare discrete run-level modes beyond tracked events | MG003c | **Checked on real runs** (V6 PASS; V4 diagnostics pass). Residual risk (unseen rare modes) accepted by the PI |
| S2 (= A2) | Runs are independent (distinct random streams) | MG003c | Adopted (by construction) |
| S3 | 12 eligible runs per decision cell are adequate | MG003b | Calibrated for decision error rates; power depends on Ω (Amendment 2 moved decisions to Ω = 10⁸) |
| S4 | Bias scales as B(τ)/Ω asymptotically | MG004c | Conditional (needs B ≠ 0); does not fix the fitted exponent over a finite ladder |

### Interpretive assumptions (not tested by gates)

| # | Assumption | Status |
|---|---|---|
| I1 | A Model G soliton corresponds to a physical particle (LaViolette's interpretation) | Adopted as a working hypothesis; not tested |
| I2 | A lower ∫G is "mass" | **Open (PI decision at Phase 5).** In standard Model G the net ∫φ_G of a stationary particle is 0 |
| I3 | Intrinsic token-count noise represents the "fluctuating ether" from which particles nucleate | Open. *Conditional on the unverified density figure and the stated assumptions:* if Ω_physical ≳ 10²⁰ (§2), this noise would be negligible in reality |
| I4 | A 1D ring captures the qualitative physics | Open. Some effects are 1D-specific (e.g. b = 29 near the window edge only in 1D) |

---

## 8. Where things live

- Frozen specs: `spec/` (current protocol `phase2_protocol_v1.3.2.md`; amendments `phase2_amendment_*.md`)
- Gates: `GATES.md`. Stage history: `STAGE_LOG.md`. Reviews: `review/`.
- Code: `src/`. Results: `results/MGxxx/`. Configs: `configs/`.

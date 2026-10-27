# Phase 2 Protocol v1.2: Stochastic Token Rules on a Fixed Lattice (MG003b)

Status: **SUPERSEDED by `phase2_protocol_v1.3.md` (MG003c)** before sign-off. Originally: PROPOSED 2026-10-27. Supersedes v1.1 (never frozen) after ChatGPT review round 2
(`review/MG003_chatgpt_r2.md`, response `review/MG003_response_r2.md`). Becomes frozen on PI sign-off.

Changes from v1.1 (MG003b):
- §3.2 wording: the exact-Euler property holds when exit probabilities are valid and caps are inactive.
- §5.3 rebuilt as a calibrated inference layer: the U-statistic is kept as the point estimate; decisions
  use t-based bounds with three-way PASS / FAIL / INCONCLUSIVE outcomes. A calibration study
  (`src/calibrate_bounds.py`, `src/calibrate_diff.py`, `results/MG003b/`) sets **12 runs per cell**.
- The α CI uses a corrected run-level bootstrap (distinct-pair U-statistic, pre-selected cells) with a
  calibrated widening factor.
- V2 is an equivalence test with predefined margins.
Rules: `spec/rules_v1.1.md`, unchanged. Gate definitions: `GATES.md` (MG2, MG3).

Changes from v1.0:
- new integrator (§3)
- bounded τ selection (§3.3)
- unbiased bias estimator and alignment validation (§5)
- tracking, survival and selection rules (§5)
- first-passage model (§6)
- fully specified MG3 adaptive design and budget (`GATES.md`)
- deterministic event confirmation (§5.4)
- RNG stream assignment (§7)

## 1. Substrate (unchanged from v1.0)

- **Primary lattice:** 1D periodic ring, h = 0.2, L = 100 (N = 500). A cycle graph, the simplest
  Wolfram-style spatial graph.
- **Secondary lattice:** h = 0.1, L = 100 (N = 1000), used for the MG2-C6 lattice-independence check.
- **Noise parameter:** Ω (tokens per node at unit concentration). Ω_len = Ω/h is reported alongside it.
- **Boundary:** periodic. All distances on the ring are circular distances.

## 2. Deterministic references (unchanged; `results/MG003/`)

| Quantity (1D, b = 29) | continuum | h = 0.1 | h = 0.2 |
|---|---|---|---|
| φ_X(0) | −8.3645 | −8.6165 | −9.0225 |
| seed threshold (scale) | 0.975 | 0.968 | 0.949 |
| Turing threshold b_c | 30.0441 | 30.0476 | 30.0458 |
| single-soliton window [b_s, b_r] (operational estimates) | b_s = 28.993 | [28.978, 29.412] | [28.929, 29.397] |
| continuous-time LNA Ω·Var (G, X, Y) | — | 15.57, 63.46, 2.29 | 15.59, 112.69, 2.73 |

**Window edges are operational estimates, not established bifurcation points.**
- b_s: seed-free continuation from the b = 29 soliton, "exists" if φ_X(0) < −1 at t = 300.
- b_r: the standard seed gives more than one deep minimum at t = 400.

Continuum lower edges: 28.993 (1D), 27.343 (2D), 26.972 (3D). b = 29 sits close to the edge only in 1D.

## 3. Stochastic integrator

### 3.1 Why v1.0's scheme was replaced

v1.0 used Lie-split, exponential-loss/explicit-gain tau-leaping. Its Ω → ∞ limit is a discrete-time map
whose fixed points are shifted from the lattice ODE's by O(τ). The homogeneous X and Y shift by +2.0% and
−3.9% at τ. Started from the stationary lattice-ODE soliton, that map **loses the soliton** at both b and
at τ/2 (E = 108% / 103% / 107% / 103%). ChatGPT found this, and MG003a reproduced it exactly
(`src/audit_mean_map.py`, `results/MG003a/audit_split_scheme.txt`). Increasing Ω cannot fix it.

### 3.2 Unsplit multinomial tau-leaping (`src/lattice.py::tau_leap_euler_1d`)

Every channel (reactions and hops) is drawn from the **start-of-step state**. Each token chooses
independently among its exits, with probability rate × τ for each:

| Species | Exits (per-token probability) |
|---|---|
| G | R2 → X: pτ; R1r → ∅: (q−p)τ; hop left / right: τ/h² each |
| X | R2r → G: gτ; R3 → Y: bτ; R5 → ∅: (1−g)τ; hop left / right: d_x τ/h² each |
| Y | R4 → X: n_X(n_X−1)τ/Ω²; hop left / right: d_Y τ/h² each |

Sources: R1 adds Poisson(Ω a τ) new G per node. R_seed (MG2-C5 only) removes Poisson(Ω|χ|τ) X tokens,
capped at the X tokens that stayed.

**Property** (confirmed independently by ChatGPT, round 2). When all exit probabilities are valid
(no renormalization) and the seed cap is inactive, the conditional expected increment is
E[Δc | n] = τ[F(c) + (0, −XY/Ω, +XY/Ω)]. The finite-Ω term comes from the intended falling-factorial R4
propensity. As Ω → ∞ the mean map is exactly forward Euler of the lattice ODE, so its fixed points are
identical to the ODE's, and stationary solitons are exact fixed points. Renormalized ("clipped") exits
change the map and are counted. Fixed-point preservation does **not** by itself establish transient
accuracy, or stability at an arbitrary τ. Those are covered by §3.3, V2 and V5.

### 3.3 Bounded τ selection (deterministic and analytic; done in MG003a)

Candidates are τ_k = 0.4h²/d_Y / 2^k, for k = 1, 2, 3. The first one meeting all three tests below at both
b = 29 and b = 29.16 is selected:
- (i) the mean map started from the lattice-ODE soliton keeps it: E < 10⁻⁶ at t = 100
- (ii) the exact discrete-time LNA (`lna_discrete_euler`, which includes O(τ²) multinomial terms) is
  within 2% of the continuous-time LNA
- (iii) the maximum per-token exit-probability sum on the deterministic soliton is ≤ 0.8.

If no k ≤ 3 passes, the protocol fails and escalates to the PI.

**Result** (`results/MG003a/tau_selection*.txt`):

| Lattice | Selected τ | Mean-map E | LNA deviation (G / X / Y) | Max exit-prob sum |
|---|---|---|---|---|
| h = 0.2 | 6.667×10⁻⁴ (k = 1) | 1×10⁻¹⁴ | 0.00 / −0.64 / −0.83% | 0.65 |
| h = 0.1 | 1.667×10⁻⁴ (k = 1) | 1×10⁻¹⁴ | 0.00 / −0.20 / −0.12% | 0.46 |

**Ω-dependent refinement (fluctuation clipping):** at small Ω, X fluctuations can push the R4 exit
probability over the limit. At Ω = 3 the clip fraction was 1.4×10⁻³ of node-steps. Rule: for each Ω,
the clip fraction (clips per node-step, from the V3 pilot) must be ≤ 10⁻⁴. If not, use τ/2, then τ/4.
If it still fails at τ/4, that Ω is excluded and the exclusion is reported.

Cost: about 0.4–0.9 s per unit time at h = 0.2 (2-core container), and about 4× that at h = 0.1.

## 4. Pre-gate validations (all must pass before any gate run)

- **V1 — exact SSA comparison.** Ring N = 30, h = 0.2, b = 29, Ω = 3 and Ω = 30, homogeneous start,
  T = 500, t < 50 discarded. Gillespie direct-method SSA vs `tau_leap_euler_1d` (at the τ that the
  clip rule selects for that Ω).
  - Measure single-node mean and variance of G, X, Y, and the nearest-neighbour covariance of X.
  - **Pass:** each statistic agrees within max(3%, 3 batch-means SE).
- **V2 — τ-halving equivalence on the MG2 statistics.** b = 29, Ω = 300, 12 runs each at τ and τ/2
  (MG2 setup). This is an equivalence test with predefined margins, not a non-rejection test.
  - **Profiles:** for each of G, X, Y, the upper 95% bound on ‖μ_τ − μ_{τ/2}‖/‖φ*‖ (two-sample bound,
    §5.3) must be < 0.25%, a quarter of MG2-C1's tolerance.
  - **Fluctuations:** for each of G, X, Y, the 90% CI (Welch t, run-level values) of the difference in
    far-field Ω·Var must lie within ±5% of the τ value.
  - PASS if all hold; FAIL if any point estimate lies outside its margin; INCONCLUSIVE otherwise. If
    INCONCLUSIVE, one extension of 12 + 12 runs is allowed, after which INCONCLUSIVE counts as FAIL.
- **V3 — clip pilot.** One run of T = 20 per Ω in each gate's Ω set, homogeneous start. This sets the
  per-Ω τ (§3.3). Only the clip counter is read from these runs.
- **V4 — alignment validation.** Synthetic ensembles: φ* plus a noise field taken from independent
  *homogeneous-state* runs at the same Ω and b, with random circular shifts. These go through the full
  detection/alignment/averaging pipeline (§5) at Ω ∈ {100, 1000, 10000}.
  - 12 synthetic runs per Ω. Report the alignment-induced bias (point estimate and bounds) and its
    Ω-scaling.
  - **Pass:** the upper 95% bound of the alignment-induced bias at Ω = 10⁴ is < 0.3% for G, X, Y (a third
    of MG2-C1's tolerance). FAIL if the lower bound is > 0.3%; INCONCLUSIVE (blocking) otherwise.

## 5. Measurement definitions

### 5.1 Detection and tracking

- **Smoothing:** Gaussian with σ = 0.4 *physical units* (2 nodes at h = 0.2, 4 at h = 0.1). It is used
  for detection and alignment only.
- **Particles:** local minima of smoothed φ_X below −0.5·|A*| (A* = the lattice-ODE φ_X(0) at that b and h).
  Minima closer than λ*/2 (circular) are merged and the deeper one is kept.
- **Tracking** (snapshots every 0.5): link each particle to the nearest particle of the previous snapshot
  within circular distance 1.5 (greedy, nearest first).
  - Unmatched for ≥ 20 time units → **loss**.
  - A new particle persisting ≥ 20 → **birth**. A birth within 2λ* of an existing particle counts as
    **replication**.
  - Two tracks linking to one particle → **merge**.
- **Position:** sub-node position from a parabolic fit to the smoothed minimum.

### 5.2 Alignment and averages

- Each snapshot is shifted by Fourier (band-limited) interpolation on the periodic lattice so the tracked
  particle sits at x = 0. This replaces v1.0's linear interpolation, which smooths.
- **MG2 averages are survival-conditioned.** Only snapshots where the original track exists and is the
  only particle are used. Survival (time to loss, replication or merge) is reported separately for every
  (b, Ω). MG2-C4 is unconditional.
- Run-level profile m_i = the time average of run i's eligible aligned snapshots (t ≥ 10). A run with
  fewer than 20 eligible snapshots is excluded from profile statistics, but counted in survival.
- **Far-field variance:** nodes at circular distance > 20 from every particle, time-averaged. Reported
  for each of G, X, Y.

### 5.3 Profile error, bias estimator and calibrated decisions (`src/stats_mg2.py`)

For field S, with run-level profiles m_i (i = 1…M), δ_i = m_i − φ*_S on |x| < 10, and run mean m̄:
- **Point estimate (reported):** the U-statistic Ê²_b,S = Σ_{i≠j}⟨δ_i, δ_j⟩ / (M(M−1)‖φ*_S‖²).
  It is unbiased for the squared relative bias of the survival-conditioned, aligned ensemble, and it
  can be negative. E_b,S = √max(Ê², 0) is reported.
- **Decision bounds:** by the triangle inequality, |‖μ − φ*‖ − ‖m̄ − φ*‖| ≤ ‖m̄ − μ‖, and ‖m̄ − μ‖ is
  bounded by k·√(tr Σ̂/M) with k = t_{0.975, M−1}. Relative bounds:
  [max(E_tot − k√(tr Σ̂/M)/‖φ*‖, 0), E_tot + k√(tr Σ̂/M)/‖φ*‖].
- **Three-way decision** against a tolerance ε:
  - **PASS** if the upper bound < ε;
  - **FAIL** if the lower bound > ε;
  - **INCONCLUSIVE** otherwise. Insufficient precision is never a pass.
- **Two-sample version** (V2, MG3-C2b): bounds on ‖μ_A − μ_B‖ with Welch-type variance and k from the
  smaller sample.

**Calibration** (synthetic ensembles with known bias; six run-distribution families: isotropic Gaussian,
rank-1 along the bias, rank-3, two-point ±, which is the reviewer's counterexample family, and
multivariate t₃; bias 0–2%, run noise 0.1–3%):

| Runs per cell | Worst coverage of [lower, upper] | Worst false-PASS rate (true bias ≥ 1%) | Worst false-FAIL rate (true bias ≤ 1%) |
|---|---|---|---|
| 4 | 0.664 | 0.322 | 0.333 |
| 6 | 0.802 | 0.194 | 0.189 |
| **12** | **0.945** | **0.034** | **0.035** |

At 4 or 6 runs no procedure can be calibrated against discrete run distributions. When every run
happens to be identical (probability (3/4)^M in the two-point family), the sample carries no
information about the spread. **All profile-bias decisions therefore use M = 12 runs per cell.**
The two-sample bound at 12 per side covers ≥ 0.950, with false-equivalence ≤ 0.026 at the 0.25% V2
margin. At 8 per side it covers ≥ 0.962, with false-equivalence ≤ 0.017 at the 10% MG3 margin
(`results/MG003b/calibration*.txt`). On the reviewer's exact counterexample (three runs at +4%, three
at −4%, true bias 2%), the bounds are [0, 4.6%]: INCONCLUSIVE, not a false PASS.

**Exponent α (MG2-C2):**
- **Cell selection** (pre-registered, decided once on the data): Ω ≥ 100 cells whose bias lower bound
  is > 0.
- **Estimate:** α̂ from a least-squares fit of log Ê²_b,X = c − 2α log Ω over the selected cells.
  It needs ≥ 3 cells, all with Ê² > 0.
- **CI:** run-level bootstrap (2000 resamples within each selected cell). Bootstrap U-statistics use
  only pairs of *distinct* original runs, since pairing a duplicated run with itself biases the
  bootstrap upward (found and fixed in MG003b). The percentile interval is widened about α̂ by the
  calibrated factor **w = 2.0**.
- **Calibration of the C2 decision itself** (`src/calibrate_alpha_decision.py`,
  `results/MG003b/calibration_alpha_decision.txt`): M = 12; six run families; run-noise/bias ratio at
  Ω = 100 of 0.1, 0.3 and 1; 300 replicates.
  - Under the correct scaling (true α = 1), worst P(FAIL) = 0.020.
  - Under noise-like scaling (true α = 0.5), worst P(PASS) = 0.003.
  - Noise-dominated cases come out mostly INCONCLUSIVE, as intended.
  - Per-scenario *conditional* interval coverage at w = 2.0 is ≥ 0.92 wherever the interval is usually
    defined (`calibration.txt`, Part B v2). The decision rule, not the interval, is what is pre-registered.
- If fewer than 3 cells are selected, or fewer than 80% of bootstrap fits are valid, the result is
  INCONCLUSIVE.

### 5.4 Deterministic confirmation of nucleation (MG3)

A detection (a birth from a homogeneous start) is a **confirmed nucleation** if the concentration field
at detection + 20 time units, used as the initial condition for the seed-free lattice ODE at the same b,
converges by t = 300 to a stationary state containing at least one soliton (φ_X min < −0.5|A*|,
max |∂φ/∂t| < 10⁻⁶).

Detections that relax to the homogeneous state are **excursions**. This is a committor-style test: it asks
whether noise carried the state across the deterministic separatrix.

## 6. First-passage model (MG3)

- **Initial state:** independent Poisson counts at the homogeneous steady state.
- **Burn-in:** t_b = 50, about 10 relaxation times of the slowest homogeneous mode (1/0.195 ≈ 5.1).
  Confirmed events with detection time < t_b are **startup events**. They are reported, the run is
  left-truncated, and it contributes no exposure.
- **Censoring:** confirmed events count if detected in [t_b, T − 20]. Runs without one are right-censored
  at T_c = T − 20.
- **Likelihood:** a constant hazard k per unit length. MLE k̂ = n / (L · Σ exposure), using only each
  run's first confirmed event. 95% CI: exact Poisson (Garwood) on n. Zero events give a one-sided 95%
  upper bound of 2.996 / (L · Σ exposure). This is a bound, **not** a demonstration that nucleation is
  impossible.
- **Stationary-hazard check:** a conditional binomial test of events in [t_b, 200) vs [200, T_c] given the
  exposures. If p < 0.01, the constant-hazard model is rejected for that cell, and the cell cannot enter
  the C3 fit.
- **Rare-event τ check (V5, inside MG3):** at the middle Stage-B cell, 16 runs at τ/2. The rate ratio's
  95% CI (conditional binomial) must include 1, and the point ratio must lie in [0.67, 1.5]. Otherwise all
  MG3 runs repeat at τ/2. At most 2 halvings, then MG3 is INCONCLUSIVE.

## 7. Random streams

Each trajectory's seed comes from
`numpy.random.SeedSequence(entropy=20261027, spawn_key=(gate, b_idx, Ω_idx, stage, replicate, segment)).generate_state(1)[0]`.
- gate: 0 = validations, 2 = MG2, 3 = MG3.
- b_idx and Ω_idx are indices into the gate config lists.
- stage: e.g. Stage A = 0, B = 1, C4 = 2, null = 3, τ-check = 4.
- replicate: globally unique within its cell.
- segment = 0 for an uninterrupted trajectory. Each continuation gets segment + 1 (numba's RNG state is
  not carried over), so streams never repeat.

All keys are written to the result files.

## 8. Compute budget (2 cores)

| Block | Cap (core-hours) |
|---|---|
| V1–V4 | 3 |
| MG2 (incl. C5, C6, nulls; 12 runs per cell) | 8 |
| MG3 (Stage A, B, extensions, C4, null, V5) | 10 (hard cap; see `GATES.md`) |

If a cap is reached, the block stops at its current step and reports INCONCLUSIVE with all the data
collected so far.

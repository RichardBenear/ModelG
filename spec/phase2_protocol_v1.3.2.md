# Phase 2 Protocol v1.3.2: Stochastic Token Rules on a Fixed Lattice (MG004d)

Status: **FROZEN 2026-10-28 (MG004d)**. This is v1.3.1 with **Amendment 2** applied
(`spec/phase2_amendment_2.md`, revision 2), approved by the PI after ChatGPT review. v1.3 and v1.3.1 are
preserved unchanged, and the validation records obtained under them stand: V2 registered FAIL (profile
equivalence unresolved), V4 INCONCLUSIVE, V1/V3/V6 PASS. Changes are marked **[A2-amend]**.

**[A2-amend] summary:**
- **Single high-Ω decision point Ω = 10⁸** for MG2-C1, MG2-C6, V4 and the V2 profile test.
- **Bounded allocation:** 12 eligible runs (cap 24) plus at most one 12-run extension. No further
  escalation of Ω or M; an INCONCLUSIVE stands.
- Tolerances unchanged, with no relaxation for G.
- **V2:** bound-based profile rule and signed variance rule (below). The Ω = 300 variance-equivalence
  pass is retained as evidence at that operating point only.
- **Time step and α:** bias(Ω, τ) = B(τ)/Ω + O(Ω⁻²). The time step changes B, and preserves the
  asymptotic exponent only if B(τ) ≠ 0. This does not guarantee the fitted exponent over the tested
  ladder.
- **New V7** (descriptive detector characterization at MG3 Ω values) with a **PI checkpoint before
  MG3**.

Changes from v1.2 (MG003c):
- The confidence claims are now explicitly **conditional on stated assumptions A1–A2** (§5.3).
- A1 is checked on every decision cell (diagnostics), on a pre-gate pilot of real runs (V6), and across
  each gate (pooled audit).
- The complete decision procedure (minimum eligible runs, survival top-up, diagnostics, one extension)
  is calibrated, and A1-violating families are reported as held-out cases.
- Two-sample tests use a Welch bound with deterministic normalization.
- The α-calibration record is corrected: w = 2.0 is calibrated for the C2 decision rule.

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
- **V2 — τ-halving equivalence on the MG2 statistics.** **[A2-amend]** b = 29.08, **Ω = 10⁸**, 12 eligible
  runs each at τ and τ/2 (cap 24, one extension). Profile rule (unsigned norm): PASS if the upper bound
  < 0.25%; FAIL if the lower bound > 0.25%; else INCONCLUSIVE. Variance rule (signed difference, 90% CI
  [lo, hi], margin m = 5%): PASS if −m < lo and hi < m; FAIL if lo > m or hi < −m; else INCONCLUSIVE.
  The v1.3.1 run at Ω = 300 (registered FAIL; profile equivalence unresolved; variance equivalence
  passed) is retained. **[A1-amend history]** b = 29.08 (non-gate; also V6's point), originally
  Ω = 300. The original
  location (b = 29) yielded 3 eligible of 24 and was aborted (`results/MG004/v2_attempt1_b29_aborted.txt`).
  A pass establishes time-step accuracy at this operating point only. This is an equivalence test with predefined margins, not a non-rejection test.
  - **Profiles:** for each of G, X, Y, the upper 95% bound on ‖μ_τ − μ_{τ/2}‖/‖φ*‖ (two-sample bound,
    §5.3) must be < 0.25%, a quarter of MG2-C1's tolerance.
  - **Fluctuations:** for each of G, X, Y, the 90% CI (Welch t, run-level values) of the difference in
    far-field Ω·Var must lie within ±5% of the τ value.
  - ~~PASS if all hold; FAIL if any point estimate lies outside its margin; INCONCLUSIVE otherwise.~~
    **[Erratum MG004g]** This v1.3/v1.3.1 sentence was carried over by mistake. Under v1.3.2 the
    Amendment 2 rules stated at the top of this V2 entry govern (profiles bound-based; variances signed),
    and the V2 run under v1.3.2 (`results/MG004/v2b.*`) applied them. If INCONCLUSIVE, one extension of
    12 + 12 runs is allowed, after which INCONCLUSIVE counts as FAIL.
- **V3 — clip pilot.** One run of T = 20 per Ω in each gate's Ω set, homogeneous start. This sets the
  per-Ω τ (§3.3). Only the clip counter is read from these runs.
- **V6 — A1 pilot on real runs** (MG003c). Two non-gate cells, b = 29.08 at Ω = 300 and Ω = 3000, 12
  eligible runs each, MG2 setup.
  - **Pass:** the per-cell A1 diagnostics pass in both cells, and the pooled 24-run studentised audit
    flags nothing.
  - Reported: the block-deviation distributions (QQ data), and the ratio of the run-level variance to
    block variance/K, which checks that blocks act as near-independent samples.
  - If V6 fails, no gate run proceeds and the protocol returns to the PI and reviewer.
- **V7 — detector characterization at MG3 Ω values [A2-amend; descriptive].** Homogeneous start, b = 29,
  every MG3 grid Ω (1–300), 4 runs, T = 100, per-Ω τ from V3. Reported:
  - the sub-threshold detection rate and detections per snapshot (tracking load);
  - the persistence-time distribution of tracks;
  - every persistent detection (≥ 20 time units) classified by §5.4 deterministic confirmation as
    **confirmed** (relaxes to a stationary soliton: genuine nucleation, i.e. signal), **excursion**
    (relaxes to homogeneous), or **unresolved**. At most 10 per run are confirmed; the rest are counted.

  No thresholds are set. **PI checkpoint:** the PI reviews V7 before any MG3 run.
- **V4 — alignment validation.** Synthetic ensembles: φ* plus a noise field taken from independent
  *homogeneous-state* runs at the same Ω and b, with random circular shifts. These go through the full
  detection/alignment/averaging pipeline (§5) at Ω ∈ {100, 1000, 10000} (v1.3/v1.3.1). **[A2-amend]** v1.3.2: Ω ∈ {10⁴, 10⁶, 10⁸}, decided at 10⁸.
  - 12 synthetic runs per Ω. Report the alignment-induced bias (point estimate and bounds) and its
    Ω-scaling.
  - **Pass:** the upper 95% bound of the alignment-induced bias at Ω = 10⁸ **[A2-amend; was 10⁴]** is < 0.3% for G, X, Y (a third
    of MG2-C1's tolerance). FAIL if the lower bound is > 0.3%; INCONCLUSIVE (blocking) otherwise.

## 5. Measurement definitions

### 5.1 Detection and tracking

- **Smoothing:** Gaussian with σ = 0.4 *physical units* (2 nodes at h = 0.2, 4 at h = 0.1). It is used
  for detection and alignment only.
- **Particles:** local minima of smoothed φ_X below −0.5·|A*_s|, where **A*_s is the minimum of the
  lattice-ODE reference φ_X (at that b and h) smoothed with exactly this smoothing operation**, computed
  in code. **[A1-amend]** This replaces −0.5·|A*| with the unsmoothed amplitude, which put the threshold
  at 81% of the smoothed depth. Values: −2.7685 (b = 29), −2.8962 (b = 29.08), −2.9908 (b = 29.16).
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
- Run-level profile m_i = the time average of run i's eligible aligned snapshots (t ≥ 10). Each run is
  also split into K = 6 consecutive blocks of 15 time units (block profiles, for A1 diagnostics).
  A run with any loss, replication or merge before t = 100, or fewer than 20 eligible snapshots in any
  block, is not eligible for profile statistics, but is counted in survival.
- **Minimum eligible runs:** each decision cell needs **12 eligible runs**. Runs are attempted in
  stream order until 12 are eligible or 24 have been attempted; fewer than 12 eligible gives an
  INCONCLUSIVE cell ("insufficient eligible runs"). The extension adds 12 more eligible runs under the
  same cap rule.
- **Far-field variance:** nodes at circular distance > 20 from every detected particle (all nodes in
  snapshots with none detected), over snapshots with t ≥ 10, about the run's own far-field mean.
  Reported for each of G, X, Y. **[A1-amend] C3 data rule:** it is measured on the first 12 *attempted*
  runs of each cell in stream order, whether eligible or not. Fewer than 12 runs with ≥ 100 such snapshots
  gives an INCONCLUSIVE C3 cell, and C3 cannot then PASS.

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

**Assumptions (the confidence claims below hold only under these):**
- **A1:** the eligible runs' aligned run-level profiles (and their K = 6 block profiles) are i.i.d.
  with approximately Gaussian fluctuations. In particular there are **no rare discrete run-level
  modes** beyond the tracked events (loss, replication, merge), which survival conditioning removes.
- **A2:** runs are independent. This holds by construction: distinct random streams, §7.

Physical basis for A1: every run starts from the same deterministic state. The run-level profile is a
time average over 90 time units, about 18 relaxation times of the slowest homogeneous mode. The
discrete changes this system can make (losing the particle, replicating, merging) are exactly the
tracked events. **Without tail assumptions, no fixed run count gives a distribution-free guarantee.**
The reviewer's example (run errors of 0.5% with p = 0.9, 6.5% with p = 0.1) falsely passes with
probability 0.9¹² = 28%. A1 is therefore checked rather than presumed:

- **Per-cell diagnostics** (`stats_mg2.a1_diagnostics`). Projections onto the run-mean bias direction
  and the first two principal components of the within-run block deviations are tested three ways:
  1. Shapiro–Wilk on the 72 within-run block deviations (flag if p < 10⁻³);
  2. block-level robust |z| > 5;
  3. run-level externally studentised leave-one-out residuals beyond the Bonferroni t_{M−2} critical
     value at overall level 1%. A lone deviating run among identical runs is flagged.

  Any flag makes the cell **INCONCLUSIVE ("A1 not supported")**. The run-level test replaced a MAD
  robust-z, whose null false-alarm rate at M = 12 was 21%.
- **V6 — pre-gate pilot on real runs** (below). A1 must hold on actual survival-conditioned run
  averages before any gate run.
- **Gate-level pooled audit.** After all of a gate's cells are run, the run-level studentised test is
  repeated over **all** runs of the gate pooled (each cell standardised by its own mean and SD;
  Bonferroni over every run and direction, overall level 1%). Pooling gives power against rare modes
  that no single cell has. Under A2, a 3% run-level mode **shared across all cells** would produce no
  minority run among the 168 MG2 runs with probability 0.97¹⁶⁸ ≈ 0.6%. That figure is the chance of the
  mode being entirely absent from the data, **not** the audit's overall failure probability, which also
  depends on the audit's detection power when a minority run is present and on whether the mode is
  shared across cells. (Clarified at freeze per review round 4.)
  Any flag downgrades every PASS in that gate to **INCONCLUSIVE pending investigation**, and the
  flagged runs are reported.

**Calibration of the complete single-cell procedure** (`src/calibrate_procedure.py`,
`results/MG003c/procedure_calibration.txt`). Synthetic runs with K = 6 blocks. The procedure includes
survival filtering (eligibility 1.0 or 0.7, at least 12 eligible from at most 24 attempted), the A1
diagnostics, and one 12-run extension. Biases 0–2% (including 0.9, 1.0 and 1.1%), run noise 0.1–1%,
600 replicates per cell:

| Family | False PASS (bias ≥ 1%) | False FAIL (bias ≤ 1%) |
|---|---|---|
| Gaussian isotropic (A1 holds) | 0.000 | 0.000 |
| Gaussian rank-3 (A1 holds) | 0.002 | 0.002 |
| t₃ blocks (stress; mostly flagged) | 0.000 | 0.000 |
| *Held out (A1 violated):* run mixture 0.9 / 0.1 | 0.312 | 0.003 |
| *Held out:* run mixture 0.97 / 0.03 | 0.698 | 0.000 |
| *Held out:* block-level jumps 0.1 | 0.002 | 0.000 |

- Power: P(PASS) ≥ 0.985 at true bias 0–0.5% with 0.1% run noise.
- Diagnostic flag rates (12 runs): Gaussian **0.015** (false alarm); t₃ blocks 0.785; run mixture 0.9
  0.384; run mixture 0.97 0.273; block jumps 1.000.
- **Reading this honestly:** single-cell diagnostics cannot catch a rare run-level mode when no minority
  run appears in the sample. That residual risk is what the gate-level pooled audit and V6 address. It
  is not eliminated, and it is stated as a limitation.

**Earlier fixed-sample calibration** (v1.2, `results/MG003b/`) is retained. At 4 / 6 / 12 runs the
worst coverage was 0.664 / 0.802 / 0.945, which is why decisions use 12 runs.

**Two-sample bound** (V2, MG3-C2b): `stats_mg2.welch_diff_bounds` (Welch–Satterthwaite df on trace
variances). The normalization is the **deterministic** lattice-ODE soliton norm ‖φ*_S‖, never a sample
norm. Calibrated with unequal noise (SD ratio 1/3, 1, 3) and unequal M (8/12, 12/8, 12/12) at the 0.25%
and 10% margins: worst coverage **0.990**, worst false-equivalence **0.001**.

**Exponent α (MG2-C2):**
- **Cell selection** (pre-registered, decided once on the data): Ω ≥ 100 cells whose bias lower bound
  is > 0.
- **Estimate:** α̂ from a least-squares fit of log Ê²_b,X = c − 2α log Ω over the selected cells.
  It needs ≥ 3 cells, all with Ê² > 0.
- **CI:** run-level bootstrap (2000 resamples within each selected cell). Bootstrap U-statistics use
  only pairs of *distinct* original runs, since pairing a duplicated run with itself biases the
  bootstrap upward (found and fixed in MG003b). The percentile interval is widened about α̂ by the
  calibrated factor **w = 2.0**.
- **Calibration of the C2 decision itself** (the pre-registered quantity; interval coverage is
  informational only. The corrected record is in `results/MG003b/calibration.txt`, Part B v2: undefined
  intervals are reported separately, never as zero coverage) (`src/calibrate_alpha_decision.py`,
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
| V1–V4, V6 | 4 |
| MG2 (incl. C5, C6, nulls; 12 runs per cell) | 8 |
| MG3 (Stage A, B, extensions, C4, null, V5) | 10 (hard cap; see `GATES.md`) |

If a cap is reached, the block stops at its current step and reports INCONCLUSIVE with all the data
collected so far.

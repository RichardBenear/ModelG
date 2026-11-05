# MG9b — How does the noise-activation barrier depend on b? Design proposal, DRAFT 1

*Stage MG034. Claude, 2026-11-05. Status: **DRAFT 1, not registered. No MG9b run has been made.** Requested by the PI
(thread message 02:32 UTC, "Let's proceed to MG9b now"), after MG9d closed. Naming: MG9b has been reserved since MG9 for
the stochastic barrier S(b). MG8 stays reserved for multiway/foliation.*

## 0. Starting point and disclosed prior knowledge

- **MG3 / MG3b (b = 29.00, h = 0.2, L = 100, τ = 6.667e-4).** Confirmed nucleations from a homogeneous Poisson start
  have a first-event rate k(Ω) that falls with the system size Ω. The rate is linear in Ω on a log scale:
  log k = α − S·Ω, with **S = 0.064 per unit Ω** (profile 95% interval 0.047–0.082, calibrated as approximate; lack of
  fit p = 0.89). The cell rates were 1.62e-5 (Ω = 130) to 1.13e-4 (Ω = 100) per unit length and time. τ against τ/2
  at Ω = 120 was EQUIVALENT within a factor of 2 (MG3b Part 1). Cost was about 280 core-s per run at τ on the PI's
  machine (9 workers).
- **MG9 (deterministic).** The single-particle window on h = 0.2 is [b_s, b_r] = [28.929, 29.396]. The local onset is
  SUBCRITICAL (L₃ > 0). The deterministic seed threshold falls with b: 0.949 at 29.00, 0.887 at 29.16, 0.763 at 29.50.
  MG9 called this a protocol-specific basin-access proxy, not a stochastic action.
- **MG9c / MG9d.** Extended patches near b_r creep slowly. MG9d found the late drift of such a patch depends on its
  position relative to the h = 0.2 lattice (MESH EFFECT SUPPORTED). So localized states on this lattice feel the grid,
  and MG9b is explicitly a **lattice-conditional** study. ChatGPT MG033 advised MG9b to include its own deterministic
  event-confirmation checks at each new b.
- **What the prior knowledge suggests (not a registered prediction).** If the falling seed threshold reflects a lower
  barrier, S should fall as b rises. Because this was known before the design, the primary question is stated two-sided
  (§1), and the expected direction is disclosed here rather than built into the decision rule.

## 1. Question and scope

**Q (primary): Does the activation slope S(b) change with b inside the particle window?** S(b) is MG3b's slope of log k
against Ω, measured with MG3's protocol at b = 29.20 and b = 29.35, and combined with MG3b's b = 29.00 data.

- The outcomes are two-sided: DECREASING, INCREASING, or NOT RESOLVED (§5).
- Descriptive only: k(b) at a common Ω, the prefactor α(b), and S(b) against the seed threshold.
- Out of scope: the noise law beyond MG3's tau-leap chemical Langevin model, the continuum limit, h = 0.1, and the
  asymptotic fate of nucleated particles beyond MG3's confirmation horizon.
- Interpretation is conditional on h = 0.2, τ, the tau-leap scheme and the MG3 event definition. S is an empirical
  slope, not a proven large-deviation action.

## 2. Design

**b values.** 29.20 and 29.35, both inside the window and at least 0.04 from b_r. b = 29.00 comes from MG3b, with the
same code, lattice, τ, event rule and estimator. 29.35 is the highest b below b_r with an existing MG9/MG9c deterministic
record (an 18-cell patch creeps there, so a nucleated particle should still be a single localized structure on the
confirmation horizon).

**Run protocol (unchanged from MG3/MG3b).**
- Homogeneous Poisson start, no seed, T = 500, τ = 6.667e-4.
- Events are confirmed nucleations (protocol §5.4), with burn-in 50, censoring at 480, and the first event per run.
- Smoothing, detection, tracker (A4) and crowding exclusion (A3) are as in MG3.
- The reference soliton for each b, `runner.reference(b, 0.2)`, is computed in the main process before the pool starts.
- Random streams come from `runner.seed_for` with a new gate index for MG9b, so they do not overlap MG3/MG3b.

**Part A — Ω placement scan (per b).** Rates shift with b, so MG3b's Ω values cannot be reused blindly.
- Descending grid Ω ∈ {250, 200, 170, 150, 140, 130, 120, 110, 100, 90, 80}, 4 runs per Ω, at τ.
- A V3 clip pilot (T = 20) and the A3 crowding screen run before each Ω cell. A crowded cell stops the descent, as in MG3.
- The scan stops after the first Ω where at least 3 of 4 runs have a confirmed event in [50, 480].
- **Top cell Ω_T** = the highest scanned Ω with at least 2 of 4 confirmed events. If none qualifies, the b is
  NOT PLACED and gets no Part B (reported; the trend uses the other b values).

**Part B — rate cells (per b).** Three cells: Ω_T, Ω_T − 15 and Ω_T − 30.
- The 30-unit lever arm equals MG3b's four cells spaced by 10, and costs one cell less.
- A cell that fails its clip pilot or crowding screen is replaced once by Ω_T + 15. If that also fails, the b has two
  cells and its S is labelled "two-cell".
- MG3b Part 2 batching: batches of 16 runs until at least 20 post-burn-in events at the end of a batch, cap 64 runs. A
  started batch is completed. Capped cells stay in the fit, flagged "below target precision".
- Part A runs are not pooled into Part B. Part A is for placement only and is reported descriptively.

**Time step.** All runs at τ. MG3b's EQUIVALENT result at b = 29.00 is carried as an assumption, not re-tested. All
rates and slopes are labelled "conditional on τ". An optional τ/2 arm is §8 decision 3.

## 3. Lattice-conditional checks (per b; deterministic; before Part B)

- **V-b1 reference stationarity.** Integrate the reference soliton deterministically to t = 2000. Record its maximum
  residual and centre drift. If it is not a single stationary particle (residual > 1e-6, or drift > λ/4), the b is
  flagged REFERENCE-UNSTABLE, and its events are still counted but reported with the flag.
- **V-b2 confirmation controls.** In Part A, each detection's §5.4 outcome (t_end = 300) is also re-run to t = 2000.
  The first 10 such pairs per b are compared. Agreement below 9/10 flags CONFIRMATION-SENSITIVE at that b.
- **V-b3 solver tolerance.** For the Part B cell Ω_T at each b, every confirmation is repeated at rtol 1e-10 and atol
  1e-12. Any change of class flags CONFIRMATION-TOLERANCE at that b. This answers MG9d's §6 note on the confirmation
  step.
- **Event shape (descriptive).** For each confirmed event: the number of persistent structures at t = 300, and whether
  it is one particle or a growing patch. At 29.35 a creeping patch is possible, so this is reported per b.

A flag does not stop the study. It is carried into the reading for that b (§5).

## 4. Estimation (MG3b methods, unchanged)

- Per cell: events, exposure, startup count, k̂ with its CI (Garwood), and risk sets at t = 50, 100 and 150.
- Per b: the Poisson GLM log k = α_b − S_b·Ω with an exposure offset. Profile-likelihood 95% CI for S_b, and deviance
  lack of fit (df = cells − 2; df = 0 for three cells, so lack of fit is reported for b = 29.00 only).
- **Trend model (primary).** All cells from the three b values in one Poisson GLM:
  log k = α_b − (S₀ + S₁·(b − 29.175))·Ω, with a free α_b per b. The quantity tested is **S₁ = dS/db**, with a
  profile-likelihood 95% CI.
- **Calibration (synthetic, before any run).** The whole procedure is simulated, including Part A placement, batching,
  caps and the trend fit. The cases are S₁ ∈ {−0.10, 0, +0.10} per unit Ω per unit b, with S₀ = 0.064 and MG3b-like
  prefactors, and 2000 replicates each. The S₁ CI is labelled "calibrated" if coverage ≥ 0.93 in all three cases,
  and "approximate" otherwise. The power to detect S₁ = −0.10 is reported.

## 5. Decision rules

**Primary (Q), ordered:**
1. **NOT RESOLVED (insufficient):** fewer than two of the three b values have an S_b estimate. That happens when a b
   is NOT PLACED or has fewer than two cells.
2. **DECREASING:** the 95% CI of S₁ lies entirely below 0.
3. **INCREASING:** the 95% CI of S₁ lies entirely above 0.
4. **NOT RESOLVED:** otherwise. This is not evidence of no dependence. The CI is reported, and the power from §4 is
   quoted beside it.

**Qualifiers (added to the label, not changing it):**
- **LATTICE-FLAGGED (b):** any V-b flag at that b.
- **CAPPED (b):** any capped cell at that b.
- **APPROXIMATE:** if the calibration label is "approximate".
- **NON-MONOTONE (descriptive):** S_b at the middle b falls outside the range of the other two by more than its CI
  half-width.

**Wording.** DECREASING reads "the activation slope falls with b over 29.00–29.35 on the h = 0.2 lattice, at τ". It
supports, but does not prove, a lower barrier nearer b_r. No continuum or Landau-exponent claim is made.

## 6. Budget, caps and records

| Item | Run cap | Projected core-h (PC, 9 workers) |
|---|---|---|
| Calibration (synthetic) | — | ≤ 1 |
| V-b1 (2 deterministic runs to t = 2000) | 2 | < 0.1 |
| Part A, per b (≤ 11 Ω × 4 runs; mostly full-length runs) | 44 | ≤ 5.5 |
| Part B, per b (3 cells × ≤ 64 runs; expected about 48 each) | 192 | about 11 (cap 15) |
| V-b2/V-b3 extra confirmations | — | < 0.5 |
| **Total (two b values)** | **474 runs** | **about 34 projected; cutoff 50** |

- About 4–5 hours of wall time on the PI's PC at 9 workers. 1D runs are light on memory; 9 workers ran MG3b.
- **Cutoff enforcement.** The ledger is persisted per run. The cutoff is checked before each batch. The deadline is
  cooperative, not a strict wall-time guarantee (ChatGPT MG033). On reaching the cutoff, work stops and the PI decides.
  Nothing is replaced, and no cell is extended beyond its cap.
- **Order.** Calibration; V-b1 for both b; Part A at 29.35, then 29.20; Part B at 29.35, then 29.20. 29.35 goes first
  because it is the newer regime and the more likely to need a PI decision.
- **Records committed.** Per-run records (event times, exposures, confirmation outcomes), Part A placement records,
  V-b check records, the ledger, calibration results and the evaluation, with an independent recompute. No field
  snapshots are committed (repo-size preference), except the V-b1 final reference fields (small).

## 7. Relevance

- A DECREASING S supports a barrier that falls toward b_r, consistent with the falling seed threshold and the
  subcritical onset. It is a measured link between the deterministic window and stochastic particle creation.
- NOT RESOLVED leaves MG3b's b = 29.00 slope as the only measured barrier. A larger b range (closer to b_s and b_r)
  or more cells would then be a separate study.
- The asymptotic pinning question from MG9d stays optional and separate.

## 8. Decisions for the PI

**PI decision (2026-11-05 02:39 UTC, thread message): all five recommendations adopted as written.**

1. **b values 29.20 and 29.35**, with MG3b's 29.00 reused. Recommended. (An alternative is 29.10 and 29.30, further
   from b_r but with a smaller expected effect.)
2. **Three cells per b, spaced 15** (lever arm 30), instead of MG3b's four spaced 10. Recommended. It saves about 7
   core-h with a similar slope precision.
3. **No τ/2 arm.** Recommended; MG3b's EQUIVALENT is carried as an assumption with the "conditional on τ" label. The
   alternative is one τ/2 arm of 64 runs at Ω_T for b = 29.35 (about 11 core-h more).
4. **Budget:** about 34 core-h projected, cutoff 50, on the PI's PC. Recommended.
5. **Review:** ChatGPT reviews this draft before registration, then the pre-run implementation before any run.
   Recommended.

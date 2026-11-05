# MG9b — How does the noise-activation barrier depend on b? Design proposal, DRAFT 2

*Stage MG034. Claude, 2026-11-05. Status: **DRAFT 2, not registered. No MG9b run has been made.** Draft 1
(`spec/MG9b_proposal_v1.md`, superseded) was reviewed by ChatGPT in MG034 (`review/MG034_chatgpt.md`); every item is
accepted (`review/MG034_response.md`). The PI's five §8 decisions (02:39 UTC) are unchanged. Naming: MG9b has been
reserved since MG9 for the stochastic barrier S(b). MG8 stays reserved for multiway/foliation.*

**Changes from draft 1, in brief.** Outcomes renamed to linear-trend labels, with per-b slopes and a linear-model check
(§1, §5). Lack-of-fit df corrected: three cells give 1 df (§4). Calibration fully specified, including the historical
data and failure branches; failing calibration makes a directional result exploratory (§4.3, §5). V-b controls given
exact definitions, a late confirmation window and a fixed sample; a bad reference stops that b (§3). Placement and
replacement made deterministic with unique candidates (§2). Pilots, screens and controls added to caps and budget;
partial analysis at a budget stop defined (§6). τ and record limitations stated directly (§2, §6).

## 0. Starting point and disclosed prior knowledge

- **MG3 / MG3b (b = 29.00, h = 0.2, L = 100, τ = 6.667e-4).** Confirmed nucleations from a homogeneous Poisson start
  have a first-event rate k(Ω) that falls with Ω: log k = α − S·Ω, with **S = 0.064 per unit Ω** (profile 95%
  interval 0.047–0.082, labelled approximate; deviance 0.23 on 2 df, p = 0.89). Cells (runs / post-burn-in events /
  startup events): Ω = 130 (48/25/0), 120 (64/44/9), 110 (48/28/19), 100 (32/21/11). The Ω = 120 cell was MG3b Part 1's
  fixed 64-run arm; the others used event-count stopping. τ against τ/2 at Ω = 120 was EQUIVALENT within a factor of 2.
- **MG9 (deterministic).** Single-particle window on h = 0.2: [b_s, b_r] = [28.929, 29.396]. The local onset is
  SUBCRITICAL. The deterministic seed threshold falls with b: 0.949 at 29.00, 0.887 at 29.16, 0.763 at 29.50. MG9
  called this a protocol-specific basin-access proxy, not a stochastic action.
- **MG9c / MG9d.** Extended patches near b_r creep slowly, and MG9d found the late drift depends on position relative
  to the h = 0.2 lattice (MESH EFFECT SUPPORTED). MG9b is therefore explicitly **lattice-conditional**.
- **Disclosed expectation (not a registered prediction).** If the falling seed threshold reflects a lower barrier, S
  should fall as b rises. The primary question is two-sided.

## 1. Question and scope

**Q (primary): Under a linear model in b, does the activation slope S change with b inside the particle window?**
S_b is MG3b's slope of log k against Ω, measured with MG3's protocol at b = 29.20 and 29.35 and combined with MG3b's
b = 29.00 data. The tested quantity is the linear trend S₁ in S(b) = S₀ + S₁·(b − 29.175).

- Outcomes: **NEGATIVE LINEAR TREND, POSITIVE LINEAR TREND, or NOT RESOLVED** (§5). A linear-trend label is evidence
  for a trend under the linear model. It is not proof that S falls (or rises) monotonically across the three b values,
  and it is not an established derivative of an unrestricted S(b).
- Always reported beside the label: the **unrestricted per-b slopes S_b** with their intervals, and the linear-model
  check (§4.2).
- Descriptive only: k(b) at a common Ω, the prefactor α_b, S_b against the MG9 seed threshold, landmark sensitivities.
- **What is measured.** The event is MG3's finite-horizon field-level confirmed event (protocol §5.4). The slope is
  therefore for **protocol-confirmed structures**. It is read as a single-particle creation barrier at a given b only if
  that b's V-b controls pass (§3).
- Out of scope: noise laws other than MG3's tau-leap chemical Langevin model, the continuum limit, h = 0.1, τ
  independence at the new b values, and the fate of nucleated structures beyond the confirmation horizon. S is an
  empirical slope, not a proven large-deviation action, and no physical energy-barrier claim follows.

## 2. Design

**b values.** 29.20 and 29.35, both inside the window and at least 0.04 below b_r. b = 29.00 is MG3b's data, used as
recorded. 29.35 is the highest b below b_r with an MG9/MG9c deterministic record.

**Run protocol (unchanged from MG3/MG3b, with b_det fixed).**
- Homogeneous Poisson start, no seed, T = 500, τ = 6.667e-4, with MG3's stop 50 time units after the first confirmed
  event. Events are confirmed nucleations (§5.4): burn-in 50, censoring at 480, first event per run. Smoothing,
  detection, tracker (A4) and crowding exclusion (A3) as in MG3.
- **`b_det = b`** in every `mg3.new_state` call (its default is 29.0). The driver asserts that each record's b_det
  equals its b. The reference `runner.reference(b, 0.2)` is computed in the main process before the pool starts.
- Random streams: `runner.seed_for` with a new gate index for MG9b, with distinct stage indices for Part A runs, clip
  pilots, Part B screens and Part B production, so no stream overlaps MG3/MG3b or another MG9b item.

**Time step.** All runs at τ, with no τ/2 arm (PI decision 3). MG3b's τ-equivalence at b = 29.00 does **not**
validate time-step independence at 29.20 or 29.35. MG9b is a study of the fixed-τ algorithm, and every rate and slope
is labelled "at τ = 6.667e-4".

### 2.1 Part A — Ω placement scan (per b)

1. Descending grid Ω ∈ {250, 200, 170, 150, 140, 130, 120, 110, 100, 90, 80}.
2. At each Ω, first a **clip pilot at τ**: one homogeneous run to T = 20 with MG3's clip measure. A clip fraction
   above 1e-4 does **not** halve τ (unlike MG3). It makes that Ω unavailable and stops the descent.
3. Then 4 runs. The A3 crowding screen uses these 4 runs (mean detections per snapshot for t < 50, threshold 2). A
   crowded Ω is excluded and stops the descent.
4. The scan stops after the first Ω where at least 3 of 4 runs have a confirmed event in [50, 480], or at a stop from
   step 2 or 3, or at the end of the grid.
5. **Top cell Ω_T** = the highest scanned, non-excluded Ω with at least 2 of 4 confirmed events. If none qualifies,
   the b is **NOT PLACED**: it gets no Part B and no S_b.
6. Part A runs are never pooled into Part B. They are reported descriptively.

### 2.2 Part B — rate cells (per b)

- **Cells.** Ω_T always (it passed its screens in Part A). Two more cells are taken, in order, from the unique
  candidate list **Ω_T − 15, Ω_T − 30, Ω_T + 15, Ω_T + 30**. A candidate is accepted if it passes a pre-production
  screen. Candidates are tried in order until two are accepted or the list is used up (at most 4 attempts). Every
  candidate is a distinct Ω, so no two cells can coincide.
- **Pre-production screen** (per candidate): a clip pilot at τ (as in §2.1) and a separate 4-run crowding screen
  integrated to t = 50, on their own seeds. These runs are ledgered and reported, never pooled.
- **Production.** MG3b Part 2 batching: batches of 16 runs until at least 20 post-burn-in events at the end of a
  batch, cap 64 runs. A started batch is completed and all its runs count. A capped cell stays in the fit, flagged
  "below target precision". Ω_T's production runs are new runs (Part A is not pooled).
- **Order of cells within a b:** Ω_T first, then accepted candidates in list order.
- **When an S_b exists.** At least 2 distinct production cells with positive exposure, and a finite maximum-likelihood
  S_b. With exactly 2 cells, S_b is labelled "two-cell" (lack of fit 0 df).
- **Screening is not replacement of data.** Screening decides which candidate cells enter production. Once a cell's
  production has started, its data are kept and it is never replaced or extended beyond its cap (§6).

## 3. Lattice-conditional controls (per b)

Constants: λ* = `runner.reference(b, 0.2)["lam"]` (twice the mean zero spacing of the reference tail) and A* its
unsmoothed depth, as used by `confirm_field`. "Deep minima" and "linking" are exactly as in `confirm.py`: local minima
of φ_X below −0.5|A*|, linked from check to check with each step's circular displacement ≤ λ*.

**V-b1 reference validity (deterministic, before Part A at that b).**
- Integrate the seed-free lattice ODE from the reference state `runner.reference(b, 0.2)["phi"]` to t = 2000 at the
  default tolerances (rtol 1e-8, atol 1e-10), with outputs every 10 units on [1900, 2000].
- **Pass** requires all of: (a) the integration completes; (b) exactly one linked deep-minimum chain over the 11
  checks 1900, 1910, …, 2000, and exactly one deep minimum at each check; (c) centre drift |c(2000) − c(0)| ≤ λ*/4,
  where c(t) is the vertex of the parabola through the deepest φ_X minimum and its two neighbours (circular distance);
  (d) the depth min φ_X at t = 2000 within 5% of A*.
- Recorded descriptively: max |RHS| at 0 and 2000, the centre track, and the width.
- **A failure stops that b.** Part A and Part B are not run there, and the b is reported REFERENCE-INVALID, with no
  S_b. The reference and the event rule are never replaced after outcomes are seen.

**V-b2 late confirmation control (from production confirmations).**
- A new function confirms the same field over a **late window 1900, 1910, …, 2000** (11 checks, same linking rule),
  with excursion (max |φ| < 1e-3) and unresolved judged at t = 2000. This is needed because `confirm_field` always
  evaluates at its fixed checks 200–300, whatever t_end is (MG034 §3).
- **Sample:** the first 10 **unique confirmation fields** at that b, in the deterministic order Part A (Ω descending,
  then run index, then confirmation time), then Part B (cell order, run index, time). A unique field is one
  (run, snapshot) pair; repeated detections that share one cached field count once. Only these 10 fields get the late
  continuation.
- **Comparison:** binary, persistent against not persistent (excursion and unresolved together), between the standard
  (200–300) and late (1900–2000) outcomes. An incomplete continuation counts as a disagreement.
- **Pass:** at least 9 of 10 agree. Below 9: CONFIRMATION-SENSITIVE. Fewer than 10 unique fields at that b:
  CONTROL-INSUFFICIENT. Neither is a pass. Three-class agreement is reported descriptively.

**V-b3 solver tolerance control.** `confirm_field` gains optional `rtol` and `atol` arguments passed to
`lattice_ode_1d` (defaults unchanged, so existing calls are bit-for-bit identical; a test checks this). Every
confirmation in the Ω_T production cell is repeated at rtol 1e-10, atol 1e-12. Paired classes, n_persistent and solver
completion are saved. **Pass:** no change between persistent and not persistent, and every tight run completes.
Otherwise CONFIRMATION-TOLERANCE.

**Event shape (descriptive).** For each confirmed event: the number of persistent structures at t = 300 and whether it
is one deep minimum or a multi-minimum patch.

**What a V-b2 or V-b3 flag does.** It does not stop the study, and the b stays in the trend fit. The result for that
b is worded "a protocol-confirmed structure rate; the particle-barrier interpretation is not validated at this b"
(qualifier LATTICE-FLAGGED). Numerical incompletion is never counted as a class or as a passed control.

## 4. Estimation

### 4.1 Per cell and per b (MG3b methods)
- Per cell: runs, events, exposure, startup count, k̂ with its Garwood CI, risk sets at t = 50, 100 and 150, and
  **landmark sensitivities**: k̂ with left truncation at 100 and at 150 (events, exposure, CI). These are descriptive;
  risk sets alone do not test a constant hazard, and no stationarity claim is made.
- Per b: Poisson GLM log k = α_b − S_b·Ω with an exposure offset, profile-likelihood 95% CI for S_b, and deviance lack
  of fit on **cells − 2 df**: 2 df at b = 29.00 (four cells), 1 df for a three-cell b, 0 for a two-cell b (not
  reported). The 1-df deviance has limited power and only an asymptotic reference distribution; this is stated with it.

### 4.2 Trend model (primary) and linear-model check
- All cells from the available b values in one Poisson GLM: log k = α_b − (S₀ + S₁·(b − 29.175))·Ω, with a free α_b
  per b. S₁ gets a profile-likelihood 95% CI.
- **Linear-model check (three b values only).** Likelihood-ratio test of the trend model against the model with a
  separate S_b per b (1 extra parameter, 1 df). With two b values the two models are identical and the check is
  unavailable (stated).
- **Range.** The fitted b range and the b values used are stated in every result.

### 4.3 Calibration (synthetic, in the cloud, before any production run)

**Data-generating process.**
- Run-level first-passage model: an independent startup event (event before t = 50) with probability p_s; otherwise a
  constant hazard k_b(Ω)·L after t = 50 (L = 100), censored at 480. A run's exposure is from 50 to its event or 480.
- Rates: k_b(Ω) = 3.41e-5 · exp(−S_b·(Ω − Ω0_b)) per unit length and time, with S_b = 0.064 + S₁·(b − 29.00), so b =
  29.00 always has MG3b's S.
- **Historical b = 29.00 regenerated under its original rule** (not held fixed): cells 130, 120, 110, 100 with Ω0 =
  120; Ω = 120 a fixed 64 runs; the others batches of 16 to ≥ 20 events, cap 64; p_s per cell 0.00, 0.14, 0.40, 0.34
  (the observed fractions).
- **New b values simulated through the full §2 algorithm:** Part A placement (descending grid, 4 runs, the 3-of-4 stop,
  the 2-of-4 top rule), Part B candidates in order, batching and caps.
- Screens (clip and crowding) cannot be generated from first-passage times. Each screen fails independently with
  probability π_f. Otherwise calibration is conditional on screens passing and on confirmed events following the
  first-passage model.

**Scenarios (2000 replicates each, fixed seed).**

| Scenario | Ω0 at 29.20 / 29.35 | p_s (new b) | π_f | S₁ values |
|---|---|---|---|---|
| O (ordinary) | 140 / 160 | 0.25 | 0 | −0.10, −0.05, 0, +0.10 |
| P (startup stress) | 140 / 160 | 0.50 | 0 | −0.10, 0 |
| Q (placement and screen stress) | 140 / 75 | 0.25 | 0.2 | −0.10, 0 |
| N (non-monotone, O settings) | 140 / 160 | 0.25 | 0 | S_b = 0.064, 0.080, 0.020 |

Q is chosen so that 29.35 sits near the grid floor (often NOT PLACED or few cells) and replacements and caps occur.

**Reported, per case and by available-b branch (three b, two b, fewer) where a branch has ≥ 100 replicates,** each
with its Monte Carlo SE (about 0.006 at 0.93 with 2000 replicates): the fraction with no estimable trend; coverage of
the true S₁; the rate of a directional call on each side; power for S₁ = −0.10 and −0.05 (scenario O); the
linear-model check's rejection rate (all cases, including N); and the distribution of cells per b, caps and runs used.

**Handling of failures.** A replicate where the design yields fewer than two b values with S_b counts as "no estimable
trend": in power and false-positive rates it counts as no directional call. Coverage is computed over replicates with
an estimable trend, and a fit that fails numerically there (non-finite estimate, or a profile bound that cannot be
found) counts as **not covering**. No replicate is dropped from any denominator.

**Calibration label.** "Calibrated" requires, in all eight linear cases, coverage ≥ 0.93 and a directional-call rate
on each side ≤ 0.035 when S₁ = 0. This is a check on a finite grid, not a guarantee of 95% coverage. Otherwise the
label is "not calibrated" and §5's exploratory rule applies.

## 5. Decision rules

**Primary (Q), applied in order:**
1. **NOT RESOLVED (insufficient):** fewer than two b values have an S_b (a b that is REFERENCE-INVALID, NOT PLACED, or
   with fewer than two production cells has none).
2. **NEGATIVE LINEAR TREND:** the 95% CI of S₁ lies entirely below 0.
3. **POSITIVE LINEAR TREND:** the 95% CI of S₁ lies entirely above 0.
4. **NOT RESOLVED:** otherwise. This means the trend is unresolved, not that S is independent of b; the CI and the
   §4.3 power are reported.

**If calibration is "not calibrated":** a directional label from rule 2 or 3 is reported as **EXPLORATORY** (for
example "NEGATIVE LINEAR TREND (EXPLORATORY)"), and the confirmatory primary outcome is NOT RESOLVED.

**Qualifiers (added to the label, not changing it):**
- **LINEAR MODEL MISFIT:** the §4.2 check has p < 0.05. The headline then leads with the per-b slopes; a trend label is
  reported only as a summary of the linear fit, and no monotone or derivative statement is made.
- **NON-MONOTONE (descriptive):** the middle S_b lies outside the range of the other two by more than its CI half-width.
- **LATTICE-FLAGGED (b):** a V-b2 or V-b3 flag (CONFIRMATION-SENSITIVE, CONTROL-INSUFFICIENT or CONFIRMATION-TOLERANCE)
  at that b.
- **CAPPED (b)** / **BUDGET-SHORT (b):** a capped cell, or a cell stopped by the budget, at that b.
- **TWO-B:** only two b values are available; the range is stated and the linear-model check is unavailable.

**Wording.** NEGATIVE LINEAR TREND reads "under a linear model, the activation slope of protocol-confirmed structures
decreases with b over [fitted range] on the h = 0.2 lattice, at τ = 6.667e-4". It is consistent with, but does not
prove, a lower barrier nearer b_r. No continuum, Landau-exponent or physical energy-barrier claim is made.

## 6. Budget, caps and records

| Item | Run cap | Projected core-h (PC, 9 workers) |
|---|---|---|
| Calibration (synthetic, cloud; separate cloud ledger) | — | ≤ 4 (cloud) |
| V-b1 (one deterministic run to t = 2000 per b) | 2 | < 0.1 |
| Part A clip pilots (T = 20) | 11 per b | < 0.2 |
| Part A runs (≤ 11 Ω × 4) | 44 per b | ≤ 5.5 per b |
| Part B screens (≤ 4 candidates × [1 clip pilot + 4 runs to t = 50]) | 20 per b | < 0.3 per b |
| Part B production (3 cells × ≤ 64; expected about 48 each) | 192 per b | about 11 per b |
| V-b2 late continuations (10 to t = 2000) and V-b3 tight confirmations (Ω_T cell) | 10 + Ω_T cell's confirmations, per b | < 0.5 per b |
| **Total on the PC (two b values)** | **534 stochastic runs + 2 reference runs + controls** | **about 35 projected; cutoff 50** |

- About 4–5 hours of wall time on the PI's PC at 9 workers. All PC work counts toward the one 50 core-h cutoff.
- **Cutoff enforcement.** The ledger is persisted per run, and the cutoff is checked before each batch, pilot, screen
  and control. The deadline is cooperative, not a strict wall-time guarantee. Runs in flight finish and count.
- **Partial analysis at a budget stop (frozen now).** All completed production runs are analysed under §4–§5 as they
  stand. A cell stopped below its event target stays in the fit flagged BUDGET-SHORT (like a capped cell). A b without
  two production cells has no S_b. Part A records and unfinished controls are reported as they stand, and a missing
  V-b2/V-b3 control counts as CONTROL-INSUFFICIENT. The PI then decides on any continuation, which would be a separate,
  disclosed decision. Nothing completed is discarded.
- **Order.** Calibration; V-b1 for both b; at 29.35 Part A, Part B screens, Part B production (V-b3 with its Ω_T cell);
  then the same at 29.20. V-b2 runs as soon as a b has 10 unique confirmation fields, at the latest after its Part B.
- **Records committed.** Code commit and platform; seeds (gate, stage and replicate indices); per-run records (event
  times, exposure, startup flag, stop reason); every confirmation (class, n_persistent, near_detection, solver
  completion, max step); clip pilots and screens; V-b1 records and its final field; the V-b2 sample's ten
  confirmation fields (small) with paired outcomes; V-b3 paired outcomes; the ledger; calibration results; the
  evaluation with an independent recompute.
- **Limitation (stated).** These records let the statistics be recomputed exactly. They do not allow independent
  field-level checking of production confirmations without replaying the simulations, because production fields are
  not archived (repo-size preference). The seeds and code make such a replay possible.

## 7. Relevance

- A NEGATIVE LINEAR TREND supports a barrier that falls toward b_r, consistent with the falling seed threshold and the
  subcritical onset: a measured link between the deterministic window and stochastic structure creation, on this
  lattice and τ.
- NOT RESOLVED means the trend is unresolved. The new per-b slopes S(29.20) and S(29.35) are still reported with their
  intervals and flags. A wider b range or more cells would be a separate study.
- The asymptotic pinning question from MG9d stays optional and separate.

## 8. Decisions for the PI

**PI decision (2026-11-05 02:39 UTC, thread message): all five recommendations adopted as written.** Draft 2 keeps
them; MG034 §5 endorsed all five.

1. **b values 29.20 and 29.35**, with MG3b's 29.00 reused.
2. **Three cells per b**, from Ω_T and the candidate list Ω_T − 15, Ω_T − 30, Ω_T + 15, Ω_T + 30 (draft 2 makes the
   replacement order deterministic; the spacing is unchanged).
3. **No τ/2 arm.** Results are labelled "at τ", and MG3b's equivalence is not carried over as evidence.
4. **Budget:** about 35 core-h projected on the PI's PC, cutoff 50, now including pilots, screens and controls. The
   synthetic calibration runs in the cloud.
5. **Review:** ChatGPT reviews this draft, then the calibration results and the implementation, before any production
   run.

No new PI decision is needed for draft 2. Its one judgment call is the exploratory rule in §5 (an uncalibrated
directional result does not count as a confirmatory trend), chosen because MG034 asked for that status to be frozen.

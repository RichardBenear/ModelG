# MG3b — Follow-up to MG3: nucleation rates, time step and uncertainty (DRAFT 2 for review)

Status: **DRAFT 2** (Claude, 2026-10-29). Revised after ChatGPT's review of draft 1 (`review/MG006b_chatgpt.md`;
response `review/MG006b_response.md`). Draft 1 is preserved as `spec/MG3b_followup_proposal_v1.md`. Nothing here is
registered: it needs PI approval and entry in `GATES.md` before any lattice run. MG3's registered FAIL, and the accepted
A1 residual risk, are unchanged by anything here. This is a separate study, designed after MG3's data were seen.

## Summary of the design
1. **Part 3 first.** Calibrate every interval and test against the exact MG3b sampling and stopping rules, by synthetic
   simulation only, and freeze the methods and the selection rule before any lattice run.
2. **Part 1** decides the production time step from a **direct comparison of post-burn-in nucleation rates** at
   Ω = 120 (τ against τ/2), using an **equivalence criterion** with three outcomes: EQUIVALENT, DIFFERENT, INCONCLUSIVE.
   The startup fractions are reported descriptively and do not decide anything.
3. **Part 2** measures rates at Ω ∈ {130, 120, 110, 100} with ≥ 20 events per cell where the cap allows. Every cell,
   capped or not, enters one joint censored-likelihood fit for the activated slope. All rates are stated as conditional
   on the production τ.
4. **Part 4 (full-trajectory confirmation) is removed** from MG3b and deferred until detection-level attribution exists.
5. **Budget** is enforced by registered run caps, plus a core-hour ceiling checked before every batch.

Motivation (MG3 record and the MG006 review):
- Nucleation is confirmed and its rate rises steeply as Ω falls, but the rates rest on 5–11 events per cell.
- Time-step independence is unresolved. MG3's V5 compared post-burn-in rates at Ω = 120 with a CI of [0.37, 2.2], and
  startup events were 5/16 at τ against 0/16 at τ/2.
- Stationarity was checked only weakly (no late exposure at Ω = 100).
- The interval procedures were not calibrated for stopped, left-truncated first-passage data.

Setting throughout: b = 29, h = 0.2, L = 100, homogeneous Poisson start; τ = 6.667e-4 as in MG3. Detection, tracking
(with the A4 fix) and §5.4 field-level confirmation are as in MG3, including the stop 50 time units after the first
confirmed event, burn-in t_b = 50 and censoring at t_c = 480 (by detection time). Fresh streams: gate id 4, with a
distinct stage index per arm (Part 1 τ: 0, τ/2: 1, τ/4: 2; Part 2: 3). Arms use independent streams, so the analysis
treats them as independent samples.

Every confirmation field (persistent, excursion and unresolved) is saved with each run record, along with the
per-snapshot detection counts. This costs about 6 kB per confirmation. It keeps a later detector study possible without
re-running.

## Part 3 — Calibration (runs first; synthetic only; no lattice simulation)

### Data-generating model
Each synthetic run reproduces MG3b's bookkeeping exactly: a startup event (a confirmed first event detected at t < 50)
with probability **p_s**, which ends the run at t_s + 50 with no exposure; otherwise a post-burn-in first event from a
hazard h(t) on t ≥ 50, censored at 480, with the run ending 50 after the event or at 500. **p_s is varied independently
of h(t)**, so the calibration does not assume any link between startup frequency and the later hazard. Batching,
event-count stopping, caps, and the reuse of Part 1's Ω = 120 arm in Part 2 are simulated exactly as registered below.

### Configurations
- p_s ∈ {0, 0.2, 0.4, 0.6}.
- Constant post-burn-in hazard per run (k·L) ∈ {5e-4, 1.5e-3, 5e-3, 1.5e-2} per time unit (Ω ≈ 150 to 100 in MG3).
- Non-constant alternatives on t ≥ 50: Weibull shape 0.7 and 1.5; a step ×2 at t = 200; a step ×0.5 at t = 200; an
  early excess (×2 on [50, 100)). Each is scaled to the same mean hazard over [50, 480] as its constant counterpart.
- Part 1: rate ratio k(τ/2)/k(τ) ∈ {0.5, 0.71, 0.8, 1, 1.25, 1.41, 2}, crossed with p_s pairs (0.35, 0.05) and
  (0.35, 0.35), at the Ω = 120 rate (k·L = 4.6e-3).
- Part 2 slope: four cells with S ∈ {0.05, 0.073, 0.10} per unit Ω and k(120) = 4.6e-5, with p_s per cell
  (0.15, 0.35, 0.33, 0.31) and the alternative (0, 0.6, 0.6, 0.6). One extra "capped" configuration sets Ω = 130's rate
  to one third of the activated value, so the 64-run cap binds there in most replicates.

### Candidates (frozen now)
- Per-cell rate CI: (a) Garwood on n with exposure treated as fixed; (b) exponential profile-likelihood; (c) parametric
  bootstrap percentile, B = 2000.
- Rate ratio CI (Part 1): (a) conditional binomial, Clopper–Pearson (MG3's `rate_ratio`); (b) conditional binomial,
  mid-p; (c) Wald on log ratio, with 0.5 added to each count.
- Activated slope CI (Part 2): (a) Poisson GLM, log k = α − S·Ω, events with exposure offset, Wald; (b) the same GLM,
  profile likelihood; (c) parametric bootstrap from the fitted GLM, percentile, B = 2000, re-running the full batching
  and stopping rule in each bootstrap replicate. MG3's weighted least squares on log k̂ is reported for continuity only.
  It is not a candidate, because it cannot take a cell with zero events.
- Stationary-hazard check: conditional binomial on the events before and after a split point, given exposures. The
  split is (a) fixed at t = 200 (as MG3) or (b) the cell's exposure median.

### Monte Carlo precision and the coverage requirement
- R = 4000 replicates per configuration. The Monte Carlo SE of an estimated coverage near 0.95 is about 0.0034.
- A CI method is **eligible** if, in **every constant-hazard configuration**, its estimated two-sided coverage is
  ≥ 0.94 and each one-sided miss rate is ≤ 0.035. For the 90% intervals used in Part 1's equivalence decision, the
  requirement is coverage ≥ 0.89, with each one-sided miss ≤ 0.06.
- **Selection rule (deterministic):** among the eligible methods, choose the one with the smallest mean log-width at the
  central configuration (p_s = 0.35, k·L = 4.6e-3, S = 0.073). Widths within 1% of each other count as a tie, and a
  tie goes to the earlier letter above (the simpler method).
- **If no candidate is eligible:** use the one with the highest minimum coverage across configurations. Report its
  intervals as "approximate (calibrated minimum coverage X)", and no claim may be stated at 95% on that interval.
  Calibration never adds runs or changes a cap, except through the one registered rule in Part 1.
- **Stationary check:** a split is eligible if its size is ≤ 0.012 at α = 0.01 in every constant-hazard configuration.
  Among the eligible splits, choose the one with the higher mean power against the ×2 and ×0.5 steps, with ties going to
  (a). Non-constant configurations also report the bias of k̂ relative to the exposure-weighted mean hazard. This shows
  what "the rate" means when the hazard is not constant.
- **Part 1 operating characteristics:** the probabilities of EQUIVALENT, DIFFERENT and INCONCLUSIVE at each ratio. The
  size requirement is P(EQUIVALENT) ≤ 0.05 at ratios 0.5 and 2.
- Outputs: `results/MG3b/calibration/` (script, seeds and JSON), committed. It is run and reviewed **before** any lattice
  run.

A preliminary check (Claude, scratch, not a registered calibration) used 64 runs per arm, Clopper–Pearson, p_s 0.35
against 0.05, and 3000 replicates. It gave P(EQUIVALENT) = 0.84, 0.65, 0.44 and 0.04 at ratios 1, 0.8, 0.71 and 0.5.
P(DIFFERENT) was 0.04 at ratio 1.

## Part 1 — Time step (decides the production τ)

### Design
- Ω = 120. **64 runs at τ and 64 at τ/2.** The runs use the full MG3 run protocol (T ≤ 500, the stop 50 after the
  first confirmed event, burn-in 50, censoring 480), so each arm yields post-burn-in first-passage data. Expected yield:
  about 36 events at τ and about 50 at τ/2, if τ/2 really has few startup events.
- **Registered calibration-dependent rule:** if Part 3 gives P(EQUIVALENT | ratio 1) < 0.80 for 64 runs per arm under
  p_s (0.35, 0.05), both arms are 96 runs instead. This is the only way calibration can change a sample size.

### Primary decision: post-burn-in rate ratio r = k(τ/2) / k(τ)
- Equivalence margin: **0.5 ≤ r ≤ 2** (see "Why this margin" below).
- **EQUIVALENT:** the 90% CI for r (the Part 3 method) lies inside [0.5, 2]. This is two one-sided tests at α = 0.05.
- **DIFFERENT:** not EQUIVALENT, and the 95% CI excludes 1.
- **INCONCLUSIVE:** otherwise.

### What each outcome does (all automatic; no decision after the data)
| Part 1 outcome | Next step | Production τ |
|---|---|---|
| EQUIVALENT | Part 2 at τ. The τ arm is Part 2's Ω = 120 cell | τ |
| INCONCLUSIVE | Part 2 at τ/2. The τ/2 arm is Part 2's Ω = 120 cell | τ/2, labelled "time-step agreement not established" |
| DIFFERENT | Part 1b: 64 runs at τ/4, then the same test for τ/4 against τ/2 | see below |

Part 1b (τ/4 against τ/2), same margin and rule:
| Part 1b outcome | Next step |
|---|---|
| EQUIVALENT | Part 2 at τ/2 (the τ/2 arm is the Ω = 120 cell) |
| DIFFERENT or INCONCLUSIVE | **MG3b stops.** Its rate results are INCONCLUSIVE on time-step grounds. Going on (τ/4 production, τ/8) would be a new, separately registered study, by PI decision |

This closes the case draft 1 left open, where τ/2 and τ/4 disagree. Going to τ/2 when the result is INCONCLUSIVE is a
deliberate choice. The tau-leap error falls as τ falls, so when agreement isn't shown, the smaller step is the safer
default. Its cost is known in advance (see Budget).

### Why this margin
- A ratio of 2 is about the width of MG3's per-cell rate CIs. It corresponds to about 10 units of Ω on the activated
  curve (S ≈ 0.073).
- A tighter margin, [0.67, 1.5], needs about 75 events per arm for 80% power. That means about 130 runs per arm, and
  Part 1 would cost about 37 core-hours. The precision of the slope gains little from that, so I recommend [0.5, 2].
- A τ effect that is a constant factor across Ω changes the prefactor but not the slope S. A τ effect that varies with
  Ω would bias S, and one comparison at Ω = 120 cannot rule that out. So **all Part 2 rates and the slope are reported
  as conditional on the production τ**, whatever Part 1 finds. EQUIVALENT is stated as "no time-step effect larger than
  a factor 2 at Ω = 120", not as time-step independence.

### Startup (descriptive, does not decide anything)
- Reported for each arm: the fraction of runs with a startup event and the full survival curve of the first confirmed
  event from t = 0. Also the difference between arms, with a Newcombe 95% CI and the Fisher exact p.
- Wording: "a difference was detected at the 0.05 level" or "no difference was detected at the 0.05 level". Neither
  wording implies agreement.

## Part 2 — Rates with more events

### Design
- Ω ∈ {130, 120, 110, 100} at the production τ. Ω = 140 is not included (ChatGPT's answer to question 3; Claude agrees).
- The Ω = 120 cell is Part 1's arm at the production τ: a fixed 64 (or 96) runs, with no event-count stopping.
- Other cells: batches of 16 runs until **≥ 20 post-burn-in events** at the end of a batch, **cap 64 runs**. A batch
  that has started is always completed, and its runs all count.

### Analysis (all through the Part 3 methods)
- Per cell: n, exposure, startup count, k̂ with its CI, and the risk set at t = 50, 100 and 150.
- **Capped cells:** a cell that reaches 64 runs with fewer than 20 events (including zero) **stays in the joint
  censored-likelihood fit** with its actual events and exposure. Nothing is excluded for being capped. Its k̂ is
  flagged "below target precision". Part 3's capped configuration calibrates this case. If the slope CI method is not
  eligible in the capped configuration, the slope is reported as "approximate" whenever a cell is capped.
- **Activated slope:** S from the Poisson GLM over the four cells, with the CI method chosen in Part 3. A deviance
  lack-of-fit statistic is reported (df 2), with its size taken from Part 3. Pre-registered sensitivity: the same fit
  without Ω = 100, because MG006 §8 called Ω = 100 a diagnostic boundary.
- **Burn-in landmarks** from the same histories: k̂ with left truncation at 100 and at 150. For each, the risk set
  (runs alive and without an event at the landmark), the events, the exposure and the CI are reported. A direct
  comparison of k on [50, 100) against [100, 480], from disjoint exposure, uses the Part 3 rate-ratio method.
- **Stationary check:** the conditional binomial with the split chosen in Part 3. It is labelled **informative** only if
  its power against a ×2 step at the split, at α = 0.01, is ≥ 0.5. That power is computed from the realised window
  exposures and the cell's k̂ (it uses exposures only, not how the events split). Otherwise it is labelled
  **uninformative**, and the p value is shown for completeness only.

### Predicted yield and cost (from MG3's rates and S = 0.073; scratch simulation of the batching rule)
| Ω | p_s assumed | mean runs | mean events | P(capped) | core-h at τ | core-h at τ/2 |
|---|---|---|---|---|---|---|
| 130 | 0.15 | 46 | 24 | ≈ 0 | 4.2 | 11.5 |
| 120 | (Part 1 arm) | 64 | 36 (τ) / 50 (τ/2) | — | (in Part 1) | (in Part 1) |
| 110 | 0.33 | 36 | 24 | ≈ 0 | 1.7 | 4.5 |
| 100 | 0.31 | 35 | 24 | ≈ 0 | 1.3 | 3.4 |

Ω = 130's rate is an extrapolation. If it is three times lower than predicted, the cell caps. The capped rule above
covers that case.

## Part 4 — removed (deferred)
Draft 1's Part 4 confirmed every birth over full trajectories. ChatGPT is right that this does not measure detector
precision. `confirm_field` classifies the whole continued field: once one particle exists, a later unrelated detection
also gets "persistent". `confirm_field` does return a crude attribution flag (`near_detection`: a persistent structure
within 2λ* of the detection), but it has not been validated, and it can be fooled by a pre-existing particle nearby or
by drift. Part 4 is therefore dropped from MG3b. A detector-precision study needs validated detection-level attribution
first, and it is a separate proposal. MG3b saves every confirmation field, so the data stay available.

This does not affect MG3's C2(a) as scored. In MG3, births were confirmed only up to and including a run's first
"persistent" outcome. So no counted birth was confirmed while an earlier confirmed particle was present. C2(a)'s stated
scope is unchanged.

## Budget and enforcement

### Measured costs (PI's machine, 9 workers, 2026-10-29)
0.89 core-s per simulated time unit per run at τ, and 2.43 at τ/2. τ/4 has not been measured; 5.0 is assumed. A
short τ/4 timing run, counted in the ledger, is done before Part 1, and the τ/4 projection and ceiling are rescaled to
it by the rule below. Each confirmation costs about 7 s (cloud figure). Projections include 15% overhead.

### Registered caps and ceilings
Run caps are binding. Core-hour ceilings are stop triggers set at about 1.5 × the projection.
| Item | Run cap | Projected core-h | Ceiling core-h |
|---|---|---|---|
| Part 3 calibration | — (synthetic) | < 1 | 3 |
| Part 1, τ arm | 64 (96 by the Part 1 rule) | 3.8 | 6 (9) |
| Part 1, τ/2 arm | 64 (96) | 13.0 | 20 (30) |
| Part 1b, τ/4 arm (only if DIFFERENT) | 64 | 26.7 (at 5.0 core-s/unit) | 40 × (measured τ/4 cost ÷ 5.0) |
| Part 2 at τ: Ω = 130, 110, 100 | 64 each | 7.2 | 12 |
| Part 2 at τ/2: Ω = 130, 110, 100 | 64 each | 19.4 | 30 |

| Path | Projected core-h | Ceiling | Wall time at 9 workers (approx.) |
|---|---|---|---|
| EQUIVALENT (τ) | ≈ 24 | 41 | ≈ 3 h |
| INCONCLUSIVE (τ/2) | ≈ 36 | 59 | ≈ 4 h |
| DIFFERENT, then τ/4 EQUIVALENT (τ/2) | ≈ 63 | 99 | ≈ 7 h |

### Enforcement rules (implemented in the driver and tested before any run)
1. **Ledger.** Every core-second on the production machine is written to `results/MG3b/budget_ledger.csv`: complete,
   partial, failed and excluded runs, references, timing runs, calibration and tests. Each line records the item, run
   key, wall time, worker count and platform. The evaluator **sums the ledger**. It does not infer compliance from the
   absence of a log marker (MG006 §8).
2. **Check before every batch, not only at cell start.** Before submitting a batch, the driver computes the used
   core-hours plus the projected cost of that batch (batch size × mean cost per completed run in that item so far, or
   the prior projection if none, × 1.25). If that exceeds the item's ceiling, it does not submit. It logs `[BUDGET]`,
   saves state, and stops. Runs already in flight finish and are counted. Part 1's arms are submitted in batches of 16
   for this purpose. The analysis still uses the full fixed sample.
3. **What a budget stop means.** The affected item is INCONCLUSIVE by the budget rule. At a budget stop the driver
   reports only run counts and costs, not rates or event counts. A PI extension can be decided on that basis and is
   recorded in `GATES.md` before any further run. It is reported as a budget exception, as for MG3.
4. **Escalation is pre-authorised.** Moving to τ/2 production or to Part 1b follows the Part 1 rule automatically. Each
   branch has the ceiling listed above. No branch beyond those listed can be entered within MG3b.
5. **Calibration-dependent changes** are limited to (a) the choice among the frozen candidates and (b) the 64 → 96 rule
   in Part 1. Both are written to `results/MG3b/calibration/decisions.json`, which is committed before the first lattice
   run. The driver reads it and refuses to start without it.
6. **Platform.** Costs and results record the machine, OS, Python and package versions, and the worker count.

## Order of work
1. ChatGPT reviews this draft; Claude responds and revises; the PI approves; the design is registered in `GATES.md` as
   MG3b, together with a pre-run implementation specification.
2. Part 3 calibration (synthetic). The results and `decisions.json` are committed and sent for review.
3. Driver with the ledger and budget checks, plus tests. The τ/4 timing run.
4. Part 1, then Part 1b if triggered, then Part 2. Evaluation, independent recomputation from the raw files, and the
   record.

## Questions for review
1. Is the equivalence margin [0.5, 2] on the post-burn-in rate ratio acceptable, given the cost of a tighter one?
2. Is going to τ/2 production on an INCONCLUSIVE Part 1 acceptable, or should INCONCLUSIVE also trigger Part 1b?
3. Is the "informative" rule for the stationary check (power ≥ 0.5 against a ×2 step, from realised exposures)
   adequate?
4. Are the Part 3 eligibility thresholds (≥ 0.94 coverage and ≤ 0.035 per tail in every constant-hazard configuration,
   R = 4000) appropriate, given the number of configurations (a method whose true coverage is 0.95 fails a given
   configuration by Monte Carlo chance about 0.2% of the time)?

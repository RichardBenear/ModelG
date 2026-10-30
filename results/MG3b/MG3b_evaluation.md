# MG3b — evaluation of the production run

*Labels corrected 2026-10-30 after ChatGPT's review (`review/MG008_chatgpt.md`, response `review/MG008_response.md`): the
three-cell sensitivity interval and the landmark and early/late intervals are nominal, not calibrated. No number changed.*

Claude, 2026-10-30. Registered design: `GATES.md` (MG3b) and `spec/MG3b_followup_proposal.md` (frozen at `202c00f`).
Methods: `results/MG3b/calibration/decisions.json`, with the addendum 1 labels. Run on the PI's machine (RB_LibraryPC,
Windows 11, Python 3.13.2, numpy 2.2.4, scipy 1.16.1, numba 0.68.0, 9 workers), 2026-10-29 local time. Evaluator:
`python src/evaluate_mg3b.py` → `evaluation.json`. Independent recomputation from the raw run records, sharing no code
with the evaluator: `python -I results/MG3b/independent/recompute.py` → `independent/recompute_output.json`. Every
number below agrees between the two.

MG3b is an estimation study. It has no pass/fail verdict and does not re-score MG3, whose FAIL stands.

## Conclusions first
1. **Time step: EQUIVALENT.** At Ω = 120, the post-burn-in rate at τ/2 is 0.95 times the rate at τ. The 90% CI is
   0.65–1.40, inside the registered margin [0.5, 2]. Production therefore ran at τ, and Part 1b was not triggered. This
   means no time-step effect larger than a factor of 2 at Ω = 120. It does not establish time-step independence, and
   all rates and the slope below are conditional on τ = 6.667e-4.
2. **Rates:** all three new cells reached the 20-event target. None was capped.

   | Ω | Runs | Events | Startup | Exposure | k̂ (per unit length and time) | 95% CI |
   |---|---|---|---|---|---|---|
   | 130 | 48 | 25 | 0 | 15 461 | 1.62e-5 | 1.05e-5 – 2.39e-5 (Garwood) |
   | 120 | 64 | 44 | 9 | 12 892 | 3.41e-5 | 2.48e-5 – 4.58e-5 (approximate: observed minimum branch-conditional coverage 0.943) |
   | 110 | 48 | 28 | 19 | 4 464 | 6.27e-5 | 4.17e-5 – 9.07e-5 (Garwood) |
   | 100 | 32 | 21 | 11 | 1 866 | 1.13e-4 | 6.97e-5 – 1.72e-4 (Garwood) |

3. **Activated slope:** S = **0.064** per unit Ω, profile-likelihood 95% interval **0.047–0.082**, labelled
   approximate (observed minimum coverage 0.937 over qualifying calibration cells; not a guaranteed floor). The
   deviance lack of fit is 0.23 on 2 df (p = 0.89), so the four cells are consistent with log k linear in Ω. Without
   Ω = 100 (the registered sensitivity fit): S = 0.068, nominal 95% profile interval 0.041–0.094. That interval is
   not calibrated: the 0.937 figure belongs to the four-cell fit only (review MG008, item 1). MG3's 0.073 (95% CI 0.052–0.094) lies inside both
   intervals. The Wald interval (0.047–0.082) is almost the same as the profile interval.
4. **Stationary check: uninformative in every cell.** None rejects at α = 0.01. The smallest p is 0.0125, at Ω = 110.
   The simulated power against a ×2 or ×0.5 change at the split is 0.04–0.26, far below the 0.5 needed for the
   informative label. As registered, this says nothing about whether the hazard is stationary.
5. **Cost:** 26.0 core-hours of production plus 0.75 for calibration, 26.75 in all. Part 2 used 9.94 core-hours, so it
   would also have fitted under the original 12-hour ceiling. The PI's ceiling extension was never exercised.

## Part 1 (time step)
- τ arm: 44 events, exposure 12 891.5. τ/2 arm: 39 events, exposure 12 017.5.
- r = k(τ/2)/k(τ) = 0.951. Clopper–Pearson (conditional binomial): 90% CI 0.646–1.397, 95% CI 0.602–1.497.
- The decision is EQUIVALENT, the same as the driver's (`plan.json`).
- **Startup (descriptive):** 9/64 (14%) at τ and 10/64 (16%) at τ/2. Difference −0.016, Newcombe 95% CI −0.142 to 0.111,
  Fisher p = 1.0: no difference was detected at the 0.05 level. MG3's 5/16 against 0/16 asymmetry did not recur. The
  first-event times of every run are in `evaluation.json` (`part1.first_event_times`), for the survival curves.

## Part 2 (rates), more detail
**Risk sets and landmarks** (k̂ with left truncation at the landmark). The landmark and early/late intervals are
**nominal 95% intervals and descriptive sensitivities**. Their sampling differs from the calibrated principal
estimates, so no coverage claim is made for them (review MG008, item 2). The landmark intervals are Garwood and the
early/late ones Clopper–Pearson.

| Ω | Risk set at 50 / 100 / 150 | k̂ from 100 | k̂ from 150 | k[100, 480] / k[50, 100), 95% CI |
|---|---|---|---|---|
| 130 | 48 / 45 / 44 | 1.67e-5 (22 events) | 1.92e-5 (21) | 1.28 (0.38–6.7) |
| 120 | 55 / 47 / 42 | 3.45e-5 (36) | 3.78e-5 (31) | 1.06 (0.48–2.6) |
| 110 | 29 / 25 / 18 | 7.70e-5 (24) | 8.55e-5 (17) | 2.60 (0.89–10.3) |
| 100 | 21 / 13 / 8 | 1.33e-4 (13) | 1.73e-4 (8) | 1.47 (0.57–4.1) |

In every cell the late-window rate is at least as high as the early one. No single comparison excludes 1, and they
are not independent tests. The landmark estimates at Ω = 100 rest on few runs (8 at the 150 landmark).

**Stationary check** (split (b), the cell's exposure median; power from 2000 simulations of the cell's own procedure):

| Ω | Split (t) | p | Power ×2 | Power ×0.5 | Label |
|---|---|---|---|---|---|
| 130 | 228 | 0.23 | 0.10 | 0.08 | uninformative |
| 120 | 195 | 0.45 | 0.26 | 0.16 | uninformative |
| 110 | 138 | 0.0125 | 0.08 | 0.04 | uninformative |
| 100 | 104 | 0.38 | 0.08 | 0.04 | uninformative |

**Startup fractions** differ a lot between cells: 0/48 at Ω = 130, 9/64 at 120, 19/48 (40%) at 110 and 11/32 (34%) at
100. At Ω = 110 and 100, a startup run gives no exposure, so these cells needed more runs per event. The fractions are
descriptive and enter no estimate.

**Prediction check (descriptive):** the registration predicted S near 0.073 and Part 1 EQUIVALENT. Both held. The
Ω = 120 rate (3.4e-5) is below MG3's 4.6e-5, but the two CIs overlap. MG3's per-cell CIs were wide (10 events).

**Confirmation fields saved:** P1_tau 79 (53 persistent, 26 excursion), P1_tau2 69 (49/20), Ω = 130 36 (25/11),
Ω = 110 65 (47/18), Ω = 100 53 (32/21). These are kept for later work and are not analysed here.

## Deviations and notes (stated plainly)
1. **PI budget decision after Part 1 was seen** (GATES.md, 2026-10-30 02:56 UTC). The decision removed the Part 2
   ceiling. It was not exercised, because Part 2 used 9.94 of the original 12 core-hours.
2. **Evaluator written after the data were in.** The analyses are those the proposal registers (Part 2 "Analysis").
   Details the proposal leaves open, chosen while writing the evaluator, are listed here so the reviewer can judge them:
   - Landmark k̂: runs alive and without an event at the landmark, with exposure counted from the landmark to the event
     or to 480.
   - The early/late comparison splits at t = 100, using the Part 1 rate-ratio method (Clopper–Pearson) at 95%.
   - Stationary power: hazard k̂ before the cell's observed exposure-median split and c·k̂ after it. The p̂_s is the
     cell's startup fraction. Ω = 120 is simulated as a fixed 64 runs, the others with stop20. 2000 replicates, seeds
     20261030 + cell index.
   - Slope fit with Ω centred at its mean. The profile interval comes from the χ²₁ cut, with α profiled out.
   - The Ω = 120 rate takes the addendum 1 label for the τ branch (0.943).
3. **Platform:** the run was on Windows with older numpy and scipy than the cloud machine, where the evaluator ran. The
   statistics use only the saved records, so they don't depend on that.

## Files
- `plan.json`, `run_log.txt`, `timing_tau4.json`, `clip_pilots.json` (both pilots clip fraction 0), `budget_ledger.csv`.
- `runs/<cell>/<rep>.pkl`: 258 run records, including every confirmation field.
- `evaluation.json`: the evaluator's output. `independent/`: the independent recomputation.

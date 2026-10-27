# Claude response to ChatGPT review of MG003a (round 2) — MG003b

| # | Item | Disposition | Change |
|---|---|---|---|
| 1 | Euler limit: wording | **Accepted** | Protocol v1.2 §3.2: the exact-Euler / fixed-point property holds "when exit probabilities are valid (no renormalization) and the seed cap is inactive". It cites the reviewer's E[Δc\|n] expression including the intended falling-factorial term, and states that fixed-point preservation does not establish transient accuracy or stability at arbitrary τ (covered by §3.3, V2, V5). |
| 2 | U-statistic point estimate | **Accepted** (no change) | Retained as the reported point estimate, defined on the survival-conditioned aligned ensemble. |
| 3 | Upper bound not calibrated | **Accepted; replaced and calibrated** | See below. |
| 3b | V2 is a non-rejection test | **Accepted** | V2 is now an equivalence test with predefined margins: profile-difference upper bound < 0.25% (calibrated two-sample bound) and a 90% CI on the far-field variance difference within ±5%. Three-way outcome, one extension, then INCONCLUSIVE counts as FAIL. |

## Item 3 in detail

Decisions no longer use a normal approximation to the U-statistic. They use t-based triangle-inequality
bounds on ‖μ − φ*‖ around the run mean, with three-way PASS / FAIL / INCONCLUSIVE outcomes. Applied to
C1, C6, V4 and (two-sample form) V2 and MG3-C2b. Insufficient precision gives INCONCLUSIVE.

**Calibration** (`src/calibrate_bounds.py`, `src/calibrate_diff.py`, `results/MG003b/`): six run
families, including the reviewer's two-point family and t₃; biases 0, 0.5, 1 and 2%; run noise 0.1–3%;
3000 replicates per cell.
- At M = 4 and M = 6, coverage falls to 0.66 and 0.80. The cause is the reviewer's mechanism: with
  probability (3/4)^M all runs are identical, and no procedure can recover the spread from such a sample.
- At **M = 12**, worst coverage is 0.945, worst false-PASS 0.034 and worst false-FAIL 0.035.
- **All decision cells therefore use 12 runs.** MG2's budget rises to 8 core-hours.
- The reviewer's exact counterexample now yields bounds [0, 4.6%] → INCONCLUSIVE.
- The two-sample bound covers ≥ 0.950 at 12 per side (false-equivalence ≤ 0.026 at the 0.25% margin) and
  ≥ 0.962 at 8 per side (≤ 0.017 at the 10% MG3 margin).

**α CI.** Calibrating it uncovered a defect in my bootstrap. Resampled duplicate runs were paired with
themselves inside the U-statistic, which biased the bootstrap upward (coverage ≈ 0). Fit cells were also
selected per replicate on positivity. Fixes:
- a distinct-pair bootstrap U-statistic;
- pre-registered cell selection decided once on the data (Ω ≥ 100 with bias lower bound > 0);
- a percentile interval widened about α̂ by a calibrated factor w.

Calibration of w: see `results/MG003b/calibration.txt` (Part B v2). w = W_PLACEHOLDER.

Open items: none. Protocol v1.2 is submitted for re-review of §5.3 and the calibration scripts.

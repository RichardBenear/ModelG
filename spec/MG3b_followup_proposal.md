# MG3b — Follow-up to MG3: nucleation rates, time step and uncertainty (DRAFT for review)

Status: **DRAFT** (Claude, 2026-10-29). For ChatGPT review, then PI approval and registration in `GATES.md` before any
run. MG3's registered FAIL is unchanged by anything here. This is a separate study, labelled as a follow-up designed
after MG3's data were seen.

Motivation (MG3 record and the MG006 review):
- Nucleation is confirmed and its rate rises steeply as Ω falls, but the rates rest on 5–11 events per cell.
- Time-step independence is unresolved: startup events at Ω = 120 were 5/16 at τ against 0/16 at τ/2.
- Stationarity was checked only weakly.
- The interval procedures have not been calibrated for stopped, left-truncated first-passage data.

Setting throughout: b = 29, h = 0.2, L = 100, homogeneous Poisson start. Detection, tracking (with the A4 fix) and
§5.4 field-level confirmation are as in MG3. Fresh streams: gate id 4.

## Part 1 — Burn-in and time-step check (decides the production τ)
- Ω = 120. **48 runs at τ and 48 at τ/2**, each to t = 70 (enough to confirm any event detected before t = 50).
- Statistic: the fraction of runs with a confirmed event detected at t < 50.
- Rule: two-sided Fisher exact test.
  - **p < 0.01** → early dynamics depend on τ. Part 2 runs at τ/2, and Part 1 adds 48 runs at τ/4 to check τ/2.
  - **Otherwise** τ is retained, and the MG3 asymmetry is reported as not reproduced.
  - The difference and its CI are reported either way.
- Cost: ≈ 3.3 core-hours (≈ 1.1 core-s per time unit at τ and 2.4 at τ/2, from MG3).

## Part 2 — Rates with more events
- Ω ∈ {130, 120, 110, 100}, at the τ set by Part 1.
- Batches of 16 runs per cell until **≥ 20 post-burn-in events** (cap 64 runs per cell). Burn-in 50 and censoring at 480
  as in MG3, with the early stop 50 after the first confirmed event.
- Reported per cell: k̂, the calibrated CI (Part 3), and the stationary-hazard check where there is late-window exposure.
  There is also a burn-in sensitivity, k̂ with burn-in 100.
- Activated-scaling fit over the 4 cells. The slope CI uses the method Part 3 calibrates.
- Predicted cost from MG3's rates: Ω = 130 ≈ 4.4, 120 ≈ 2.3, 110 ≈ 1.5, 100 ≈ 1.5, **total ≈ 10 core-hours at τ** (about
  double at τ/2). Ω = 130 is the expensive cell; dropping it saves 4.4 h but leaves 3 cells, as in MG3.

## Part 3 — Uncertainty calibration for the actual design (cheap; done first)
- Synthetic first-passage data with MG3/MG3b's exact sampling: left truncation at 50, censoring at 480, early stop,
  batch-wise event-count stopping. Constant hazards span the observed range; a non-constant (Weibull) hazard is the
  alternative.
- Coverage is measured for the Garwood-based rate CI, the conditional-binomial rate-ratio CI, the stationary-hazard
  test (size and power), and three activated-fit slope CIs: normal with estimated variances, t, and parametric
  bootstrap. The method used in Part 2 is the one with coverage ≥ 0.94.
- Cost: minutes, no lattice simulation.

## Part 4 (optional, descriptive) — full-trajectory confirmation fraction
- 16 runs at Ω = 120 to T = 500 **without** the early stop, with every birth confirmed. This gives the detector's
  confirmation fraction over full trajectories, the quantity MG3's C2(a) did not measure.
- It is reported as a scientific question about the 80% threshold. It does not re-score MG3.
- Cost ≈ 2.8 core-hours.

## Budget and options for the PI
| Option | Contents | Core-hours |
|---|---|---|
| A | Parts 1 + 2 + 3 at τ | ≈ 13 |
| A + 4 | as A, plus the full-trajectory study | ≈ 16 |
| B | Part 1 says τ/2: Parts 1 (with τ/4) + 2 at τ/2 + 3 | ≈ 26 |

Part 1 runs first, so the B case would be known before most of the budget is committed. If it occurs, the driver stops
for a PI decision.

## Questions for review
1. Is the Part 1 decision rule (Fisher p < 0.01) appropriate, or should it be an equivalence test with a margin?
2. Is event-count stopping (≥ 20 events, in batches of 16, cap 64) acceptable given Part 3's calibration?
3. Should a cell at Ω = 140 be included, despite costing ≈ 9 core-hours for 20 events?

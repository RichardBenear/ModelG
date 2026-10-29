# Phase 2 Protocol v1.3.3 — Amendment 4 (PROPOSED, MG005a): corrections before MG3

Status: **PROPOSED, revised after ChatGPT review** (`review/MG005_amendment4_chatgpt.md`: approve with two wording
qualifications, incorporated and marked [r]). Awaiting PI approval. If approved it is applied as **protocol v1.3.4**.
Origin: the MG005 review of the MG2 result (`review/MG005_chatgpt.md`, findings P1 tracking, P1 audit, and the MG3
stopping-rule cleanup). No MG3 gate run has been made. Nothing here re-scores MG2, whose registered FAIL stands.

## 1. Tracker: the first homogeneous birth was omitted (implementation correction)

**Defect.** `measure.track(..., x0=None)` gave its first new track id 0, and event generation treated id 0 as the
seeded original particle. In a homogeneous-start run the first persistent track therefore produced **no birth** (and
could produce a spurious "loss"). This contradicts the intended definition; it is not a new criterion.

**Fix.** Only the track created from a supplied `x0` is the original (`orig=True`). All other tracks, including the
first one in a homogeneous run, are ordinary tracks whose persistence ≥ 20 makes a birth.

**Regression tests** (`tests/test_track.py`, 8 tests, all pass on the fix): ChatGPT's input (one detection at x = 0,
t = 0.5–30 → one birth at 0.5); two homogeneous births; a short track (no birth); seeded start with no events; seeded
loss; seeded start plus a birth; `analyze_homogeneous` on synthetic fields built from the reference soliton (first
birth, and burn-in exclusion with `t_start`); a null field. On the unfixed code the `analyze_homogeneous` test fails,
returning births [70] instead of [10, 70]. Censoring at t = 480 is applied by the MG3 driver (not yet written), which
will carry its own test.

**Impact on earlier records: none.**
- MG2 seeded runs always pass `x0`, so their event lists are unchanged.
- MG2 C5 uses detection plus field confirmation, not the event list.
- MG2 N1 had no detections.
- V7 counts persistent tracks from track durations, not from the event list.

## 2. Gate-level A1 audit: the registered pooled audit cannot flag (replacement for MG3 onward)

**Defect.** Protocol §5.3 standardises each cell by its own mean and SD before pooling. For n = 12 a standardised
value is bounded by (n − 1)/√n = 3.18, so the pooled leave-one-out statistic over 13 cells is bounded by ≈ 3.43, while
the critical value is 4.65. The audit cannot flag any data set of this shape. Calibration confirms it (table below:
0 flags even with one run displaced by 20 run-level SDs). **V6's pooled 24-run audit has the same defect**
(bound ≈ 4.49 against 5.05). V6's per-cell A1 diagnostics, which are valid, remain the substantive part of its pass.
MG2's "audit ok" is recorded as non-informative.

**Replacement (for MG3 and later gates).** Within-cell audit (`stats_mg2.audit_within_cell`): each run's projection
is externally studentised against the other runs of its own cell (exactly t_{M−2} under i.i.d. Gaussian runs).
Bonferroni is applied over every run × field × direction test at overall level 1%.
- Applied descriptively to MG2 (`results/MG005/audit_corrected.txt`): 1404 tests, max |t| 4.2175 against 8.47; no
  run flagged. This agrees with ChatGPT's independent value (4.21751). It is not a retroactive registered replacement.
- **Calibration** (`src/calibrate_audit.py` → `results/MG005/audit_calibration.txt`; 13 cells × 12 runs × 9
  projections, 2000 simulations per row; one run in one cell displaced by δ run-level SDs). A separate 200 000-cell
  null check gives a gate-level false-alarm rate of 0.87% (nominal ≤ 1%).

| δ (run-level SDs) | registered pooled audit flag rate | within-cell audit flag rate |
|---|---|---|
| 0 (null) | 0.000 | 0.0145 (SE 0.003; the 200 000-cell check gives 0.0087) |
| 4 | 0.000 | 0.018 |
| 6 | 0.000 | 0.112 |
| 8 | 0.000 | 0.384 |
| 10 | 0.000 | 0.763 |
| 15 | 0.000 | 0.998 |
| 20 | 0.000 | 1.000 |

- **Stated limitation.** With 12 runs per cell this audit is the per-cell A1 run-level test at a stricter, gate-wide
  error rate (critical value 8.47 against the per-cell ≈ 5.9). It **adds no power** over the per-cell diagnostics:
  anything it flags would already fail the per-cell run-level check. In the calibration it detects a single outlying
  run reliably only at about 12–15 run-level SDs (76% at 10, 99.8% at 15). [r] That power curve applies to **one
  displaced run in one scalar projection under the simulated Gaussian model**. It does not measure sensitivity to every
  kind of rare mode.
- **Residual risk (A1).** The rare-mode risk that pooling was meant to address (protocol §5.3) is **not covered by any
  test**. Absence of rare modes remains an assumption. Passing these diagnostics does not establish it, and does not
  validate confidence bounds against such modes. This limitation is carried into the interpretation of MG3 and later
  gates. [r] The route considered here to a stronger pooled audit, borrowing within-run block variance, would require
  a separately justified scale estimate, because block variance is not automatically the variance of a run mean
  (temporal correlation, variation between runs). Other statistical designs are possible. None is proposed here.

## 3. MG3 Stage A: terminal outcomes on the amended grid

The gate text still said "down to Ω = 1", left over from the v1.3 grid. Replace it with the following. The
*screened set* is the grid cells actually run in Stage A that were not excluded by the crowding rule.
- **(a) No confirmed event in any screened cell.** This covers the grid exhausted at Ω = 50, or descent stopped by an
  excluded cell. **C1 FAIL**, reported as an exposure-limited bound over the screened Ω range, with 95% upper bounds
  per Ω. No Stage B. The pre-stated interpretation reads "no events at any screened Ω" (was "any Ω ≥ 1").
- **(b) Confirmed events, but no 4/4 cell before the screened set ends.**
  - S = screened cells with 1–3 of 4 confirmed.
  - Neighbours are added only from screened, non-excluded grid values. If |S| < 3 once these are exhausted, Stage B
    runs on the available cells. C3 is then INCONCLUSIVE under its existing rule (fewer than 3 cells with ≥ 5 events).
  - No midpoint is inserted, because there is no 0/4 → 4/4 transition to bracket.
- **(c) Midpoints** are inserted only between two adjacent screened, non-excluded grid values that bracket a
  0/4 → 4/4 transition. Each inserted point gets a V3 clip pilot and the crowding screen before its Stage-A cell.
- **(d) 4/4 already at Ω = 300.** The activated window lies above the registered grid. Stage B runs at Ω = 300 only,
  where C1 and C2 can be evaluated, and C3 is INCONCLUSIVE.
- **(e) Ω = 300 excluded by crowding.** There is no screened set, and MG3 is INCONCLUSIVE.

The crowding exclusion, hazard checks, field-level confirmation (v1.3.3 §5.4), V5 rare-event τ check and the
10 core-hour cap are unchanged.

## 4. MG2-informed expectations for MG3 (notes, no criterion change)

- At b = 29 and Ω ≤ 300, seeded particles are often lost within T = 100 (MG2: 5 of 24 runs eligible at Ω = 300). The
  C2(b) same-Ω reference ensemble therefore may not reach 8 eligible runs. If not, the pre-stated replacement
  ((a) plus (c)) applies, and is reported.
- At b = 29.16 and 29.3 (C4), particles that nucleate may later replicate (MG2 C4). MG3 uses only the first event per
  run, so this does not affect the rate estimate.

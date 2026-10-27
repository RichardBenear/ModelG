# Phase 2 Protocol v1.3 — Amendment 1 (PROPOSED, MG004)

Status: PROPOSED 2026-10-27. It requires a PI ruling, and ChatGPT review is recommended, before any
further validation or gate run. Protocol v1.3 stays frozen. If approved, this amendment is applied as
v1.3.1, and every affected step is re-run from scratch.

## What happened

- V3 (clip pilot) PASS, and V1 (exact SSA) PASS, at the protocol-specified settings
  (`results/MG004/v1.*`, `v3.*`).
- V2 started at its registered cell (b = 29, Ω = 300, τ). **Only 3 of 24 attempted runs were eligible.**
  The V2 batch was stopped. The log is kept unedited: `results/MG004/v2_attempt1_b29_aborted.log`.
- Diagnosis (`results/MG004/amendment1_diagnostics.txt`):
  - **Real noise-induced annihilation.** At b = 29, Ω = 300 the raw (unsmoothed) φ_X minimum rises from
    −9 to about −2.5 within 10–20 time units, and the particle is gone. This is not a pipeline artifact:
    - At Ω = 10⁵ the same particle holds at −9.0 for the whole run.
    - At b = 29.08, Ω = 300 it persists in 2 of 2 test runs.

    b = 29 sits only 0.07 above the operational lattice window edge b_s = 28.929, and at Ω = 300 noise
    carries it over the saddle.
  - **Detection-threshold specification flaw.** §5.1 sets the threshold at −0.5·|A*| with A* the
    *unsmoothed* lattice-ODE amplitude (−9.02 → −4.51), but applies it to the *smoothed* field. The
    smoothed reference particle is only −5.54 deep, so the threshold sits at **81%** of the actual depth
    instead of the intended half, which makes surviving particles flicker in and out of detection.

## Proposed changes

1. **Detection threshold (§5.1; affects MG2, MG3, V2, V4, V6).** Use −0.5·|A*_s|, where A*_s is the
   minimum of the *smoothed* lattice-ODE reference φ_X at that b and h: −2.77 at b = 29, −2.90 at
   b = 29.08. This matches the stated intent ("half depth"). MG3's deterministic confirmation step
   (§5.4) already guards against noise excursions, which the lower threshold makes more frequent at
   small Ω.
2. **V2 location.** Move V2 from (b = 29, Ω = 300) to **(b = 29.08, Ω = 300)**, the V6 parameter
   point (non-gate), with distinct stream keys. V2 tests time-step accuracy, which needs surviving
   particles. b = 29 at Ω = 300 cannot supply 12 eligible runs. All other V2 settings are unchanged.
3. **No change to the MG2/MG3 criteria.** Expect many b = 29 MG2 cells at low Ω to be INCONCLUSIVE
   ("insufficient eligible runs"), and C2 at b = 29 may be INCONCLUSIVE. The protocol already handles
   both outcomes. b = 29 survival is a reported measurement (MG2-C4), and the b = 29.16 criteria are
   unaffected.

## Disclosure

The aborted V2 attempt was a pre-registered validation, but its cell coincides with the MG2 cell
(b = 29, Ω = 300). It therefore revealed one MG2 *measurement* in advance: survival at that cell, about
3/24 to t = 100. That quantity is reported-only (not a pass/fail criterion), and nothing in this
amendment was chosen to change any criterion's outcome. The gate's own runs will use their own streams.

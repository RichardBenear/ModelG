# Phase 2 Protocol v1.3 — Amendment 1 (PROPOSED, MG004)

Status: **APPROVED 2026-10-27** by the PI (Richard Benear) after ChatGPT review
(`review/MG004_amendment1_chatgpt.md`), with that review's corrections incorporated (marked [rev]).
Applied as **protocol v1.3.1** (`spec/phase2_protocol_v1.3.1.md`). Protocol v1.3 is preserved unchanged.
All validations and gate runs from here on are identified as obtained under v1.3.1, and the aborted
V2 attempt is preserved as an aborted v1.3 attempt.

## What happened

- V3 (clip pilot) PASS, and V1 (exact SSA) PASS, at the protocol-specified settings
  (`results/MG004/v1.*`, `v3.*`).
- V2 started at its registered cell (b = 29, Ω = 300, τ). **Only 3 of 24 attempted runs were eligible.**
  The V2 batch was stopped. The log is kept unedited: `results/MG004/v2_attempt1_b29_aborted.txt`.
  [rev: renamed from .log, which the repo's .gitignore had excluded, so it was missing from GitHub.
  Content unchanged. The aborted attempt saved no run-level JSON; this log is the complete record.]
- Diagnosis (`results/MG004/amendment1_diagnostics.txt`):
  - **Raw-amplitude decay, consistent with noise-induced annihilation near the operational lower
    edge** [rev: wording]. In the illustrated runs at b = 29, Ω = 300, the raw (unsmoothed) φ_X minimum
    rises from −9 to about −2.5 within 10–20 time units. This is not a detector artifact:
    - At Ω = 10⁵ the same particle holds at −9.0 to t = 30.
    - At b = 29.08, Ω = 300 it persists in 2 of 2 test runs **to t = 30** [rev: horizon stated]. These
      runs support attempting V2 there; they do not predict that 12 runs will remain eligible to t = 100.

    b = 29 sits only 0.07 above the operational lattice window edge b_s = 28.929. A saddle crossing is
    the natural interpretation, but these traces do not directly demonstrate the mechanism.
  - [rev] "3 eligible of 24" is **not** an established survival rate. Eligibility combines tracking,
    detection and snapshot requirements, and the detector flaw below affects those classifications.
  - **Detection-threshold specification flaw.** §5.1 sets the threshold at −0.5·|A*| with A* the
    *unsmoothed* lattice-ODE amplitude (−9.02 → −4.51), but applies it to the *smoothed* field. The
    smoothed reference particle is only −5.54 deep, so the threshold sits at **81%** of the actual depth
    instead of the intended half, which makes surviving particles flicker in and out of detection.

## Proposed changes

1. **Detection threshold (§5.1; affects MG2, MG3, V2, V4, V6).** Use −0.5·|A*_s|, where A*_s is the
   minimum of the *smoothed* lattice-ODE reference φ_X at that b and h. [rev] A*_s is **computed in
   code from each reference with exactly the detector's smoothing operation** (`runner.reference` →
   `measure.smooth_periodic`, σ = 0.4), never hard-coded. Computed values: −2.7685 (b = 29),
   −2.8962 (b = 29.08), −2.9908 (b = 29.16).
   [rev] This is a substantive measurement change. It affects eligibility, survival classification and
   MG3 detection times; deterministic confirmation (§5.4) rejects excursions but does not make those
   measurements threshold-independent. All affected validations and gate runs are run fresh under
   v1.3.1. This matches the stated intent ("half depth"). MG3's deterministic confirmation step
   (§5.4) already guards against noise excursions, which the lower threshold makes more frequent at
   small Ω.
2. **V2 location.** Move V2 from (b = 29, Ω = 300) to **(b = 29.08, Ω = 300)**, the V6 parameter
   point (non-gate), with distinct stream keys. V2 tests time-step accuracy, which needs surviving
   particles; at b = 29, Ω = 300 the attempt yielded 3 eligible of 24 runs. b = 29.08 was already
   designated for V6, which gives an independent reason for choosing it. All other V2 settings are
   unchanged, **including the 12-eligible / 24-attempt cap and the stopping rule** [rev]. If V2 is
   infeasible there too, it is reported as such under the existing rules. [rev] A V2 pass establishes
   time-step accuracy **at b = 29.08, Ω = 300 only**. It does not establish the accuracy of annihilation
   probabilities near b = 29, which remain separately reported measurements.
3. **No change to the MG2/MG3 numerical tolerances.** [rev] Outcome handling when data are missing:
   - **MG2-C3 data rule:** far-field Ω·Var is measured on the first 12 *attempted* runs of each cell, in
     stream order and whether or not they are eligible. It uses all snapshots with t ≥ 10 and all nodes at
     circular distance > 20 from every detected particle (all nodes when none is detected). A cell where
     fewer than 12 runs supply ≥ 100 such snapshots gives a C3 value of **INCONCLUSIVE for that cell**.
     C3 then cannot PASS: it is INCONCLUSIVE overall unless some cell FAILs. Missing data never becomes
     a pass or a silently omitted cell.
   - C2 is INCONCLUSIVE with fewer than 3 selected cells (unchanged).

   Further expectations: Expect many b = 29 MG2 cells at low Ω to be INCONCLUSIVE
   ("insufficient eligible runs"), and C2 at b = 29 may be INCONCLUSIVE. The protocol already handles
   both outcomes. b = 29 survival is a reported measurement (MG2-C4), and the b = 29.16 criteria are
   unaffected.

## Disclosure

The aborted V2 attempt was a pre-registered validation, but its cell coincides with the MG2 cell
(b = 29, Ω = 300). It therefore gave advance information about one MG2 *measurement* at that cell:
3 of 24 runs eligible under the old detector. That is not a survival rate [rev], and survival at b = 29
is reported-only. The amendment was not designed to steer any criterion. [rev] Changing the detector
**can** change criterion outcomes, since eligibility, survival classification and detection times all
depend on it. That is why every affected step is re-run under v1.3.1. The gate's own runs use their own
streams.

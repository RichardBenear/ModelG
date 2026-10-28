# Phase 2 Protocol v1.3.2 — Amendment 3 (PROPOSED, MG004e): V7 checkpoint findings

Status: **APPROVED 2026-10-28** by the PI (Richard Benear), after ChatGPT review (`review/MG004_amendment3_chatgpt.md`),
with that review's corrections incorporated (marked [r]). Applied as **protocol v1.3.3**. Originally: PROPOSED for
the PI checkpoint required before MG3.
No gate run has been made. Validation status under v1.3.2: V1, V2, V3, V4, V6 PASS; V7 (descriptive)
done, `results/MG004/v7.*`.

## V7 findings (b = 29, homogeneous start, L = 100, T = 100, 4 runs per Ω)

| Ω | detections / snapshot | snapshots with any | tracks | persistent (≥ 20) | classification of up to 10 per run |
|---|---|---|---|---|---|
| 1 | 25.3 | 100% | 304 | 286 | 40 unresolved |
| 5 | 21.1 | 100% | 264 | 208 | 40 unresolved |
| 10 | 19.3 | 100% | 216 | 160 | 40 unresolved |
| 30 | 13.0 | 98.9% | 200 | 141 | 40 unresolved |
| 50 | 8.3 | 96.6% | 192 | 91 | 40 unresolved |
| 100 | 1.10 | 53.2% | 56 | 11 | 2 excursion, 9 unresolved |
| 300 | 0.00 | 0% | 0 | 0 | — |

1. **Three detector regimes over the observed exposure** [r: wording]: **crowded** (Ω ≤ 50),
   **transitional** (Ω ≈ 100), **quiet** (Ω = 300).
   - **Ω ≤ 50: crowded.** The field carries 8–25 persistent minima, i.e. up
     to about one per Turing wavelength (100/λ* ≈ 33). This is **consistent with**
     noise-sustained ("stochastic") Turing patterning below the deterministic threshold. That phenomenon
     is established in the Brusselator (Biancalani, Fanelli & Di Patti, arXiv:0910.4984), Model G's X–Y
     core, which is a **precedent, not verification** of Model G's mechanism [r]. It is a **candidate
     interpretation pending E1**.
   - **Ω ≈ 100: transitional.** Sparse, mostly short-lived detections; 11 persistent tracks in 400
     time units of ring time.
   - **Ω = 300: quiet.** No detection in 400 time units.
2. **The confirmation criterion (§5.4) cannot confirm what it should.** Every classified persistent
   detection at Ω ≤ 50 was "unresolved". Yet the seed-free continuations **reach soliton depth** (median
   φ_X minimum −9.4 to −9.7 against A* = −9.02); they are simply still evolving at t = 300 (median
   max |∂φ/∂t| 2e-4 to 3e-2), as multi-soliton states that keep rearranging slowly. The requirement
   max |∂φ/∂t| < 10⁻⁶ by t = 300 is the wrong test of "crossed the separatrix". It would also affect
   MG2-C5 (seeded nucleation is judged by the same rule).
3. **Consequence for MG3 as registered.** Its first-passage model assumes rare single detections out
   of a quiet homogeneous state, and its Stage-A stop rule needs confirmed events.
   - **V7 confirmed zero events under the registered rule**, and Ω = 100's persistent detections remain
     mostly unresolved. Whether, and where, isolated rare nucleation occurs is **not yet established**
     [r: "nucleation is everywhere" and "only around Ω = 100–300" withdrawn as premature]. What V7 does
     show is that the crowded regime makes isolated-event bookkeeping impossible at Ω ≤ 50 for this
     detector.
   - The registered grid has only 300 and 100 in the transitional range.
   - With the registered confirmation rule, Stage A would scan down to Ω = 1 without triggering its
     stop rule.

## Proposed changes

1. **Confirmation criterion (§5.4; MG2-C5 and MG3), as finalized** [r]. Implemented in `src/confirm.py`.
   - **Field-level.** The test applies to the **whole field** at the confirmation time (detection + 20).
     Several detections from the same snapshot share one outcome and are **not** independent
     confirmations; each field is counted once [r]. A detection-level attribution flag records whether a
     persistent structure lies within 2λ* of the detection.
   - **Outcome name:** "**deterministically persistent particle-depth structure over continuation times
     200–300**" (class `persistent`) [r]. The seed-free lattice ODE runs to t = 300. A structure is
     persistent if deep minima (unsmoothed φ_X < −0.5|A*|) can be **linked as the same structure** across
     all 11 checks t = 200, 210, …, 300, with each step's circular displacement ≤ λ*. This identity
     requirement answers the reviewer's point that a minimum need not be the same structure each time.
     `excursion`: max |φ| < 10⁻³ at t = 300. `unresolved`: otherwise.
   - **Interpretation:** basin membership (the state ends in particle territory), **with a stated
     finite-horizon limitation** [r]. A sufficiently long transient could be misclassified.
   - **Controls** (`src/confirm_controls.py`, `results/MG004g/confirm_controls.txt`). Ground truth for
     each control is the deterministic fate after continuing to t = 2000. Results: **15 of 15 agree.**
     They cover a stationary single soliton; two-soliton states at separation 20, 9 and 6.6 (the last
     spawns a third particle, matching Pulver & LaViolette's observation at 6.628); a two-soliton overlap
     at separation 5 that annihilates; four irregularly spaced solitons, slowly evolving; seeded states at
     0.5, 0.9, 0.99 and 0.999 of the lattice threshold, all of which decay; 1.01 and 1.1, which persist;
     and φ* continued just below the operational window edge (b = 28.928, 28.92, 28.90), which decays.
     No control showed a long transient misclassified at the t = 300 horizon. That supports, but does
     not prove, the finite-horizon interpretation.
2. **MG3 Ω grid, focused on the activated window.** Stage A grid {300, 250, 200, 150, 120, 100, 80, 65,
   50}, descending, with the same stop rule and Stage B rules.
   - **Crowding exclusion (new, pre-declared)** [r: renamed; it is a sampling rule, not a diagnosis of
     patterning]: a cell whose mean detections per snapshot exceed 2 in the window t < 50 is **excluded**
     from the isolated-event study.
     - It is applied to **every** Stage-A and Stage-B candidate, including inserted midpoints.
     - An excluded cell **cannot re-enter** the rate fit through the neighbour-selection rule.
     - Excluded cells are reported with their detection statistics.
     - Reaching an excluded cell **stops the descent**. This is a pre-declared budget choice; it does
       not assert that every lower Ω is more crowded.
   - **Clip pilots** (V3 procedure) were completed for the new grid values; see
     `results/MG004g/v3_newgrid.txt`.
   - **Expected consequence:** the C3 fit, which needs ≥ 3 Stage-B cells with ≥ 5 events, may be
     INCONCLUSIVE if the activated window is narrower than three grid points. That outcome is already
     defined, and is acceptable.
3. **New exploratory study (not a gate): E1, stochastic Turing patterns at Ω ≤ 50.** Structure factor
   S(k) of φ_X on homogeneous-start runs. Compare the peak position and height with the LNA spectrum
   (`lattice.lna_*` machinery, mode-resolved). Report the Ω at which patterning sets in.
4. **No other changes.** MG2's criteria are unaffected except C5's confirmation rule. MG2 cells at
   Ω = 10 and 30 are expected to be INCONCLUSIVE (no isolatable particle in the patterned regime), and
   are reported as such.

## Physics note [r: corrected]

LaViolette's picture has particles nucleating out of a fluctuating subquantum ether. In this discrete
token version the fluctuation strength is set by Ω. What the evidence supports so far:
- **Increasing Ω suppresses intrinsic fluctuations.** That makes spontaneous nucleation rarer, and
  noise-induced loss of a stable seeded soliton less likely. At b = 29.08, Ω = 3000, 12 of 12 runs kept
  their particle (V6), and the deterministic fixed point is independently verified. Finite runs cannot
  establish that spontaneous nucleation never occurs at large Ω.
- **At small Ω** (≤ 50 at b = 29), the detector is crowded. Noise-sustained patterning is the candidate
  interpretation, pending E1.
- **Whether isolated, rare nucleation exists in an intermediate window, and where, is what MG3 will
  measure.** It is not yet established.
- At b = 29 specifically, noise-induced loss of seeded particles at Ω = 300 is substantial (survival to
  t = 100 well below 1 at b = 29.08; see MG004 records). It is weaker at larger Ω.

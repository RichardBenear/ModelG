# Phase 2 Protocol v1.3.2 — Amendment 3 (PROPOSED, MG004e): V7 checkpoint findings

Status: PROPOSED 2026-10-28, for the PI checkpoint required before MG3. ChatGPT review recommended.
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

1. **Three regimes.**
   - **Ω ≤ 50: there is no quiet homogeneous state.** The field carries 8–25 persistent minima, i.e. up
     to about one per Turing wavelength (100/λ* ≈ 33). This is consistent with **noise-sustained
     ("stochastic") Turing patterns** below the deterministic threshold, a known phenomenon in the
     Brusselator (Biancalani, Fanelli & Di Patti, arXiv:0910.4984), which is Model G's X–Y core. It is
     not yet verified here by structure-factor analysis.
   - **Ω ≈ 100: transitional.** Sparse, mostly short-lived detections; 11 persistent tracks in 400
     time units of ring time.
   - **Ω = 300: quiet.** No detection in 400 time units.
2. **The confirmation criterion (§5.4) cannot confirm what it should.** Every classified persistent
   detection at Ω ≤ 50 was "unresolved". Yet the seed-free continuations **reach soliton depth** (median
   φ_X minimum −9.4 to −9.7 against A* = −9.02); they are simply still evolving at t = 300 (median
   max |∂φ/∂t| 2e-4 to 3e-2), as multi-soliton states that keep rearranging slowly. The requirement
   max |∂φ/∂t| < 10⁻⁶ by t = 300 is the wrong test of "crossed the separatrix". It would also affect
   MG2-C5 (seeded nucleation is judged by the same rule).
3. **Consequence for MG3 as registered.** Its first-passage model assumes rare single nucleations out
   of a metastable homogeneous state, and its Stage-A stop rule needs confirmed events.
   - At b = 29, that regime exists only in a narrow band of roughly Ω ≈ 100–300.
   - The registered grid has only 300 and 100 in that band.
   - With the current confirmation rule, Stage A would scan down to Ω = 1 without ever triggering its
     stop rule, and C1 would FAIL for "no confirmed events" while nucleation is in fact everywhere.

## Proposed changes

1. **Confirmation criterion (§5.4; MG2-C5 and MG3).** Seed-free lattice-ODE continuation from the field
   at detection + 20:
   - **confirmed** if a soliton-depth minimum (φ_X < −0.5|A*|) is present at all 11 check times
     t = 200, 210, …, 300;
   - **excursion** if max |φ| < 10⁻³ at t = 300;
   - **unresolved** otherwise.

   The number of soliton-depth minima at t = 300 is recorded, so single- and multi-particle outcomes
   are distinguished. Rationale: this tests basin membership (the state ends in particle territory)
   rather than full relaxation of slow multi-particle dynamics.
2. **MG3 Ω grid, focused on the activated window.** Stage A grid {300, 250, 200, 150, 120, 100, 80, 65,
   50}, descending, with the same stop rule and Stage B rules.
   - **Regime screen (new, pre-declared):** a Stage-A cell whose mean detections per snapshot exceed 2
     in the burn-in window t < 50 is classed as **noise-sustained patterning**. It is reported, never
     enters the activated-rate fit (C3), and stops the descent, since lower Ω is more patterned.
   - **Expected consequence:** the C3 fit, which needs ≥ 3 Stage-B cells with ≥ 5 events, may be
     INCONCLUSIVE if the activated window is narrower than three grid points. That outcome is already
     defined, and is acceptable.
3. **New exploratory study (not a gate): E1, stochastic Turing patterns at Ω ≤ 50.** Structure factor
   S(k) of φ_X on homogeneous-start runs. Compare the peak position and height with the LNA spectrum
   (`lattice.lna_*` machinery, mode-resolved). Report the Ω at which patterning sets in.
4. **No other changes.** MG2's criteria are unaffected except C5's confirmation rule. MG2 cells at
   Ω = 10 and 30 are expected to be INCONCLUSIVE (no isolatable particle in the patterned regime), and
   are reported as such.

## Physics note

LaViolette's picture has particles nucleating out of a fluctuating subquantum ether. In this discrete
token version the fluctuation strength is set by Ω.
- At **small Ω**, the "ether" is not a calm medium at all: it is a noise-sustained pattern carrying
  particle-depth structures everywhere.
- **Isolated, rare particle nucleation**, the SQK scenario, exists only in a window of intermediate Ω
  (about 100–300 at b = 29, on this lattice).
- At **large Ω**, particles neither nucleate spontaneously nor (at b = 29) survive long without decaying.

Mapping this window, and how it depends on b, is now MG3's real content.

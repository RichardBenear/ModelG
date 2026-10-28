# Phase 2 Protocol v1.3.1 — Amendment 2 (PROPOSED, MG004b)

Status: PROPOSED 2026-10-28. It needs a PI ruling, and ChatGPT review is recommended. No gate run has
been made. v1.3.1 and all validation records stay as registered.

## Validation record under v1.3.1

| Validation | Verdict | Record |
|---|---|---|
| V3 clip pilot | **PASS**: τ/2 for Ω ≤ 5; τ for Ω ≥ 7; nothing excluded | `results/MG004/v3.*` |
| V1 exact SSA | **PASS**, 14/14 statistics. The X-variance test at Ω = 30 is weak (±45%) | `results/MG004/v1.*` |
| V2 τ-halving (b = 29.08, Ω = 300) | **FAIL as registered** | `results/MG004/v2.*`, `v2_v1.3.1_run.txt` |
| V4 alignment | **INCONCLUSIVE (blocking)** | `results/MG004/v4.*`, `v4_v1.3.1_run.txt` |
| V6 A1 pilot (b = 29.08, Ω = 300 / 3000) | **PASS**: per-cell diagnostics pass for G, X, Y in both cells; pooled audit max \|t\| 3.1 vs critical 5.05 | `results/MG004/v6.*`, `v6_v1.3.1_run.txt` |

### V2
- Eligible runs: 22 of 48 attempted in each arm, after the one extension.
- Far-field variance equivalence **passed** for all fields. The 90% CIs of the τ vs τ/2 difference are
  well inside ±5%, e.g. X: [−2.34, +2.02] against a margin of ±6.17.
- Profile differences: point estimates 11.9% (G), 1.7% (X), 1.9% (Y); upper bounds 37.8%, 6.6%, 6.8%;
  margin 0.25%.
- The registered FAIL comes from the rule "FAIL if a point estimate lies outside its margin". The point
  estimate ‖m̄_τ − m̄_τ/2‖ is inflated by sampling noise, since its expected square is ‖Δμ‖² +
  tr(Σ_A/M_A + Σ_B/M_B). The bound-based rule used everywhere else in §5.3 (FAIL only if the lower
  bound > margin) would give **INCONCLUSIVE**: the lower bounds are 0. In both readings the experiment
  cannot resolve a 0.25% margin at Ω = 300, with half-widths of 5–26%. **Time-step accuracy of the
  profile was neither demonstrated nor refuted.**

### V4
Alignment-induced bias bounds with 12 runs:

| Ω | G upper | X upper | Y upper |
|---|---|---|---|
| 10³ | 22.4% | 3.1% | 3.2% |
| 10⁴ | 7.4% | 0.98% | 1.0% |

The tolerance is 0.3%. The half-widths scale as Ω^(−1/2) (G: 15.2% → 4.9%). Point estimates are small
(X and Y 0.00% at 10⁴), and the A1 diagnostics pass everywhere.

**At Ω = 100 no synthetic run was eligible.** With the amended threshold (−2.77), noise minima in the
homogeneous field fall below the threshold, producing multiple detections per snapshot.

## Root cause (a design-level power failure, found before any gate run)

The calibration studies (MG003b–c) assumed run-level noise of 0.1–3% but never measured the actual
noise of aligned run profiles. Two things push it high:
- **G has a small amplitude** (smoothed core about −0.17), so its relative noise is large;
- particle wandering and slow near-threshold modes limit how much averaging over a T = 100 run helps.

As registered:
- **MG2-C1 and C6 at Ω = 10⁴ with 12 runs are very likely INCONCLUSIVE** for G (projected half-width
  about 5%). MG2 could therefore not PASS as designed.
- **Cells with Ω ≤ 100 will mostly be INCONCLUSIVE** (no isolatable particle).
- **C2 (α) is likely INCONCLUSIVE:** bias ∝ 1/Ω while noise ∝ Ω^(−1/2), so bias is only significant at
  low Ω, where detection fails.

These are permitted outcomes under the protocol, but running a gate that is predictably INCONCLUSIVE
wastes about 8 core-hours.

## Proposed changes

1. **Move MG2's decision Ω to 10⁶ (C1, C6) and V4's to 10⁶.** Add 10⁵ and 10⁶ to the MG2 ladder; V4
   ladder {10⁴, 10⁵, 10⁶}.
   - Kernel cost does not grow with Ω (measured 0.85 / 0.66 / 0.64 s per unit time at 10⁴ / 10⁵ / 10⁶).
   - Projected 12-run half-widths at 10⁶ from the measured Ω^(−1/2) scaling: G ≈ 0.49%, X ≈ 0.07%,
     Y ≈ 0.07% (V4 at 0.3% resolvable for X and Y; G marginal).
   - Tolerances unchanged: C1/C6 ε = 1%; V4 ε = 0.3%.
   - **Alternative 1b, if the PI prefers:** apply V4's 0.3% to X and Y, and give G a tolerance of 1%
     (G is reported, never dropped).
2. **V2 redesign.**
   - Replace "FAIL if point estimate outside margin" with the §5.3 bound rule everywhere it appears
     (V2 profiles and variances): FAIL only if the lower bound exceeds the margin.
   - Move the profile-equivalence test to **b = 29.08, Ω = 10⁶**, where C1/C6 are decided (margin
     0.25%; 12 + 12 runs; one extension).
   - Keep the Ω = 300 variance-equivalence result, which already passed, as the fluctuation-accuracy
     evidence.
   - Rationale for C2: in this scheme the τ-induced change of the stationary mean enters at O(τ/Ω),
     the same 1/Ω order as the finite-Ω bias. It rescales the bias amplitude by 1 + O(τ·rate) but
     cannot change the exponent α that C2 tests.
3. **New pre-gate check V7 before MG3 (detector characterization at small Ω).** On homogeneous runs
   at b = 29 for each MG3 grid Ω (1–300), 4 runs, T = 100, report:
   - the rate of sub-threshold detections, their persistence-time distribution, and the fraction
     persisting ≥ 20;
   - the tracking load (detections per snapshot).

   V7 is descriptive; it sets no thresholds. If persistent noise minima are common, MG3's design
   (smoothing scale, or detection on deterministic confirmation alone) returns to the PI before any
   MG3 run.
4. **No change** to the reaction model, the integrator, A1–A2, or the C3/C4/C5 criteria.

## What this does not fix

- MG2 cells at Ω ≤ 100 will still mostly be INCONCLUSIVE because no particle can be isolated. They are
  reported as such.
- C2 may still be INCONCLUSIVE. Under the protocol that does not block the gate.

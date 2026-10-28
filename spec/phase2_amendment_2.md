# Phase 2 Protocol v1.3.1 — Amendment 2 (PROPOSED, MG004b)

Status: PROPOSED, revision 2, 2026-10-28. Revised after ChatGPT review (`review/MG004_amendment2_chatgpt.md`);
changes are marked [r2]. It needs a PI ruling. No gate run has been made. v1.3.1 and all validation records
stay as registered: V2's registered FAIL and V4's INCONCLUSIVE remain on record, and the scientific reading
of V2 is that profile equivalence was **unresolved**.

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

## Proposed changes (revision 2)

### 1. Bounded high-Ω allocation for the precision tests [r2: power recomputed]

Projection basis: the half-widths and upper bounds measured in V2 and V4, scaled as (Ω·M)^(−1/2).
These are approximate (the t critical values change slightly with M).

Measured at Ω = 10⁴ with 12 runs (V4): G upper bound 7.36%, half-width 4.87%; X and Y about 1.0% / 0.71%.

**Revision 1 was wrong for G.** At Ω = 10⁶ the projected G upper bound is ≈ 0.74% for V4 (tolerance
0.3%) and the V2 half-width ≈ 0.61% (margin 0.25%). Both are unresolvable.

**Proposal: one decision point, Ω = 10⁸**, for MG2-C1, MG2-C6, V4 and the V2 profile test.

| Test | Projected upper bound at Ω = 10⁸ (12 runs per arm) | Tolerance |
|---|---|---|
| V4 alignment, G | ≈ 0.074% | 0.30% |
| V4 alignment, X / Y | ≈ 0.01% | 0.30% |
| V2 profile difference, G (two-sample, upper ≈ 2× half-width) | ≈ 0.09% | 0.25% |
| MG2-C1 / C6, G (noise part; intrinsic bias ∝ 1/Ω is negligible) | ≈ 0.07% | 1% |

- **Cost check** (`results/MG004/amendment2_highomega_check.txt`): 0.77 s per unit time at Ω = 10⁸
  (1.8 s at 10⁷, probably process contention); no clipping; counts ≤ 2×10⁹ fit in int64. The
  single-snapshot X deviation from φ* falls as Ω^(−1/2) (0.35% → 0.11%, 10⁷ → 10⁸).
- **Bounded allocation and stopping rule:** each decision cell gets exactly 12 eligible runs (cap 24
  attempted) plus at most **one** 12-run extension, as in §5.3. There is **no further increase of Ω or
  M** after an INCONCLUSIVE result: it stands.
- **Tolerances unchanged:** C1/C6 ε = 1%; V4 0.3% for **all** fields, G included [r2: no relaxation
  for G, per the reviewer's recommendation]; V2 margin 0.25%.
- The registered Ω = 10⁴ cells remain in the MG2 ladder and are reported. Ω = 10⁵, 10⁶ and 10⁸ are
  added. V4's ladder becomes {10⁴, 10⁶, 10⁸}, decided at 10⁸.
- **What 10⁸ establishes:** convergence of the token rules to the lattice mean field (C1) and the
  absence of alignment artifacts at that noise level (V4). The *approach* to the limit is C2's job.

### 2. V2 decision rules [r2: signed rule for variances]

- **Profiles** (unsigned norm): PASS if the upper bound < margin; FAIL if the lower bound > margin;
  INCONCLUSIVE otherwise. This replaces the noise-inflated point-estimate rule, consistent with §5.3.
- **Variances** (signed difference d = Var_τ/2 − Var_τ, 90% CI [lo, hi], margin m = 5% of Var_τ):
  - PASS if −m < lo and hi < m (CI entirely inside the margins);
  - FAIL if lo > m or hi < −m (CI entirely outside);
  - INCONCLUSIVE otherwise.
- **Location:** the profile test moves to b = 29.08, Ω = 10⁸ (12 + 12 runs; one extension). The
  Ω = 300 variance equivalence, which passed for G, X and Y, is retained as evidence **at that operating
  point only** [r2].

### 3. Asymptotic wording on τ and α [r2]

Under a regular weak-noise expansion, bias(Ω, τ) = B(τ)/Ω + O(Ω⁻²). The time step changes the leading
coefficient B(τ) but preserves the asymptotic exponent, **provided B(τ) ≠ 0**.
- This does **not** guarantee an unchanged *fitted* exponent over the actual finite-Ω ladder, given
  survival conditioning, near-threshold dynamics and higher-order terms.
- A high-Ω V2 pass establishes the stated absolute time-step accuracy at Ω = 10⁸. It does **not**
  validate the low-Ω bias coefficient.
- C2's α is therefore reported as the exponent of the τ-scheme's bias over the tested ladder.

### 4. Expected MG2 outcomes, stated accurately [r2]

Raising the decision Ω makes C1, C6 and V4 resolvable. **It does not make MG2 capable of passing by
itself:**
- **C3 is mandatory.** At Ω = 100 (and possibly 300), detector crowding (V4: no isolatable particle at
  Ω = 100) may leave too little far-field data. The C3 cell is then INCONCLUSIVE, and C3 cannot PASS.
- **C4 is a hard criterion.** Every b_mid run at every Ω ≥ 100 must keep its original particle as the
  only particle. **A confirmed loss or replication is a FAIL**, not an insufficient-data INCONCLUSIVE.
  Validation data at b = 29.08, Ω = 300 (eligibility 22/48, 22/48, 12/18 under the amended detector)
  make losses at b = 29.16 for Ω ≤ 300 likely.
- **C4 is not changed.** Revising it after seeing validation data would be post hoc. A C4 FAIL would be
  a legitimate, informative outcome: particles at the window centre are not stable against intrinsic
  noise at Ω ≤ about 300.
- C1, C2, C3, C5 and C6 are evaluated and reported regardless. The gate verdict follows the registered
  rule.

### 5. V7 detector characterization, then a PI checkpoint before MG3 [r2]

- Homogeneous runs at b = 29 for each MG3 grid Ω (1–300): 4 runs, T = 100.
- Every persistent detection (≥ 20 time units) is classified by §5.4 deterministic confirmation as
  **confirmed** (relaxes to a soliton; this is genuine spontaneous nucleation, i.e. signal), an
  **excursion** (relaxes to homogeneous), or **unresolved** (neither by t = 300).
- Also reported: the sub-threshold detection rate, the persistence-time distribution, and the tracking
  load (detections per snapshot).
- V7 is descriptive and sets no thresholds. **Explicit PI checkpoint:** the PI reviews V7 before any
  MG3 run. A noisy detection workload is a measurement-design issue, to be resolved by then; confirmed
  nucleations are early signal for MG3 and are reported as such.

### 6. No change

The reaction model, the integrator, assumptions A1–A2, and the C3, C4 and C5 criteria are unchanged.

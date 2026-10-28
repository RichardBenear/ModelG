# Stage Log

| Stage | Date | Phase | Summary | Gate | Commit |
|---|---|---|---|---|---|
| MG000 | 2026-10-26 | 0 | Repo scaffold; charter (README), rules spec v1.0 frozen, gate table created | — | 1bc0015 |
| MG001 | 2026-10-27 | 0 | Published Model G pinned (Pulver & LaViolette 2013); linear check confirms subcritical, 3.5% below Turing threshold; rules v1.1 adds G⇌X back-reaction; MG1 criteria proposed | MG1 proposed | 110bf8b |
| MG001a | 2026-10-27 | 0 | ChatGPT review of MG001: 6 items + 2 minor, all accepted; w scaling fixed; integral identity derived (∫φ_G = ∫φ_X = 0 in infinite space); gradient benchmark and measurement defs added | MG1 proposal amended | 5d74cfa |
| MG001b | 2026-10-27 | 0 | Review round 2: identity confirmed (det q−gp); three conclusions narrowed; MG1 diagnostic adds time-derivative balances; MG1 pre-registered with PI sign-off | MG1 pre-registered | d56026a |
| MG002 | 2026-10-27 | 1 | Continuum solver; MG1 run: 14/15 criteria pass; 1D not stationary at t=100 (registered FAIL); paper's full integrals (row e) not reproduced; nucleation thresholds measured | MG1 FAIL as registered, PI ruling requested | 7c4b194 |
| MG002a | 2026-10-27 | 1 | ChatGPT independently replicated 3D (≤0.0012%) and the 1D failure; seeded balance fixed; wording softened; wavelength-window sensitivity; PI ruling: MG1 PASS by exception | MG1 PASS (exception) | ba7176e |
| MG003 | 2026-10-27 | 2 | Phase 2 protocol and MG2/MG3 proposals; lattice references; narrow single-soliton window found (published b=29 is 0.02% above its lower edge in 1D) | MG2, MG3 proposed | dd9e9be |
| MG003a | 2026-10-27 | 2 | ChatGPT review: 7 items accepted; integrator replaced (deterministic limit now exactly Euler of lattice ODE; soliton exact fixed point); bounded tau selection; unbiased bias estimator; full MG3 design; protocol v1.1 | MG2, MG3 re-proposed | 6ecd0bd |
| MG003b | 2026-10-27 | 2 | ChatGPT review r2: Euler limit confirmed; inference layer rebuilt and calibrated (12 runs per cell; three-way decisions; alpha bootstrap bug found and fixed; C2 decision calibrated); protocol v1.2 | MG2, MG3 re-proposed | d5b381c |
| MG003c | 2026-10-27 | 2 | ChatGPT review r3: confidence claims made conditional on stated assumptions A1–A2, checked per cell, per gate (pooled audit) and by a real-run pilot (V6); complete procedure calibrated with held-out A1-violating families; Part B record corrected; Welch two-sample with deterministic normalization; protocol v1.3 | MG2, MG3 re-proposed | 232ee36 |
| MG003d | 2026-10-27 | 2 | ChatGPT closed spec review (r4, no errors); PI sign-off; protocol v1.3 FROZEN; MG2/MG3 PRE-REGISTERED; configs written | MG2, MG3 pre-registered | 0f5615e |
| MG004 | 2026-10-27 | 2 | Pipeline + SSA built; V3 PASS, V1 PASS; V2 halted (3/24 eligible at b=29, Omega=300): real noise-induced annihilation near window edge + detection-threshold flaw; Amendment 1 proposed | validations paused | e5bed86 |
| MG004a | 2026-10-27 | 2 | Amendment 1 approved (PI) with ChatGPT corrections; protocol v1.3.1; threshold from smoothed reference in code; V2 at b=29.08; C3 data rule; aborted-V2 log restored | v1.3.1 frozen | e8f6877 |
| MG004b | 2026-10-28 | 2 | Validations under v1.3.1: V1, V3, V6 PASS; V2 FAIL as registered; V4 INCONCLUSIVE (blocking). Root cause: run-level profile noise far above tolerances (design power failure). Amendment 2 proposed | validations blocked | e4cb48e |
| MG004c | 2026-10-28 | 2 | ChatGPT review of Amendment 2: G power understated; Amendment 2 revised (decision Omega 1e8, bounded allocation, signed variance rule, conditional alpha wording, accurate C3/C4 expectations, V7 with PI checkpoint) | awaiting PI ruling | (this commit) |

## MG000 — Project scaffold
- Created charter, `spec/rules_v1.0.md`, `GATES.md`, folder layout.
- Next (MG001): pin LaViolette's published Model G equations, reverse-step variants and parameter sets from primary sources into `spec/model_g_published.md`.


## MG001 — Published Model G pinned
- Source: Pulver & LaViolette, Int. J. Gen. Syst. 42 (2013) 519–541. Transcribed in two independent passes → `spec/model_g_published.md`.
- Reference set (eq. 17): d_x=1, d_y=12, a=14, b=29, g=1/10, p=q=1, s=u=w=0. Only reverse step used: G ⇌ X.
- Independent linear check (`src/linear_stability.py`, `results/MG001/linear_stability.txt`): steady state G₀=X₀=15.556, Y₀=1.864; homogeneous state linearly stable (max Re σ = −0.195); Turing threshold b_c = 30.044, so b = 29 is 3.5% subcritical; least-stable λ = 3.149 vs published nonlinear λ₀ = 3.08 (2.3% apart).
- Findings affecting the plan:
  1. Rules v1.0 lacked the G ⇌ X back-reaction → `spec/rules_v1.1.md` (v1.0 superseded).
  2. Full-space φ_G integral is slightly positive in standard Model G; the G-well is a core feature. MG6's "G-deficit" must be defined on the core or use the 2Dδ variant (PI decision at Phase 5).
  3. X core is ~94% depleted in 3D → largest relative noise at the core; relevant to Ω choice in MG2/MG3.
  4. Paper gives no grid/step/method → MG1 includes an explicit grid-convergence criterion.
- MG1 pass/fail criteria proposed in `GATES.md`, awaiting PI sign-off.
- Next (MG002): implement the continuum solver and run MG1 after sign-off.

## MG001a — Review round 1 (ChatGPT) on MG001
- Review: `review/MG001_chatgpt.md`. Response: `review/MG001_response.md`. All items accepted; one (item 4) extended.
- New result: steady-state integral identity. With conservative diffusion, a stationary localized soliton has ∫φ_G = ∫φ_X = 0 exactly in infinite space. Table 1's full integrals require boundary fluxes larger than the integrals themselves (`results/MG001/integral_identity.txt`). The 2Dδ "positive mass" exists only because δ∇² is non-conservative.
- Consequence for Phase 5: graph hopping is conservative, so standard Model G particles carry no net G monopole. A 1/r far field is not expected without a modification. This is flagged as a PI decision.
- MG1 proposal amended with an integral-identity diagnostic.

## MG001b — Review round 2; MG1 pre-registered
- ChatGPT confirmed the integral identity (determinant q − gp, nonzero for the reference set) and asked for three conclusions to be narrowed. All accepted using its wording: required fluxes are inferred, not measured; δ∇² permits the 2Dδ result but is not shown to be its sole cause; the Phase 5 implication is a conditional constraint on coupling models.
- The MG1 diagnostic now checks the time-dependent integral balances with a full error budget.
- ChatGPT: no further mathematical review round needed. PI signed off on the MG1 criteria → PRE-REGISTERED.
- Next (MG002): continuum solver, convergence study, MG1 runs.

## MG002 — Continuum solver and MG1
- Solver `src/mg_continuum.py`: conservative finite-volume radial method of lines (n = 1, 2, 3), BDF with analytic sparse Jacobian. The discrete integral balances hold to round-off by construction.
- MG1 (results/MG002/MG1_evaluation.md, profiles.png): λ₀, amplitudes, core radii, RMS radii and core integrals reproduce Table 1 to ≤ 0.5% (one 1.4%) in 1D/2D/3D and are grid- and tolerance-converged. Registered verdict FAIL on one cell: the 1D soliton is still relaxing at t = 100 (rate ≈ 0.047), and the paper's 1D column matches that snapshot.
- Paper's full integrals of φ_G, φ_X (row e) are not reproduced; ours are ~0 as the identity requires, and measured boundary fluxes are too small to explain the paper's values. Likely a numerical artifact in the paper.
- Nucleation threshold (seed scale): 1D 0.975, 2D 0.572, 3D 0.486.
- Review request for ChatGPT: `src/mg_continuum.py`, `results/MG002/MG1_evaluation.md`, the MG1 Result block in `GATES.md`. Independent replication of at least the 3D amplitudes and λ₀ is requested.
- Next (pending PI ruling on MG1): Phase 2 spec — stochastic token rules on a lattice (MG2/MG3 pre-registration).

## MG002a — Review of MG002; PI ruling on MG1
- ChatGPT independent replication (different discretization and integrator) matches 3D amplitudes and λ₀ within 0.0012% and confirms the 1D stationarity failure (3.49e-4).
- Accepted: seed term added to the summed integral balance (`analyze`); "pre-stationary snapshot" and "tail quadrature" inferences softened; wavelength-window sensitivity reported (0.2–0.4%; far-field windows give ≈ 3.085).
- PI ruling: MG1 PASS with documented exception (1D criterion 4). The original FAIL record is preserved in GATES.md.
- Next (MG003): Phase 2 spec — stochastic token rules on a lattice; MG2/MG3 pre-registration for PI sign-off.

## MG003 — Phase 2 protocol; MG2/MG3 proposed
- `spec/phase2_protocol_v1.0.md`: 1D periodic ring (h = 0.2, L = 100), rules v1.1, binomial/Poisson tau-leaping (`src/lattice.py`, numba), validations V1 (vs exact SSA) and V2 (τ-halving), measurement definitions, compute budget.
- Deterministic lattice references (`results/MG003/`): soliton offset vs continuum 3% (h = 0.1) and 8% (h = 0.2); seed thresholds; lattice b_c ≈ 30.046; LNA variances. X fluctuations are strongly amplified (Ω·Var ≈ 113 vs Poisson 15.6 at h = 0.2).
- **Finding:** the stationary single soliton exists only in a narrow window of b. Lattice h = 0.2: [28.929, 29.397]. Below the window it decays; above it, it self-replicates and fills the domain. Continuum lower edge b_s = 28.993 (1D), 27.343 (2D), 26.972 (3D), so the published b = 29 is 0.02% above the edge in 1D but well inside the window in 2D/3D. This accounts for MG1's slow 1D settling and the 0.975 seed threshold.
- Only stochastic run so far: a throughput benchmark on the homogeneous state. No gate outcome has been observed.
- Review request for ChatGPT: `spec/phase2_protocol_v1.0.md`, `src/lattice.py` (kernel and LNA), the MG2/MG3 blocks in `GATES.md`.
- Next: PI sign-off → V1/V2 → MG2 → MG3.

## MG003a — Review of MG003; Phase 2 protocol v1.1
- Review item 1 reproduced exactly (`results/MG003a/audit_split_scheme.txt`): the v1.0 tau-leap's Ω → ∞ limit loses the soliton (E = 108/103/107/103%). New kernel `tau_leap_euler_1d`: unsplit multinomial exits with probability rate × τ. Its mean map is exactly forward Euler of the lattice ODE, so stationary solitons are exact fixed points (E = 1e-14 after t = 100).
- Bounded τ selection (`results/MG003a/tau_selection*.txt`): h = 0.2 → τ = 6.667e-4; h = 0.1 → τ = 1.667e-4. The scheme's exact discrete-time LNA is within 0.9% of the continuous LNA. A clip rule halves τ at small Ω (Ω = 3 pilot: clip fraction 1.4e-3).
- Items 2–7 addressed in `spec/phase2_protocol_v1.1.md` and the MG2/MG3 blocks of `GATES.md`: U-statistic bias estimator with jackknife SEs and bootstrap CI on α; Fourier alignment plus validation V4; tracking/survival rules; first-passage model (burn-in, censoring, Garwood CIs, stationary-hazard test, rare-event τ check V5); deterministic (committor-style) confirmation of nucleation; fully specified MG3 adaptive design with fallbacks and a hard budget; SeedSequence stream keys.
- Next: ChatGPT re-review of v1.1 → PI sign-off → V1–V4 → MG2 → MG3.

## MG003b — Review round 2; calibrated inference layer; protocol v1.2
- ChatGPT confirmed the new kernel's Ω → ∞ limit is exactly forward Euler (when exits are valid and caps inactive) and accepted the U-statistic. It rejected the jackknife-normal upper bound, with a counterexample.
- New inference layer (`src/stats_mg2.py`): t-based triangle bounds and three-way PASS/FAIL/INCONCLUSIVE decisions, plus a two-sample form. Calibration (`results/MG003b/`) shows no procedure is reliable at 4–6 runs (worst coverage 0.66 / 0.80), because identical-run samples carry no spread information. At **12 runs**: coverage ≥ 0.945, false-PASS ≤ 0.034, false-FAIL ≤ 0.035. The two-sample bound is also calibrated (0.25% and 10% margins).
- Calibration exposed a bug in my α bootstrap (duplicate runs paired with themselves inside the U-statistic) → distinct-pair bootstrap, pre-selected cells, widening w = 2.0. C2 decision calibrated: false-FAIL ≤ 0.020, false-PASS ≤ 0.003.
- V2 is now an equivalence test with predefined margins. MG2 cells → 12 runs; MG2 budget 8 core-hours.
- Next: ChatGPT re-review of protocol v1.2 §5.3 and the calibration scripts → PI sign-off → V1–V4 → MG2 → MG3.

## MG003c — Review round 3; assumption-conditional inference; protocol v1.3
- Accepted ChatGPT's point that no fixed run count gives a distribution-free guarantee. Its rare-mixture example falsely passes with probability 0.9¹² = 28%. Confidence claims are now explicitly conditional on A1 (i.i.d., approximately Gaussian run-level profiles, with no rare discrete modes beyond tracked events) and A2 (independent streams).
- A1 is checked three ways: per-cell diagnostics (block Shapiro, block robust z, run-level studentized leave-one-out with Bonferroni; Gaussian false alarm 1.5%); a pre-gate pilot on real runs (V6); and a gate-level pooled audit across all runs (power against rare modes that single cells lack).
- Complete single-cell procedure calibrated (`results/MG003c/procedure_calibration.txt`), including survival top-up, minimum 12 eligible runs, diagnostics and one extension. A1-satisfying families: false-PASS/FAIL ≤ 0.002, power ≥ 0.985. Held-out A1-violating run mixtures: false-PASS up to 0.70. This residual risk is reported, not hidden.
- Two-sample: Welch bound with deterministic normalization; unequal noise and M calibrated (coverage ≥ 0.990, false-equivalence ≤ 0.001).
- Part B record corrected: undefined intervals reported separately; w = 2.0 labelled as calibrated for the C2 decision rule.
- Next: ChatGPT re-review → PI sign-off (or a PI ruling on any remaining disagreement) → V1–V4, V6 → MG2 → MG3.

## MG003d — Spec frozen; MG2/MG3 pre-registered
- ChatGPT round 4: no remaining blocking error. It noted that the residual A1 risk is a PI decision and that 0.97¹⁶⁸ is the probability of a shared mode being absent from the data, not the audit's failure rate. That wording was clarified at freeze.
- PI (Richard Benear) signed off; protocol v1.3 frozen; MG2 and MG3 pre-registered; `configs/` written.
- Next (MG004): implement the measurement pipeline (tracking, alignment, blocks) and exact SSA; run V1–V4 and V6.

## MG004 — Validations (in progress); Amendment 1 proposed
- Built: `src/measure.py` (§5 pipeline), `src/runner.py` (streams, references, cell runner), `src/ssa.py` (exact SSA), `src/validations.py`.
- V3 clip pilot PASS: τ/2 for Ω ≤ 5, τ for Ω ≥ 7; no Ω excluded.
- V1 exact-SSA comparison PASS (14/14 statistics). At Ω = 30 the X-variance tolerance is wide (±45%) because batch noise is large there, so V1 is a weak test of that statistic.
- V2 halted: at its registered cell (b = 29, Ω = 300) only 3/24 runs were eligible. Diagnosis: real noise-induced annihilation of the near-edge b = 29 particle (it holds at Ω = 1e5 and at b = 29.08). Separately, the detection threshold was specified against the unsmoothed amplitude but applied to the smoothed field (81% of depth instead of 50%).
- `spec/phase2_amendment_1.md` proposed: smoothed-depth threshold; V2 moved to b = 29.08. Awaiting PI ruling.

## MG004a — Amendment 1 approved; protocol v1.3.1
- PI approved Amendment 1. ChatGPT's corrections are incorporated: threshold computed in code from each reference with the detector's own smoothing; V2 claims limited to its operating point, with the cap and stopping rule kept; "3 eligible of 24" not presented as a survival rate; mechanism wording softened; explicit C3 data rule for missing data; disclosure that changing the detector can change outcomes; aborted-V2 log restored (it had been excluded by .gitignore).
- Code: `runner.reference` returns A_s and the threshold; `measure.analyze_run` far-field variance includes particle-free snapshots; validations save every attempted run's record.
- Next: V2, V4, V6 under v1.3.1. V1 and V3 are unaffected (no detection involved).

## MG004b — Validations under v1.3.1; Amendment 2 proposed
- V6 PASS: A1 holds on real survival-conditioned runs at b = 29.08, Ω = 300 / 3000 (all fields, both cells; pooled audit clean).
- V2 FAIL as registered. The profile difference is noise-dominated (bounds 6.6–37.8% vs a 0.25% margin), and the point-estimate FAIL rule is inflated by noise. Variance equivalence passed.
- V4 INCONCLUSIVE (blocking): alignment-bias bounds at Ω = 1e4 are 7.4% (G) and about 1% (X, Y) vs a 0.3% tolerance; half-widths scale as Ω^(−1/2). At Ω = 100 no synthetic run was eligible (noise minima below the threshold).
- Root cause: actual run-level noise of aligned profiles was never measured during calibration. As registered, MG2 would predictably be INCONCLUSIVE.
- Amendment 2 proposed: decision Ω raised to 1e6 (kernel cost flat in Ω); V2 FAIL rule made bound-based and profile test moved to Ω = 1e6; new descriptive detector check V7 before MG3.
- b = 29.08, Ω = 300 survival to t = 100 (eligibility): V2 22/48 + 22/48; V6 12/18.

## MG004c — Amendment 2, revision 2
- ChatGPT showed that Ω = 1e6 leaves G unresolvable for V4 (≈ 0.74% vs 0.3%) and V2 (half-width ≈ 0.61% vs 0.25%). Revision 2 puts the decision point at Ω = 1e8 (projected G upper bounds ≈ 0.07–0.09%), with a bounded allocation and no further escalation. Kernel checked at 1e7/1e8 (cost flat; int64 safe).
- Signed variance-equivalence rule; τ/α claim made conditional on B(τ) ≠ 0 and asymptotic; C3/C4 expectations stated accurately (C4 likely FAIL at b_mid for Ω ≤ 300 based on validation survival data; C4 deliberately not changed); V7 classifies persistent detections as confirmed / excursion / unresolved, with a PI checkpoint before MG3.

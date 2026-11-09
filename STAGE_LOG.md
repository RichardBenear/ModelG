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
| MG004c | 2026-10-28 | 2 | ChatGPT review of Amendment 2: G power understated; Amendment 2 revised (decision Omega 1e8, bounded allocation, signed variance rule, conditional alpha wording, accurate C3/C4 expectations, V7 with PI checkpoint) | awaiting PI ruling | a34f15e |
| MG004d | 2026-10-28 | 2 | Amendment 2 (r2) approved by PI; protocol v1.3.2 frozen; V2b/V4b at Omega=1e8 and V7 implemented; validations rerun | v1.3.2 frozen | 6011830 |
| MG004e | 2026-10-28 | 2 | V2, V4 PASS (Omega=1e8); V7 done: three noise regimes at b=29 (patterned Omega<=50, transitional ~100, quiet 300); confirmation criterion too strict; Amendment 3 proposed (PI checkpoint before MG3) | PI checkpoint | b5c944a |
| MG004f | 2026-10-28 | — | Glossary and assumptions register added (PI request): Omega definitions and values, token = etheron, notation clash with LaViolette's sink Omega, assumptions M1–M10, N1–N5, S1–S4, I1–I4 | — | 8b2f742 |

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

## MG004d — Amendment 2 approved; protocol v1.3.2
- PI approved Amendment 2 revision 2. Protocol v1.3.2: decision point Ω = 1e8 for MG2-C1, C6, V4 and the V2 profile test; bounded allocation; bound-based and signed V2 rules; new descriptive V7 with a PI checkpoint before MG3. The MG2 ladder adds 1e5, 1e6 and 1e8.
- Code: `validations.py` v2b, v4b, v7 (the v1.3.1 functions are kept for the record). Stream keys are distinct from all earlier validations.
- V2 (v1.3.2) **PASS** at b = 29.08, Ω = 1e8 (12/12 eligible each arm): profile differences ≤ 0.027% point, ≤ 0.088% upper (G); variance equivalence PASS for G, X, Y.
- V4 (v1.3.2) **PASS** at Ω = 1e8: alignment-bias upper bounds 0.137% (G), 0.154% (X), 0.143% (Y) vs 0.3%. Observation: the X/Y upper bounds plateau (0.20% at 1e6 → 0.15% at 1e8) instead of scaling as Ω^(−1/2). This suggests an alignment floor of about 0.1% (sub-node position estimate); it is far below MG2's 1% tolerance.
- Validation status under v1.3.2: V1, V2, V3, V4, V6 PASS. V7 (descriptive) next, then the PI checkpoint before MG3. MG2 may proceed after V7 is reviewed or in parallel (V7 gates MG3 only).

## MG004e — V7 and the PI checkpoint
- V7 (descriptive): at Ω ≤ 50 there are 8–25 persistent minima per snapshot (noise-sustained Turing-like patterning; cf. Biancalani, Fanelli & Di Patti, arXiv:0910.4984). Ω = 100 is transitional; Ω = 300 is quiet. Seed-free continuations from Ω ≤ 50 fields reach soliton depth but are still evolving at t = 300, so the §5.4 stationarity requirement classifies them as "unresolved".
- Amendment 3 proposed: basin-membership confirmation rule (affects MG3 and MG2-C5); MG3 grid focused on Ω 50–300 with a pre-declared regime screen; exploratory E1 (stochastic Turing patterns vs LNA).

## MG004g — Amendment 3 finalized and approved; protocol v1.3.3
- ChatGPT review of Amendment 3 and the glossary accepted in full: corrected the large-Ω statement; regime wording ("crowded / transitional / quiet detector regimes"; stochastic Turing patterning as a candidate interpretation pending E1); confirmation renamed and made field-level with identity linking; crowding exclusion applied to all candidates without re-entry; clip pilots for the new grid; glossary limits separated (Euler map vs lattice ODE vs PDE); noise formula qualified (Ω^(−1/2), Poisson only a benchmark).
- Confirmation controls against t = 2000 deterministic truth: **15/15 agree**, including slowly evolving multi-soliton states, near-threshold seeds (0.999 decays; 1.01 persists) and near-edge decay. Seeds at separation 6.6 spawn a third particle, matching Pulver & LaViolette's 6.628 observation.
- Clip pilots for the new grid values: no clipping; τ standard.
- Erratum: stale V2 sentence in v1.3.2 struck (the code had applied the amended rule).
- PI approval (conditional on considering the review) → protocol v1.3.3 frozen. Next: MG2.

## MG005 — MG2 gate run (protocol v1.3.3)
- Driver `src/run_mg2.py`: resumable cell-level checkpoints (raw cells in `results/MG2/cells/`, committed after a container restart, 4.3 MB); budget tracked, raised from 8 to 12 core-hours by PI decision; decisive cells first (C1 at Ω = 1e8, C6, ladder, C5, N1; low Ω deferred as D3).
- `measure.analyze_run` now also reports the far-field snapshot count (C3 data rule) and per-snapshot particle positions (exploratory D_p).

- **Result: MG2 FAIL as registered** (C3 FAIL at b = 29, Ω = 100; C4 FAIL from noise-induced replication at b_mid, Ω ≤ 300). C1, C5, C6, N1, N2 PASS; C2 INCONCLUSIVE (bias unresolved at 12 runs); pooled audit ok. 11.6 core-hours. Details in GATES.md (MG2 result) and `results/MG2/MG2_evaluation.md`.
- Post-gate diagnostic `src/diag_mg2_replication.py`: replication events are physical (two persistent particles under deterministic continuation).
- Awaiting ChatGPT review (`review/MG005_review_request.md`) and PI decision on next steps.

## MG005a — MG2 review; corrections before MG3; Amendment 4 proposed
- ChatGPT review (`review/MG005_chatgpt.md`): independent re-scoring confirms every MG2 decision. Registered FAIL retained; option (a) recommended (MG3 after corrections, by PI decision). Response: `review/MG005_response.md`.
- Record corrections: C4 is 12 replication events in 11 of 24 runs (not 12 of 24); C3 causal attribution withdrawn; b_mid design statement softened.
- Tracker fix (`measure.track`): the first homogeneous birth was omitted. Regression tests in `tests/test_track.py`. Earlier records unaffected.
- Audit defect: the registered pooled audit (and V6's) cannot flag data of these shapes. A corrected within-cell audit flags no MG2 run (descriptive); calibration in `results/MG005/`. It adds no power over per-cell diagnostics, and the residual risk is stated.
- `spec/phase2_amendment_4.md` PROPOSED: tracker fix, audit replacement, MG3 Stage-A terminal outcomes. Awaiting ChatGPT review and PI approval before MG3.
- ChatGPT reviewed Amendment 4 (`review/MG005_amendment4_chatgpt.md`): approve with two wording qualifications (incorporated). **PI approved Amendment 4 and authorized MG3 (option a), 2026-10-28.** Applied as `spec/phase2_protocol_v1.3.4.md`; `configs/MG3_main.json` updated to the A3 grid (it still listed the v1.3 grid).

## MG006 — MG3 gate run (protocol v1.3.4)
- PI decisions before the run: MG3 cap 18 core-hours (from 10), and each run stops 50 time units after its first confirmed event.
- `src/mg3.py` (segmented simulation with online detection, tracking, confirmation; first-passage statistics; Stage-B selection), `src/run_mg3.py` (resumable driver), `tests/test_mg3.py` (14 tests; "15" in the first record was a miscount). Pre-run implementation specification recorded in GATES.md before any MG3 run.

- Driver change during the run (before N1 restarted): the budget is checked only at cell start, so a forced 2-hour restart does not re-check a cell already started.
- **Result: MG3 FAIL as registered (C2(a))**; C1, C3, C4, V5 PASS; N1 had no confirmed event; N2 PASS. Rates 4.05e−6, 4.63e−5, 1.58e−4 at Ω = 150, 120, 100; S = 0.0726 per unit Ω. 19.03 core-hours (the overrun is inside the last cell, N1). Details in GATES.md (MG3 result) and `results/MG3/MG3_evaluation.md`. Review request: `review/MG006_review_request.md`.

## MG006a — MG3 review; follow-up drafted
- ChatGPT review (summary in `review/MG006_chatgpt.md`; full review on the PI's machine): MG3's FAIL preserved. Numbers reproduce (168 records; 36 first-event fields re-confirmed independently). Qualifications adopted in GATES.md: C2(a) scope, C3 estimated variances, weak stationarity, C4's three short waits, V5 not resolving time-step independence, interval exactness unjustified for this design, Stage-A trigger ≤ 2.2%.
- Response: `review/MG006_response.md`. Follow-up `spec/MG3b_followup_proposal.md` DRAFT: burn-in/time-step check first, rates with ≥ 20 events per cell, simulation-calibrated intervals, and an optional full-trajectory confirmation study. Awaiting review and PI approval.
- Full review committed by the PI (`review/MG006_replication/`). Record corrected per its §§2–8: wording on excursions (18/54, a minority), attribution (`near`), profile eligibility ≠ survival, C3 and V5 wording, interval exactness, the C4 contrast, the Stage-A power figure (2.2%), the test count (14), and the budget recorded as an explicit exception with a ledger (`results/MG3/budget_ledger.txt`; Ω = 80 checkpoints committed). The follow-up proposal is to be revised to the review's §8 before registration.

## MG006b — MG3b follow-up proposal revised after review
- ChatGPT reviewed draft 1 of the follow-up (`review/MG006b_chatgpt.md`, verbatim). Response: `review/MG006b_response.md`. Draft 1 preserved as `spec/MG3b_followup_proposal_v1.md` (superseded).
- `spec/MG3b_followup_proposal.md` DRAFT 2: Part 3 calibration first, with startup probability varied independently of the hazard and the methods, grid, Monte Carlo precision and selection rule frozen; Part 1 decides τ by equivalence (margin [0.5, 2]) on the post-burn-in rate ratio at Ω = 120, 64 runs per arm, outcomes EQUIVALENT / DIFFERENT / INCONCLUSIVE with τ/4 escalation closed; startup comparison descriptive only; Part 2 capped cells kept in a joint Poisson-GLM likelihood, rates conditional on τ; Part 4 removed (field-level confirmation is not detection-level attribution); budget enforcement by run caps, per-batch ceiling checks and a summed ledger. Projected 24 / 36 / 63 core-hours by path.
- Not registered. Awaiting PI discussion, ChatGPT review of draft 2, and PI approval.

- ChatGPT reviewed draft 2 (`review/MG006b_chatgpt_r2.md`): design supported; three calibration corrections (capped-cell configuration preserving the activated model; safe INCONCLUSIVE fallback for Part 1 and a rule if no stationary split qualifies; full decision-tree calibration with coverage by production branch) and two smaller ones (tighter-margin power was overstated in draft 2; benchmark calibration runtime). All accepted: draft 3 (`spec/MG3b_followup_proposal.md`; draft 2 preserved as `_v2.md`), response `review/MG006b_response_r2.md`. Stationary-check power now simulated from the actual procedure. Part 3 ceiling 8 core-hours with a registered R = 1000 fallback for bootstrap candidates. Not registered; awaiting PI approval.

## MG007 — MG3b registered
- **PI approved draft 3 of the MG3b follow-up, 2026-10-29** ("I'm satisfied. go ahead and register it"), after two ChatGPT review rounds. Registered in `GATES.md` (MG3b) with a pre-run implementation specification (streams, run records, ledger columns, τ/4 timing runs, driver order). Design frozen as `spec/MG3b_followup_proposal.md` at commit `202c00f`. Estimation study, no pass/fail; MG3's FAIL unchanged. No MG3b run of any kind has been made.
- Next: Part 3 synthetic calibration (`results/MG3b/calibration/`, with `decisions.json`), sent for review before the driver and any lattice run.

## MG007a — MG3b Part 3 calibration (synthetic)
- `src/calib_mg3b.py` and `tests/test_calib_mg3b.py` (10 pass). The benchmark projected 0.60 core-hours, so R = 4000 throughout. The main run took 0.62 core-hours on 4 cloud workers; the total with benchmarks and diagnostics is 0.74.
- Decisions (`results/MG3b/calibration/decisions.json`; report `results/MG3b/calibration/README.md`):
  - per-cell rate CI: Garwood;
  - Part 1 ratio CI: Clopper–Pearson; 64 runs per arm;
  - slope CI: **no candidate eligible**. The cause is branch-conditional tail asymmetry, with an upper miss up to 0.047 in a few selected branches. Mean coverage is 0.952, and Monte Carlo noise contributes. The registered fallback is the profile likelihood, labelled approximate (minimum 0.937);
  - stationary split: exposure median. The fixed split is slightly anticonservative. Power is low, so the check will mostly be uninformative.
- Deviations: the implementation choices were recorded in GATES.md after the run rather than before; the GLM fit was fixed before the main run; the bootstrap branch threshold was corrected after the run, with the selection unchanged.
- Next: ChatGPT review of the calibration, then the MG3b driver.

## MG007b — Calibration review and addendum 1
- ChatGPT reviewed the Part 3 calibration (`review/MG007_chatgpt.md`). The method decisions are kept. Two post-selection gaps were found: the reused Ω = 120 rate's coverage within branches, and Part 1b's selected comparison.
- Addendum 1 was registered in GATES.md (commit `bccf62a`) before it was run. The τ/4 arm is fixed at 64 runs. Results:
  - the Ω = 120 rate is labelled approximate in every branch (observed minima 0.943 / 0.875 / 0.841);
  - the Part 1b whole-study error is ≤ 0.040 and its conditional error is ≤ 0.161; both are reported.
- These reproduce ChatGPT's figures. The 0.937 slope figure is reworded as an observed minimum over qualifying cells. Response: `review/MG007_response.md`.

## MG007c — MG3b driver
- `src/mg3b.py` and `src/run_mg3b.py` reuse MG3's simulation, detection, tracking and confirmation. They add:
  - storage of every confirmation field;
  - the run ledger and the per-batch ceiling check;
  - the stop20 rule and the Part 1 → 1b → production-branch logic;
  - τ/4 timing runs, clip pilots for Ω = 130 and 110, and the `decisions.json` gate.
- `tests/test_mg3b.py`: 7 pass.
- A plumbing smoke run of `main()` was made in a scratch folder on the cloud machine, with runs cut to t = 60, batches of 2 and a cap of 4. Nothing was committed. It exercised the pilots, Part 1, the branch, Part 2, the ledger and resume (a re-run did nothing and added no ledger rows).
- The driver's implementation choices are recorded in GATES.md before any production run.
- Next: the production run on the PI's machine (`MODELG_WORKERS=9 python src/run_mg3b.py`). Projected 24–63 core-hours, about 3–7 h of wall time.

## MG007d — MG3b production run and evaluation
- Ran on the PI's machine (RB_LibraryPC, 9 workers) through Remote Control, 2026-10-29 local time, about 4.0 h wall,
  with no errors and no budget stop. Run records committed from the PC (258 records, about 6 MB).
- PI decisions during the run: τ/4 ceiling back to the unscaled 40 core-hours (before Part 1 was known; not used) and
  no core-hour ceiling for Part 2 at τ (after Part 1 was known; not exercised). Both are in GATES.md.
- Part 1 EQUIVALENT (r = 0.95, 90% CI 0.65–1.40), so production at τ. Part 2: Ω = 130, 110, 100 reached 25, 28 and 21
  events in 48, 48 and 32 runs. S = 0.064 (approximate profile 95% interval 0.047–0.082). The stationary check is
  uninformative in every cell. 26.75 core-hours in all.
- `src/evaluate_mg3b.py` (written after the run, implementing the registered analyses; open details listed in the
  record). Independent recomputation in `results/MG3b/independent/`. Record: `results/MG3b/MG3b_evaluation.md`.
- ChatGPT review (`review/MG008_chatgpt.md`): supports accepting MG3b as an estimation result and reproduced every
  number, including all 302 confirmation classifications. Two label clarifications were applied
  (`review/MG008_response.md`): the three-cell sensitivity interval and the landmark and early/late intervals are
  nominal, not calibrated. No number changed.
- PI sign-off of the result: 2026-10-30 12:50 UTC ("I sign off on MG3b"). MG3b is closed.

## Deferred work (PI decisions)
Items deliberately postponed so the project reaches its main questions first. They can be revisited later.
- **D1 (2026-10-28):** Noisy existence-window edge study: survival S(δ, Ω, T) with δ = b − b_s(h); noisy edge at S = 1/2; test of the saddle-node scaling δ ∝ Ω^(−2/3). Deferred by PI.
- **D2 (2026-10-28):** Repeating key MG3 cells at h = 0.1, and reporting b relative to each lattice's window. Deferred by PI. Until then, near-edge results at b = 29 (e.g. survival) are stated as h = 0.2-specific (assumption M11).
- **D3 (2026-10-28):** MG2 ladder cells at Ω = 10 and 30 (both b). Deferred by PI when the MG2 cap was raised from 8 to 12 core-hours. Not used by any MG2 criterion.

## MG009 — MG4 design, draft 1
- `spec/MG4_proposal.md` (draft 1, not registered). MG4 tests particles on frozen irregular 2D graphs: a triangular
  comparator, Delaunay (primary) and random geometric graphs, and a random regular graph as a negative control.
  Part 0 calibrates each graph's spacing from its homogenized diffusivity. Part A (deterministic graph ODE) carries
  survival, profile and spectral (Weyl) criteria. Part B (stochastic, Ω·ρ matched, 16 runs) carries survival under noise.
- No MG4 dynamics was run. Only graph statistics (degrees, edge-length moments) were computed to size τ and the cost.
- Wolfram hypergraph snapshots are proposed as a separate later study (MG4b).
- Draft 2 (same day): the PI chose "2D graphs, noise-free only" (13:41 UTC). The stochastic Part B is removed and
  deferred; the primary h returns to 0.2; projected cost about 14 core-hours, all deterministic.
- Next: PI choices on the questions in §12, ChatGPT design review, then registration in GATES.md.

## MG009a — MG4 design review and draft 3
- ChatGPT reviewed draft 2 (`review/MG009_chatgpt.md`). All five items were accepted (`review/MG009_response.md`).
- Draft 3 changes:
  - C1 becomes a fixed benchmark panel with explicit counts. Draft 2's Clopper–Pearson claim was wrong.
  - C3 splits into a spectral count against the exact torus and a coarse-grained operator test on the soliton.
  - All families use one commensurate rectangular torus.
  - A measurement specification is added, with a precedence rule for outcome classes.
  - Stationarity is now reported separately from survival.
- Drafts 1 and 2 are preserved as `spec/MG4_proposal_v1.md` and `_v2.md`. No MG4 dynamics has been run.
- Next: ChatGPT re-review of draft 3, then PI approval and registration.

## MG009b — MG4 review round 2 and draft 4
- ChatGPT reviewed draft 3 (`review/MG009_chatgpt_r2.md`). It supports the design and registration once four definitions
  are added. All four were accepted (`review/MG009_response_r2.md`): C2 componentwise medians, the eigenvalue sign
  convention, an absolute-plus-relative balance tolerance, and gate precedence. The draft also states two interpretation
  limits.
- Draft 3 is preserved as `spec/MG4_proposal_v3.md`. No MG4 dynamics has been run.
- Next: PI approval, then registration in GATES.md.

## MG009c — MG4 registered
- The PI approved draft 4 on 2026-10-30 at 14:19 UTC ("I approve MG4 draft 4."). MG4 is entered in GATES.md, with the
  design frozen at `773255d`. No MG4 run preceded the registration.
- Next: Part 0 code with tests first (`src/graphs.py`, `src/graph_ode.py`, `src/measure_graph.py`, `tests/test_mg4.py`).
  Every implementation choice is recorded in GATES.md before the run it affects.

## MG010 — MG4 Part 0 code
- `src/graphs.py` covers the commensurate torus, the T, D, R and N3 builders, the Laplacian, the homogenized corrector,
  shift-invert spectra and the exact torus count. `src/graph_ode.py` covers the graph ODE, the continuum 2D initial
  condition and the node-sum balances. `src/measure_graph.py` covers smoothing, particles, classes, radial profiles and λ.
- `tests/test_mg4.py`: 9 tests pass. They cover:
  - T's corrector against 1.5a²;
  - T's spectrum against the analytic one;
  - Delaunay E = 3N on the torus;
  - the Jacobian against finite differences;
  - the balances, including at the exact homogeneous state;
  - detection and class precedence on synthetic fields.
- Part 0 implementation choices are recorded in GATES.md (MG4) before any Part 0 run. A small-torus smoke run (code
  validation only) gave 80 s for 4 340 nodes, and the particle was not yet stationary at t = 300.
- Next: run Part 0.
- **Part 0 run (2026-10-30, cloud):** complete, with no stop. The 20 × 20 torus is used. T's window is b ∈ (27.33,
  29.20), and b = 29 is inside it. T is stationary (max|dφ/dt| 6.7e-12) with depth −13.58. Box check L2 1.8e-4.
  - Cost 10.57 core-h against the 8 ceiling (2.04 of it estimated lost work). The driver did not check the ceiling.
    The PI recorded it as an exception and continued (GATES.md).
- **Solver speed-up:** nested-dissection LU ordering, with the same steps and rounding-level differences. T is 4.7×
  faster and R 2.7×. A 10th test was added. Recorded in GATES.md before use.
- **Part A scope (PI, 2026-10-31):** gate runs only (T 4, D 32, R 32) plus the h = 0.3 check, with the ceiling raised to
  80 core-h. Seeded runs are deferred and the N3 dynamics dropped (the LU is nearly dense). Part A runs on the PI's PC.

## MG010 (result) — MG4 Part A evaluated: FAIL
- Part A ran on the PI's PC (2026-10-31, 00:30–12:00 UTC, 9 workers): 80 runs, 91.6 core-h. The gate group overran its 80
  ceiling while runs were in flight; the PI recorded it as an exception.
- Every D and R run multiplied: the particle copies itself until 33–37 fill the 20 × 20 torus. T is single at all 4
  placements. C1 FAILs (0/32) and C3 FAILs (disclosed before the runs). C2 is INCONCLUSIVE (no single runs).
  **MG4: FAIL.**
- Recomputed independently with separate detection code: same classes. The balance-check flags on 11 R runs are float64
  rounding (about 1e-15 in extended precision).
- Next: PI sign-off, then ChatGPT review. A follow-up (D/R existence window in b) would be a new registered study.
- ChatGPT review MG010 (2026-10-31): accepts FAIL. Corrections applied:
  - The evaluator now keeps the balance flags and reads their recorded resolution, so it scores FAIL.
  - All 11 flagged runs are in the balance evidence.
  - The mechanism is worded as a hypothesis, and the scope of the FAIL is stated.
- PI signed off MG4 (FAIL) 2026-10-31 13:52 UTC. MG4 CLOSED.

## MG011 — MG4c design (draft 1)
- `spec/MG4c_proposal.md` draft 1 (2026-10-31), PI request after MG4 closed: does unweighted Delaunay have a single-particle window in b? Plus one cotangent-Laplacian diagnostic at b = 29. To run on the PI's PC. Not registered; no run made. Next: ChatGPT review.
- MG011a: ChatGPT reviewed draft 1 (review/MG011_chatgpt.md) and accepted the study but asked for six fixes. Draft 2 addresses all six (review/MG011_response.md); draft 1 is preserved as _v1.
- MG011b: ChatGPT round 2 on draft 2 (review/MG011_replication/review_r2.md). Draft 3 adds a fixed T-reference amplitude eligibility rule and three wording/branch fixes (review/MG011_response_r2.md).
- MG011c: ChatGPT round 3 supports draft 3 unchanged. The PI approved it (16:25 UTC), and MG4c was REGISTERED in GATES.md (draft 3 frozen at a946550). Next: driver, operator code and tests; pre-run implementation specification in GATES.md.

## MG012 — MG4c implementation (before any run)
- Code: `src/run_mg4c.py`, `src/evaluate_mg4c.py`, `graphs.cotan`, operator and area-weight arguments in `graph_ode`. Tests: `tests/test_mg4.py` (11, incl. cotangent construction) and `tests/test_mg4c.py` (12: decision rules on constructed tables, eligibility, mpmath recheck, synthetic end-to-end passes). All pass on the cloud container (Python 3.13, mpmath 1.4.1).
- The synthetic end-to-end test found a driver bug before any run: a finished but uncollected run was counted twice against the run cap, which would have refused the 102nd run. Fixed.
- Pre-run implementation specification appended to GATES.md MG4c (MG012). It discloses a cloud trial of the Part 3 T control that was stopped before finishing, with no output read.
- Next: PI merges; run on the PI's PC with MODELG_WORKERS=9.

## MG013 — MG4c result
- Run on the PI's PC, 2026-10-31 (started 10:42 local), from main at 494403a. All 102 runs, 11.33 core-h, no balance flags. Raw results were committed from the PC at edd989d (branch mg4c-results).
- Evaluator: W1 + O1. The four D runs are all single from b 27.2656 to 28.7266. Edge brackets are [28.7266, 28.7305] (upper; D1 p1 multiplies) and [27.2617, 27.2656] (lower; D1 p0 decays, unsettled). T is not an eligible reference below 27.375, so the lower edge is detector-conditional. The cotangent operator keeps D0 p0 single at b = 29; the T control matches to 4.5e-13.
- Independent recompute (`src/recompute_mg4c.py`) agrees with every class and threshold. Next: ChatGPT review, then PI sign-off.
- MG013a: ChatGPT review MG012 reproduces W1 + O1 and both brackets. It raised one deviation, recorded in GATES.md: the T-operator check was scaled (1e-12 × max|L|), not absolute as worded. The measured difference is 3.5e-12; the operators are identical in exact arithmetic, and the T control passed at 4.5e-13. Interpretation limits added (stationarity 63/72; the weighted single is not stationary; the edges are finite-time detector boundaries). Response: review/MG012_response.md. The PI decides whether O1 stands.
- MG013b: PI accepted O1 with the deviation recorded and signed off W1 + O1 (2026-10-31 20:24 UTC). MG4c CLOSED. ChatGPT replication files (review/MG012_replication/) still to be committed from the PI's PC.

## MG014 — MG5 design (draft 1)
- `spec/MG5_proposal.md` draft 1 (2026-10-31), PI request after MG4c closed. Phase 4: a settled single particle (MG4c's
  b = 28.0 final states on D0/D1) on a Delaunay graph rewritten by a field-blind birth–death rule (nodes die at rate ν,
  are born uniformly at rate νρ; the snapshot ensemble stays MG4's D ensemble). Concentrations interpolated at births.
  Part A: survival (PASS/PARTIAL/FAIL) and pinned/mobile at ν = 0.002, 0.02, 0.2 (12 runs, t = 500). Part B: annealing at
  b = 29, frozen vs ν = 0.2 (8 runs). New solver: fixed-step ROS2 (BDF restarts cost about 75 s per epoch).
- Design-stage computation: MG4c particle positions re-read (pinning shifts 0.23–0.76, unchanged t = 200–300); timing
  and accuracy pilot on an unregistered graph only (no outcome read). One pilot script that started from MG4c's
  registered D0 p0 field was stopped unread and its output deleted (disclosed in proposal §9).
- Not registered; no MG5 run made. Next: PI decisions (proposal §10), then ChatGPT review.
- MG014a: PI accepted the three recommendations (2026-10-31 21:14 UTC): interpolating transfer rule, Part B included, no
  weighted operator. Recorded in proposal §10. Next: ChatGPT review of draft 1.
- MG014b: ChatGPT review MG014 of draft 1 (review/MG014_chatgpt.md; full checks on the PI's PC). All 5 items accepted
  (review/MG014_response.md). Draft 2: tick-synchronous rule (births interpolate on the tick-start mesh; one rebuild per
  tick; tick check Δ 0.2 vs 0.1); each panel's own h (the stored states are not equilibria at h = 0.2, confirmed); matched
  frozen controls in Part A; Part B keeps absolute concentrations (φ_Y shift −0.0642857); step and tick replay checks;
  A1 reworded as recovery under rewriting plus transfer; full check records, jump accounting, run statuses. 24 runs,
  about 20 core-h. Draft 1 preserved as spec/MG5_proposal_v1.md.
- MG014c: PI decision on storage (2026-10-31 22:44 UTC): commit only check records, accounting and run records; full
  fields stay on the PI's PC with a committed SHA-256 manifest. Proposal §7 and §11 updated. General PI preference: keep
  the repo as small as reasonable.
- MG014d: ChatGPT review MG014 round 2 (review/MG014_replication/review_r2.md, committed by the PI at c7534ca) supports
  implementation and Part 0 after two corrections. Both are accepted in draft 3 (review/MG014_response_r2.md):
  - Event lists are not committed. A replay specification plus a hash check replaces them; the lists stay on the PI's
    PC with the fields.
  - An unresolved interpolation-audit flag makes a run FAILED-NUMERICAL, with an exact-rational recheck procedure.
  Implementation details are fixed now: boundary timeline, ephemeral nodes, row ordering, strict Part 0 counts, and a
  full-implementation timing restatement. The interpretation wording for MOBILE and A3 is tightened. Draft 2 is
  preserved as spec/MG5_proposal_v2.md. Next: PI approval to register.
- MG014e: PI approved draft 3 (2026-10-31 23:17 UTC, decision card "Approve and register"). MG5 REGISTERED in GATES.md;
  draft 3 frozen at e2a25e3. No MG5 run preceded registration. Next: implementation (tests first) and the pre-run
  implementation specification.

## MG015 — MG5 implementation (before any run)
- Code: `src/evolve.py`, `src/run_mg5.py`, `src/run_mg5_part0.py`, `src/evaluate_mg5.py`, `src/recompute_mg5.py`;
  `tests/test_mg5.py` (11 tests pass). `results/MG5/local/` is gitignored (PC-only event lists and fields).
- Found while testing: ROS2's convergence ratio is 2.5–3.3 during the stiff initial transient from the raw continuum
  profile, and 3.7–3.8 (second order) once it has passed. The test checks order after the transient. The registered
  Part 0 solver check (ROS2 vs BDF, 50 units, ≤ 1e-3) covers the transient.
- Pre-run implementation specification appended to GATES.md MG5. Next: PI merges; Part 0 on the PI's PC
  (`$env:MODELG_WORKERS=9; python src/run_mg5.py part0`), then Parts A and B.
- MG015a (2026-11-01): on the PI's PC (Windows), test 10 of `tests/test_mg5.py` crashed at temporary-directory cleanup,
  after all its assertions had passed. `run_mg5.load_checks` left the `.npz` open, and Windows will not delete an open
  file. Rerun with cleanup errors ignored, tests 10 and 11 pass there (88 s, 2 s). Part 0 was started at 17:52 local
  from main 7ae167a. Part 0 never reloads check files, so it does not use the affected code. Fix:
  `load_checks` and `recompute_mg5` now copy the arrays out and close the file. All 11 tests pass in the cloud. The fix
  does not touch any code Part 0 runs, and it is needed before the evaluator runs on Windows.

## MG016 — MG5 Part 0 (2026-10-31, PI's PC)
- Part 0 ran from main `7ae167a`, 17:52–20:11 local (2.3 h wall). All three checks pass: solver 1.5e-4 (≤ 1e-3); step
  check centre 5.6e-4, depth 0.09%; tick check centre 0.060, depth 1.34% (≤ 0.1 and 2%). Counts were 1 at every check in
  all four runs. No flags. Decision: proceed with dt = 0.05, Δ = 0.2 (no changes). Recomputed independently from the
  committed check files.
- Restated budget: 33.3 core-h (design 20), under the 40 core-h launch cutoff; run cap 24. Recorded in GATES.md before
  Parts A and B. Next: PI merges, then `$env:MODELG_WORKERS=9; python src/run_mg5.py` on the PC.

## MG017 — MG5 first launch interrupted; budget and execution amendment (2026-11-01)
- First launch (PI's PC, 9 workers, 06:23 local) stopped at about 13:55 on a console control signal, cause unknown. 2 of
  24 runs COMPLETE and kept (`A_nu0.02_D0_p0`, `_D0_p1`; 6.23 core-h on the ledger). Nine in-flight runs were lost
  (about 65 core-h off the ledger); they restart from t = 0 with the same seeds. Nobody read their fields or checks.
- Runs took about 3.5 times the Part 0 rate because of memory pressure with 9 workers (RAM nearly full, page file in
  use, about 0.4 core per worker, machine about 70% idle), plus unrestricted BLAS threads.
- Claude's MG016 statement on run order was wrong: the driver ran ν = 0.2 first and the frozen controls last, and the
  40 core-h cutoff would have kept the controls from starting.
- Amendment (PI approved, 20:58 UTC): frozen controls first, 5 workers, one BLAS thread per worker, cutoff 80 core-h,
  run cap 24. No criterion, rule, seed or run changes. Tests pass (11). Records from the interrupted launch are merged
  from branch `mg5-interrupted-1` (commit `9c678e0`).

## MG018 — MG5 results (2026-11-02, DRAFT pending PI and ChatGPT review)
- All 24 runs COMPLETE on the PI's PC (39.96 core-h on the ledger, about 65 more off-ledger from the lost first launch).
  Part A frozen controls valid. **Gate MG5: PASS** (single 4/4 at every rate). Mobility: PINNED, PINNED, MOBILE for
  ν = 0.002, 0.02, 0.2. **Part B: A1** (frozen 4/4 multiplied; ν = 0.2 4/4 single), read as recovery under rewriting
  plus transfer, not annealing.
- Independent recompute agrees on all 24 classes. Details and caveats in GATES.md MG5 "Result (MG018)".
- Correction (Claude): I first said the run finished 0.04 core-h under the 40 core-h cutoff. The cutoff was raised to 80
  by MG017, so there was no margin issue.
- Open: PI sign-off; ChatGPT review of the result; MG6 design (directed response to a G gradient) is not started.

## MG018a — ChatGPT review of the MG5 results (2026-11-02)
- ChatGPT (MG018, `review/MG018_chatgpt.md`, replication in `review/MG018_replication/`) independently reproduced every
  class, track, D_p, the gate (PASS) and Part B (A1). It also reproduced the Part 0 comparisons, the jump accounting and all
  16 event lists. It recommends acceptance, with no blocking error.
- All three nonblocking items are accepted (`review/MG018_response.md`):
  - `recompute_mg5.py` now also checks tracks, D_p and the aggregate outcomes, and exits 1 on a difference. Everything
    agrees.
  - Mobility labels are descriptive.
  - Survival is not invariance of shape.
- Total compute, made explicit: 39.96 ledger + about 65 lost-launch core-h.
- Next: PI sign-off to close MG5; MG6 draft 1 (`spec/MG6_proposal.md`) awaits PI decisions (§9) and ChatGPT review.
- MG5 CLOSED 2026-11-02 13:39 UTC: Gate PASS, Part B A1 (PI sign-off in the project thread).

## MG019a — MG6 design draft 2 after ChatGPT review MG019 (2026-11-02)
- ChatGPT (MG019) supports the study but asked for six changes before registration:
  - a valid coupled birth–death process;
  - a density proxy, not curvature, with two baselines;
  - corrected coupling arithmetic and a quasistatic prediction;
  - a calibrated run-level paired rule in place of the 3σ rule;
  - falloff on annular density;
  - complete numerical and audit rules.
- All six are accepted in draft 2 (`spec/MG6_proposal.md`; response `review/MG019_response.md`). Two errors of mine are
  corrected: F2 clipping is active at κ = 4, and the κ4/κ2 ratio is about 1.7, not 2.
- PI (14:13 UTC) adopted ChatGPT's recommended decisions:
  - F1 with bound e^0.4;
  - κ panel {2, 4};
  - shifted null;
  - M₋ well integral, with the b-scan deferred;
  - DEC rerun;
  - ν = 0.2 only;
  - R = 1.5 and window 100–200;
  - 20 runs, cap 24, provisional 40 core-h cutoff.
- Draft 1 is preserved as `spec/MG6_proposal_v1.md`. Not registered. No MG6 run. Next: ChatGPT reviews draft 2; then PI
  registration.

## MG019b — MG6 design draft 3 after ChatGPT review MG019 round 2 (2026-11-02)
- ChatGPT r2 (`review/MG019_chatgpt_r2.md`) keeps all eight PI decisions and asks for targeted corrections. All accepted
  in draft 3 (`spec/MG6_proposal.md`; response `review/MG019_response_r2.md`):
  - a full-field nodal-maximum bound check at every boundary;
  - a 398 × 402 quadrature that tiles the actual torus, with exact or sub-lattice disc weights and four quadrature checks;
  - the audit scope (offline recompute covers the primary gate; quadrature cross-checked in-run);
  - exact Part 0 formulas (tolerance reference, blocks, ephemeral exclusion, Σp(1 − p) variance, immigration–death
    covariance);
  - falloff stated as expected to be UNRESOLVED (start-state wells cross zero at 1.24–1.39);
  - the conditional working model;
  - the absolute-deficit sign reported separately.
- Draft 2 is preserved as `spec/MG6_proposal_v2.md`. Not registered; no MG6 run. Next: the PI decides registration
  (ChatGPT r2 supports proceeding to implementation and Part 0 after these corrections).
- MG019c: PI approved draft 3 (2026-11-02 14:45 UTC, decision card "Register draft 3"). MG6 REGISTERED in GATES.md;
  draft 3 frozen at `106b3e2`. No MG6 run preceded registration. Next: implementation (tests first), a pre-run
  implementation specification in GATES.md, then Part 0 on the PI's PC.

## MG020 — MG6 implementation (before any Part 0 or registered run), 2026-11-02
- Code:
  - `src/mg6.py`, `src/quad6_check.py`, `src/run_mg6.py`, `src/run_mg6_part0.py`, `src/evaluate_mg6.py` and
    `src/recompute_mg6.py`;
  - `tests/test_mg6.py`, with 10 tests, all passing in the cloud.
- The pre-run implementation specification is appended to GATES.md MG6.
- Disclosure: a smoke test ran on registered start state D0_p0 with an unregistered candidate seed (κ = 4, 20 units).
  Its values are recorded in item 12 of the specification. Nothing was changed in response.
- MG020a (PI decision 15:36 UTC, before Part 0): 8 workers instead of 5; the Part 0 timing now runs 8 concurrent runs.
  Speed setting only (GATES MG6 specification item 13).
- Next: PI merges; then Part 0 on the PI's PC (`$env:MODELG_WORKERS=8; python src/run_mg6.py part0`, launched with no
  console attached).

## MG020b/MG021 — MG6 Part 0 STOP and amendment 1, 2026-11-02
- Part 0 ran on the PI's PC (08:43–11:58 local, 8 workers, code `410696d`). Checks 1–3 passed. The tick check failed
  part (a) at Δ 0.2 vs 0.1 and again at Δ 0.1 vs 0.05; part (b) passed both times. Registered outcome: STOP. Records
  committed (`5d351dc`). The overwritten `check_d0.1_dt0.05.npz` is identical to the run it replaced (same seed and
  configuration; D̄₀ equal to every digit).
- Diagnosis (unregistered): the centroid centre has estimator noise of about 0.1 (quadratic-fit comparison). At Δ 0.1
  the part (a) failure falls within it (0.36 → 0.11); at Δ 0.2 a real separation remains (0.26–0.38).
- PI decision (21:18 UTC, card): amend to Δ 0.1 with part (b) governing. Amendment MG021 drafted
  (`spec/MG6_amendment_1.md`), with `part0-mg021` code to complete checks 5–8 and the 8-worker timing at Δ 0.1;
  provisional cutoff 80 core-h. Test added (11 tests).
- Disclosure: checks 5 and 6 were computed in the cloud from committed files while testing (bound ok, no drift).
- Next: merge; ChatGPT reviews MG021; then `$env:MODELG_WORKERS=8; python src/run_mg6.py part0-mg021` on the PI's PC.
- ChatGPT review MG021 (committed by the PI): supports Δ 0.1 as a disclosed, narrow exception; independently reproduced
  the Part 0 records, the fitted-centre diagnostic and checks 5–6 from records; asked for four guard fixes and narrower
  wording. All accepted (`review/MG021_response.md`); amendment draft 2 (draft 1 kept as `spec/MG6_amendment_1_v1.md`);
  `run_ticks6` now records φ_G maxima including the terminal boundary; 12 tests.
- Next: merge; then `$env:MODELG_WORKERS=8; python src/run_mg6.py part0-mg021` on the PI's PC.

## MG021c/d — MG6 Part 0 completed; budget raised to 110 core-h, 2026-11-02/03
- `part0-mg021` on the PI's PC: checks 5–8 all pass (details in GATES.md); projected 90.73 core-h at 81.7 s/unit with 8
  concurrent runs, above the 80 cutoff, so `proceed = false` as designed. Records committed (`0248e78`); numbers
  re-checked in the cloud from the committed JSON.
- PI decision (00:14 UTC): cutoff 110 core-h (`results/MG6/part0/budget_decision_mg021.json`); driver reads it. 13 tests.
- Next: merge; launch the 20 registered runs on the PI's PC (`$env:MODELG_WORKERS=8; python src/run_mg6.py`).

## MG022 — MG6 registered runs and result, 2026-11-03
- 20/20 COMPLETE on the PI's PC (17:39–05:08 local, 8 workers), 83.09 core-h; records `09b72c2`.
- Gate DEMONSTRATED (CP2 − SH2 −67.1, t −36.2; CP2 − DEC −66.9, t −34.1); feedback SURVIVES (CP2, CP4); absolute
  deficit negative; falloff UNRESOLVED. Evaluator and recompute agree on the PC and in the cloud.
- `recompute_mg6.py --cross-platform` added: the `ell` candidate array differs in the last bits between Windows and Linux
  (numpy exponential uses libm log on a few draws); nothing else differs that matters.
- Next: merge; ChatGPT review; PI sign-off.
- ChatGPT review MG022: strict Windows replay reproduces everything; recommends accepting DEMONSTRATED and SURVIVES, no
  rerun. Three record corrections made (κ = 4 falloff and Dc CP4 − DEC added, both falloffs UNRESOLVED; recompute compares
  stored statistics; rounding and cost labels). Response `review/MG022_response.md`. Awaiting PI sign-off.
- **MG6 CLOSED 2026-11-03 13:24 UTC** (PI sign-off in the project thread): DEMONSTRATED, SURVIVES, falloff UNRESOLVED.

## MG023 — MG7 design (mass scaling b-scan), draft 1, 2026-11-03
- PI asked for MG7 (b-scan, falloff, shuffled null; ν = 0.02 dropped). Draft 1 `spec/MG7_proposal.md`; static design
  table `spec/MG7_design_calc.py` (MG4c settled states; no dynamics): P₀(1.5) at κ = 2 grows from −35 to −45 nodes and
  M₋ by 46% from b = 27.5 to 28.5; the φ_G zero crossing stays at r ≈ 1.3–1.4, so falloff stays UNRESOLVED.
- Correction: I had said the shuffled null needs no new runs; it does (4 runs at b = 28 if included).
- Next: PI decisions (§10), ChatGPT review.


## MG023a — MG7 design draft 2 after ChatGPT review MG023, 2026-11-03
- ChatGPT MG023 reproduced the static design table with the production detector and quadrature and asked for six changes
  before registration; all adopted (`review/MG023_response.md`). Key catch: MG6's start loader always loads b = 28, so MG7
  gets its own endpoint-field loader with hashes and b in every run key.
- PI adopted ChatGPT's §10 recommendations (13:50 UTC): endpoints 27.5/28.5 with MG6's b = 28 midpoint, CP2 and SH2, no
  shuffled null, endpoint step and tick checks in Part 0, window fixed, 90 core-h cutoff, cap 22, 8 workers.
- Draft 2 `spec/MG7_proposal.md` (draft 1 kept as `spec/MG7_proposal_v1.md`): gate LARGER AT HIGHER b (directional),
  16 runs, about 74–75 core-h plus about 25 for Part 0.
- Next: merge; ChatGPT round 2 or registration (PI decides); then implementation.

## MG023b — MG7 registered, 2026-11-03
- PI (13:55 UTC, project thread): "register you plan as is". Draft 2 registered without a second review round; content
  frozen at `5b44cde`; GATES.md MG7 entry added. The multiway/foliation gate is renumbered MG8.
- Next: implementation (endpoint-field loader and manifest, driver, evaluator, recompute, tests) and the pre-run
  implementation specification in GATES.md; then Part 0 on the PI's PC.

## MG024 — MG7 implementation and pre-run specification, 2026-11-03
- Code: `run_mg7.py` (own start-field loader with the committed endpoint manifest; refuses any b but 27.5/28.5; no
  shift), `run_mg7_part0.py`, `evaluate_mg7.py`, `recompute_mg7.py`; MG6's engine reused byte-identical to `c795770`
  (tested). `tests/test_mg7.py`: 8 tests pass. Pre-run specification in GATES.md (MG7, items 1–10).
- Check 1 computed statically in the cloud: all starts single, max φ_G 0.0564, mean P₀^start −35.00 / −44.78 (matches
  ChatGPT's replication). A 0.3-unit smoke run of the driver and of Part 0 in a temporary folder completed (not kept).
- Next: merge; Part 0 on the PI's PC (`$env:MODELG_WORKERS=8; python src/run_mg7.py part0`, about 22–25 core-h, several
  hours); then the 16 runs.

## MG024a — MG7 Part 0, 2026-11-03
- Ran on the PI's PC 07:29–10:42 local (code `ca2c513`, 8 workers), 8.63 core-h: **PROCEED**. Step check ΔD̄₀ 0.17 /
  0.17 nodes (limit 2), tick check 0.5 / 1.0 (limit 6), bound 0.062 (limit 0.14); projection 74.2 core-h (cutoff 90).
  Records `9b20219`; outcome in GATES.md.
- Next: merge; launch the 16 registered runs (`$env:MODELG_WORKERS=8; python src/run_mg7.py`, about 10–11 h).

## MG025 — MG7 registered runs and result, 2026-11-04
- 16/16 COMPLETE on the PI's PC (2026-11-03 12:44–21:53 local, 8 workers), 73.01 core-h; records `282726f`.
- **Gate NOT SCORED (feedback):** at b = 27.5 the particle dissolved in 7/8 runs (SH2 4/4, CP2 3/4); b = 28.5 SURVIVES,
  8/8 single and pinned. The rewriting, not the coupling, removes the b = 27.5 particle.
- Post hoc descriptive (not registered): CP2 − SH2 = −72.0 at 28.5 vs −67.1 at 28.0 (MG6), different streams.
- Recompute bug (NaN at no-particle checks) fixed after the data; cloud cross-platform recompute agrees; strict PC rerun
  (2026-11-04) agrees, hashes included. Next: merge; strict recompute on PC; ChatGPT review; PI sign-off.

## MG025b — ChatGPT review MG025 of the MG7 result, 2026-11-04
- Strict Windows replay and an independent audit reproduce everything; ChatGPT recommends accepting NOT SCORED
  (feedback), 27.5 CHANGED, 28.5 SURVIVES, with no new runs.
- Corrections adopted (`review/MG025_response.md`): no survival edge is claimed (b = 27.5 failed the finite-panel
  criterion; no DEC arm, so uniform rewriting is not isolated); threshold sensitivity was computed (unchanged verdicts);
  early/late means at 28.5 and the list of unavailable quantities added. Awaiting PI sign-off.
- **MG7 CLOSED 2026-11-04 14:38 UTC** (PI sign-off in the project thread): NOT SCORED (feedback); 27.5 CHANGED, 28.5
  SURVIVES. Mass scaling remains open.

## MG026 — MG9 design proposal, draft 1, 2026-11-04
- PI asked (project thread, 16:03 UTC) for a draft study of whether b acts as the control parameter of a subcritical
  Landau-type instability, after discussing Higgs-like readings of b with Claude and ChatGPT.
- `spec/MG9_proposal.md` (draft 1, not registered): deterministic only. A1 linear coefficient r(b); A2 cubic Landau
  coefficient L₃ (P1: subcritical), with a fitted-growth validation; A3 front pinning against the particle window (P2,
  primary: COINCIDES / OVERLAPS / DISJOINT); descriptive seed-threshold barrier proxy. Stochastic S(b) deferred to a
  conditional MG9b. About 1–2 core-h, cloud. Prior knowledge that bears on P1 and D2 is disclosed in §2.
- No MG9 computation run. One timing check of the existing b = 29 reference ODE (about 8 s per 400 time units).
- Next: PI decisions in §10; ChatGPT design review; registration.

## MG026a — MG9 draft 2 after ChatGPT review MG026, 2026-11-04
- ChatGPT review `review/MG026_chatgpt.md` (static checks in `review/MG026_replication/`): recommended draft 2 before
  registration; seven items. The PI adopted its recommended decisions (16:23 UTC, project thread).
- All seven items accepted (`review/MG026_response.md`). Draft 2 `spec/MG9_proposal.md` (draft 1 kept as
  `spec/MG9_proposal_v1.md`): P1 = local onset type (either sign valid); P2 = finite-time front arrest against the
  operational window, never "snaking"; frozen domains, normalization, L₃ formula, V2 fit rules; persistence bracket
  (not a fold); interface-envelope observable, UNRESOLVED and robustness repeats; benchmark recheck; budget 193 + 27
  reserve, cap 220, cutoff 8 core-h. Seed-threshold extension deferred.
- No MG9 computation run. Next: PI decides on a ChatGPT round 2 or registration.

## MG026b — MG9 draft 3 after ChatGPT round 2, 2026-11-04
- ChatGPT r2 (`review/MG026_chatgpt_r2.md`; synthetic checks `review/MG026_replication/draft2_checks.*`): supports the
  scope; three items (physical-width envelope, ordered classification, failed preparation ≠ evidence) and seven
  implementation points.
- All accepted into draft 3's text (`review/MG026_response_r2.md`; draft 2 kept as `spec/MG9_proposal_v2.md`):
  λ-width boxcar envelope with synthetic validation and patch tracking; one ordered classification procedure;
  PREP-UNRESOLVED → U, PREP-LOST/CHANGED → N with censored edges; robustness follows final edges; frozen benchmark
  brackets (b_s keeps MG003's potential carry, b_r its seed); `mg9_ode` wrapper exposing solver success; P1 outcome order
  incl. DEGENERATE; budget 210 + 20 reserve, cap 230, cutoff 10 core-h, six-stage stop priority.
- No MG9 computation run (a synthetic kernel check only). Next: PI decides registration or a round 3.

## MG026c — MG9 registered, 2026-11-04
- PI chose "Register draft 3" on the decision card (16:40 UTC), without a ChatGPT round 3 (ChatGPT r2 had no remaining
  scope objection once its items were fixed). Draft 3 content frozen at `7810b8b`; GATES.md MG9 entry added.
- Next: implementation (code, tests with the synthetic envelope/classification validations) and the pre-run
  implementation specification in GATES.md; ChatGPT reviews those before the first run.

## MG027 — MG9 implementation, pre-run specification and amendment 1, 2026-11-04
- Code: `src/landau_mg9.py`, `src/run_mg9.py`, `src/evaluate_mg9.py`; `tests/test_mg9.py` (16 pass). The only Model G
  integration in the tests is MG003's b = 29 reference (wrapper agreement with `lattice_ode_1d`, ≤ 1e-12).
- The registered §6.2 envelope failed its own synthetic validation (edge error 0.27 λ vs λ/10). Per the registered rule
  the study is blocked until fixed: **amendment 1** (`spec/MG9_amendment_1.md`) proposes complex demodulation with a
  double λ-boxcar (0.036 λ). Awaiting PI sign-off.
- Pre-run implementation specification appended to GATES.md (15 items). V1 thresholds and α reproduce ChatGPT's static
  values. No P1, sweep, front or benchmark run made.
- Next: PI sign-off on amendment 1; ChatGPT review of the code and pre-run specification; then stage p1.
- PI approved amendment 1 (decision card, 16:59 UTC). Next: ChatGPT review of the code and pre-run specification; then stage p1.
- ChatGPT review MG027 (`review/MG027_chatgpt.md`): supports amendment 1 (envelope reproduced independently; worst edge error 0.037 λ). One wording fix accepted: cancellation at 2k is exact only for the continuous boxcar (discrete single-pass gain 1.47e-3 / 1.99e-4). Response `review/MG027_response.md`. The review covered the envelope only.

## MG028 — MG9 runs, 2026-11-04
- PI chose "Start the runs now" on the decision card (17:10 UTC), with ChatGPT's review limited to the envelope (MG027). Production code: main at f7ef9e9 plus the documentation-only 54b2a08. Stage p1 launched 17:11 UTC in the cloud (Python 3.13.18, numpy 2.5.3, scipy 1.18.1).
- **Stage p1 (17:11–17:17 UTC, 4 integrations, 310 core-s): P1 = SUBCRITICAL.** L₃ (real cosine amplitude, original time) = 0.013653 (i, D02 mode 34), 0.014162 (ii), 0.013128 (iii), 0.012602 (iv); all positive, no near-degeneracy (terms 0.0156, 0.0131, −0.0151 for (i)). Code check passed. V2 passed: fitted 0.013650 (A₀ = 0.1) and 0.013641 (0.2), within 0.1%; A₀ = 0.05 fit 0.013653 but below the 25% identifiability change (20%); tight-tolerance repeat 0.0136500. Independent recompute of (i) with finite-difference Taylor coefficients of the raw reaction terms: 0.013653.
- Sweep, bench, grid, refine, d01 ran 17:17–17:42 UTC (150 integrations in all, 0.49 core-h). **P2 = UNRESOLVED:** no
  arrested point; 29.25–29.35 creep slowly (UNRESOLVED-SLOW), 29.40–30.00 invade, below 29.22 the 34-cell pattern is
  unavailable (30-cell wavelength change, then lost below 28.62). Benchmarks reproduce MG003 to 0.001. Result in
  GATES.md. Descriptive: creep→invasion between 29.35 and 29.40 brackets b_r. Next: ChatGPT result review; PI sign-off.
- ChatGPT result review MG028 (`review/MG028_chatgpt.md`): recommends accepting P1 SUBCRITICAL and P2 UNRESOLVED; independent recompute of L₃, V2 fits, all 16 front labels, 110 state hashes. One wording fix accepted ("no near-degeneracy"). Response `review/MG028_response.md`. Awaiting PI sign-off.
- **MG9 CLOSED** (PI decision card 19:21 UTC, "Sign off and close MG9"): P1 SUBCRITICAL; P2 UNRESOLVED. S(b) deferred to MG9b.

## MG029 — MG9c design (slow-front follow-up), 2026-11-04
- ChatGPT recommended (PI relay, 19:33 UTC) a follow-up on MG9's three slow fronts; the PI asked for a draft (19:36 UTC).
- Prior look at committed MG9 tracks only (no integration): the 29.25–29.35 interfaces relax exponentially (time constant about 715; about 0.22 left to move; exponential rms 2.5e-4 vs linear 5.8e-2); cell count constant at 18. Disclosed in the draft §0.
- `spec/MG9c_proposal.md` DRAFT 1: late-time windows from t = 1800 (T = 12,000, extension 30,000); labels SETTLED / ADVANCES; robustness D01, D02L, ±5% perturbation; Q2 bisection of the settling edge between 29.35 and 29.40 vs b_r (tolerance 0.01, PI decision). Cap 60 integrations, cutoff 2 core-h, cloud. Not registered.
- ChatGPT review MG029 (`review/MG029_chatgpt.md`): proceed after a short revision. All items accepted into **draft 2** (draft 1 kept as `_v1.md`; response `review/MG029_response.md`). Corrections: 29.35's detector count went 18 → 19 at t = 45 (Claude had generalised from 29.30); Q2 uses interval arithmetic and NEAR b_r wording (draft 1's "cannot be found above b_r" was wrong); opposite-sign rule, MG9's separate terminal and structural windows, fixed-horizon perturbation runs, PASS/FAIL/UNRESOLVED robustness, replay-mismatch rule; D02L criterion uses half-width 25; snapshots at every window boundary. Not registered.
- ChatGPT round 2 (`review/MG029_chatgpt_r2.md`): draft 2 ready for registration; one nonblocking §10 MIXED wording correction, applied. Response `review/MG029_response_r2.md`. Awaiting PI registration.
- **MG9c REGISTERED** (MG029c; PI decision card 19:56 UTC, "Register draft 2", §9 defaults). Draft 2 frozen at `372336b`. GATES.md entry added. Next: code, synthetic tests, pre-run specification; then the cloud run.
- MG029d: code (`src/run_mg9c.py`, `src/evaluate_mg9c.py`), 20 tests (`tests/test_mg9c.py`, all pass; MG9's 16 still pass) and a 10-item pre-run specification in GATES.md. No MG9c run yet. Next: stage q1 in the cloud.

## MG030 — MG9c runs, 2026-11-04
- PI chose "Start the runs now" (decision card 20:40 UTC). Code `cee5498`, cloud (versions as MG9). 11 integrations, 0.23 core-h.
- Replays of MG9's [0, 1800] tracks exact (max diff 0). **Q1 = UNRESOLVED:** 29.25/29.30/29.35 all UNRESOLVED-SLOW at T = 30,000, failing only the residual criterion (1.084–1.092e-5 vs 1e-5); edges otherwise stationary (|D₃| ≤ 0.009, |v| ≤ 4.5e-7, counts constant). Q2 not eligible.
- Descriptive: 29.22 SETTLED (residual 9.45e-6), robustness ROBUSTNESS-UNRESOLVED (D01 preparation PREP-CHANGED; D02L and ±5% PASS); width trials settled, D02 hw 25 + λ/4 locks to the hw 25 edges. Late motion at all four b is a slow rigid drift (~5e-7 per unit time, decelerating) with the width frozen. Result in GATES.md. Next: ChatGPT result review; PI sign-off.
- ChatGPT result review MG030 (`review/MG030_chatgpt.md`): accept Q1 UNRESOLVED and Q2 NOT ELIGIBLE; all labels, hashes and residuals reproduced. Two descriptive corrections accepted (width rate: "nearly constant" over [15,000, 30,000], not "< 4e-8 after 6000"; width-trial edges within 1.6e-4, not 1e-4) and one documentation slip in pre-run item 8. Response `review/MG030_response.md`. Awaiting PI sign-off.
- MG030b (exploratory, post-result; no integration): centre diagnostic requested in `review/MG030_amplitude_center_prompt.md`. `src/diag_mg9c_centre.py`, `results/MG9c/diagnostics/`. Validated on synthetic patches (translation recovered to 1e-6; an asymmetric amplitude change moves weighted centroids by 0.035 with minima and alignment fixed). Late interval [15,000, 30,000]: edge midpoint, E²- and φ_X²-weighted centres, interior minima and profile alignment all agree (−0.0081 to −0.0093); minima move together (sd ~1e-4); translation explains 98–99% of the change. Supports coherent translation; cause not identified. No classification changes.
- **MG9c CLOSED** (PI decision card 22:29 UTC, "Sign off and close MG9c"): Q1 UNRESOLVED; Q2 NOT ELIGIBLE.

## MG031 — MG9d design (cause of MG9c's late drift), 2026-11-05
- PI request (00:30 UTC, after discussion with ChatGPT): a small preregistered diagnostic starting from MG9c's b = 29.30 t = 30,000 field. `spec/MG9d_proposal.md` DRAFT 1: seven runs (baseline; tighter solver; shifts h/4, h/2, 3h/4; node reflection; same patch interpolated to h = 0.1), T = 10,000, diagnostics every 50 units, ordered decision rules (INCONCLUSIVE allowed), cap 10 integrations, cutoff 1 core-h, MG9b-impact section. Static checks only (spectral content, shift/roll, reflection map). Not registered; no integration.
- Side fix (MG030d): centre-diagnostic corrections per ChatGPT's supplement review (`review/MG030_centre_response.md`).
- MG031 review (ChatGPT, `review/MG031_chatgpt.md`): keep the seven-run design, revise rules before registration. Five items, all accepted (`review/MG031_response.md`). One was an error of mine: fractional Fourier shifts do not round-trip at round-off (7.8e-6 to 1.6e-5, Nyquist coefficient).
- MG031a: `spec/MG9d_proposal.md` DRAFT 2 (draft 1 → `spec/MG9d_proposal_v1.md`, superseded). C0 removed (state-transfer checks instead); signed velocities; Nyquist-filtered common start for all arms; motion floor and validity rules for every run; complete ordered C1–C5 rules; reflection sign applied once; MG9b wording revised. Not registered; no integration.
- MG031 round 2 (ChatGPT, `review/MG031_chatgpt_r2.md`): design sound; three decision-rule corrections before registration. All accepted (`review/MG031_response_r2.md`).
- MG031b: `spec/MG9d_proposal.md` DRAFT 3 (draft 2 → `spec/MG9d_proposal_v2.md`, superseded). Changes: NO RESOLVED DRIFT needs a, b and v unresolved on R0 and R1; C3 NOT EVALUATED routes to INCONCLUSIVE; new C4 label DECREASED BELOW RESOLUTION; E_pat frozen at the exact recorded value; the filtered start's residual is disclosed. Design unchanged. Not registered; no integration.

## MG032 — MG9d implementation and pre-run specification, 2026-11-05
- PI approved draft 3 for implementation (card, 01:42 UTC: "Approve draft 3 as written").
- Code `src/run_mg9d.py`, `src/evaluate_mg9d.py`; tests `tests/test_mg9d.py` (15 pass). Static only: state-transfer checks pass; measurement floor δ = 1.1e-9 → v_f = 5e-9.
- Pre-run spec `spec/MG9d_prerun.md`. Findings: the Nyquist filter raises the initial max|RHS| 150× (1.1e-5 → 1.65e-3; disclosure fix P0); shifted and D01 starts need early shape adjustments (8.5e-4 to 1.8e-2 relative), so validity items 5–6 measured from t′ = 0 would invalidate them by construction (proposal P1: measure from t′ = 2000); C_edge's static error 2.3e-4 exceeds item 6's 2e-4 floor (proposal P2: 5e-4). Awaiting ChatGPT implementation review and PI decisions. No integration.
- PI adopted P1 and P2 (card, 01:50 UTC); written into draft 3 §3 items 5–6, with P0 in §0. Next: ChatGPT implementation review, then registration.
- MG031 round 3 (ChatGPT, `review/MG031_chatgpt_r3.md`): design supported, P0–P2 accepted; asked to freeze the velocity estimator. Accepted (`review/MG031_response_r3.md`): endpoint averages written into draft 3 §3 and disclosed as a change from draft 1's least-squares v; least squares descriptive only; transient wording qualified; 16th test added. No integration.
- MG9d REGISTERED 2026-11-05 (PI card 02:00 UTC, "Register; runs after code review"); draft 3 frozen at 69ad83e. Runs wait for ChatGPT's implementation review and a PI go.
- ChatGPT implementation review MG032 (`review/MG032_chatgpt.md`): science and calibration pass; two operational fixes (failure provenance, strict runtime cutoff) plus cleanup, all made (`review/MG032_response.md`, prerun §E; 22 tests). PI approved starting the runs once fixed (02:12 UTC).

## MG033 — MG9d runs and result, 2026-11-05
- PI go 02:12 UTC (after the MG032 fixes). The 7 arms ran in the cloud: 0.084 core-h, all complete.
- RESULT: MESH EFFECT SUPPORTED (on D02). C1 NUMERICS-AGREE, C3 SYMMETRIC, C2 ALIGNMENT-DEPENDENT, C4 DECREASING (ρ = 0.845), C5 R4 INVALID (item 5: shape kept changing until about t′ = 6000).
- Descriptive: shifted arms drift at −1.41e-6, +6.0e-7 and +1.18e-6 against R0's −3.45e-7, which looks like lattice pinning. On D01 the centre is static and the residual is 2.1e-9.
- Independent recompute from the snapshots matches. Next: ChatGPT result review, PI sign-off.
- ChatGPT MG033 accepted the result (independent replication); wording corrections made (`review/MG033_response.md`). MG9d CLOSED 2026-11-05 (PI message 02:32 UTC: proceed to MG9b).

## MG034 — MG9b design (stochastic barrier S(b)), 2026-11-05
- PI request (02:32 UTC): proceed to MG9b. `spec/MG9b_proposal.md` DRAFT 1:
  - Design: MG3 protocol at b = 29.20 and 29.35. A Part A Ω-placement scan per b, then three Part B cells (Ω_T, −15, −30) with MG3b batching (≥ 20 events, cap 64).
  - Primary question: the trend S₁ = dS/db from a joint Poisson GLM with MG3b's b = 29.00 data. Labels DECREASING / INCREASING / NOT RESOLVED, with a synthetic calibration.
  - Lattice-conditional checks V-b1 to V-b3 (reference stationarity, t = 2000 confirmation controls, confirmation solver tolerance), per MG9d/MG033.
  - Budget about 34 core-h projected, cutoff 50, on the PI's PC.
  - Not registered; no run.
- PI adopted all §8 recommendations as written (02:39 UTC). Next: ChatGPT design review.
- ChatGPT MG034 (`review/MG034_chatgpt.md`): revise before registration; PI choices kept. All points accepted (`review/MG034_response.md`). Draft 2 (`spec/MG9b_proposal.md`; draft 1 kept as `_v1`, superseded):
  - Outcomes NEGATIVE/POSITIVE LINEAR TREND / NOT RESOLVED, with per-b slopes and a 1-df linear-model check. Lack-of-fit df corrected (three cells give 1 df; my error).
  - Calibration fully specified (historical b = 29.00 regenerated under MG3b's rule; ordinary, startup-stress, placement-stress and non-monotone scenarios). Not calibrated → a directional result is EXPLORATORY.
  - V-b1 to V-b3 defined exactly (late window 1900–2000, ten unique fields, rtol/atol in the confirmation path, b_det = b). A V-b1 failure stops that b.
  - Unique replacement candidates; pilots and screens inside the 50 core-h cutoff (about 35 projected); partial analysis at a budget stop frozen.
  - Not registered; no run. Next: ChatGPT review of draft 2.
- ChatGPT MG034 round 2 (`review/MG034_chatgpt_r2.md`): substantial improvement, short correction pass. All accepted (`review/MG034_response_r2.md`). Draft 3 (draft 2 kept as `_v2`):
  - V-b2 requires all ten pairs to complete (else CONTROL-NUMERICAL-FAILURE) before the 9-of-10 rule; controls read as confirmation robustness, not single-particle identity.
  - Calibration simulates MG3b's Part 1 τ/τ2 selection and keeps EQUIVALENT replicates; 4000 replicates per case; per-branch labels (≥ 400 retained replicates).
  - BUDGET-SHORT directional results are EXPLORATORY; structural identifiability separated from fit failure (NOT RESOLVED (NUMERICAL)).
  - Not registered; no run.
- MG9b REGISTERED 2026-11-05 (PI card 03:01 UTC, "Register draft 3 after merge"); draft 3 frozen at 185aceb. Next: code, pre-run spec and synthetic calibration, then ChatGPT implementation review and a PI go. No MG9b run.
- MG9b implementation (after registration; no production run): `src/mg9b.py`, `src/run_mg9b.py`, `src/evaluate_mg9b.py`, `src/calib_mg9b.py`, 20 tests; `confirm_field` gains checks/rtol/atol and a completion flag. Synthetic calibration (cloud, 0.36 core-h): B3 and B2-20 calibrated, B2-35 not supported (exploratory); B3 coverage 0.946–0.955, false-positive rates 0.021–0.028 per side; power 0.73 at S₁ = −0.10. Cloud preview of V-b1: both references pass. Pre-run spec `spec/MG9b_prerun.md`. Next: ChatGPT implementation review, then PI go.
- ChatGPT MG034 round 3 (`review/MG034_chatgpt_r3.md`): design accepted, one correction (study-level budget truncation) plus wording. Draft 3 was already registered, so the correction is amendment 1 (`spec/MG9b_amendment_1.md`), implemented and tested (24 tests); response `review/MG034_response_r3.md`. Awaiting PI sign-off on amendment 1; no production run.
- PI adopted amendment 1 (card, 03:27 UTC). Next: ChatGPT implementation review, then the PI go.
- ChatGPT amendment-1 implementation review (`review/MG034_amendment1_implementation.md`): keep amendment 1; four implementation fixes before the go. All made (`review/MG034_amendment1_response.md`); 27 tests pass. No production run.

## MG035 — MG9b runs and result, 2026-11-05
- PI go 04:08 UTC. Run on the PI's PC (Windows 11, Python 3.13.2, numpy 2.2.4; code `4699a6b`; 9 workers), 04:10–12:15 UTC. 280 runs, ledger 26.28 core-h (cutoff 50, no budget stop). Records committed as `5a54f7a` (fields/ kept on the PC).
- RESULT: NOT RESOLVED (branch B3, calibrated, no qualifiers). S₁ = −0.023, 95% CI [−0.094, +0.048]; per-b slopes 0.064 / 0.048 / 0.057; LR check p = 0.30. All controls PASS. `spec/MG9b_result.md`.
- Deviation: two V-b3 tight-tolerance confirmations at 29.20 ran 5.4 h and 4.4 h (median 18 s), so V-b3 at 29.20 cost 10.2 core-h and the run took about 8 h; both completed and agreed. No per-job time limit was registered; recommend one for future studies.
- Independent recomputation from raw run JSONs (`src/check_mg9b.py`) matches the evaluator. Next: ChatGPT result review, PI sign-off.
- Post-hoc (PI request after ChatGPT MG035): event-shape mix by cell, `review/MG035_event_shapes.md`. Single-minimum share falls with lower Ω at both b (29.20: 61/55/44%, trend p = 0.23; 29.35: 50/43/28%, p = 0.09); totals 40/75 and 30/76 reproduced. Descriptive only; registered result unchanged; follow-up at overlapping Ω not recommended as a priority. Addendum (same day, independent recheck): pooled stratified trend z = 2.06, p = 0.04 (post hoc); overlap at Ω ≈ 200 estimated at 15–30 core-h; recommendation unchanged (no follow-up now; a registered two-cell study at Ω ≈ 200 if shape becomes a question).
- ChatGPT MG035 result review (`review/MG035_chatgpt.md`): accept and close as NOT RESOLVED; independently reproduced the fit, controls (V-b2 replayed through t = 2000) and ledger. Two wording corrections applied to `spec/MG9b_result.md`. Event-shape review (`review/MG035_event_shapes_chatgpt.md`): accepted; pattern archived as a hypothesis. Response `review/MG035_response.md`. Awaiting PI sign-off to close MG9b.
- MG9b CLOSED 2026-11-05 (PI decision card, 14:03 UTC): NOT RESOLVED. Shape pattern archived as a hypothesis.

## MG036 — MG8 design (multiway / update order), 2026-11-05
- PI request (14:09 UTC): plan MG8, Central Question 5, whether different admissible orders of applying the rewrite rules give measurably different outcomes. Confirmed against README, InitialPlan Phase 6, rules_v1.1 §4 and GATES.md (renumbered from MG7 at MG023b).
- Draft 1, `spec/MG8_proposal.md`: Part A structural arguments (the continuum equations already have a medium rest frame; token rules are state-confluent except resource conflicts; the registered spatial tick is order-free; directional sweeps advect each species at about D²τ/h³, Y at 15 per time unit at τ₀). Part B gate on the 1D lattice: schedulers J (registered), R (sweep), P (random sequential), L (mirror check); deterministic τ ladder τ₀ to τ₀/64 and stochastic runs at Ω = 1000. Gates MG8-F (order matters at the production step) and MG8-C (VANISHES / PERSISTS as τ → 0); Q5 mapping fixed in advance. About 8–13 core-h, cutoff 20.
- Not registered; no run. Next: PI decisions (§9), then ChatGPT design review.
- ChatGPT MG036 (`review/MG036_chatgpt.md`): revise before registration. All items accepted (`review/MG036_response.md`). Draft 1's A4 was wrong: a sweep does not advect tokens (the expected displacement is zero; I omitted backward hops on re-processing). Draft 2 (`spec/MG8_proposal.md`; draft 1 kept as `_v1`, superseded):
  - A4 corrected: under a sweep, effective diffusion D/(1 − r) (Y +25% at τ₀), an odd dispersive term κ₃ ≈ D²τ/h, and reaction re-exposure; all O(τ). Exact small-ring benchmark in Part 0.
  - Primary conclusion framed as scheduler sensitivity of this lattice approximation; Q5 stays open; PERSISTS label and the "R bounds milder orders" claim removed.
  - Separate FATE / SHAPE / DRIFT labels with margins (deterministic, primary); stochastic arm descriptive with intervals; convergence rule with floors, span, slope interval, residuals and a horizon continuation.
  - Executable scheduler and stream semantics; synthetic validation; clipping diagnostics; committed snapshots; Part 0 priced outside the 20 core-h cutoff.
  - Not registered; no run. Next: PI decisions (§9), ChatGPT round 2.
- ChatGPT MG036 round 2 (`review/MG036_chatgpt_r2.md`): substantially improved; one more revision before registration. All accepted (`review/MG036_response_r2.md`). Draft 3 (`spec/MG8_proposal.md`; draft 2 kept as `_v2`):
  - B1 defined as random-seam expected-count trajectories, three seam replicates per level; SEAM-DEPENDENT label; seam spread gates refinement eligibility.
  - Mutually exclusive SHAPE/DRIFT labels with INSUFFICIENT RESOLUTION; BELOW MEASUREMENT RESOLUTION (not convergence); horizon check on both branches against J at t = 600.
  - Part 0 benchmark: conditional maps for every order, explicit averaging, basis-vector matrices, augmented cohort benchmark, frozen refinement check.
  - Independent streams per arm, Welch drift-difference interval, specified eligibility and bootstrap; stored decision inputs; A4 claims narrowed.
  - 57 runs (cap 64), about 5–18 core-h, cutoff 20. Not registered; no run. Next: PI decisions (§9), short ChatGPT round 3 check.
- ChatGPT MG036 round 3 (`review/MG036_chatgpt_r3.md`): main objections resolved, no redesign; one correction before registration. ChatGPT independently reproduced the refinement-benchmark ratios (R 0.4967, 0.4984; P 0.5039, 0.5018). All accepted (`review/MG036_response_r3.md`). Draft 4 (`spec/MG8_proposal.md`; draft 3 kept as `_v3`): Part 0 stochastic benchmark tests each replicate against its conditional-order distribution and uses across-replicate SEs against the average; common-shift alignment objective (min of max species error) with tie-breaking and brute-force synthetic reference; seam spread and split-half edge cases. Not registered; no run. Next: PI decisions (§9) and registration.
- MG8 REGISTERED 2026-11-05 (MG036d; PI decision card 20:02 UTC, "Adopt all, register"): draft 4 frozen at 3ab4ec8 with all §9 recommendations adopted. GATES.md MG8 entry added. No MG8 run. Next: code, tests and Part 0, then a ChatGPT implementation review and the PI go before production.

## MG037 — MG8 implementation and Part 0, 2026-11-05
- Code: `src/mg8.py` (schedulers, deterministic maps, Part 0 operators, measurements, labels), `src/run_mg8.py` (driver), `src/evaluate_mg8.py`, `src/part0_mg8.py`; 18 tests in `tests/test_mg8.py`. Pre-run specification `spec/MG8_prerun.md` (continuations as k = 6 records run to 600; Nyquist-free alignment; segmented RNG; disclosure of short development runs).
- Part 0 (cloud): attempt 1 STOP on check 1(b), caused by a fault in the test's tail pooling (one pooled cell with expected count 0.3), not the simulator. Pooling fixed; attempt 2 passes every check; projected production 3.9 core-h. Awaiting the PI's ruling on attempt 2 and a ChatGPT implementation review. No production run.
- 2026-11-05 22:00 UTC: the PI accepted Part 0 attempt 2 (`decisions.json` now `proceed = true`, with the decision recorded). Production still waits for ChatGPT's implementation review and the PI's go.

## MG037a — MG8 implementation review (ChatGPT MG037), 2026-11-05
- ChatGPT reviewed the code and Part 0 at `e050fbe` (`review/MG037_chatgpt.md`): keep the design and the accepted Part 0;
  fix six implementation items before production. All six accepted (`review/MG037_response.md`).
- Fixed: decision windows by unit index with endpoints included and per-unit finiteness; implementation checks gate all
  labels and stop the driver; TRACK-FAIL rescans every end of the track; atomic writes, ledger reconciliation and a
  persisted run state used for the §6 BUDGET-SHORT rules; seam-sequence hashes and B2 block profiles stored;
  diffusion-only cross-check of the production kernels (supplementary Part 0 check, passes). Tests 24, all pass.
- Recorded deviation D1 (`spec/MG8_prerun.md` §7): the scored shape error uses Nyquist-free fields; the literal E_raw
  and the removed-mode E_nyq are reported beside it. Awaiting the PI's sign-off. No production run.
- 2026-11-05 23:37 UTC: the PI accepted D1 and gave the production go (PR #81 merged). Production launched in the cloud
  (Linux, 4 cores, 4 workers).

## MG038 — MG8 production and result, 2026-11-06
- Production in the cloud (Linux, 4 workers), 2026-11-05 23:37 to 2026-11-06 00:31 UTC: 55 records, 3.09 core-h, state
  complete. All implementation checks pass (mirrors bit-identical).
- RESULT: production step **SENSITIVE** (R loses the particle at t ≈ 8 at τ₀, all three seams); refinement **NOT RESOLVED
  BY THIS LADDER** (eligible k = 3–6, E ∝ τ with slope 1.000, but E_min/E_max = 0.125 > 0.1). B2: losses J 7/8, R 8/8,
  P 8/8; drift and shape NOT ESTIMATED. Independent recompute (`src/recompute_mg8.py`) agrees on all 106 quantities.
- Descriptive (unregistered): the sweep shifts the homogeneous steady state by an amount proportional to τ (−0.259 in G
  and X at τ₀); that offset dominates E and, at τ₀, dissolves the particle. Details in `results/MG8/MG8_results.md`.

## MG038a — MG8 result review (ChatGPT MG038), 2026-11-06
- ChatGPT accepted the registered results (`review/MG038_chatgpt.md`) and asked for four reporting corrections, all made
  (`review/MG038_response.md`): J-control losses (all after t = 100; the Ω = 1000 planning assumption was insufficient for
  T = 500); "unchanged particle" withdrawn; causal reading marked provisional and the homogeneous-start diagnostic
  archived (`src/diag_mg8_homogeneous.py`, `results/MG8/diag_homogeneous.json`); seam spread 6.5e-8 to 2.8e-6, not zero;
  full B2 counts discarded (retention deviation). Labels unchanged.
- **MG8 CLOSED 2026-11-06 02:11 UTC** (PI decision card): SENSITIVE at τ₀; refinement NOT RESOLVED BY THIS LADDER; Central
  Question 5 stays open.

## MG039 — Role/protocol transition and MG10 explicit-reservoir draft, 2026-11-06

- PI instruction in the project chat: ChatGPT / Codex takes planning, code, validation, execution, analysis and documentation; Richard retains PI decisions; Claude or another external reviewer is optional. The PI explicitly authorized updating the role/review protocol across documents and preparing MG10.
- Added RESEARCH_PROTOCOL.md as the active prospective authority. Updated README.md and CLAUDE.md; added tracked AGENTS.md for Codex; added a current HANDOFF.md section and a supersession notice to historical InitialPlan.md. Updated living glossary process terms and GATES.md process guidance. Historical reviews, frozen criteria and completed scientific records are unchanged.
- Added spec/MG10_proposal.md Draft 1: explicit A/B stocks and Z/W outputs, closed-count accounting, chemostatted baseline and large-capacity limit, numerical semantics and SSA checks, staged qualification/recovery/depletion/interventions, provenance and retention.
- Same-author static design check: review/MG039_internal_check.md. Verified reaction count balances, local documentation links, and preservation of historical handoff/initial-plan bodies; no simulation was run.
- Archive source: RichardBenear/StochasticVacuum at 2161211ccc9d8427af98ebeefcabf7f379b83a37. Reuse remains selective and needs ModelG-specific validation; no old control PASS is a new physics result.
- Status: MG10 DRAFT, NOT REGISTERED. No MG10 trajectory, calibration or production run. Next: static design check and numerical/budget specification for PI consideration. Role transition itself does not authorize production.

## MG039b — MG10 decisions and local-compute preference, 2026-11-06

- PI approved MG10 Draft 1 decisions 1, 2, 3 and 4. External review remains Richard's per-run/batch decision.
- Prefer the local PC for runs projected to pressure assistant/service usage limits. Authorized maximum: eight worker processes, leaving capacity for background activity; numerical-library threads should be one per worker.
- Updated active protocol, AGENTS.md, CLAUDE.md, handoff, README and GATES guidance. Preserved MG10 Draft 1 as spec/MG10_proposal_v1.md; Draft 2 records the approval without inventing unspecified numerical caps.
- No MG10 trajectory or production run; registration awaits completion of numerical criteria, run caps, cost projection and run-level review disposition. Do not request approval again for unchanged decisions already approved.

## MG040 — MG10 numerical candidate and static capacity audit, 2026-11-06

- PI authorized finalizing numerical criteria and budgets. Static calculation (spec/MG10_design_calc.py; results/MG10/design/static_design.json), not dynamics: deterministic A availability is exp(-t/(C*T_ref)); a homogeneous quasistatic G-response estimate gives 8.28% pattern-normalized discrepancy at C=10000. Approximate C for 1% is 82775; this is not a rigorous soliton bound.
- PI explicitly approved extending the capacity grid to 100000 and 1000000. Capacity changes stored values, not lattice size or individual-token loops.
- Added spec/MG10_run_plan.md: primary unshifted error normalization; physical checkpoints; 1% recovery margin; 0.1% numerical-reference discrepancy limit; full fixed deterministic matrix of 64 trajectories; local <=8 workers, 4 core-hour production ceiling, 0.25 qualification/timing cap and 600-second job timeout; no automatic extensions.
- Draft 3 preserves Draft 2. No model trajectory, production, qualification or timing run; only static arithmetic on an existing committed reference. Firm runtime projection awaits implementation benchmarking. External review remains Richard's per-run/batch decision.

## MG041 — Initial MG10 qualification stopped (2026-11-06)

PI: “No review needed. Go ahead.” External review is waived for this batch. Implemented the decoupled A exact-case qualification first; four of 21 Euler paths fail the drafted continuous-exponential tolerance. Two unit tests pass, separate Decimal recomputation agrees on every label. This identifies a qualification-specification error, not a scientific recovery result. No full model, timing matrix or production dispatched. Evidence in results/MG10/qualification; proposed correction in spec/MG10_qualification_correction.md. The unchanged draft remains unregistered pending PI adoption of the correction and the remaining work.


## MG042 — Corrected qualification and prospective registration (2026-11-06)

PI adopted MG041 correction and authorized proceeding. Deterministic reservoirs, independent sparse BDF reference, bounded owned-process driver, evaluator and separate checker implemented. Seven qualification tests PASS after eliminating identically constant ODE coordinates; failed development attempt preserved. All-capacity one-node exact cases and actual stop/resume checks PASS. Two declared T10 timing cases independently audited. Projection 1.802273 core-hours for 64 jobs, within 4 core-hour ceiling with reserve. Raw local retention verified. Registered before scientific production, with source/spec/reference hashes and unchanged recovery thresholds. See review/MG042_internal_check.md and results/MG10/registration.json.


## MG043 — MG10 first batch complete: UNRESOLVED (2026-11-06)

All 64 registered runs complete, with 0.21418268 local production core-hours and no timeout or budget stop. All ledger/availability checks pass; separate raw recomputation agrees. Both top capacities recover for both starts (soliton normalized errors about 0.825% and 0.0825%), but C=1 homogeneous numerical disagreement forces overall UNRESOLVED. Posthoc existing-data diagnostics also show FFT-roundoff detections on exactly uniform depleted backgrounds; no frozen detector or verdict changed. Full raw fields/fluxes, hashes, evaluation, recomputation and ledger retained. Report: results/MG10/MG10_results.md. Follow-up recommended on homogeneous invariance, controlled perturbations and depleted-background measurement. No additional production authorized here; PI closure pending.


## MG044 — MG10-C1 qualification and prospective registration (2026-11-06)

PI authorized the recommended next study. Independent symmetry-preserving references, all-mode time-dependent linear response, planned mode selection, controlled nonlinear amplitudes and a prominence-qualified minimum detector are specified in spec/MG10C1_proposal.md. Six tests and separate raw checking pass; 23 detector fixtures retained. Failed development interpolation/checker attempts preserved and charged. Six T5 timings project 0.1285662 core-hours for ten fixed jobs, below new 0.5 production core-hour ceiling with reserve; <=8 local workers. Registered before production with full hashes and durable raw retention. Original MG10 is not rescored. No external review requested with execution instruction.


## MG045 — MG10-C1 results (2026-11-06)

All ten jobs complete after prospective registration1451074, separate recomputation agrees. Homogeneous reference resolved, controlled amplification supported, contrast-minimum qualification passes. Uniform E2 follows reduced ODE; seeded E2/O+ both match linear mode32 response in frozen windows (worst relative error0.03295%, bound5%). All23 fixtures pass; uniform E2 yields no resolved minima at all301 checkpoints. Production0.03238218 core-hours, no cap/timeout/replacement. Complete raw fields, dense references, all250-mode responses, ledgers, fixtures and environment retained. Exact homogeneity differs from perturbed stability; no physical noise law or self-sustaining-particle claim. Original MG10 remains UNRESOLVED; PI closure pending.


## MG046 — PI closes MG10-C1 and authorizes next study (2026-11-06)

Richard explicitly signed off MG10-C1 and instructed proceeding. Preserve all MG045 outcomes and original MG10 UNRESOLVED result. Next: MG10-S1, bounded finite-reactant stochastic approximation and exact small-system qualification before any full-lattice finite-supply particle claims. Scope/criteria/budgets will be registered prospectively.


## MG047 — Finite-reactant sampler qualification registered (2026-11-06)

Bounded A and joint B/X leap, independent exact SSA, integer reaction/hop ledgers, confidence-band evaluator/checker and generic owned-process controller implemented. Five sampler plus three apparatus tests PASS; analytic and synthetic raw checks and all16 timing groups independently verified. Registered16 groups x32768 short paths, 96 primary marginal comparisons with99% simultaneous bounds/.05 margin; all1024 chunk seeds fixed and unique. Forecast0.2867493 core-hours fits new1 core-hour cap with reserve. Full fields, ledgers and seeds retained; no scientific group before registration. This qualifies a prerequisite, not physical particles or all parameter regimes.


## MG048 — MG10-S1 numerical qualification result (2026-11-06)

QUALIFIED ON TESTED SMALL SYSTEMS. All524288 short trajectories in16 fixed groups complete after aeb0b39 registration. Largest simultaneous99% CDF upper .0361553 versus margin .05; all96 primary comparisons qualify. Exact integer population/reaction/hop ledgers and probability bounds pass; all1024 seeds verified. Independent raw arithmetic/CDF checker agrees. Production0.04894908 local core-hours, no timeout/cap stop/replacement. Complete raw records, environment and report retained. Full joint law, actual lattice/high-Omega noise and particle persistence are outside qualification; future stage must freeze those checks. MG10-S1 closure pending; MG10-C1 PI closure recorded MG046; original MG10 remains UNRESOLVED.


## MG049 — PI closes S1 and authorizes S2 (2026-11-06)

Richard signed off MG10-S1 and instructed proceeding. S1 CLOSED as qualified on tested small systems. Next is actual lattice-spacing/density validation including spatial-noise observables; original MG10 not rescored. New scope/criteria/budget will be frozen before runs.


## MG050 — S2 qualified and prospectively registered (2026-11-06)

Actualh.2/Omega100,1000 full500-node H/S snapshot-law oracle and short16-node SSA comparisons implemented with spatial powers/cross-products. Eight tests and independent raw development checker PASS. Before scientific data, static confidence-width review revises candidateR8192 to16384 while preserving old draft/calibration and unchanged margins/physics. All32 class timing measurements reused to forecast0.850748 core-hours, fitting2-hour ceiling with25% reserve; qualifications102.715845 seconds under900. Fixed88 jobs,32 logical groups,524288 paths,1280 unique chunk seeds. Source/spec hashes and durable reconstructible raw retention registered before production. S1 PI closure recorded MG049; original MG10 unchanged.


## MG051 — Actual-lattice spatial-noise qualification result (2026-11-06)

QUALIFIED WITHIN REGISTERED S2 SCOPE. All88 fixed jobs/32 groups/524288 trajectories complete after77eafbf. Full500-node H/S snapshot laws and short16-node SSA spatial laws at h.2/Omega100,1000 qualify under frozen margins; maximum upper bounds.0573459/.0507281 versus.08. Separate checker verifies every one of1280 raw shards and288 comparisons, exact integer states/hashes, conditional means/noise matrices and seeds. Production0.362493 core-hours, no timeout/cap/replacement; observed monitor gap<=.251052. Full1.73GB archive committed/uploaded in raw-data batches. Long-horizon particle behavior remains for a separately designed study; original MG10 not rescored. S2 PI closure pending.


## MG052 — S2 closure and finite-supply pilot preparation (2026-11-06)

PI signed off S2 and instructed proceeding. S2 CLOSED within its registered scope, without rescoring. Preparing MG10-D1: long-horizon seeded/homogeneous finite-supply controls, explicit chemostat input ledgers, and numerical refinement. Development timing is disclosed separately; production remains unregistered pending a concrete numerical plan and budget.

MG052 preparation complete: Draft2, qualified controls and cost forecast ready for PI budget/scope decision. Twenty implementation checks and three synthetic measurement checks pass. Flat-field RMS roundoff found/fixed prospectively without changing sampled dynamics; failed fixture evidence retained. Independent raw reconstruction agrees. Forecast41.57core-hours with25% allowance; candidate60-hour cap,8workers,3600seconds/job. Original48-hour candidate fails reserve. No production.


## MG053 — PI approves and D1 is prospectively registered

Richard approved the run and60-core-hour cap. Completed checkpoint/controller and end-to-end synthetic reporting qualification. Eighty fixed trajectories, T60, local8workers with one library thread; measured forecast41.682623core-hours with slowdown/startup allowance. Registration and source/spec/reference/seed freeze written before production. No scientific outcome inspected. Full checkpoint retention and automatic incomplete-safe evaluation/independent audit prepared.

MG053 PRODUCTION LAUNCH: started2026-11-07T00:46:02.469149+00:00; registration b21fb85 pushed before launch. Eight owned workers active; all first-wave paths saved their first physical-time checkpoint with accounting/safety checks. Full80-path denominator unchanged. No result/verdict yet. Automatic report and independent check run after batch stops; full raw upload and PI closure remain pending.


## MG054 — D1 completed; posthoc source discrepancy retained

All80 trajectories complete in4.40981elapsedhours,33.333427worker-core-hours versus60 approved. Integer accounting and separate480-summary audit PASS. A-source law check exposes approximately0.3% HIGH fine-step deficit; suspect numerical sampler cancellation, not resource physics. Preserve all original results and caution interpretation. PI authorizes next steps. Full production archive is being committed/uploaded; subsequent diagnostic/correction must be prospectively specified.


## MG055 — Source-law check, stable correction and registered comparison

Confirmed legacy fine-step zero-event probability differs from stable finite-binomial law. Revised geometric-spacing finite source qualified; independent80-digit histogram/CDF and source/state reconstruction agree. Original and NumPy diagnostics/threshold outcomes retained without relabeling. New32-job CHEMO/HIGH comparison registered with fresh streams and full checkpoint retention; estimated~1.94hours on8workers. Uses24-core-hour subcap within remaining60-hour approval; no larger budget assumed. Production registration commit/push must precede launch.

MG055 CORRECTED PRODUCTION START: 2026-11-07T13:37:34.306685+00:00, registrationff12250 pushed before launch. Eight active owned workers and checkpoint saves verified. Thirty-two fresh trajectories,24-hour subcap within prior60-hour approval; no corrected result yet. Process live: do not start another driver. Automatic report/independent audit after stop.


## MG056 — Corrected D2 complete and inspected

32/32 paths, no failures/timeouts/cap; independent192-summary/raw audit PASS. Stable finite-binomial source avoids old deficit; small remaining mean offsets and seeded contrast variation retained with explicit limits, no stochastic/particle equivalence declaration. Full raw archives and report prepared for publication. No next production registered; PI closure separate.

MG056 RETENTION: all32 D2 raw archives (45614898bytes) published at3c1573b68b5f51e002edd346fb15c8fadf227b68; localSHA256 matches raw manifest and committed blob sizes, origin/main points to archive commit. Original D1 preserved. PI closure pending.


## MG057 — PI closes MG10-D2 (2026-11-07)

Richard: “I approve closure for D2.” D2 CLOSED as COMPLETE DESCRIPTIVE PILOT. Registered result, independent checks and verified full retention accepted within their scope; all numerical/statistical/particle-interpretation limitations preserved. No rescoring of earlier studies and no new production registration or execution.


## MG058 — Register single-motif measurement qualification

PI instructed next step after D2 closure. P1 qualifies localization/reference matching and conservative sampled continuity on synthetic controls, before any particle/supply-intervention claim. Scope:68 frames,15 sequences,7 invalid-input assertions,300seconds total, one local worker. Constants/source/expectations frozen before generation/evaluation; independent dense/direct oracle. No new model trajectories, old-result rescoring or automatic next production.


## MG059 — Single-motif measurement qualification result

All90 registered fixture/API assertions pass, independent dense/direct oracle agrees on511 raw frame assessments and continuity classifications. QUALIFIED ON REGISTERED SYNTHETIC FIXTURES, 2.016197000ownedseconds versus300 cap. No expensive simulation or old-result rescoring. X-only/fixed-reference/noise/lineage limitations preserved; actual model-informed generalization remains a separate study. Raw arrays and decisions prepared for upload/verification; PI closure pending.

MG059 RETENTION: all16 raw P1 input archives (280039bytes) uploaded at111936587ce404f0e87c5ca55f7ad89dddf1aff3; source matrix/expectations/seeds and classifications retained. LocalSHA256 matches and committed blob sizes verified against pushed archive. P1 PI closure pending.


## MG060 — PI closes MG10-P1 (2026-11-07)

Richard approved P1 closure. CLOSED as QUALIFIED ON REGISTERED SYNTHETIC FIXTURES, preserving90 passed assertions,511 separately recomputed frame assessments and verified full archive. All scope/identity/noise limitations retained. No old-result rescoring or new production authorization inferred from closure.


## MG061 — Prepare model-noise observability calibration

P2 draft and candidate1024-path/128-job matrix prepared. Frozen P1 matcher plus per-field G/X/Y noise/projections; four primary H match-frequency caps, seeded retention descriptive. First preparation metadata collision preserved/corrected before repeat, unchanged model/statistical criteria. Independent preparation/state/metric/CP checks pass, no scientific production. Reviewable approval packet and measured forecast ready; PI protocol/budget decision required before final production registration.


## MG062 — PI approves and P2 is prospectively registered

Approved plan:1024 fresh individually seeded short model paths,128 jobs, four primary H endpoint match-cap tests and descriptive S/GXY observability. Preparation and full synthetic complete/incomplete pipeline passed independent verification. New4-core-hour cap, expected~16minutes/8workers; source/spec/reference/seeds and actual initial-state definitions frozen before production. Prior P1/D2 remain closed and unchanged. Full raw archives/progress and automatic report/recomputation retained; PI closure later.

MG062 P2 PRODUCTION LAUNCH: 2026-11-07T20:18:48.488238+00:00; registrationeaef80d pushed before scientific samples. Eight local owned workers active; saved complete paths/integer checks verified. Exact1024/128 denominator unchanged,4-hour cap/180seconds job; no result yet. Automatic report and separate audit follow any stop; full production upload and PI closure later.


## MG063 — P2 completed and independently recomputed

All1024 paths and128jobs complete;4 H match-frequency caps supported with0/128 each, upper3.3655%; seeded128/128 each descriptive. Independent4096 instrument frames, raw integer/source/state/seed and metric/exact-bound audit PASS. G global contrast is weaker than model spatial noise but planted known-phase ensemble projections remain informative; no blanket G-unobservable statement. No old-result rescoring/threshold adaptation. Full raw archive prepared; PI closure separate.

MG063 RETENTION: all128 raw P2 chunks/1024 paths (54412800bytes) published at9ce8e610538581276e285a15601fef861aca3596; localSHA256 matches and committed blob sizes verified against origin/main. Complete raw/event/input/metadata/evaluation/audit archive retained. PI closure pending.


## MG064 — PI closes MG10-P2 (2026-11-07)

Richard: “I sign off of P2.” P2 CLOSED with its registered short-horizon unseeded-match-cap outcome. Complete1024-path data, independent audit and verified full archive accepted within scope; no physical survival/stationarity/coupled-field equivalence implied. Earlier labels preserved; next scientific production remains unregistered.


## MG065 — Prepare temporal persistence and coupled averaging

P3 draft/preparation complete:64 H/S/two-step full-maintained paths toT25,20-unit sampled visibility and conditional/unconditional coupled-field averages. Episode/alignment/noise-correlation/no-reacquisition/partial-path/stream checks pass with separate direct interpolation/episode/raw reconstruction. Timing/initial-state provenance retained; new24-core-hour candidate requires PI approval. No scientific run or prior-study rescoring.


## MG066 — PI approves and P3 is prospectively registered

Approved P3 protocol/24-hour cap. Full and partial synthetic pipeline, actual quarter-stream equivalence to held development data and separate direct temporal/average/statistical audits pass.64 fresh individually seeded baseline paths, H/S/two steps,R16,T25/101frames; fixed0..20 averages all completed paths, conditional aligned S with initial20 continuity only. Exact registration/source/reference/streams/budgets/retention frozen; push before science. Automatic report/audit after stop and full upload/PI closure later. Earlier results unchanged.

MG066 P3 PRODUCTION START: 2026-11-07T21:37:04.563159+00:00; registrationd0e9fcc pushed before science. Eight owned local workers and quarter checkpoints verified, no controller stop.64 fixed paths,T25,24-hour aggregate cap/1800seconds path, no duplicate/restart. Automatic report/audit after stop; full archive/PI closure later.


## MG067 — P3 complete and independently inspected

64/64 complete, initial S20-unit events32/32 total and Hany20 events0/32 total, inference per-step R16 only. One fine S censors at22.75 after reaching20. Independent6464 frames/96 temporal averages/raw/statistical audit PASS. Temporal residual averaging improves G visibility~4.4-4.5, without sqrt81 or physical survival/lineage claims; all aligned eligibility/uncertainties retained. Full archive prepared for publication/verification. PI closure separate; no next production registered.

MG067 RETENTION: all64 complete P3 raw archives (146550116bytes),101 frames/path and full decision profiles/series published atef3b43aecefb2f244da456b47aaec7c0f3386d4b; SHA256 matches and committed blob sizes verified against origin/main. Original studies preserved; PI closure pending.


## MG068 — Add approved project summary and close P3 (2026-11-07)

Inserted the exact one-paragraph summary Richard approved into P3 results. His conditional P3 approval is satisfied; P3 CLOSED as COMPLETE TEMPORAL ESTIMATION. Recorded the standing one-paragraph project-contribution/limits summary requirement in protocol, agent instructions and handoff. No scientific criteria/data or old outcomes changed; no next production registered.


## MG069 — Prepare controlled resource-location intervention (2026-11-07)

Prepared MG10-E1 draft and froze/pushed preparation source at26a5627 before qualification. Source-rate gate, not deleted stock/edge: local and remote51-node A shutoffs, uniform.898 and unchanged control, with matched H backgrounds.32 fixed-T5 parents,128 forks,T20,8 balanced blocks/step and two steps; fixed-origin X primary interaction. Qualification/separate raw/statistical/fork audits PASS, all preparation data retained. Measured22.642613core-hour forecast/~2.83hours on8 fits proposed NEW36-core-hour cap; PI decision required. No scientific production or favorable-effect estimate. Final production implementation/pipeline audit and registration follow approval; P3 remains closed and unchanged.


## MG070 — PI approves E1 and authorizes production

Richard: “I approve the E1 study. Proceed running it.” Approved MG069 intervention design and36-core-hour cap unchanged. Final streamed driver, full/prefix/missing balanced-block analysis and separate raw/projection/statistical checker prepared; held-development replay and algebraic pipeline qualification precedes exact frozen registration/push/production. No new scientific samples yet. External review not requested.

MG070 FINAL QUALIFICATION/REGISTRATION: all20 held replays bit-identical; full and incomplete/prefix fixtures and separate audits pass, total preparation335.932328seconds<900. Exact32-parent/128-branch matrix and hashes registered, prior seeds unchanged; push registration before launch. Parent/branch aggregate charge, whole-batch failure stops and raw prefix retention qualified by proven owned controller and explicit pipeline checks.

MG070 E1 PRODUCTION START: 2026-11-08T03:04:22.395260+00:00; registration1249d0a pushed before samples. DriverPID24540, persistent execution session83731; eight active parent workers, all eight have advancing atomic quarter checkpoints and no stop. Automatic branch phase after all32 parents, then evaluation/separate audit. Initial hidden launcher did not create run state/ledger; retained empty launcher logs; foreground persistent driver is the verified active owner. No scientific path retry or duplicate launch after ledger creation. Use existing session/state to inspect, never rerun production command.


## MG071 — E1 complete, resource-location interaction supported (2026-11-08)

32 parents/128 branches complete, separate160 raw/field/fork/statistical checks pass. Both steps show negative registered seed-by-location primary intervals; all16 block contrasts negative. Report begins with required project-level contribution/limits paragraph. Local A delivery matters for seeded reference-field response beyond unseeded forcing, without particle death/physical-source/entanglement claims. Production19.808474core-hours under36; no errors/retry/criteria amendment. All full raw decision inputs retained; publish/verify archive next. PI closure separate.

MG071 PUBLICATION BLOCK: local commit334caca contains complete results and160 raw archives (290504596bytes). Pull was current. Auto-review rejected push to public origin/main, stating run approval did not clearly authorize public exposure of full possibly unpublished dataset. No push retry/workaround; PI permission requested. Local audit/data complete and preserved, remote publication/verification pending.

MG071 PUBLICATION AUTHORIZATION (2026-11-08): Richard: “i approve publishing the results of E1,” answering the explicit request to publish E1 results and full raw archive to the public ModelG repository. Approved public publication of local result/archive commit334caca and supporting records; previous automatic-review block resolved by explicit PI permission. Scientific closure remains pending. Publish then verify all160 raw files against origin/main.

MG071 RETENTION VERIFIED: explicit PI publication approval applied; results/full160 raw archives (290504596bytes) published at18adea684d2025e2f7433044ba9ec67d165cc5d1. Every origin/main raw blob matches audited SHA256 and local bytes. Earlier approval rejection preserved, resolved by explicit authorization. Scientific closure pending.


## MG072 — PI closes MG10-E1 (2026-11-08)

MG072 PI CLOSURE (2026-11-08): Richard: “I sign off”, following verified publication of E1 results/full archive. MG10-E1 CLOSED as COMPLETE INTERVENTION ESTIMATION.32/32 parents and128/128 branches complete, separate160-file raw/fork/projection/statistical audit PASS; all160 published archives match audited hashes. Registered direction-consistent model-pilot resource-location interaction retained; small-n approximate intervals, fixed-origin measurement, finite horizon and physical-mortality/reservoir-origin/entanglement limits unchanged. No criteria/data rescore, no new scientific production authorized.


## MG073 — Develop E2 reversible supply study (2026-11-08)

MG073 E2 DEVELOPMENT: Richard authorized developing E2. DRAFT reversible local A shutoff/restoration design,32 fresh fixed-T5 parents/128 T25 forks, pulse ends post5, four arms/H-S/two steps/n8 balanced blocks. Primary local-restore versus continuous-off seed interaction; secondary movement/shape/history diagnostics, no physical revival/quantum memory claim. Frozen preparation048d705 preceded8 short qualification/16 timing paths.11 assertions/7 field fixtures and24 raw/dense/statistical audits PASS, no failure/retry/tuning. Controlled157.545613seconds, total allowance162.907576<900. Forecast29.558218core-hours/~3h42 on8, proposed NEW40core-hour cap below30 reserve; PI final design/budget/public production retention pending. No science or old rescoring; full production implementation/raw-pipeline qualification and registration remains after approval.


## MG074 — PI approves E2 protocol, budget and public retention (2026-11-08)

Richard: “I approve”, answering the request for E2 design,40-core-hour local cap and public publication of results/full raw archive. Approved numerical design unchanged. Final atomic-writer held replays and complete/missing/prefix/empty raw-file audits prepared before exact registration/push/production. No E2 scientific sample yet; external review not requested.

MG074 FINAL QUALIFICATION/REGISTRATION: MG074 E2 REGISTRATION: Approved32-parent/128-branch/T25 recovery protocol and40-core-hour cap/public result+full raw retention registered prospectively.20 held writer streams match every snapshot/statistic; complete128/missing127/prefix127 and empty-parent-stop raw-file pipelines pass separate integer/dense movement/profile/paired/Student audits. Missing primary counts7/8, empty0/0; no fabrication. Total preparation336.239390seconds<900.320 candidate/warm values unchanged, unique/disjoint with40 development and previous registrations. Forecast29.558218core-hours<=30 reserve.144000aggregate worker-second cap,8workers/threads1,600parent/2400branch timeouts, all parents before branches, prior phase charge carried over, whole-batch stop/no retry. Commit/push exact registration before science; no old result rescored.

MG074 E2 PRODUCTION START: 2026-11-08T14:50:31.580704+00:00 (Nov8~07:50 Mountain), registratione2fc7f4 published before science. Persistent execution session97858, driverPID35028. Eight active parent workers all saving advancing quarter checkpoints; no stop. Fixed32-parent phase then128 branches under aggregate40core-hours, automatic evaluator/separate audit after stop. No duplicate/retry; inspect existing session/state rather than rerun. Public results/full raw retention already PI-approved; archive verification and scientific closure after completion.


## MG075 — Preserve incomplete E2 and diagnose forecast stop

Richard authorized diagnosis/continuation planning (“go ahead”). Original32/1/7prefix/120unstarted record preserved,5.452379core-hours,40-file audit PASS, zero complete blocks. Stop was registered forecast policy, not exhaustion/crash; no adaptive-step reductions and~95.5% completed-branch time in simulation. New held-seed performance probes capped300worker-seconds prospectively frozen; no production restart or new ceiling authorized.

MG075 DIAGNOSIS COMPLETE: MG075 DIAGNOSIS/AMENDMENT DRAFT: Original incomplete E2 preserved and all40 raw published hashes verified.9 held old/new timing pairs bit-identical, E2/old ratio.932..1.017; current CPU-bound probes do not establish past sustained-slowdown hardware cause. New admission-only forecast stop drains active jobs; hard cap/timeout/failure still terminate, synthetic tests PASS. Exact-prefix guard replays held pulse and rejects altered input.137.006255diagnosticseconds<300,9.494611controlseconds<30. Original19628.563985worker-seconds retained; proposed same-seed7replays+120unstarted, reuse32parents/1completed. Newforecast95.763055aggregate core-hours/~11h17remaining8workers; proposed144aggregate cap/5400seconds job DRAFT pending PI. No scientific restart, no rescore/new seed/denominator change. Final continuation orchestration/mixed-provenance audit must qualify before registration/science.


## MG076 — PI approves E2-C1 continuation amendment

Richard: “yes, i approve”, responding to144core-hour aggregate amendment including original5.452379charge,5400/job and forecast admission-drain policy. Same32-parent/128-branch scientific matrix/endpoints/seeds; reuse32parents+1complete, exact-replay7prefixes,120unstarted. New namespace preserves original incomplete report/data. Final mixed-provenance/prefix/aggregate-ledger synthetic pipeline must pass and registration publish before scientific replay/new branch. Public E2 retention approval persists; external review not requested.

MG076 QUALIFIED/REGISTERED: MG076 E2-C1 REGISTRATION: PI approved144total core-hour continuation cap inclusive of original19628.563984799664worker-seconds,5400/job,forecast admission-drain.32reused parents/1reused complete branch,7exact replay+120unstarted original seeds; unchanged128matrix/endpoints/arms/windows/n8 inference. Final mixed-provenance full/missing/prefix160-file/direct-statistical/source/fork/replay/aggregate-charge pipeline PASS,77.777412seconds<300; original fake tamper rejected. Original real E2 source/data/report preserved. Exact127job registration/source/reference/original-file manifest committed/pushed before science.518400cap,8workers/threads1,.25poll,hard cap/timeout/failure still terminate,no retry. Public E2 results/full-raw publication already approved. Forecast95.763055aggregate/~11h17remaining, uncertain fine throughput. No original rescoring/new seed/outcome selection.

MG076 E2-C1 START: 2026-11-08T16:07:54.373455+00:00 (Nov8~09:07 Mountain), registrationeacd1ac published before scientific replay/new branch. DriverPID17180, persistent session46751.32parents+1complete branch copied byte-exact, original files unchanged. Eight active workers (7replays+1unstarted fine), all advancing quarter checkpoints; prefix guards matching available histories so far, not yet full old-prefix completion. No failure/admission/hard stop. Ledger prior19628.563984799664 preserved under518400total. Automatic combined evaluator/audit after stop; original E2 remains incomplete. Public full-data retention already approved; no automatic restart.


## MG077 — E2 continuation complete and audited

MG077 E2-C1 RESULT: COMPLETE RESTORATION ESTIMATION,32parents/128branches;127continuation jobs complete with no admission/hard stop/timeout/failure.7exact prefix replay proofs,160raw/fork/source/direct-movement/statistical audits PASS. Primary coarse+.102771 nominal97.5%CI[.025840,.179701], fine+.140728 CI[.063955,.217502], n8/step: registered direction-consistent restoration response. Restored fixed-origin seed excess.11/.15 vs sham.95/.93; raw restoration largely unseeded background, not full recovery/physical identity/quantum memory. Aggregate28.363830core-hours incloriginal5.452379, continuation174.610150minutes, finishedNov8~12:02Mountain. Conservative11h17forecast exceeded actual runtime; earlier throughput variation cause remains unproven. Original stopped E2 not rescored/overwritten. Public full archive publication/verification follows, already PI-approved; scientific closure pending.

MG077 RETENTION: all160 combined E2-C1 raw archives (356072369bytes) published at7fabb8c51b6781a7be33f5ce50eb9c332304b78b; every published blob matches audited SHA256/local bytes.7replay proofs retained and all87 original-manifest files unchanged. Public publication explicitly approved; PI scientific closure pending.


## MG078 — PI closes MG10-E2-C1 (2026-11-08)

MG078 PI CLOSURE (2026-11-08): Richard: “i sign off”, following completed E2-C1 result and verified public full archive. MG10-E2-C1 CLOSED as COMPLETE RESTORATION ESTIMATION.32parents/128branches,127continuation jobs,7exact original-prefix replay proofs and160-file raw/field/statistical audit PASS; full published raw hashes verified. Registered modest direction-consistent restoration interaction retained; full-recovery/equivalence/physical identity/quantum-memory limitations unchanged. Aggregate28.363830core-hours including original stopped attempt preserved. Original E2 remains INCOMPLETE/immutable; no frozen criteria/data rescored, no next scientific production authorized.


## MG079 — Build fixed Environmental Network qualification

MG079 ENVIRONMENTAL NETWORK DEVELOPMENT: PI authorized building test. MG11-N1 DRAFT fixed symmetric/equal-volume A/B graph and finite closed contacts built, no active-lattice change/replenishment. Exact small-system SSA and continuum/linear limits;8grouped assertions/15physical raw plus deterministic direct checks PASS, source2cf1358 frozen before development, no failure/retry/tuning. Planned32jobs/131072paths,3analytic endpoint laws+independent full-tiny SSA population laws, frozen DKW margins/families.64prospective+15dev streams distinct/disjoint. Development19.686674worker-seconds, total allowance24.906960<900. Forecast.090118core-hours/~1minute8workers, proposed NEW1hour cap and public results/full raw retention pending PI;195035136byte array bound. No scientific production/old rescore or particle/full-lattice stochastic/native-vacuum/entanglement claim; final production/pooling/incomplete audit registration remains after approval.


## MG080 — PI approves Environmental Network qualification

Richard: “yes, i approve of all those items”, confirming fixed network model,criteria,1-core-hour local cap and public results/full raw archive. Scientific matrix/thresholds/seeds unchanged from MG079. Final atomic128pathstream/writer, pooling, sortedCDF/joint-generator audit and pass/fail/missing/prefix/empty algebraic pipeline qualification prepared before registration/sampling. No scientific N1 sample drawn yet. External review not requested.

MG080 FINAL QUALIFICATION/REGISTRATION: MG080 N1 APPROVAL/REGISTRATION: Richard approved all model/criteria/1hour/public retention items.32canonical jobs/131072paths plus32warmups,64candidate streams unchanged; 6analytic bounds D+sqrt(log24000/(2*32768))<=.03 and28two-sample bounds D+2sqrt(log224000/(2*16384))<=.08, combined>=99.9%nominal DKW family coverage. Formal evidence/source unchanged. Five held streaming/writer replays exact; full131072 passing, deliberate valid-accounting distribution failure, missing/prefix127104 and empty algebraic pipelines/independent audits PASS,18.769082seconds<300. No production seed drawn during qualification.3600worker-second cap,180/job,8threads1/.25poll,forecast324.424602seconds,2GBfree. Frozen registration/code/criteria/seeds/formal/reference hashes published before samples; no old rescore or broad particle/vacuum/quantum claim.


## MG081 — Fixed Environmental Network qualification complete

MG081 N1 RESULT: QUALIFIED WITHIN REGISTERED N1 SCOPE,32/32jobs,131072scientific paths+32excluded recorded warmups.6/6analytic bounds max.014669808<.03;28/28tiny two-implementation population bounds max.048666983<.08. Formal8grouped checks/source preserved. Separate64raw accounting/graph/seed/joint-generator/sortedCDF/statistical/owned-charge audit PASS. No failure/timeout/cap/admission pause/replacement/criteria amendment.105.096124worker-seconds=.029193core-hours<1cap,18.769seconds elapsed through audit. Full raw/config/counters/probability/audit archives prepared for approved public publication/hash verification. Small endpoint/formal/linear scope only, no full nonlinear particle/long-lattice/vacuum/energy/quantum claim. Old MG10 records unchanged; PI closure pending, no next production automatic.

MG081 RETENTION: all64main/warm raw archives (1932265bytes) published atf0e14d440a53ffb79a752a4f74493a043413dec5; each origin/main blob matches audited SHA256/local bytes. Public publication explicitly approved, all scientific samples/config/graph/edge histories/probability/audit records retained. PI scientific closure pending.


## MG082 — PI closes MG11-N1 (2026-11-08)

MG082 PI CLOSURE (2026-11-08): Richard: “i sign off on MG11-N1”. MG11-N1 CLOSED as QUALIFIED WITHIN REGISTERED N1 SCOPE.32/32jobs,131072scientific paths+32excluded recorded warmups,6/6analytic and28/28two-sample criteria supported,64raw separate audit and full published hash-verified archive accepted. Formal/linear/small-system endpoint limits remain unchanged; no nonlinear full-lattice particle, vacuum/energy/entanglement claim. Source/data/registered criteria/earlier outcomes unchanged, no next scientific production authorized.


## MG083 — Build matched-connectivity geometry-supply tests

MG083 GEOMETRY-SUPPLY TEST DEVELOPMENT: PI authorized writing tests after N1closure. MG11-G1DRAFT matched6node prism/K3,3 graphs,9edges/degree3/escape1/source0contact3/distance1 and stock/contact inputs equal; connectivity differs. Frozen N1 exact SSA reused. Primary BcaptureT4 difference,65536/cell/4cells/64jobs/262144paths, conservative exact-binomial interval;16endpoint DKWmodel checks separately. Analytic predictions/rationale recorded pre-MC, no outcome-selected R/graph/time.7formal groups/20raw and independent strict-state/ledger/BFS/generator-action/beta audits PASS;128candidate/18devstream distinct. No failure/tune/retry. Frozen0f52ea7 before development;13.759981controlledseconds,20.405549total allowance<900. Forecast935.620429seconds=.259895core-hours/~2minutes8workers with IOquadratic scaling, proposed1hourcap/public full archive pending PI.1.635GBmain raw bound/chunkwise analysis. No scientific run, no particle/fullgeometry/quantum claim. Final production/evaluator/incomplete pipeline and registration after approval.


## MG084 — MG11-G1 approval and prospective execution (2026-11-08)

PI: “i approve all of them. Go ahead witht he run.” Approved MG083 design/criteria, one core-hour production cap and public full raw/results publication. Frozen spec/MG11G1_proposal.md plus spec/MG11G1_prerun.md govern. Final pipeline qualification (300-second separate allowance, held streams/algebraic fixtures only) precedes registration and scientific draws. 64jobs/262144paths, max8workers/thread1, no retry. Primary B T4 contrast and sixteen endpoint-law criteria unchanged. Commit/push source and qualification plan before tests; exact registration pushed before production. Scientific closure remains PI decision.


MG084 final qualification PASS: four held writer replays and complete/missing-B/missing-A/prefix/empty/valid-accounting law-failure audits passed in 66.623406 seconds (<300), no failures/retries. Full synthetic262144/prefix258176 algebraic paths are not scientific data. Exact approved registration saved in results/MG11G1/registration.json; source/job hashes frozen. Publish registration before production.


## MG085 — MG11-G1 results (2026-11-08), PI closure pending

COMPLETE CONNECTIVITY-SUPPLY ESTIMATION.64/64jobs,262144scientific paths+64excluded warmups. Primary B T4 bipartite-prism+.015182495 (1.518250percentage points), registered interval[.007453630,.022909648] above0; analytical prediction+.011529941.16/16endpoint CDF criteria supported, maxupper.015339591<.03. Separate128raw/state/source/seed/generator/CDF/beta/ledger audit PASS. No failures/timeouts/stops/retries/amendments. 353.770201worker-seconds=0.098270core-hours<1cap; end-to-end61.825802seconds. No active workers. Supports selected classical connectivity-dependent supply, not particles/native vacuum/energy/entanglement. Full raw/results publication approved, post-push verification follows; report results/MG11G1/MG11G1_results.md. PI closure separate.


MG085 retention verified: all128 production main/warm raw archives (9205383bytes) published at6df002421f1e607c39f92e4d58553b268630da63; every origin/main blob matches audited SHA256 and local bytes. All159 preparation/qualification binary inputs/fixtures also match published blobs. See results/MG11G1/retention_verified.json. No active scientific workers; PI scientific closure pending.


## MG086 — MG11-G1 PI closure (2026-11-08)

Richard: “i sign off on MG11-G1.” CLOSED as COMPLETE CONNECTIVITY-SUPPLY ESTIMATION, with the registered positive B contrast supported. Accepted MG085 results,16/16 endpoint-law checks,128-archive separate audit and verified full public retention. Frozen design, counts, estimates and criteria unchanged. Scope remains selected fixed classical network delivery probes; particle persistence, emergent vacuum, physical energy and entanglement not established. No next production run authorized by this closure.


## MG087 — MG11-G2 Step1 build and deterministic qualification (2026-11-08)

PI authorized: “yes, I like that approach. Go ahead and build the Step 1.” Reuse frozen N1 full nonlinear continuum law, add explicit reaction/active/environment/import ledgers, qualify zero-transport legacy mapping, maintained supplies, nonlinear shared-finite-bath fast mixing, and finite-capacity recovery. Fixed small nonuniform fields,9coupled trajectories plus separately assembled Radau references; no stochastic draws, particle starts or geometry-effect production. Criteria/configurations and300-second single-process qualification ceiling in spec/MG11G2_step1.md. Commit/push source/plan before integration. No retries/criterion tuning; retain failures. Historical kernels frozen. External review not requested; same-author separate recomputation mandatory. Downstream pattern-response design/budget remains separate.


### MG087 qualification limitation and MG087a separate extension (2026-11-08)

Original Step1 qualification FAILED its frozen mixing criterion: kappa1000 error.005213697869783119>.005. All9coupled trajectories retained; formal drift/ledger/zero-resource/boundary/config checks pass, capacity errors(.00333532853,.000334077976,.0000334133652) satisfy frozen recovery test. The initial driver stopped before separate audit; original data/source/criteria frozen, no retry or tolerance change. Complete the planned separate audit using a disclosed diagnostic path, preserving failure. New MG087a prospectively fixes ONE additional kappa10000 point, unchanged fixture/solver/error target, to check whether stronger mixing reaches the same limit. No further automatic extension if it fails. Scoped to authorized deterministic Step1 build qualification, no stochastic/particle production. Separate spec/hash snapshot and all failure evidence pushed before follow-up integration.


### MG087/MG087a Step1 implementation outcome (2026-11-08)

Build complete, coupling qualified on tested small deterministic systems with original mixing-range failure retained. All9original trajectories pass separate equation/Radau/ledger/config audit; pointwise legacy/maintained drift, zero-resource/boundary and5invalid-input tests pass. Capacity C1/10/100 error decreases to3.34133652e-5<.001. Original kappa1000error.00521369787>.005 remains FAILED. Prospectively frozen additional kappa10000 gives0.000524351935444<same.005 target and passes numerical/accounting audit; no original rerun or tolerance change. Ten total coupled trajectories, no random sampling. Max normalized population/reference discrepancy2.10711192494e-07<=2e-6, ledgers/balances<1e-8. Follow-up audit4.6499945seconds/extra point.68508seconds, one thread. Full raw/source/failure evidence published and verification follows. Report results/MG11G2/MG11G2_step1_results.md. PI acceptance pending; no pattern/stochastic production authorized.


Retention verified: all22 original/follow-up raw/reference archives (50258bytes) published at1e10a7d6f46f5dfc494477ee7f4669fb02715bed; every blob matches its manifest SHA256/local bytes. Original MG087 source/data hashes unchanged, all pre-existing source files unchanged from signed-off G1 baseline. See retention_verified.json.


## MG088 — deterministic Pattern-Response preparation (2026-11-08)

Richard: “Go ahead with the Pattern-Response study.” Prepare concrete matched-environment design/implementation/qualification/timing now; final numerical criteria/budget/public scientific retention approval under RESEARCH_PROTOCOL.md before production. Step1 combined MG087/MG087a supports progression without rescoring its failed original mixing range. Fixed N500 active ring; equal-stock/degree/hop-rate LOCAL vs CROSS_RING environmental graphs; isolated/maintained controls; C1/C100; H/S; two BDF settings (28jobs). T25scientific primary C1late20..25seed-excess X contrast, deterministic numerical refinement guard, no CI/particle-death claim. Development28actual-writer T.25paths plus14separate population references/algebraic complete-missing-prefix-empty/wrong-primary tests,900elapsed ceiling/thread1, source/plan push BEFORE integrations. Candidate production cap2core-hours, forecast must fit1.5core-hour reserve; max8workers. Full spec spec/MG11G2_pattern_response.md. No scientific production registered yet.


### MG088 original preparation failure and MG088b numerical correction

Original28short actual-writer paths retained; separate Radau comparison failed ISOLATED_C1 S base accumulated Z2.36016656e-6>2e-6 (fine2.15533333e-8). No T25scientific run. Original source/plan/spec/raw/failure frozen. Versioned all six new pattern files as_v2, outputs pattern_response_v2; tighten base rtol/atol1e-9/1e-11 andfine1e-11/1e-13, no model/observable/threshold changes. Reserve30seconds for original failed development, remaining870seconds <=900aggregate allowance. All original hashes frozen in v2plan; publish before corrected qualification. No third automatic revision. Candidate production remains unregistered pending concrete measured criteria/budget/retention approval.


### MG088b Pattern-Response readiness (2026-11-08)

Corrected development PASS:28T.25actual-writer paths,14separate short Radau references,28raw ledgers,56Jacobian columns,graph/reference/legacy/rotation/measurement and complete/missing/prefix/empty/wrong-primary/invariance pipelines. Max normalized population discrepancy1.10796506192e-07<=2e-6; no threshold changes. Original failure/source/raw preserved and original hashes unchanged. New preparation29.717779seconds+30reserved<900. Forecast0.462538core-hours/~3.47minutes8workers, within75%reserve of proposed2core-hour cap. Concrete28job/T25design/criteria/cap/full public scientific raw retention pending PI approval per RESEARCH_PROTOCOL workflow3. No production registration/draw. Review packet results/MG11G2/pattern_response_v2/MG11G2_pattern_preparation.md; spec_v2.


Readiness retention verified: all103 original/corrected development binary inputs/fixtures (12448415bytes) published at7d8c3b98459cb1f3b6b7651ed5f828110acd347b, every blob matches manifest SHA256/local bytes. Original failed-version hashes unchanged; corrected source/pipeline hashes valid. Full T25scientific run remains unregistered/unstarted pending final concrete PI execution approval.


## MG089 — Pattern-Response registered for production (2026-11-08)

Richard: “Looks good. Approved and continue.” Approves the concrete MG088b 28-job design, frozen numerical/response criteria,2core-hour/7200worker-second cap, max8localworkers/thread1 and full public scientific raw/results archive. No external review requested. registration.json records exact jobs/source/reference/preparation/pipeline hashes and environment. Forecast.462538core-hours/~3.47minutes8workers plus analysis, within75%reserve. Timeout600/job, admission-drain controller,.25poll,no retries; publish registration before first T25integration. Original coarse development failure and Step1 mixing-range failure remain frozen. Result PI closure separate; no stochastic follow-on authorized.


## MG090 — deterministic Pattern-Response results (2026-11-08 local), PI closure pending

COMPLETE DETERMINISTIC PATTERN-RESPONSE ESTIMATION; RESOLVED CONNECTIVITY CONTRAST.28/28jobs,101frames each, separate28raw audit PASS. Primary C1X CROSS_RING-LOCAL+.054506911595 fine/+.054506880909 base; refinement3.06858225e-8, threshold.01; max active refinement1.46383350e-7<1e-4. Maintained continuity25passes. C1local reference match19vs cross16.75, both lost before primary20..25window: larger projection NOT better preservation. C100secondary contrast+.000676948011, all reference continuity25. A uniform graph-independent; Bredistribution mechanism. No failures/stops/timeouts/retries/amendments. Cost0.039478305core-hours<2, elapsed29.463478seconds through audit, finishedNov8~17:30Mountain. Full scientific raw/results publication approved; post-push verification follows. Report results/MG11G2/pattern_response_v2/MG11G2_pattern_results.md. No next production authorized; PI scientific closure separate.


MG090 scientific retention verified: all28raw archives (42096175bytes) published atcbecb8ac013b9cfd6f6548f4496110b27cab2b04; every origin/main binary blob matches its audited SHA256/local bytes, including byte-identical figure. Text metadata matches after Git line-ending normalization and JSON parsing (initial raw byte comparison of text correctly exposed CRLF/LF difference, no scientific raw mismatch). Registered sources/preparation/pipeline and original failed-version input hashes unchanged. See production_retention_verified.json, including text-hash reproduction note. No active scientific workers; PI closure pending.


## MG091 — PI closure of MG11-G2 Pattern-Response (2026-11-08 local)

Richard: “I sign off. Go ahead with the next stage.” CLOSED as COMPLETE DETERMINISTIC PATTERN-RESPONSE ESTIMATION, RESOLVED CONNECTIVITY CONTRAST under registered criteria. Accepts MG090 numerical/raw audit and verified full public archive, with the explicit limit that higher projection did not establish longer reference-pattern persistence. Original preparation and Step1 failures remain preserved, not rescored. Next-stage preparation authorized; its concrete modeling choices/criteria/measured budget will be recorded before production. No stochastic result inferred.


## MG092 — B-only transport switch preparation (2026-11-08 local)

PI authorized next stage after MG091 closure. Prepare conditional mechanism experiment: exact C1G2parents atT10 (8donors, both origins/H/S/settings) each fork LOCAL/CROSS_RING/OFF Btransport only; Agraph/reactions/stocks/output stores unchanged, post counters reset explicitly.24postT15branches, same finalabsolute25, primary equal-origin average late20..25cross-local seed-excessX; deterministic refinement/sham guards, no CI/quantum-memory claim. Development24shortpost.25writer branches+24separate Radau references, graph/fork/ledger/Jacobian/sham/OFF and synthetic full/failed-sham/missing/prefix/empty/wrong-primary tests,900elapsed/thread1. Source/spec/plan/donor hashes pushed before any integrations. Candidate1core-hour cap/max8workers, concrete measured design/criteria/budget/public retention return to PI before science under RESEARCH_PROTOCOL workflow3. No previous outcomes rescored.


### MG092 B-switch readiness

All development checks PASS, no revision/retry:24post.25writer paths/24separate Radau refs,8short shams, exact forks/store preservation,56Jacobian columns, B-OFF/A-invariance/two-graph ledgers and full/failed-sham/missing/prefix/empty/wrong-primary/estimand fixtures. Donor8hashes unchanged. Max population discrepancy4.00834599193e-09, short sham3.20905702011e-10<2e-6. Elapsed36.283009seconds<900. Forecast0.307928core-hours/~2.31minutes8workers, within reserve of proposed1core-hour cap. Concrete24post15branch design/criteria/cap/public scientific raw retention awaiting PI execution approval per workflow3. No full scientific branch registered/run. Review results/MG11G3/MG11G3_preparation.md.


MG092 readiness retention verified: all81binary decision inputs (73newdevelopment/reference/fixture files+8reused closed donors;28442222bytes) published at070b30d94f629674e17e63837e992cfd73f53ba2; every published blob matches manifest SHA256/local bytes. Plans/results/pipeline/environment JSON semantically match after Git line-ending normalization; source/donor hashes unchanged. No full scientific branch run. See retention_verified.json.


## MG093 — B-transport switch study registered (2026-11-08 local)

PI: “approved”, in reply to concrete24branch design,1core-hour cap and full public results/raw archive approval question. Fixed MG092 source/model/criteria/donor hashes;24postT15branches from8savedT10parents, BLOCAL/CROSS/OFF only, A/reactions/stores held fixed within donor. Max8workers/thread1,3600worker-second ceiling,300/job,.25poll, forecast admission/drain,no retries. Forecast.307928core-hours within75%reserve. Publication authorized; push registration before first scientific branch. External review not requested. Scientific closure/any next study separate.


## MG094 — MG11-G3 B-switch results (2026-11-08 local), PI closure pending

COMPLETE B-TRANSPORT SWITCH ESTIMATION; RESOLVED ONGOING B-TRANSPORT CONTRAST.24/24branches/61frames,24raw audit,12refinements/8shams PASS. Primary equal-origin X+.033337305924fine/+.033337293927base, refinement1.19968889e-8, fixed threshold.01. Origin-specificX+.011288177799/+.055386434048descriptive. Maximum active refinement8.91622542e-8 and sham8.46935690e-8<1e-4. OFF increases Xprojection+.153782145729/+.090054386698vs sham but not preserved reference. Cross-prepared→LOCALmatches19vs crosssham16.75; all finite records lose reference before late20..25, no physical lifetime/quantum claim. No failures/stops/timeouts/retries/amendments. 0.027857858core-hours<1cap,24.617782seconds through audit, finishedNov8~18:29Mountain. Full public scientific raw/results approved, publication/verification follows;8donors unchanged. Report results/MG11G3/MG11G3_results.md. PI closure separate; no next science authorized.


MG094 scientific retention verified: all24new branch archives (31326153bytes) and8registered closed donor archives published at40ee106f6a6b6c51c8d591ceea07a5a7334ecedc, every binary blob matches audited/registered SHA256 and local bytes. Figure byte-identical; metadata JSON equal after Git line-ending normalization; source/qualification/donor hashes unchanged. See production_retention_verified.json. No active scientific workers; PI closure pending.


## MG095 — MG11-G3 PI closure (2026-11-08 local)

Richard: “I sign off on G3. What's the next step?” CLOSED as COMPLETE B-TRANSPORT SWITCH ESTIMATION; RESOLVED ONGOING B-TRANSPORT CONTRAST under the frozen criteria. Accepts MG09424branch/12refinement/8sham/separate raw audit and verified scientific/donor retention. Conditional deterministic scope only; no better physical-particle persistence, stochastic robustness, hidden-state memory or quantum claim. Source/data/criteria unchanged. Next-stage recommendation requested, not production authorization.


## MG096 — MG11-G4 Step1 stochastic coupling build/development

PI: “yes, go ahead with step 1.” New bounded common-snapshot finite-resource network leap plus two-graph exact Numba/Python references, frozen old sources untouched. Joint B/X reservations and conditional transport exits avoid duplicated consumption/hopping; arrivals/products wait next numerical step, no fundamental-tick or new R4catalyst law. CSRtransport O(N+E), int64 domain/count/step guards. Prospective source/spec/seed/development plan pushed BEFORE any draws.900elapsed single-process/thread1 development: zero-hop/same-graph byte regressions, one-step arrival/conflict and1e8binomial probes,15actual-writer128path timings,3held stream replays, four.002N500benchmarks, synthetic complete/missing/prefix/empty/wrong-CDF/bad-ledger audits.

Candidate full numerical qualification:3fixtures*5methods*R16384=245760paths,120jobs+120excluded warmups.196formal finest-vs-exact/exact-vs-independent endpoint CDFbounds margin.08/alpha.0005; coarse/medium descriptive. Candidate2core-hour/8worker cap subject to measured75%reserve, concrete PI design/budget/public retention approval before that larger batch. No full pattern pilot. Protocol spec/MG11G4_step1.md. No source/threshold tuning/retry; preserve/version failures.


### MG096 stochastic coupling readiness

Development PASS without revisions/retries in39.310529seconds<900:8zero-hop legacy regressions/2common N1exact regressions;4096arrival/conflict and4096large-count binomial samples (upper.048663222494<.08),12invalid inputs, conditional-mean algebra,15writer128path timings,3held prefix identities,4N500microbenchmarks, full/missing/prefix/empty/wrong-CDF/bad-ledger pipeline.302candidate/devseeds distinct/disjoint from4901prior. Formal R16384qualification not run, no pattern pilot. Forecast120jobs/245760paths0.522849core-hours/~3.92minutes8workers within reserve of proposed2core-hour cap. N500microbenchmark linearT25extrapolation~10..17minutes/path is not a long-run guarantee. Concrete design/criteria/cap/public scientific retention awaits PI per workflow3. Review results/MG11G4/MG11G4_preparation.md.


MG096 preparation retention verified: all321binary inputs (320newdevelopment/reference/fixture archives+1existing cached seed profile,3730358bytes) published atcdc73e5ff23d7a64b932228ac09072cb41ebda13, every blob SHA256/local-byte match verified. One intentionally corrupt ledger is labeled and retained as a rejection fixture, not a failed scientific path. Metadata JSON equal after Git line-ending normalization, qualified sources unchanged. Formal R16384qualification and full pattern pilot not run. See retention_verified.json.


## MG097 — G4 Step1 formal stochastic qualification registered

PI: “i approve all. proceed.” Approves concrete120job/245760main-path matrix,196fixed marginal CDFcriteria,2core-hour/7200worker-second cap and full public scientific results/raw retention. Exact frozen source/spec/seeds/preparation/pipeline hashes recorded in results/MG11G4/registration.json, with versions/environment.120separate warmups excluded. Max8localworkers/librarythreads1,180seconds/job,.25poll, forecast admission/drain/no retry; measured.522849core-hours within5400reserve. Publish registration before any9603scientific draws. Optional external review not requested. No pattern pilot or fundamental-tick/physical/quantum claim authorized. PI scientific closure later.


## MG098 — MG11-G4 Step1 stochastic qualification results (2026-11-08 local)

QUALIFIED WITHIN REGISTERED SMALL-SYSTEM STOCHASTIC SCOPE.120/120jobs,245760main paths+120excluded warmups;196/196formal marginal CDFbounds and240raw separate audit PASS. Largest observedD.015197753906, upper.056927491867<.08, epsilon_each.020864868980; nominal99.95%for stated196bounds (not full joint/temporal/large-lattice law). Finest maximumtau.00025fixed before sampling, coarse/medium descriptive only. No failure/timeout/forecast/capstop/retry/seedreplacement/amendment. 0.186389931core-hours<2cap, elapsed94.206208seconds through audit, finishedNov8~19:27Mountain, no workers. No full pattern pilot. Full public scientific raw/results approved, publication/hash verification follows. Report results/MG11G4/MG11G4_results.md. PI scientific closure separate.


MG098 scientific retention verified: all240main/warm archives (15146501bytes) published at8054121ac74b70bbb5a782bedc3c8a8c70e943de; every blob matches audited SHA256/local bytes. Figure byte-identical; metadata semantically equal after Git line-ending normalization, registered source/preparation/pipeline/reference hashes unchanged. See production_retention_verified.json and its text-hash reproduction note. No active scientific workers; PI closure pending, no pattern pilot authorized.


## MG099 — MG11-G4 Step1 PI closure

Richard: “I sign off on G4. Proceed to next step”. CLOSED as QUALIFIED WITHIN REGISTERED SMALL-SYSTEM STOCHASTIC SCOPE. Accepts MG098120jobs/245760main paths,196marginal CDFbounds,240raw audit and verified public archive. Scope stays fixed small short endpoint laws, not full500node pattern accuracy or physical/quantum claims. Source/data/criteria unchanged. Next-stage pattern-pilot preparation authorized; its concrete model/statistics/measured budget/public retention approval remains before production under protocol.

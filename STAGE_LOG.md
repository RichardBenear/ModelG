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

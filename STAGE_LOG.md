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

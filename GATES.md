# Gates

Each gate's prediction, pass criterion and null control are filled in and
committed **before** the run that tests it. Status: PENDING / PASS / FAIL / REVISED.

| Gate | Phase | Question | Status |
|---|---|---|---|
| MG1 | 1 | Continuum PDE reproduces a published stationary Model G soliton | PASS by PI exception (registered FAIL on 1D stationarity preserved) |
| MG2 | 2 | Stochastic token rules converge to the same-lattice mean-field solution as Ω → ∞ | **FAIL** as registered (C3, C4); convergence criteria C1, C6 PASS — see MG2 result |
| MG3 | 2 | Particles nucleate spontaneously from intrinsic fluctuations at finite Ω | **FAIL** as registered (C2); C1, C3, C4, V5 PASS — see MG3 result |
| MG3b | 2 | Follow-up to MG3: time-step check, nucleation rates with more events, calibrated uncertainty (estimation study; no pass/fail) | **ACCEPTED 2026-10-30** (estimation study, no verdict; PI sign-off): Part 1 EQUIVALENT; S = 0.064 (approximate 95% 0.047–0.082) — see MG3b result |
| MG4 | 3 | Particles survive irregular connectivity; graph Laplacian ≈ ∇² at soliton scale (2D frozen graphs, noise-free) | **FAIL — CLOSED 2026-10-31** (every D and R run multiplied; T control passed; ChatGPT accepted; PI signed off) |
| MG4c | 3 | Follow-up to MG4: single-particle window in b on unweighted Delaunay graphs; cotangent-operator diagnostic (estimation study; no pass/fail) | **RESULT W1 + O1** (2026-10-31; 102 runs, 11.33 core-h); ChatGPT reproduced (MG012); PI signed off 2026-10-31; **CLOSED** |
| MG5 | 4 | Behavior of solitons on an evolving hypergraph (persist / drift / dissolve): a settled particle on a rewritten 2D Delaunay graph | **CLOSED 2026-11-02: PASS** (Part B: A1); ChatGPT MG018 accepted, PI signed off |
| MG6 | 5 | G-coupled growth produces a volume-growth deficit; falloff and mass scaling | **REGISTERED 2026-11-02** (draft 3 frozen at `106b3e2`). Part 0 **STOP** (tick check, 2026-11-02); amendment MG021 (Δ = 0.1, narrow waiver) after ChatGPT review MG021; Part 0 completed; cutoff 110 core-h (PI). **CLOSED 2026-11-03: DEMONSTRATED** (CP2 − SH2 = −67.1 nodes, t = −36.2; CP2 − DEC = −66.9, t = −34.1); feedback SURVIVES (CP2, CP4); falloff UNRESOLVED; ChatGPT MG022 accepted; PI signed off |
| MG7 | 5 | Mass scaling: is the coupled density deficit larger for a stronger particle (b-scan, 27.5 vs 28.5)? | **REGISTERED 2026-11-03** (draft 2 frozen at `5b44cde`); Part 0 PROCEED. **RESULT 2026-11-04: NOT SCORED (feedback)**: at b = 27.5 the particle dissolved in 7 of 8 runs (both arms); b = 28.5 SURVIVES (8/8 single). ChatGPT MG025 accepted; **CLOSED 2026-11-04** (PI signed off) |
| MG9 | 2 | Local onset type of the Turing instability (P1) and finite-time front arrest against the single-particle window (P2); deterministic, 1D lattice | **CLOSED 2026-11-04: P1 SUBCRITICAL; P2 UNRESOLVED** (registered at `7810b8b`; amendment 1; PI sign-off 19:21 UTC) |
| MG9c | 2 | Late-time fate of MG9's slow fronts (b = 29.25–29.35): settle into bounded patches or keep advancing; robustness; upper settling edge vs b_r | **CLOSED 2026-11-04: Q1 UNRESOLVED; Q2 NOT ELIGIBLE** (registered at `372336b`; ChatGPT MG030 accepted; PI sign-off 22:29 UTC) |
| MG9d | 2 | Cause of MG9c's late drift (b = 29.30 from the t = 30,000 field): numerics, mesh alignment, reflection, finer mesh; bounded diagnostic | **CLOSED 2026-11-05: MESH EFFECT SUPPORTED** (on D02; R4 INVALID so finer mesh unavailable; C4 DECREASING alongside). Registered at `69ad83e`; 7 runs, 0.08 core-h; ChatGPT MG033 accepted; PI sign-off 02:32 UTC |
| MG9b | 2 | Stochastic barrier: under a linear model, does the activation slope S (log k per unit Ω) change with b inside the particle window? b = 29.20, 29.35 with MG3b's 29.00; lattice-conditional | **CLOSED 2026-11-05: NOT RESOLVED** (PI sign-off card 14:03 UTC, after ChatGPT MG035) (S₁ = −0.023, 95% CI −0.094 to +0.048; branch B3 calibrated; all controls PASS; 280 runs, 26.3 core-h on PC) |
| MG8 | 6 | Scheduler sensitivity: does the within-step node update order (synchronous J, sweep R, random sequential P) change the particle at τ₀, and does the deterministic difference shrink as τ → 0? (multiway/foliation gate, numbered MG7 in `InitialPlan.md`) | **CLOSED 2026-11-06: SENSITIVE at τ₀; refinement NOT RESOLVED BY THIS LADDER** (registered 2026-11-05, draft 4 at `3ab4ec8`) |
| MG10 | Extension | Can explicit finite A/B reservoirs recover the fixed-reservoir Model G approximation, and what changes under depletion? | **DRAFT 3 — NOT REGISTERED**; PI approved decisions 1-4; external review per run/batch; [proposal](spec/MG10_proposal.md); no MG10 run |

## Active process update (MG039, 2026-11-06)

[RESEARCH_PROTOCOL.md](RESEARCH_PROTOCOL.md) supersedes the mandatory two-AI role/review arrangement prospectively. Historical gate-specific records and their approvals are preserved; this change does not waive scientific validation, registration, blocking issues, or PI sign-off.

## Template

```
### MGn — <title>
Pre-registered: <date>, commit <hash>
Prediction:
Pass criterion (quantitative):
Fail criterion:
Null control(s):
Configs:
Result:
Replication (ChatGPT):
PI sign-off:
```

---

### MG1 — Continuum reproduction of the published Model G soliton
Status: PRE-REGISTERED 2026-10-27 (MG001b), with PI sign-off on the criteria below; reviewed by ChatGPT (round 2). No MG1 run preceded this commit.

Reference: `spec/model_g_published.md` (Pulver & LaViolette 2013, eq. 17 parameters, Table 1).

Setup: dimensionless eqs. 7a–7c; eq. 17 parameters; seed eq. 18 added to the X equation;
φ = 0 initially; t ∈ [0, 100]. 1D Cartesian x ∈ [−50, 50]; 2D and 3D radial r ∈ [0, 50];
Dirichlet φ = 0 at the outer edge, Neumann at the center. Method-of-lines finite differences
with a stiff implicit integrator (our choice; the paper does not report its method).

Prediction: a stationary soliton forms in 1D, 2D and 3D, matching Table 1.

Pass criteria (all must hold):
1. Turing wavelength λ₀ = 3.08 ± 2% in 1D, 2D and 3D (measured as in the paper: twice the
   mean zero spacing).
2. Core amplitudes (Table 1 row a) for φ_G, φ_X, φ_Y within ±3% in 1D, 2D, 3D.
3. Core radii (row b) within ±3%.
4. Stationarity: max |∂φ/∂t| < 10⁻⁴ at t = 100.
5. Grid convergence: criteria 1–3 change by < 0.5% when the grid spacing is halved.

Fail criterion: any pass criterion fails at the converged resolution. Partial reproduction is
reported row by row, not rounded to a pass.

Null controls:
- N1: no seed → max |φ| < 10⁻⁶ at t = 100 (homogeneous state stable; consistent with the
  MG001 linear check, max Re σ = −0.195).
- N2 (exploratory, not pass/fail): seed amplitude scaled by 0.25, 0.5, 0.75 → report whether a
  soliton nucleates and estimate the nucleation threshold amplitude. This feeds MG3.

Diagnostic (not pass/fail; added MG001a, corrected MG001b): integral balances. At R = 50, 100
and 200, report I_S, the measured boundary fluxes F_S, **and the integrated time derivatives**
dI_G/dt and d(I_G + I_X + I_Y)/dt. Check both time-dependent balances in
`model_g_published.md` §9.2 before dropping any term. A pointwise tolerance |∂φ/∂t| < 10⁻⁴ does
not guarantee small integrated residuals: the 3D ball of radius 50 has volume ≈ 5.24 × 10⁵.
The balance-error budget includes spatial discretization, time integration, boundary
differentiation and quadrature. Each is estimated from the grid- and tolerance-refinement runs.
Expected: balances close within that budget, and I_G, I_X → 0 as R grows. A Table 1 row (e)
mismatch at R = 50 counts as consistent only if the *measured* fluxes and time derivatives
account for it.

Out of scope for MG1 (later stages): the 2Dδ variant, multi-particle bonding (3.303 / 3.314),
gradient drift (v = −0.0599 m).

Configs: `configs/MG1_main.json`; drivers `src/run_mg1.py`, `src/run_mg1_controls.py`, `src/evaluate_mg1.py`.
Result (MG002, 2026-10-27; full tables in `results/MG002/MG1_evaluation.md`):

| Criterion | 1D | 2D | 3D |
|---|---|---|---|
| 1 λ₀ ±2% | PASS (3.0756, −0.14%) | PASS (3.0766, −0.11%) | PASS (3.0793, −0.02%) |
| 2 amplitudes ±3% | PASS (≤0.15%) | PASS (≤0.30%) | PASS (≤0.23%) |
| 3 core radii ±3% | PASS (≤0.05%) | PASS (≤0.30%) | PASS (≤0.11%) |
| 4 max \|∂φ/∂t\| < 1e-4 at t=100 | **FAIL (3.5e-4)** | PASS (1.4e-10) | PASS (1.3e-8) |
| 5 grid halving < 0.5% | PASS (0.20%) | PASS (0.01%) | PASS (0.01%) |

- **Verdict as registered: FAIL.** 14 of 15 cells pass. The single failure is 1D stationarity at
  t = 100. The 1D soliton is still relaxing there (an exponential mode with rate ≈ 0.047, converged
  under grid and tolerance refinement). It is stationary by t = 300 (3e-8), where φ_X(0) has moved
  0.1% (−8.372 → −8.365). Our reproduction of the published protocol remains nonstationary at
  t = 100 while matching the tabulated values (within 0.15%). (Wording softened MG002a per review.)
- Unregistered rows also reproduce: RMS radii and core integrals agree within 0.5%, except 2D
  core ∫φ_G (−1.4%).
- **Table 1 row (e) (full integrals) is NOT reproduced for G and X.** Stationary I_G and I_X are
  ≤ 5e-8 (2D/3D, R = 50) against the paper's 2.88e-3 and 2.87e-2 (I_G), 0.128 and 1.93 (I_X).
  The measured boundary fluxes at R = 50 are ≤ 1.2e-7, so finite-boundary flux does **not**
  explain the paper's values. Nonstationarity does not either (2D/3D are stationary to ≤ 1e-8).
  The remaining candidates are numerical/quadrature error in the paper's integration or a
  tabulation inconsistency. I_Y differs by +1.5% (2D) and +5.6% (3D) while core integrals agree,
  so tail quadrature is a plausible explanation, not an established cause.
- Integral identity confirmed numerically: I_G, I_X are at round-off × volume level for
  R = 50, 100, 200 once stationary. Discrete balance residuals ≤ 4e-7 scale with domain volume
  (round-off). In 1D at t = 100, the nonzero I_G, I_X (≈ 2e-6, 2e-5) are carried by
  dI/dt ≠ 0 (nonstationarity), as the review anticipated.
- Null controls: N1 PASS (max |φ| ≤ 4e-14, all dims; a 1e-3 seed decays to ≤ 3e-11).
  N2 nucleation thresholds (seed-amplitude scale, judged at t = 100, sub-threshold confirmed to
  t = 400): **1D 0.975, 2D 0.572, 3D 0.486**. The published seed is only 2.5% above threshold in 1D.
- Solver note: scipy BDF step control failed on the 3D zero-seed run (state at round-off). It is
  retried automatically with an explicit first step, and the retry is recorded in `controls.json`.
  There was no effect on any seeded run.
- Wavelength window sensitivity (MG002a, `results/MG002/wavelength_windows.txt`): the pooled
  zero-spacing wavelength depends on the radial window. The registered estimator gives 3.0756 /
  3.0766 / 3.0793 (1D/2D/3D); the window 10 < r < 25 gives 3.0864 / 3.0853 / 3.0854. All windows
  pass criterion 1. The spread across windows is 0.2–0.4%, so the quoted four decimal places
  describe a particular window, not a uniquely determined wavelength.
- Seeded integral balance corrected (MG002a): the summed balance includes +∫χ dV. Residual during
  seeding (t = 10) ≤ 4e-9. The reported t = 100/300 diagnostics are unaffected (∫χ ≈ 0 there).
Replication (ChatGPT): independent implementation (centered finite differences, Radau, own
expansion of the equations), h = 0.0125, t = 100. 3D φ_G(0) −0.41148105, φ_X(0) −14.63365708,
φ_Y(0) 1.69934939, λ₀ 3.07930594, max |∂φ/∂t| 1.3090e-8: agreement with ours within 0.0012%.
1D max |∂φ/∂t| = 3.49318e-4, confirming the registered failure. Integrals and null controls not
independently replicated. Review: `review/MG002_chatgpt.md`; script and data: `review/MG002_replication/`.
PI sign-off (criteria): Richard Benear, 2026-10-27
PI sign-off (result): Richard Benear, 2026-10-27 — MG1 accepted as PASS with documented exception:
criterion 4 failed in 1D at t = 100 (registered failure record preserved above); the 1D soliton is
stationary by t = 300. Proceed to Phase 2.

---

### Phase 2 common setup (both gates)
Protocol: `spec/phase2_protocol_v1.3.3.md` = v1.3 + Amendments 1–3 (A3 PI-approved 2026-10-28) [v1.3.2 = v1.3 + Amendments 1–2 (A2 PI-approved 2026-10-28) [v1.3.1 = v1.3 + Amendment 1 (PI-approved 2026-10-27; v1.0–v1.2 superseded after ChatGPT reviews MG003 r1–r3; see
`review/MG003_response*.md`). 1D periodic ring, h = 0.2, L = 100. Rules v1.1. Unsplit multinomial
tau-leaping (its Ω → ∞ mean map is exactly forward Euler of the lattice ODE). τ = 6.667e-4, with
Ω-dependent halving under the clip rule. Validations V1–V4 and V6 (A1 pilot on real runs) must pass
before any gate run [A2: V2/V4 decided at Ω = 10⁸; V7 descriptive, with a PI checkpoint before MG3]. All statistical decisions are conditional on assumptions A1–A2 (protocol §5.3),
which are checked per cell (diagnostics → INCONCLUSIVE) and per gate (pooled audit → every PASS
downgraded to INCONCLUSIVE pending investigation).
Operating points: **b = 29** (published) and **b_mid = 29.16** (centre of the operational h = 0.2
single-soliton window [28.929, 29.397]).

### MG2 — Ω → ∞ convergence of the token rules to the lattice mean field
Status: PRE-REGISTERED 2026-10-27 (MG003d), with PI sign-off on protocol v1.3 and the criteria below. No MG2 gate run preceded this commit.

Setup: each run starts from the lattice-ODE soliton at that b (counts = round(Ω·(S₀ + φ*))), T = 100,
snapshots every 0.5, t ≥ 10 used. Ω ∈ {10, 30, 100, 300, 1000, 3000, 10⁴, **10⁵, 10⁶, 10⁸**}
[A2: 10⁵, 10⁶ and 10⁸ added], **12 eligible runs per (b, Ω)** (at most 24 attempted; protocol §5.2–5.3).
Profiles are survival-conditioned (protocol §5.2), with survival reported separately.

Prediction: the aligned mean profile converges to the same-lattice ODE soliton with O(1/Ω) bias.
Far-field fluctuations match the LNA.

Pass criteria:
- C1: at **Ω = 10⁸** [A2; was 10⁴, which remains in the ladder and is reported] and both b, the complete single-cell procedure (protocol §5.3: eligibility and
  top-up, A1 diagnostics, three-way decision with ε = 1%, one 12-run extension) for each of
  S = G, X, Y. PASS if all fields PASS; FAIL if any field FAILs; otherwise INCONCLUSIVE.
- C2: bias scaling. α̂ and its calibrated CI at each b (protocol §5.3: pre-selected cells, distinct-pair
  bootstrap, widening w).
  - **PASS** if α̂ ∈ [0.7, 1.3] and the CI excludes 0.5.
  - **FAIL** if the CI lies entirely outside [0.7, 1.3].
  - **INCONCLUSIVE** otherwise, including fewer than 3 selected cells.
- C3: far-field Ω·Var(φ_S) for each of S = G, X, Y within 10% of the continuous-time LNA at every
  Ω ≥ 100. [Amendment 1] Data rule: the first 12 attempted runs per cell, whether eligible or not
  (protocol v1.3.1 §5.2). A cell lacking data is INCONCLUSIVE, and C3 cannot then PASS. Targets: b = 29: 15.59 / 112.69 / 2.73; b = 29.16: 15.59 / 121.87 / 2.84. The scheme's own
  discrete LNA differs from these by < 1%.
- C4: persistence at b_mid. For Ω ≥ 100, every run keeps its original particle as the only particle for
  all t ∈ [0, 100] (no loss, replication or merge). At b = 29, survival is reported as a measurement
  (near-edge noise-induced decay), not a criterion.
- C5: seeded nucleation (b = 29, homogeneous start, R_seed rule, 20 runs per cell). Nucleation = a
  confirmed particle (protocol §5.4) by t = 60. [A3: the §5.4 field-level persistence rule; the field
  at min(first detection + 20, 60)]
  - Seed scale 1.5 (1.58× the lattice threshold): P ≥ 0.9 at each of Ω = 300, 1000, 3000.
  - Seed scale 0.5: P ≤ 0.1 at Ω = 3000.
  - Exploratory: P(Ω) at scale 1.0.
- C6: lattice independence. h = 0.1 lattice (τ = 1.667e-4), b = 29, **Ω = 10⁸** [A2; was 10⁴]: the complete single-cell
  procedure with ε = 1% against the h = 0.1 lattice ODE.

Gate outcome:
- **PASS** if C1, C3, C4, C5 and C6 pass and C2 is not FAIL.
- **FAIL** if any of them FAILs.
- **INCONCLUSIVE** if C1 or C6 remains INCONCLUSIVE after its extension, or if the gate-level pooled A1
  audit flags any run (every PASS is then downgraded pending investigation).

Controls:
- N1: homogeneous start, b = 29, Ω = 10⁴, 4 runs, T = 100 → no particle detected.
- N2 (the metric's discriminating power): the b = 29 runs evaluated against the b_mid reference must
  give C1 = FAIL (lower bound > 1%; the deterministic difference is 9.46%).
- V4's alignment-induced bias is reported next to every E_b.

Exploratory (reported): particle diffusion coefficient D_p(Ω) (expected ∝ 1/Ω); survival curves at b = 29.

Budget note (PI decision, 2026-10-28, recorded before any C2–C5 evaluation): the MG2 compute cap is raised from 8 to
12 core-hours. The Ω = 10 and 30 ladder cells are deferred (STAGE_LOG D3); they enter no criterion (C3 and C4 apply
at Ω ≥ 100, and C2 uses pre-selected cells), and they are reported as not run.

#### MG2 result (2026-10-28, protocol v1.3.3; full record `results/MG2/MG2_evaluation.md`)
Registered verdict: **FAIL** (C3 and C4 FAIL). Compute used 11.6 of 12 core-hours; Ω = 10, 30 not run (D3).

| Criterion | Result | Note |
|---|---|---|
| C1 (Ω = 10⁸, b = 29 and 29.16) | PASS | bias upper bounds ≤ 0.070 % (G), ≤ 0.012 % (X, Y) vs ε = 1 % |
| C2 (bias scaling α) | INCONCLUSIVE | no ladder cell has a bias lower bound > 0, so no cell is selected; 12 runs cannot resolve the O(1/Ω) bias |
| C3 (far-field Ω·Var vs LNA, Ω ≥ 100) | FAIL | single failing cell: b = 29, Ω = 100, X +17.1 % (run-level SE 1.6 on 131.9; G −0.0 %, Y +7.3 %); every cell Ω ≥ 300 within ±10 %; b = 29.16, Ω = 100 INCONCLUSIVE (data rule) |
| C4 (b_mid persistence, Ω ≥ 100) | FAIL | Ω = 300: 12 replication events in 11 of 24 runs [corrected per MG005 review; first record said "12 of 24 runs"]; Ω = 100: all 24 runs violate (121 loss/replication events; 72 births); Ω ≥ 1000: no events |
| C5 (seeded nucleation) | PASS | scale 1.5: 19/20, 20/20, 20/20 at Ω = 300/1000/3000; scale 0.5: 0/20 at Ω = 3000 |
| C6 (h = 0.1, Ω = 10⁸) | PASS | bias upper bounds ≤ 0.069 % (G), ≤ 0.010 % (X, Y) |
| N1, N2, pooled A1 audit | PASS / PASS / ok* | *The registered pooled audit is non-informative: per-cell standardisation caps its statistic at ≈ 3.43 < 4.65, so it cannot flag any data of this shape (MG005 review P1). A corrected within-cell audit (descriptive, not a retroactive replacement) flags no run: max \|t\| 4.22 vs 8.47 (`results/MG005/audit_corrected.txt`). |

Post-gate diagnostic (descriptive, `results/MG2/diag_replication.txt`): the C4 replication events are physical. Three
replicating runs re-run with identical streams reproduce their events; their fields 20 time units after the event
continue deterministically (b = 29.16, to t = 300) to **two** persistent particles. A no-event run continues to one.

The diagnostic shows genuine two-particle states in the numerical Model G system for the checked runs; it does not validate every
replication label (nor the crowded Ω = 100 labels), and does not distinguish a split from nucleation nearby.

Exploratory: Ω·D_p ≈ 0.08–0.13 from Ω = 300 to 10⁸ (particle diffusion ∝ 1/Ω). Seeded nucleation at scale 1.0:
13/20, 14/20, 18/20 at Ω = 300/1000/3000.

Interpretation (revised after the MG005 review, `review/MG005_chatgpt.md`): the high-Ω agreement (C1, C6, and C3 for
Ω ≥ 300) is a success; the low-Ω failures are limitations of that agreement. C4 is not a counterexample to mean-field
convergence. Noise at Ω ≤ 300 produced replication at b_mid and loss at b = 29 (17 of 24 runs at Ω = 300); a narrow
noisy persistence region is a plausible reading, **not a mapped result**. The deterministic distance b_r − b_mid is not a
stochastic escape barrier, so choosing the deterministic window midpoint did not guarantee noise-robust persistence,
and it is not established that another b would pass at Ω = 100. The C3 cell at Ω = 100 contains 16 births and 23
replications in its 12 measured runs, but finite-Ω corrections to the leading-order LNA (Ω·Var = C0 + O(1/Ω)), detection
and masking effects can coexist with that contamination; the saved aggregates cannot separate them, and no cause is
asserted. The registered FAIL stands; MG2 will not be re-scored by re-selecting b.

### MG3 — Spontaneous nucleation from intrinsic noise
Status: PRE-REGISTERED 2026-10-27 (MG003d), with PI sign-off on protocol v1.3 and the criteria below. No MG3 gate run preceded this commit.
Amended by A3 (v1.3.3) and **A4 (v1.3.4, approved 2026-10-28)**: Stage-A terminal outcomes, tracker fix, corrected audit. Interpretation
carries the A1 residual risk (absence of rare run-level modes is assumed, not tested).

Setup: b = 29, homogeneous Poisson start, no seed, T = 500. Events are **confirmed nucleations**
(protocol §5.4), handled by the first-passage model of protocol §6 (burn-in 50, censoring at 480,
first event per run, constant-hazard MLE, Garwood CIs, stationary-hazard check).

Adaptive design (fully specified):
- **Ω grid:** [A3] {300, 250, 200, 150, 120, 100, 80, 65, 50} (v1.3–v1.3.2: {300, 100, 50, 30, 20, 15, 10,
  7, 5, 3, 2, 1}).
- **Crowding exclusion [A3]:** a cell with mean detections per snapshot > 2 in t < 50 is excluded from
  the isolated-event study. It is applied to every Stage-A and Stage-B candidate, including midpoints.
  Excluded cells never re-enter via neighbour selection, are reported, and stop the descent (a budget
  choice). Events are confirmed by the A3 field-level rule.
- **Stage A (scan):** 4 runs per Ω, descending. Stop after the first Ω at which 4/4 runs have a
  confirmed event in [50, 150]. At most 12 cells (48 runs).
  - [A4] *Screened set* = the grid cells run in Stage A and not excluded by crowding. Terminal outcomes:
    - (a) No confirmed event in any screened cell (grid exhausted at Ω = 50, or descent stopped by an excluded
      cell): **C1 FAIL**, an exposure-limited bound over the screened Ω range. Report 95% upper bounds per Ω
      and stop. No Stage B. (v1.3–v1.3.3 text: "at any Ω down to 1".)
    - (b) Confirmed events but no 4/4 cell before the screened set ends: S = screened cells with 1–3 of 4;
      neighbours only from screened, non-excluded grid values; if |S| < 3 when these are exhausted, Stage B runs
      on the available cells (C3 then INCONCLUSIVE by its own rule). No midpoints.
    - (c) Midpoints only between two adjacent screened, non-excluded grid values bracketing a 0/4 → 4/4
      transition; each gets a V3 clip pilot and the crowding screen before its Stage-A cell.
    - (d) 4/4 already at Ω = 300: Stage B at Ω = 300 only (C1 and C2 evaluable); C3 INCONCLUSIVE.
    - (e) Ω = 300 excluded by crowding: no screened set; MG3 INCONCLUSIVE.
- **Stage B cell selection:**
  - Start from S = the grid values with 1–3 of 4 runs confirmed in Stage A.
  - If |S| < 3, add the nearest grid neighbours of S, or of the 0/4 → 4/4 transition, until |S| = 3.
  - If the transition sits between two adjacent grid values with no partial cell, insert up to 2
    geometric midpoints (2 significant figures). Each inserted point gets its own 4-run Stage-A cell
    before selection.
  - At most 4 Stage-B cells.
- **Stage B runs:** 16 per cell, T = 500. A cell with < 5 confirmed events gets **one** extension of 16
  runs. Cells still below 5 events are excluded from the C3 fit and reported with their Poisson CI or
  upper bound.
- **Hard budget:** 10 core-hours for all of MG3. If it would be exceeded, stop at the current step:
  MG3 = INCONCLUSIVE, with all data reported.

**Pre-run implementation specification (2026-10-28; recorded before any MG3 run; Claude's choices for details the
registration left open, plus two PI decisions; flagged for ChatGPT review).** Code: `src/mg3.py`, `src/run_mg3.py`;
tests `tests/test_mg3.py` (15 pass [corrected MG006a: 14 test functions, all pass]), `tests/test_track.py` (8).
- **Budget (PI decision):** 18 core-hours (was 10). From MG2's measured cost (~1.07 core-s per time unit per run plus
  ~7 s per confirmation), the registered design needs about 23 core-hours with full-length runs.
- **Early stop (PI decision):** each run stops at the first 10-unit segment boundary at least 50 time units after its first
  confirmed event's detection. Runs without one go to T = 500. Only the first event per run enters any criterion.
  The exploratory post-event studies are limited to that 50-unit window.
- **Event bookkeeping:** births (persistent tracks, §5.1 with the A4 fix) are confirmed (§5.4) in detection order at
  detection + 20. After the run's first confirmed event, later births are not confirmed. A first event at t < 50 is a
  startup event: the run is left-truncated, gives no exposure, and does not count in Stage A.
- **Stage-A statistic:** the number of the 4 runs whose first confirmed event is detected in [50, 150]. Crowding screen:
  the mean over the 4 runs of mean detections per snapshot in t < 50, using only those first 50 time units.
- **Rates:** Stage-B runs only. Stage-A runs are excluded from k̂ because they were used for selection; they count for C1.
- **Stage-B selection details:** neighbour ties go to more Stage-A events, then lower Ω. More than 4 partial cells: the 4
  lowest Ω. Events only after t = 150 (no partial cell, no 4/4): anchor at the lowest-Ω cell with an event.
  Midpoints: the geometric mean to 2 significant figures, up to 2, bisecting the remaining 0/4 → 4/4 bracket.
- **Middle cell** (V5 and C4): the median Ω of the Stage-B cells (for an even count, the lower of the two middle values).
  If V5 fails, the driver stops and returns to the PI, since repeating every run at τ/2 needs a budget decision.
- **C2(a):** persistent / (persistent + excursion + unresolved) field outcomes over the cell's runs, up to and including
  each run's first confirmed event. A shared outcome is counted once.
- **C2(b):** for the nucleated side, each Stage-B run with a first confirmed event contributes its aligned profile over
  [t_e + 30, t_e + 50] in 4 blocks of 5 time units. An eligible snapshot has exactly one detection, belonging to the
  nucleated track; a run is eligible if every block has ≥ 7 eligible snapshots (of 10). The reference side is 8
  soliton-start runs at the same Ω (MG2 protocol, T = 100, MG2 eligibility). The test is on φ_X, |x| < 10: Welch bound
  on run means, normalised by the deterministic ‖φ*_X‖ on that window, PASS if the upper bound < 10%. A1 diagnostics
  are applied to both sides. Fewer than 8 eligible on either side → the registered replacement, (a) plus (c).
- **C2(c):** λ by a pre-stated core estimator: 2 × the mean spacing of zero crossings in 1 < x < 7 of the x-symmetrised
  aligned mean φ_X, pooled over the cell's eligible nucleated runs. It is compared with the same estimator on the
  lattice-ODE soliton (3.0296 at b = 29; the tail estimator gives λ* = 3.0619), PASS within 5%. The tail-based λ* is
  not used because the soliton's tail falls below the Ω ~ 10² noise beyond x ≈ 7.
- **C3:** weighted LS of ln k̂ on Ω with weights n_i (inverse Poisson variance of ln k̂). The 95% slope CI uses those known
  variances (normal); a t-based CI (df m − 2) is reported as a sensitivity check. PASS needs the slope CI below 0 and
  weighted R² ≥ 0.9.
- **C4:** the detector at each b uses that b's reference (§5.1 thresholds −2.7685, −2.9908, −3.1225); b = 29 reuses the
  Stage-B cell. Monotonic point rates, and k(29.3)/k(29.0) with its exact conditional-binomial two-sided 95% CI lower
  bound > 1.
- **N1:** b = 28.5 dynamics and confirmation, b = 29 detector (no soliton exists at 28.5 to define one), smallest
  Stage-B Ω, 16 runs. **N2:** the Euler mean map from the exact homogeneous state at b = 29 to T = 500; max |φ| reported.
- **Clip pilots** (V3 procedure) for each new (b, Ω). **Streams:** stage A 0, B 1, C4 2, N1 3, V5 4, C2(b) reference 5,
  pilot 6; b_idx in [29, 29.16, 29.3, 28.5]; midpoints om_idx 9, 10; segment = t₀/10.

Prediction: nucleation is an activated (large-deviation) process, ln k ≈ c − S·Ω.

Pass criteria:
- C1: spontaneous nucleation exists: ≥ 10 confirmed first-nucleation events at some Ω
  (Stage A + B, post burn-in).
- C2: what nucleates is a Model G soliton. In every Stage-B cell:
  - (a) ≥ 80% of detections are confirmed (not excursions);
  - (b) the survival-conditioned aligned mean profile, 30 time units after detection, matches a
    **same-Ω reference ensemble**: 8 runs started from the deterministic soliton, MG2 protocol, T = 100.
    Test: the calibrated Welch two-sample bound (protocol §5.3) on ‖μ_nuc − μ_ref‖/‖φ*_X‖, normalized
    by the **deterministic** lattice-ODE soliton norm, with PASS if the upper bound < 10%. A1 diagnostics
    apply to both samples. It needs ≥ 8 eligible runs on each side; otherwise (b) is
    INCONCLUSIVE and the stated replacement below applies.
  - (c) λ of the nucleated mean profile is within 5% of the lattice value.
  - Survival fractions are reported.
  - If either side has < 8 eligible runs, (b) is replaced by (a) plus (c), and this is stated.
- C3: activated scaling. Over ≥ 3 Stage-B cells with ≥ 5 events that pass the stationary-hazard check,
  a weighted fit of ln k̂ vs Ω (weights from Poisson counts) gives slope −S with the 95% CI excluding 0,
  and R² ≥ 0.9. If fewer than 3 cells qualify, C3 is INCONCLUSIVE.
- C4: distance from the window edge. At the median Stage-B Ω, add b = 29.16 and b = 29.3 cells
  (16 runs each, one extension allowed), reusing b = 29 from Stage B. The point rates must increase
  monotonically in b, and k(29.3)/k(29.0) > 1 at 95% (conditional binomial).
- V5 (rare-event τ check, protocol §6) must pass.

Gate outcome:
- **PASS** if C1–C4 and V5 pass.
- **INCONCLUSIVE** if C3 is inconclusive or the budget stops the gate.
- **FAIL** otherwise.

#### MG3 result (2026-10-29, protocol v1.3.4; full record `results/MG3/MG3_evaluation.md`)
Registered verdict: **FAIL** (C2). Compute: 19.03 core-hours against the 18-hour cap. The overrun is inside N1, the last
cell, which started at 15.6 h. The budget is checked at cell start (as in MG2), so no step was cut.

| Criterion | Result | Note |
|---|---|---|
| Stage A | — | Ω = 300, 250, 200, 150: 0/4, no events in 4 × 500; Ω = 120: 0/4 in [50,150] (3 startup events, 1 at t = 330); Ω = 100: 2/4; Ω = 80 excluded by crowding (2.26 detections/snapshot). Stage B: Ω = 150, 120, 100 |
| C1 (≥ 10 confirmed post-burn-in events at some Ω) | PASS | 15 at Ω = 100, 11 at Ω = 120 (Stage A + B) |
| C2 (what nucleates is a soliton) | **FAIL** | (a) persistent fraction 5/7 = 0.71, 15/21 = 0.71, 16/26 = 0.62 against ≥ 0.80, so it fails in every cell. (b) replaced by (a)+(c), as pre-stated: the soliton-start reference ensembles kept their particle in 1, 0 and 0 of 8 runs. (c) λ_core −2.9% (Ω = 150, 1 run), −2.3% (Ω = 120, 3 runs), PASS; Ω = 100 INCONCLUSIVE (no profile-eligible run) |
| C3 (activated scaling) | PASS | k̂ = 4.05e−6, 4.63e−5, 1.58e−4 at Ω = 150, 120, 100 (5, 10, 11 events). Weighted fit S = 0.0726 per unit Ω, 95% CI [0.0516, 0.0936]; R² = 0.993. With 3 cells the t-based sensitivity interval (1 df) includes 0. [MG006 review] The weights are estimated, not known, Poisson variances; three cells (5 events in the weakest) support an activated trend over Ω = 100–150 but do not establish the mechanism |
| C4 (distance from the edge, Ω = 120) | PASS | k̂ = 4.63e−5, 1.52e−4, 2.73e−3 at b = 29, 29.16, 29.3 (monotone); k(29.3)/k(29) = 58.9, 95% CI [10.4, 229]. Most b = 29.3 runs (29/32) nucleated during burn-in; the b = 29.3 rate rests on the 3 surviving runs, whose waits after burn-in were 2, 4 and 5 time units [MG006 review] |
| V5 (τ/2, Ω = 120) | PASS | ratio 0.883, 95% CI [0.365, 2.221] (14 vs 10 events). [MG006 review] Registered comparison passed; time-step independence is **not resolved**, given the broad interval and the startup asymmetry |
| N1 (b = 28.5, Ω = 100) | no confirmed event | 16 runs to T = 500, 0 persistent detections |
| N2 | PASS | mean map stays homogeneous (exact start 0; 10⁻³ perturbation decays to 1e−11) |

Independent recomputation from the raw run files (rates, the C3 fit, the C4 interval, the C2(a) counts) agrees.
ChatGPT's independent replication (`review/MG006_chatgpt.md`) checked 168 run records with no bookkeeping mismatch (exposure,
censoring, confirmation time, early stop). It re-derived C1–C4 and V5, re-confirmed all 36 Stage-B first-event fields with a
separately implemented solver, and replayed two excursion records exactly (both decay).

Qualifications from that review, adopted:
- **C2(a)** describes outcomes through each run's first confirmed event, as pre-specified. It is not a measure of detector quality
  over full trajectories.
- **Stationarity is weakly checked.** Ω = 100 and b = 29.3 have no late-window exposure, so their p = 1 values carry no information.
- The CIs (Garwood, conditional binomial) are applied to stopped, left-truncated first-passage data. Their exactness in this
  design has not been justified; no undercoverage has been shown either.
- The Stage-A 4/4 trigger had at most 0.385⁴ ≈ 2.2% probability per cell under the exponential model.
- The budget overrun (19.03 against 18 core-hours) stands as recorded. **Budget exception, recorded explicitly:** the pre-run
  text raised the cap to 18 h but did not replace the registered hard-stop wording with the driver's cell-start check
  (which was in the code before any data). Under a literal hard-cap reading MG3 would be INCONCLUSIVE, with the observed
  C2 failure retained. Following the reviewer's recommendation, the registered verdict is kept as FAIL with this overrun
  documented as an exception. Ledger: `results/MG3/budget_ledger.txt` (18.958 h of completed work + 0.075 h in the excluded
  Ω = 80 checkpoints, now committed in `results/MG3/excluded_A_80_checkpoints/`).

Further wording corrections from the full review (`review/MG006_replication/README.md`):
- 18 of 54 sampled field outcomes (33%) are excursions: a minority, not "mostly". Preferred wording: some fields sampled
  after persistent detections relax to homogeneous under noise-free continuation; others keep particle-depth structures.
  An operational excursion is a finite-horizon result, not proof that the earlier noisy trajectory never crossed a
  basin boundary.
- Field-level confirmation does not attribute the persistent structure to the candidate detection. B_100 replicates 13
  and 14 have `near = False`. This is allowed by the protocol and is not grounds for rejection, but it limits per-track
  wording.
- The profile-eligible fractions (1/5, 3/15, 0/16) combine isolation, tracking and coverage. They are **profile
  eligibility, not survival probabilities**.
- C3: "PASS under the specified normal approximation; consistent with activated scaling over Ω = 100–150". It does not
  establish a large-deviation mechanism. At Ω = 150 the early rate is about 6.7 times the late rate (4 vs 1 events);
  not rejecting at 0.01 is weak evidence. The τ/2 cell's p = 0.060 (4 early, 10 late) is not a convincing validation either.
- Garwood and conditional-binomial "exactness" does not transfer automatically to stopped first-passage data. In a
  one-subject example Garwood covers about 97% (conservative); MG3 undercoverage has not been shown.
- C4's factor 58.9 is an imprecise, model-dependent contrast resting on 3 survivors (waits 5, 4, 2), not a
  well-measured physical factor.
- Stage A: the 4/4 trigger is poorly powered (≤ 2.2% per cell under a constant hazard from t = 0), not literally
  unreachable.

Interpretation (Claude; for review). Spontaneous nucleation from intrinsic noise is **observed and confirmed** (C1).
Its rate is consistent with an activated law ln k ≈ c − S·Ω over Ω = 100–150 (C3). It rises toward the window interior
(C4). V5 passed its registered post-burn-in comparison; startup sensitivity to τ is unresolved [wording corrected per MG006 review]. The FAIL is C2(a): about 30% of persistent detections (tracks lasting ≥ 20
time units) relax to the homogeneous state when the field at detection + 20 is continued without noise. That fits the
pre-stated reading that noise also drives transient structures that do not cross the separatrix. The 80% threshold was
set before any such data existed. C2(b) could not be tested: at b = 29 and these Ω, seeded solitons do not survive
(MG2's finding). C2(c) rests on 1 and 3 runs. After nucleating, the particle is usually not the only detection for the
next 50 time units (profile-eligible 1/5, 3/15, 0/16); whether that is from extra nucleation, replication or loss is
not established. Design limitation: under a constant hazard the Stage-A statistic (first event in [50,150], with startup
events excluded) is at most 0.385 per run, so a 4/4 stop was effectively unreachable. The descent ended at the crowding
exclusion. No criterion is re-scored. Observation for follow-up: startup events at Ω = 120 were 5/16 at τ against 0/16
at τ/2 (Fisher p = 0.043), outside any criterion.

Controls:
- N1 (detector false-positive rate): b = 28.5, below the operational lattice b_s, at the smallest
  Stage-B Ω, 16 runs, T = 500. Report the detection rate (excursions per unit length per unit time) and
  the number of confirmed events. A confirmed event at b = 28.5 would contradict the operational window
  estimate, and would be reported and investigated. It is not a gate failure, because the window edges
  are estimates.
- N2: noise-free control. The Ω → ∞ mean map from the homogeneous state at b = 29 stays homogeneous
  (deterministic).

Interpretation of a FAIL or INCONCLUSIVE (pre-stated):
- No events at any screened Ω [A4; was "any Ω ≥ 1"] → the barrier is too high for intrinsic noise at this lattice scale. This is
  an exposure-limited bound, not impossibility.
- Detections that are mostly excursions (C2a fail) → noise drives transient structures that don't cross
  the separatrix.
- Non-activated scaling or non-constant hazard → not a single-barrier, memoryless process.

Exploratory (reported): annihilation rate (soliton → homogeneous) at Stage-B Ω; structure factor before
first nucleation vs the LNA; multi-particle growth after the first event; excursion statistics.

Configs (both gates): `configs/MG2_main.json`, `configs/MG3_main.json`, `configs/validations.json`.
Result:
Replication (ChatGPT):
PI sign-off (criteria): Richard Benear, 2026-10-27 (both MG2 and MG3; protocol v1.3 frozen; ChatGPT spec review closed, round 4)
PI sign-off (result):

---

### MG3b — Follow-up to MG3: time step, nucleation rates and calibrated uncertainty
Status: **PRE-REGISTERED 2026-10-29 (MG007)**, with PI approval. No MG3b run of any kind (calibration or lattice) preceded
this commit. The registered design is `spec/MG3b_followup_proposal.md` **draft 3, as merged at commit `202c00f`**
(content frozen from that commit). Reviews: `review/MG006b_chatgpt.md` and `_r2.md`, with responses
`review/MG006b_response.md` and `_r2.md`. Drafts 1 and 2 are preserved as superseded (`_v1.md`, `_v2.md`).

Nature of the study: a **separate follow-up designed after MG3's data were seen**. It is an estimation study with
registered decision rules. It has **no pass/fail verdict**, and it does not re-score MG3 (MG3's FAIL and the A1 residual
risk are unchanged). It runs under protocol v1.3.4, with no protocol amendment. The study-specific methods (Part 3
candidates, the GLM slope, the equivalence rule) are those of the registered proposal.

Summary of the registered design (the proposal governs where this summary is shorter):
- **Order:** Part 3 (synthetic calibration), then the driver and tests, then Part 1, Part 1b if triggered, then Part 2.
  The Part 3 results and `results/MG3b/calibration/decisions.json` are committed and reviewed before the first lattice run.
- **Setting:** b = 29, h = 0.2, L = 100, homogeneous Poisson start, τ = 6.667e-4. MG3 detection, tracking (A4) and
  §5.4 confirmation; stop 50 time units after the first confirmed event; burn-in 50; censoring 480 by detection time.
  Every confirmation field is saved.
- **Part 1 (time step):** Ω = 120, 64 runs at τ and 64 at τ/2 (96 each if Part 3 gives P(EQUIVALENT | r = 1) < 0.80).
  The decision is on r = k(τ/2)/k(τ) after burn-in, with margin [0.5, 2]:
  - EQUIVALENT if the 90% CI lies inside [0.5, 2] → production at τ;
  - DIFFERENT if not EQUIVALENT and the 95% CI excludes 1 → Part 1b (64 runs at τ/4, same rule against τ/2);
  - INCONCLUSIVE otherwise → production at τ/2, labelled "time-step agreement not established".
  - Part 1b: EQUIVALENT → production at τ/2; otherwise **MG3b stops**, and its rates are INCONCLUSIVE on time-step
    grounds.
  - If no rate-ratio method is eligible in Part 3, Part 1 is INCONCLUSIVE.
  - Startup fractions are descriptive only.
- **Part 2 (rates):** Ω ∈ {130, 120, 110, 100} at the production τ. The Ω = 120 cell is Part 1's arm at that τ. The
  other cells run in batches of 16 until ≥ 20 post-burn-in events, cap 64.
  - Every cell, capped or not, enters the joint Poisson-GLM fit of log k = α − S·Ω.
  - Sensitivity fit without Ω = 100.
  - Landmarks at 100 and 150, with risk sets.
  - Stationary check labelled informative only if its simulated power against both ×2 and ×0.5 is ≥ 0.5.
  - **All rates and the slope are conditional on the production τ.** EQUIVALENT means agreement within a factor 2 at
    Ω = 120, not time-step independence.
- **Part 3 (calibration):** as frozen in the proposal.
  - The configuration grid, including the activated-preserving capped configurations and the full decision tree.
  - The candidates, and the eligibility thresholds (≥ 0.94 coverage, ≤ 0.035 per tail, in every constant-hazard
    configuration and every production branch with ≥ 1000 replicates; R = 4000).
  - The deterministic selection rule, and the fallbacks.
  - The runtime benchmark rule, with R = 1000 for the bootstrap candidates if the projection exceeds the ceiling.
- **Removed:** the full-trajectory confirmation study (draft 1, Part 4). It is deferred until detection-level
  attribution is validated.

Prediction (descriptive; not a criterion): S near MG3's 0.073 per unit Ω; Part 1 EQUIVALENT is the expected outcome,
given MG3's V5 ratio of 0.88. The startup asymmetry observed in MG3 (5/16 against 0/16) may or may not recur.

Null controls: none new. MG3's N1 (b = 28.5: no confirmed event in 16 × 500) and N2 apply to the same detector and
confirmation rule. The internal checks are the Part 1 comparison, the landmarks and the stationary check.

Budget (binding run caps; core-hour ceilings are stop triggers, checked before every batch):
| Item | Run cap | Ceiling (core-h) |
|---|---|---|
| Part 3 calibration | — | 8 |
| Part 1, τ arm | 64 (96) | 6 (9) |
| Part 1, τ/2 arm | 64 (96) | 20 (30) |
| Part 1b, τ/4 arm | 64 | 40 × (measured τ/4 cost ÷ 5.0 core-s per time unit) |
| Part 2 at τ: Ω = 130, 110, 100 | 64 each | 12 |
| Part 2 at τ/2: Ω = 130, 110, 100 | 64 each | 30 |

Enforcement rules (from the proposal):
- Every core-second on the production machine goes into `results/MG3b/budget_ledger.csv`, and the evaluator sums it.
- A batch is not submitted if used + projected cost would exceed the item's ceiling.
- A budget stop reports only counts and costs, and the affected item is INCONCLUSIVE unless the PI records an
  extension here first.
- The escalation branches are pre-authorised within these ceilings.
- Calibration can change only the method choice and the 64 → 96 rule, both through the committed `decisions.json`.

Projected cost: about 24, 36 or 63 core-hours plus calibration, depending on the branch.

**Pre-run implementation specification (2026-10-29; Claude's choices for details left open; recorded before any MG3b
run).**
- **Streams:** `runner.seed_for` with gate id 4. Stage index: Part 1 τ = 0, τ/2 = 1, τ/4 = 2, Part 2 = 3, timing runs
  = 6. b_idx = 0 (b = 29). om_idx: Ω = 130 → 0, 120 → 1, 110 → 2, 100 → 3. replicate = run number within the arm or
  cell; segment = t₀/10. Part 3 uses numpy `default_rng` with entropy 20261029, and a per-configuration child from
  `SeedSequence.spawn`.
- **Run records:** the MG3 `summarize` fields, plus every confirmation field (float32, 3 × 500) with its class, the
  detection position and `near`, and the per-snapshot detection counts. Batches are 16 runs. Run-level resume is as in
  MG3.
- **Ledger columns:** item, key, kind (run / partial / failed / excluded / reference / timing / calibration / test),
  wall_s, workers, machine, os, python, numpy, scipy, numba, utc_time.
- **Timing:** before Part 1, 2 runs of 100 time units at τ/4, Ω = 120, timing stream 6. They are counted in the ledger
  and not used in any estimate. Their mean core-s per time unit rescales the Part 1b projection and ceiling.
- **Driver order:** it refuses to start without `decisions.json`; it runs Part 1 arms interleaved by batch (τ, τ/2);
  evaluates the Part 1 rule only when both arms are complete; then Part 1b or Part 2 automatically.
- Anything further that the code needs and this text leaves open is appended here, dated, before the run it affects.
- **Part 3 implementation choices (2026-10-29; made in `src/calib_mg3b.py` before the calibration run but recorded
  here only after it — a procedural deviation, stated in `results/MG3b/calibration/README.md`):**
  - Startup times are not simulated, because they do not enter any statistic; startup runs give no exposure.
  - Per-cell configurations cover the designs fixed64, fixed96 and stop20 (16 × 4 p_s × 4 rates × 6 shapes each).
    The central configuration (p_s 0.35, k·L = 4.6e-3) is added to each design. Constant hazard gives 51
    configurations.
  - Widths: rates and ratios use the mean log(hi/lo) over replicates with lo > 0; the slope uses the mean
    log(hi − lo). The central configuration for the ratio is r = 1, p_s (0.35, 0.05), 64 runs, 95% interval. For the
    slope it is S = 0.073, 64 runs, r = 1.
  - "Highest minimum coverage" in the no-eligible fallback uses the 95% coverage. A stationary-split tie means power
    within 0.01 (absolute).
  - Decision tree: all Part 1 and Part 1b arms use the Ω = 120 cell's p_s. Part 1b ratio r2 = 1 except for r = 0.5,
    which uses r2 ∈ {1, 0.71}. The capped configurations are at S = 0.073. Each slope configuration runs at 64 and 96
    runs per arm. Replicates ending in "stop" have no slope and are counted.
  - Slope bootstrap: the Ω = 120 cell is re-simulated as a fixed sample of its size. The other cells use stop20. Startup
    loss uses each cell's own p̂_s. A per-cell bootstrap with k̂ = 0 gives the degenerate interval [0, 0].
  - Non-constant hazards: the bias target is Σ H(exposure) / (L Σ exposure), i.e. expected events per unit exposure.
  - Seeds: SeedSequence(20261029).spawn(2) gives one root for the benchmark (index 0) and one for the main run (index 1).
    Each root spawns one child per configuration, in the fixed order cell, Part 1, tree.
  - GLM fit: Newton on the profile likelihood, with step cap 0.5 and step halving (added after the first benchmark,
    before the main run).
  - Branch threshold for the bootstrap: 1000 × R_boot / 4000. The first version used the live-replicate count; this was
    corrected after the run, and the selection is unchanged (`eligibility_recheck.json`).
- **Part 3 result (2026-10-29), from `results/MG3b/calibration/decisions.json`:**
  - per-cell rate CI: Garwood;
  - Part 1 ratio CI: Clopper–Pearson;
  - 64 runs per arm (P(EQUIVALENT | r = 1) = 0.85);
  - slope CI: no candidate eligible, so the profile likelihood is used, labelled **approximate** (calibrated minimum
    coverage 0.937), with no 95% claim on the slope;
  - stationary split: exposure median (low power; usually uninformative).
  - Cost: 0.74 core-hours. Awaiting ChatGPT review before the driver and any lattice run.
- **Part 3 addendum 1 (2026-10-29; registered here before it is run; responds to `review/MG007_chatgpt.md`).**
  Two post-selection diagnostics. They change only labels, not methods, sample sizes or the decision tree. There is no
  lattice simulation, and the existing Part 3 decisions stand.
  - **Simulation:** Part 1 and Part 1b only, with the selected method (Clopper–Pearson) and 64 runs per arm. The τ/4
    arm always has 64 runs. The grid is r = k(τ/2)/k(τ) ∈ {0.5, 0.71, 0.8, 1, 1.25, 1.41, 2} × r2 = k(τ/4)/k(τ/2) ∈
    {0.5, 0.71, 1, 1.41, 2} × startup pairs (τ, τ/2 and τ/4) ∈ {(0.35, 0.35, 0.35), (0.35, 0.05, 0.05)}, at
    k·L = 4.6e-3 for τ. That is 70 configurations, each with **R = 20 000** replicates, and seeds from
    SeedSequence(20261029, spawn_key=(2,)).
  - **(A) Reused Ω = 120 rate.** Garwood 95% coverage of the selected arm's own true rate is computed within each
    production branch (τ after EQUIVALENT; τ/2 after INCONCLUSIVE; τ/2 after Part 1b EQUIVALENT), together with each
    one-sided miss. A branch qualifies in a configuration if it has ≥ 1000 replicates. **Rule:** in any branch where some
    qualifying configuration has coverage < 0.94 or a tail > 0.035, the Ω = 120 rate interval reported in that branch is
    labelled "approximate (observed minimum branch-conditional coverage X over qualifying configurations)", and no 95%
    claim is made for it. It still enters the slope fit unchanged.
  - **(B) Part 1b equivalence error.** At the boundaries r2 ∈ {0.5, 2}, for every r and startup pair, the diagnostic
    reports P(Part 1b triggered), the **conditional** error P(EQUIVALENT | triggered) and the **whole-study** error
    P(triggered and EQUIVALENT). **Rule:** a Part 1b EQUIVALENT is reported with both observed maxima. If the
    whole-study maximum exceeds 0.05, Part 1b's EQUIVALENT is labelled "approximate (selected comparison)", and the τ/2
    production rates carry that label as well as "conditional on τ/2".
  - **Wording:** 0.937 for the slope is an *observed minimum over qualifying calibration cells*, not a guaranteed
    coverage floor.
  - **Implementation correction:** in the main run's hypothetical 96-run tree configurations, the τ/4 arm had 96 runs
    where Part 1b specifies 64. The code now always uses 64. The main run's results are left as they are, and the
    selected 64-run design is unaffected.
  - **Addendum 1 result (2026-10-29; `results/MG3b/calibration/addendum1/`; 0.005 core-hours):**
    - (A) The Ω = 120 rate interval is labelled **approximate in every branch**. The observed minimum
      branch-conditional coverage is 0.943 (τ), 0.875 (τ/2 after INCONCLUSIVE) and 0.841 (τ/2 after Part 1b).
    - (B) The whole-study Part 1b equivalence error at the margin is at most 0.040 (≤ 0.05, so no label). The
      conditional error is up to 0.161, and both figures are reported with any Part 1b EQUIVALENT.
    - The Part 3 method decisions are unchanged.
- **Driver implementation choices (2026-10-29; recorded before any MG3b lattice run; code `src/mg3b.py`,
  `src/run_mg3b.py`, tests `tests/test_mg3b.py`):**
  - **Clip pilots:** Ω = 130 and 110 are not on MG3's grid. Each gets a V3-procedure clip pilot at τ: T = 20, stage 7,
    replicate 0, threshold 1e-4. If a pilot fails, MG3b stops for a PI decision (exit 4), because τ is fixed by the
    registration. Smaller time steps clip less, so a pilot at τ covers τ/2.
  - **Timing:** the two τ/4 timing runs (100 time units each) run concurrently in the pool. Their mean core-s per time
    unit divided by 5.0 rescales the Part 1b ceiling.
  - **Budget accounting:** an item's used core-hours are its ledger rows plus the walls of its partial checkpoints. Part 2
    at τ (or τ/2) is one item covering the three cells. The projected cost per run is the mean over the item's completed
    runs. The ledger records each finished run (kind `run`), the timing runs and the pilots. Walls of runs left
    unfinished at the end are added by the evaluator from `partial/`.
  - **Order:** Part 1 arms run as batches of 16, alternating τ and τ/2. The Part 1 decision is made only when both arms
    are complete. Part 1b compares τ/4 against τ/2. Part 2 cells run in the order 130, 110, 100, one batch at a time.
  - **Records:** each record is the MG3 `summarize` record plus `conf_fields`, keyed by snapshot index. Each entry holds
    the float32 field (3 × 500), its class, the confirmation time, and MG3's confirmation summary.
  - **Exit codes:** 2 = budget stop (counts and costs only), 4 = clip pilot failed, 5 = `decisions.json` missing.
  - **PI budget decision (2026-10-30, 01:51 UTC; recorded before the Part 1 result is known):** the τ/4 timing ran
    only 2 runs on the machine and measured 2.71 core-s per time unit. That scaled the Part 1b ceiling to
    40 × 0.54 = 21.7 core-hours. Under the full 9-worker load, τ runs cost about 1.0 and τ/2 about 2.1 core-s per time
    unit, so τ/4 probably needs about 4. On Claude's recommendation, the PI restored the **unscaled 40 core-hour
    ceiling** for Part 1b. The measured factor still scales the prior per-run projection. The running driver loaded the
    old rule, so if Part 1b triggers a budget stop under it, the driver is restarted (resume) with the corrected code.
    This counts as the recorded PI extension, not a new decision.
  - **PI budget decision (2026-10-30, 02:56 UTC; after the Part 1 result):** Part 1 gave EQUIVALENT, so Part 2 runs at
    τ (item `P2_tau`, registered ceiling 12 core-hours). The first 16 runs at Ω = 130 cost 1.9 core-hours (8 events).
    At that cost, the three cells could need up to about 23 core-hours at their caps, so a budget stop was likely before
    Ω = 110 and 100 finished. The PI removed the core-hour ceiling for `P2_tau`: "I will go ahead and approve a higher
    ceiling. Complete the run, regardless of time." The run caps are unchanged and still bind: stop20 with batches of
    16, at most 64 runs per cell. No estimand, stopping rule or criterion changes. **Disclosure:** the decision was made
    after the Part 1 outcome and the first Ω = 130 batch had been seen. It affects only whether the registered design
    runs to completion, not how any result is computed. The `P2_tau2` ceiling (30 h) is unchanged and not used. If the
    running driver (old rule) stops with exit 2, it is restarted from the corrected code and resumes.
  - **Production machine:** the PI's machine, with `MODELG_WORKERS=9`. The platform and worker count are written into
    `plan.json`, `timing_tau4.json` and the ledger.

Result (2026-10-30; full record `results/MG3b/MG3b_evaluation.md`; evaluator `src/evaluate_mg3b.py`, recomputed
independently by `results/MG3b/independent/recompute.py` with identical numbers). No pass/fail verdict (estimation study).
- **Part 1: EQUIVALENT.** Ω = 120: τ 44 events / exposure 12 891.5, τ/2 39 / 12 017.5. r = k(τ/2)/k(τ) = 0.951,
  Clopper–Pearson 90% CI 0.646–1.397 (inside [0.5, 2]), 95% CI 0.602–1.497. Production at τ; Part 1b not run. Startup
  9/64 vs 10/64 (Fisher p = 1.0; no difference was detected at the 0.05 level).
- **Part 2 at τ** (all rates conditional on τ = 6.667e-4; no cell capped):
  | Ω | Runs | Events | Startup | k̂ | 95% CI |
  |---|---|---|---|---|---|
  | 130 | 48 | 25 | 0 | 1.62e-5 | 1.05e-5 – 2.39e-5 |
  | 120 | 64 | 44 | 9 | 3.41e-5 | 2.48e-5 – 4.58e-5, approximate (observed minimum branch-conditional coverage 0.943) |
  | 110 | 48 | 28 | 19 | 6.27e-5 | 4.17e-5 – 9.07e-5 |
  | 100 | 32 | 21 | 11 | 1.13e-4 | 6.97e-5 – 1.72e-4 |
- **Slope:** S = 0.064 per unit Ω, profile 95% interval 0.047–0.082, **approximate** (observed minimum coverage 0.937;
  not a guaranteed floor). Deviance 0.23 on 2 df (p = 0.89). Without Ω = 100: S = 0.068, nominal (uncalibrated) 95%
  interval 0.041–0.094; the 0.937 figure applies to the four-cell fit only. MG3's 0.073 lies inside both.
- **Stationary check:** no rejection at 0.01 (smallest p = 0.0125, Ω = 110). Power against ×2 / ×0.5 is 0.04–0.26 in
  every cell, so the check is **uninformative** everywhere. Landmark rates (100, 150) and the late/early ratios are in
  the record as nominal, descriptive sensitivities (not calibrated); none of the ratios excludes 1.
- **Cost:** 26.0 core-hours of production plus 0.75 calibration. Part 2 used 9.94 core-hours, inside the original
  12-hour ceiling, so the PI's extension was not exercised.
- **Deviation:** the evaluator was written after the data were in. It implements the registered Part 2 analyses. Its
  open details (landmark exposure, early/late split at 100, the stationary power simulation, the slope fit) are listed
  in the record for review.
Replication (ChatGPT): `review/MG008_chatgpt.md` and `review/MG008_replication/` (2026-10-30). Supports accepting MG3b as an estimation result.
Independently reproduced all rates, exposures, stopping decisions and slope estimates, all 302 confirmation
classifications (206 persistent, 96 excursion) with a separate deterministic solver, the stationary-check powers, and
the budget replay (every Part 2 batch fits under the original 12-hour ceiling). It asked for two label clarifications,
both applied (`review/MG008_response.md`). No number changed, and no lattice run is needed.
PI sign-off (design): Richard Benear, 2026-10-29, project thread ("I'm satisfied. go ahead and register it and modify
GATES.md etc."), approving draft 3 after ChatGPT's two review rounds.
PI sign-off (result): Richard Benear, 2026-10-30 12:50 UTC, project thread ("I sign off on MG3b"), after ChatGPT's
review (MG008) and the two label clarifications. MG3b is accepted as an estimation result; MG3's FAIL is unchanged.

---

### MG4 — Particles on frozen irregular 2D graphs (noise-free)
Status: **CLOSED 2026-10-31 — FAIL** (result below; PI sign-off). Registered 2026-10-30 (MG009c), with PI approval. **No MG4 run of any kind preceded this commit.** The only
computation behind the design was graph statistics (degrees and edge-length moments of random point sets, with no
reaction or diffusion). The registered design is `spec/MG4_proposal.md` **draft 4**, with its content frozen as committed
at `773255d`. Any change after this point needs a versioned amendment. Reviews and responses: `review/MG009_chatgpt.md` and
`_r2.md`, with `review/MG009_response.md` and `_r2.md`. Drafts 1–3 are preserved as superseded (`spec/MG4_proposal_v1.md`,
`_v2.md`, `_v3.md`).

Scope (PI decision, 2026-10-30 13:41 UTC): 2D graphs, noise-free only. Survival under token noise is deferred, and so are
Wolfram hypergraph snapshots (MG4b). MG4's verdict applies to the Ω → ∞ mean field of rules v1.1 on these graphs.

Summary of the registered design (the proposal governs where this summary is shorter):
- **Setting:** b = 29, eq. 17 parameters, rules v1.1 unchanged. Every family uses one rectangular torus commensurate with
  the triangular lattice: 122 × 142 cells, 19.92 × 20.08 at h_t = 0.2. Each graph's h is calibrated from its
  homogenized diffusivity tensor, h² = (e₁ + e₂)/2, and accepted within ±3% of h_t.
- **Families:** T is the triangular comparator and positive control. D is Delaunay of a periodic Poisson set, the
  primary family. R is a random geometric graph with mean degree 12, largest component only. N3 is random 6-regular,
  the negative control for C3(a). D and R each use a fixed panel of 4 graphs × 8 placements.
- **Part 0:** graph construction and calibration, spectra, T's existence window and depth check, and the box check
  (20 against 30). It writes `results/MG4/part0/decisions.json`.
- **Part A:** the deterministic graph ODE to t = 300, started from the MG1 continuum 2D soliton. Measurement follows
  proposal §5.1 and outcome classes follow §5.2, with precedence multiplied > decayed > single > unresolved.
- **Criteria:**
  - C1 (survival panel): PASS if ≥ 28/32 and every graph ≥ 6/8. FAIL if ≤ 20/32 or any graph ≤ 3/8.
  - C2 (shape against T): componentwise medians, with (a) L2 < 10%, both depth measures ≤ 10%, and λ ≤ 5%, plus ≥ 80% of
    runs meeting (a).
  - C3(a) (spectral count of −L_G against the exact finite-torus count): the band mean over [0.9, 1.1] k* within ±5% on
    every graph, with all 21 ratios saved.
  - C3(b) (coarse-grained action on the soliton profiles): panel median ≤ 10% for G, X and Y; FAIL if any exceeds 20%.
  - N3 must fail C3(a).
- **Gate precedence:** proposal §7. In order: Part 0 stop; budget stop; implementation failure; a failed positive control
  makes the gate INCONCLUSIVE; FAIL dominates INCONCLUSIVE; otherwise INCONCLUSIVE; otherwise PASS.
- **Implementation checks:** the node-sum balances with their derivative terms, with |r| ≤ 10⁻¹⁰ + 10⁻⁹·M. Zero sums are
  required only for stationary runs.

Prediction (descriptive; not a criterion): D and R pass C1–C3. Disorder pins the particle (drift under one λ₀). N3 fails
C3(a). The least certain parts are C2 and C3(b) on D.

Null and positive controls: N3, as above. T is the positive control (window, depth, C3 on T, and a single particle at all
4 placements). The homogeneous state must be an exact fixed point, and a 10⁻³ perturbation must decay.

Budget: run caps as in proposal §9. Core-hour ceilings, checked before every batch: Part 0 8, Part A 25, the h = 0.3
check 3. The projection is about 14 core-hours. The cost of the 2D graph ODE is timed in Part 0, and the projection is
re-stated here before Part A starts.

Pre-run implementation specification: anything the code needs that the proposal leaves open is appended here, dated,
before the run it affects. (The MG3b lesson: record choices before the run, not after.)

**Part 0 implementation specification (2026-10-30, MG010; recorded before any Part 0 run; code `src/graphs.py`,
`src/graph_ode.py`, `src/measure_graph.py`; tests `tests/test_mg4.py`, 9 pass).**
- **Streams:** `runner.seed_for(5, 0, family, stage, replicate)` seeds a numpy `default_rng` for the point sets (stage 0)
  and the placements (stage 1). For N3 it is passed to `networkx.random_regular_graph`; networkx is version 3.x, recorded in
  `decisions.json`. Replicates: panel graphs 0–3, pilot graphs 90–93, replacements 50 onward in order, N3 graphs 0–1.
  Placements: graph × 100 + i.
- **Point counts:** pilots use D: N = round(2.39·A/h_t²) and R: N = round(c·A/h_t²) with c = 12²/(8π) (from h² = k̄r²/8;
  design-stage values). Each pilot gives c_F = (point density)·h². The panel uses the mean c_F over the 4 pilots.
  - R's connection radius always keeps the mean degree at 12 for the realised point density: r = √(12/(π ρ)). Its h and
    acceptance refer to the largest component.
- **Corrector:** the minimiser of Σ_edges (e·Δx + χ_J − χ_I)², by sparse LU with one node pinned, which gives the tensor
  D = FᵀF/N. On T this reproduces 1.5a² exactly in the tests.
- **Continuum profile** (initial condition and C3(b)): `mg_continuum.run(2, R=50, h=0.025, t_end=300)`, seed scale 1.
  The stationary state at t = 300 is cached in `results/MG4/part0/continuum2d.npz`: φ_X(0) = −13.594, max |∂φ/∂t| = 5.5e-12.
  - It is interpolated onto nodes by a cubic spline with zero slope at r = 0, and set to 0 beyond r = 50 (never reached on
    a 20 × 20 torus).
  - C3(b)'s continuum ∇²f is that solver's own conservative radial Laplacian, interpolated the same way.
- **ODE:** BDF with a sparse analytic Jacobian (checked against finite differences in the tests), rtol 1e-8, atol 1e-10,
  max step 0.5, and the 1D first-step retries. Fields are saved at the 11 checks. Stationarity is judged from max |RHS|
  at t = 300.
- **Measurement:** the smoother is a row-normalised sparse matrix, self weight counted once. Candidate, merge, centroid
  and link radii are 1.0, 1.5, 0.6 and 1.0, as registered. Positions are wrapped into the torus for neighbour searches.
- **T's threshold:** A*_s(T) is the mean, over T's 4 placements, of the minimum smoothed φ_X at t = 300, and
  θ = −0.5·|A*_s(T)|. T placements are the node nearest the torus centre, the midpoint of its +x edge, the centroid of
  its upper-right triangle, and one uniform random point (stage 1, replicate 0).
- **Positive control:**
  - Depth check: the unsmoothed minimum φ_X of the T particle at t = 300 (node placement), against −13.6, within 20%.
  - Window: start from that state and run seed-free to t = 300 at each b.
    - b_s: bisection on [26.0, 29.0], 10 steps (±0.003). A particle "exists" if at least one particle is detected at
      t = 300.
    - b_r: bisection on [29.0, 30.04], 9 steps (±0.002). "Replicating" means two or more particles at t = 300.
    - Both edges are operational finite-time estimates, as in 1D (N5).
  - Box check: the same on a 184 × 212 torus (30.05 × 29.98). Relative L2 of the radial φ_X profile over r < 6.
- **Spectra:** shift-invert Lanczos at σ = −10⁻³ (the shift keeps the factorisation non-singular), growing k until the
  cutoff 1.5 k*² is passed. Eigenvectors are kept for the plane-wave diagnostic.
- **Smoke run (code validation, not a registered run):** one T particle on a small 62 × 70 torus (10.1 × 9.9) at h = 0.2,
  on the 2-core cloud machine.
  - It took 80 s for 4 340 nodes. Its centre φ_X was −13.55 at t = 300.
  - It was not stationary at t = 300 (max |∂φ/∂t| = 6e-6), probably slow drift towards a lattice pinning site.
  - The balances held to 4e-12.
  - It is recorded here because it bears on the *stationary* label (§5.2), which is reported only.
- **Timing:** Part 0's first full-size T run is timed, and the Part A projection is re-stated here before Part A.
- **Driver (`src/run_mg4_part0.py`), added 2026-10-30 before the Part 0 run:**
  - It is resumable by step: one JSON per step in `results/MG4/part0/steps/`, plus a ledger `budget_ledger.csv`.
  - Replacements for a panel slot come from the family's unused reps 50, 51, … in order.
  - N3 has the D panel's mean node count (rounded up to even) and the D panel's mean h.
  - Spectra are computed for T, the 8 panel graphs and both N3 graphs. The window bisection uses the threshold from
    T's node-placement reference, θ_ref = −0.5·|smoothed depth|. Part A uses the registered 4-placement mean.
  - It applies the §3 stop rule automatically (`stop_for_PI`). If the box check fails, it records `torus_for_partA =
    30x30`, and the panel would then need rebuilding on that torus. That would be stated here before Part A.
- **Driver smoke test (disclosed; scratch, not committed, not a registered run).** The full driver was run once with
  `MG4_SMOKE=1`: a 62 × 70 torus (10.1 × 9.9), box 74 × 84, and only 2 bisection steps. It exercised every step.
  - Its numbers do not apply to the registered 122 × 142 torus. They are disclosed because one of them previews a
    criterion: the C3(a) band-mean count ratio was 1.058–1.065 on the four small D graphs and 1.066–1.096 on the four
    small R graphs, against 1.015 on T and 0.031 on N3.
  - It also gave: the T depth −13.59 (smoothed −10.30); the box difference 3.6% between those two small tori; T's
    max |∂φ/∂t| at t = 300 of 6e-6; c_F of 1.89 (D) and 3.85 (R); a T ODE run of 86 s for 4 340 nodes on 2 cores.
  - **No criterion, threshold or rule is changed in response.** The registered Part 0 and Part A run as written.
- **Part A driver (`src/run_mg4_partA.py`), added 2026-10-30 before any Part A run:**
  - **Order:** T's 4 survival runs come first. They fix θ = −0.5·|mean minimum smoothed φ_X at t = 300|, written to
    `threshold.json`. Then come the D and R survival runs, the seeded runs (T, D, R), N3, and the h = 0.3 check.
  - **Placements:** D/R slot s, placement i: one uniform point from `seed_for(5, 0, fam, 1, 100·s + i)`.
  - **h = 0.3 check:** its D graph is drawn from stream rep 200 at the Part 0 density constant c_F(D), with placements
    from replicate 200·100 + i. That graph's own calibrated h is used and reported, with no acceptance test.
  - **N3 (descriptive):** 2 runs per graph. The initial profile is the continuum profile at r = h × hop distance from
    node 0 or 1. Its final minimum φ_X is reported. There is no detector, since N3 has no embedding.
  - **Records:** one JSON per run, plus a compressed field file holding the final G, X, Y, and φ_X and smoothed φ_X at the
    11 checks (float32).
  - **Analysis:** done after the runs (`analyse`, `c3b`), from the saved fields with no re-integration.
    - The core depth is the [0, 0.5) bin mean of φ_X around the final centre, unsmoothed and smoothed.
    - The T reference profile is the mean over T's single-particle placements, at the same bins.
  - **Field files:** an R run's file is a few MB. Whether all of them are committed (or only the summaries and final
    fields) is decided with the PI before the results PR.
- **Part A driver smoke test (disclosed; scratch, not a registered run):** one T run and one D run on the small
  62 × 70 smoke torus, with the smoke-test graphs.
  - Small-torus outcomes, again only a preview: T *single* (drift 0.03, λ 3.14); the D run *multiplied*. C3(b)
    smoothed e_X was 0.032 on T and 0.276 on the small D graph (unsmoothed 0.06 and 2.6).
  - **No criterion, threshold or rule is changed in response.** No further dynamics will be run on non-registered
    graphs before Part A.
- **Solver speed-up (MG010, recorded 2026-10-30 before any use in a registered run):** the cost of a run is mostly
  SciPy's sparse LU of the BDF Newton matrix, whose default ordering (COLAMD) gives 50–90 times the Jacobian's nonzeros
  on these graphs.
  - From now on `graph_ode.integrate` (default `order="nd"`) renumbers the unknowns node by node in a geometric
    nested-dissection order (`graphs.nd_order`) and keeps that order in the LU (`permc_spec="NATURAL"`).
  - The equations, the solver (BDF), the tolerances, max_step and the first-step retries are all unchanged. Only the
    elimination order and therefore the rounding differ. N3 has no embedding and keeps the original path.
  - **Checks:**
    - The full Part 0 T reference run (122 × 142, t = 0–300) was repeated with the new path. It took the same steps
      (nfev 1634 in both; nlu 229) and its final state differs from the saved `Tstate_122x142.npy` by at most 1.4e-13
      (field values up to 13.6). Wall time went from 1332 s to 283 s on the cloud machine
      (`results/MG4/part0/nd_check_T122x142.json`).
    - On T for t = 0–20 the steps are again identical (nfev 510, nlu 47) and the largest difference is 5.5e-14.
    - A new unit test (`test_nd_order_matches_default_solver`) compares the two paths on a 900-point scratch Delaunay
      graph for t = 0–2 (identical step counts, difference < 1e-10). This is a solver-equivalence check, not a dynamics
      outcome. On R panel graph 0 (38 595 points), t = 0–5: identical steps (nfev 800, nlu 68), largest difference 1.9e-14, wall 1967 s → 740 s (2.7 times; cloud machine, shared with the Part 0 run).
  - Part 0 evaluations already finished keep their records. The remaining window-bisection steps and all of Part A use
    the new path. Ledger rows carry the platform, as before.

- **Part 0 budget overrun (deviation, recorded 2026-10-30 after Part 0 finished):** Part 0 used 10.57 core-hours
  against its 8 core-hour ceiling, which the registration makes a stop trigger checked before every batch.
  - The total includes three estimated rows (2.04 core-h) for work lost when the cloud's 2-hour background limit
    stopped the driver. The measured total is 8.53 core-h. Including the estimates, the ledger passed 8 at about
    19:20 UTC, during the window bisection.
  - Cause: the Part 0 driver never checked the ceiling (Claude's implementation error), and the cost of the 2D graph
    ODE was underestimated in the design. Nothing in the outcome depended on it, and Part 0 scores no criterion.
  - The Part A driver now checks its ceiling before every job is started (`run_mg4_partA.py`, `CEILING`).
  - **PI decision, 2026-10-30 23:00 UTC (decision card in the project thread):** "Record as exception, continue." MG4 is
    not stopped. The overrun stands as a recorded exception, as MG3's did, and goes to ChatGPT's review with the result.

- **Part A budget re-statement and scope (PI decision, 2026-10-31 00:27 UTC, decision card in the project thread:
  "Gate runs only, 80 core-h"; recorded before any Part A run):**
  - **Measured cost** (t = 0–300, one run, new solver path unless stated):
    - Cloud: T 283 s, D panel graph 0 1073 s, R panel graph 0 about 7 900 s. The R run was timed in segments, and the
      cloud's 2-hour limit stopped it at t = 200 (7 194 s).
    - PI's PC, original solver path: T 990 s, D 1 944 s, R 13 188 s (`results/MG4/timing/timing_RB_LibraryPC.json`).
    - Projected on the PC with the new path: T about 4 min, D about 13 min, R about 1.6 h per run.
  - **Runs (cap 80):** the gate-deciding runs only, i.e. T 4, D 32 and R 32 with the soliton initial condition. Then the
    h = 0.3 check, which is cheap: T 4 and D 8, with its own ceiling of 3 core-h unchanged.
  - **Ceiling:** Part A at h = 0.2 is raised from 25 to **80 core-h** (projection about 60). The driver checks it before
    every job is started.
  - **Not run:**
    - The 68 seeded-formation runs (descriptive; no criterion uses them) are deferred, at a projected cost of about 60
      core-h, mostly on R. They may be registered later as a separate study.
    - The 4 N3 dynamics runs (descriptive) are dropped. On a random 6-regular graph the sparse LU of the Newton matrix
      is nearly dense: measured fill at n = 1 000–4 000 is about 56% of dense, with time growing as n³. At N = 18 902
      that is about 20 GB and hours per factorisation. N3's registered role, failing C3(a), is met by its spectra
      alone; the outcome is scored at evaluation.
  - **Machine:** the PI's PC (Windows), `MODELG_WORKERS=9`. The platform is recorded in every ledger row.
  - No criterion, threshold or precedence rule changes.
- **Evaluator test and disclosure (2026-10-31, while Part A runs on the PI's PC):** `src/evaluate_mg4.py` was tested
  end to end on fabricated fields: the continuum soliton placed at each registered centre, with one run zeroed so
  that it reads as decayed. No dynamics were involved.
  - C3(a) and C3(b) do not depend on the dynamics: they use only the graphs, the placements and the continuum profile. So
    this test computed their registered values, which are disclosed here before any Part A result is seen.
    - C3(a) (Part 0 spectra): band means D 1.051–1.052 (panel 1.051) and R 1.059–1.068 (panel 1.064). These are outside
      1 ± 0.05, so FAIL for both. T 1.008 passes. N3 0.008 fails, as the negative control requires.
    - C3(b) panel medians of e_S (G, X, Y): D 0.27, 0.31, 0.30 and R 0.51, 0.52, 0.52. These are above 0.20, so FAIL for
      both. T 0.007, 0.010, 0.009 passes.
  - Reading (not a re-score): C3(b) compares the graph operator's local action on a smooth profile with ∇² at the
    calibrated (homogenized) h. That local action scales with the naive second moment of the edges, not with the
    homogenized one. So e_S is expected to be about (h_naive/h)² − 1: 0.27 on D and 0.47 on R, close to what was
    measured. The homogenization corrector acts on the solution, not on a fixed profile. C3(b) as registered therefore
    could not pass on graphs where h_naive differs from h. This is a design flaw, recorded as such. The verdict stands
    as registered.
  - Under §7, unless T fails as a positive control or an implementation check fails, MG4's verdict will be **FAIL**
    (C3 for D and R). Part A still decides C1 and C2, which are reported in full.
- **Part A budget overrun (deviation, recorded 2026-10-31 after the runs finished and before evaluation):** all 80
  runs completed on the PI's PC (RB_LibraryPC, 9 workers), from 00:30 to about 12:00 UTC.
  - The gate group used 90.52 core-h over 68 runs, against the 80 ceiling. The h = 0.3 group used 1.04 core-h over 12
    runs (ceiling 3).
  - The driver checks the ledger before it starts each job, and every job started below 80. The ledger passed 80 at
    10:09 UTC (after R3_p2), when the last R3 runs were already running. Those runs added about 10.5 core-h.
  - No run is missing, so this is a cost overrun, not a budget stop.
  - **PI decision, 2026-10-31 12:23 UTC (decision card): "Record as exception, evaluate."**
  - Data: the run records, ledger and log are committed. The field files are committed trimmed to the arrays the
    evaluator uses (`fields_eval/`: final and checks_sX, 96 MB, identical arrays). The full files, which also hold the
    unsmoothed checks_X, stay on the PI's machine (138 MB).

Result: **FAIL** (evaluated 2026-10-31, `results/MG4/evaluation.json`). The evaluator first reported NOT SCORED, because
11 runs had balance flags. Those flags are kept in the output (`balance_failures`). After the rule-3 investigation
below, their resolution was recorded in `results/MG4/balance_resolution.json`, which the evaluator reads; it then
scores FAIL (MG010 response to review, 2026-10-31). ChatGPT accepted FAIL (`review/MG010_chatgpt.md`).
- **Every D and every R run multiplied.** C1 is 0/32 single on both families (0/8 on every graph), so C1 FAILs for D and R.
  - Separately written detection code (`src/recompute_mg4.py`, `results/MG4/recompute.json`) gives the same classes for
    all 72 runs.
  - The particle copies itself until the torus is full: 33–37 particles at t = 300 (median 35) on both families. The
    count is already constant from t = 200 in the runs checked (D0_p0 37, D2_p5 36, R1_p3 36).
  - The h = 0.3 D check (reported only) behaves the same: 8/8 multiplied, 34–37 particles.
- **T, the positive control, passes:** single particle at all 4 placements, stationary (max|dφ/dt| ≤ 2.3e-9), C3(a)
  1.008, C3(b) medians 0.007/0.010/0.009. The h = 0.3 T runs are single at 4/4.
- C2 is INCONCLUSIVE for D and R: there are no single-particle runs to compare.
- C3 FAILs for D and R in both parts (values disclosed above, before the Part A results). N3 fails C3(a), as required.
- **Implementation check (precedence rule 3) — investigated before scoring.**
  - The float64 node-sum balance r2 exceeded 10⁻¹⁰ + 10⁻⁹·M in 11 R runs, with |r2| between 1.0e-10 and 2.6e-10.
    Every r1 passed.
  - Recomputed in extended precision (np.longdouble) on the saved final states, r1 is about 1e-18 and r2 about 1e-15 in
    all 11 flagged runs (`src/check_mg4_balances.py`, `results/MG4/balance_investigation.txt`, which lists all 11).
    In float64 the same states give r2 anywhere from 3e-12 to 3e-10, varying from state to state. ChatGPT's independent
    40-digit check of all 11 gives residuals below 2e-35.
  - Finding: the balance identities hold algebraically for the implemented right-hand side. The flagged residuals are
    float64 roundoff in the node sums, not a conservation error. The registered absolute floor of 10⁻¹⁰ lies below
    that roundoff level on the degree-12 R graphs. No run is excluded, and scoring proceeds.
- Precedence: no Part 0 stop. Both budget overruns are recorded exceptions (PI). The rule-3 investigation is cleared as
  above. The positive control passes. FAIL dominates (C1 and C3 for both D and R). **MG4: FAIL.**
- **Prediction refuted:** D and R were predicted to pass C1–C3.
- **Reading (descriptive, not a criterion; hypothesis, not established):** one possible explanation is that the
  irregular graphs shift the replication boundary below b = 29.
  - On T, b = 29 sits 0.2 below that boundary (29.20). Above it, T fills the torus in the same way (31–35 particles,
    Part 0 window records).
  - MG4 did not measure where the boundary lies on D or R. The C3(a) excess of modes (+5% on D, +6% on R) shows that the
    graph spectrum differs from the continuum's near the particle's scale. On its own it does not show slower
    particle-scale diffusion, and it does not show the mechanism.
  - Whether D or R have a single-particle window at other b would need a new registered study.
- **Scope of the FAIL:** MG4 rejects single-particle preservation on unweighted D and R graphs under the registered
  conditions (b = 29, eq. 17 parameters, h-calibration as registered). It does not rule out localized solutions at other
  parameters or on weighted graphs.
- Cost: Part 0 10.57 core-h, Part A 90.52 + 1.04 core-h (PI's PC); both overruns recorded as exceptions.
Replication (ChatGPT): MG010 review 2026-10-31 (`review/MG010_chatgpt.md`). It checked all 80 runs at all 11 times with no
classification disagreement, reproduced C3(b) and checked the 11 balance flags at 40 digits. It accepts FAIL, with three
corrections, all applied (`review/MG010_response.md`).
PI sign-off (design): Richard Benear, 2026-10-30 14:19 UTC, project thread ("I approve MG4 draft 4."), after ChatGPT's
two design reviews (MG009 r1 and r2; the reviewer was comfortable with registration once draft 4's definitions were in).
PI sign-off (result): Richard Benear, 2026-10-31 13:52 UTC, project thread ("I will sign off on MG4"), after ChatGPT's
MG010 review and the three corrections (PR #17). ChatGPT's replication files are in `review/MG010_replication/`.

### MG4c — Single-particle window on unweighted Delaunay graphs (follow-up to MG4; estimation study)
Status: **CLOSED 2026-10-31: W1 + O1** (PI sign-off 20:24 UTC). Registered 2026-10-31 (MG011c), with PI approval. **No MG4c run of any kind preceded the registration commit.** The
registered design is `spec/MG4c_proposal.md` **draft 3**, with its content frozen as committed at `a946550`. Any change
after this point needs a versioned amendment.
- Reviews and responses: `review/MG011_chatgpt.md`, `review/MG011_replication/` (r1–r3), `review/MG011_response.md`
  and `review/MG011_response_r2.md`.
- Drafts 1 and 2 are preserved as superseded (`spec/MG4c_proposal_v1.md`, `_v2.md`).
- MG4c does not re-score MG4. The name MG4b stays reserved for Wolfram hypergraph snapshots.

Summary of the registered design (the proposal governs where this summary is shorter):
- **Graphs and runs:** MG4's D0 and D1 panel graphs, placements p0 and p1 ("the four D runs"), with MG4's solver,
  initial condition (b = 29 continuum soliton), detector and classes. At each tested b, one T run supplies the
  threshold.
- **Threshold:** a T run is an eligible reference only if it is single under the fixed θ_elig = −1.0 and its smoothed
  minimum at t = 300 is ≤ −2.0.
  - If eligible: θ_b = ½ × that minimum. T must also be single under θ_b.
  - Otherwise: the fixed θ_29 = −5.117354, with the T reference marked invalid.
  - Each D run's class under θ_29 and its smoothed minimum are reported as well.
- **Part 1:** b ∈ {28.75, …, 27.0} in steps of 0.25, with T × 1 and D × 4 at each b. MG4's Part A runs are the b = 29
  point.
- **Part 2:** at most 6 bisection steps per edge, each step T × 1 + D × 4, with these rules:
  - Case A: only the highest all-single block's edges, on the predicate "all-single". Unbracketed edges are reported as
    such.
  - Case B: on "at least one multiplied"; an all-single discovery is W1, with its edges "not estimated".
- **Part 3:** a diagnostic amendment to rules v1.1. The cotangent (lumped-mass P1) Laplacian on D0 p0 at b = 29, plus a T
  control that must reproduce MG4's T p0 field to 1e-6. Construction checks run first, and the area-weighted balances
  apply.
- **Outcomes (no PASS/FAIL):** W1/W2/W3, INCOMPLETE, NOT SCORED; O1/O2/O3, O-invalid. Readings are limited as written
  in proposal §5. Programme policy: if the result is W2 with O2, stop and reconsider before MG5.
- **Budget:**
  - Run cap 102 (40 + 60 + 2).
  - 20 core-h is a launch cutoff: in-flight jobs finish. About 22.5 core-h is the projected completion allowance, an
    estimate, not an enforced maximum.
  - Runs on the PI's PC with `MODELG_WORKERS=9`.
- **Balance flags:** resolved only by a recheck in genuine higher precision (mpmath at 40 digits, or exact fractions)
  of the run's own b and operator, recorded before outcomes are read.

Pre-run implementation specification: appended here, dated, before the run it affects.

**Pre-run implementation specification (MG012, 2026-10-31; committed before any MG4c run on a registered graph).**
Code: `src/run_mg4c.py` (driver), `src/evaluate_mg4c.py` (evaluator), `graphs.cotan`, the `op` and `w` arguments of
`graph_ode.make_system`/`integrate`/`balances`; tests `tests/test_mg4.py` (cotangent construction) and
`tests/test_mg4c.py` (decision rules on constructed class tables, threshold eligibility, the precision recheck, and
synthetic end-to-end passes of driver plus evaluator through Case A, Case B and "multiplied at every tested b").
1. **Order.** The Part 3 construction checks (proposal §4) run first, in the main process, on D0 and on T, and are written
   to `results/MG4c/part3_checks.json`. If any fails, the two Part 3 runs are not started and Part 3 is O-invalid. The 40
   Part 1 runs and the 2 Part 3 runs then run as one batch. Then the balance recheck, then the grid table, then Part 2.
2. **Run keys and b values.** `b<b to 6 decimals>_<T|D0|D1>_p<0|1>`; Part 3 keys carry the prefix `w_`. The run itself
   uses the exact b, which its record stores. A bisection midpoint is ½(inside + outside) rounded to 9 decimals. Six
   steps from a 0.25 bracket keep tested values at least 0.0039 apart, so 6-decimal keys are unique.
3. **b = 29.** The four D classes are MG4's registered classes (`results/MG4/evaluation.json`, keys `D{0,1}_p{0,1}_sol`,
   under θ_29); T is the eligible reference there (MG4 Part A). No b = 29 run is repeated, except the Part 3 runs.
4. **Classification.** Exactly MG4's: `measure_graph.particles` and `classify` on the smoothed φ_X at the 11 checks,
   stored as float32 as in MG4 Part A. The smoothed minimum at t = 300 used for eligibility and θ_b is taken from the
   float64 smoothed field before storage. The final state (G, X, Y) is stored in float64 (MG4 stored float32).
5. **Part 3.** Both Part 3 runs are classified with θ_29 (b = 29, where MG4's T reference is eligible). The T control
   compares its final φ_X with `results/MG4/part0/steps/Tstate_122x142.npy`, MG4 Part 0's T reference run at b = 29 from
   the same centre node and initial condition (MG4's T p0 placement). The criterion is max |Δφ_X| ≤ 1e-6, plus class
   single; the maximum over all three fields is reported too. That reference was integrated before the
   nested-dissection ordering was introduced, so rounding-level differences are expected. On T, h² = 0.04 = 1.5 a² (a =
   0.1633), so L_w equals MG4's L_G.
6. **Balance recheck.** The float64 check is unchanged (node sums, A-weighted for L_w). A flagged run is rechecked from
   its stored float64 final state, its own b and its own operator, rebuilt from off-diagonal weights so that its rows sum
   to zero exactly (1/h² per edge, or the cotangent weights over A), in `mpmath` at 40 digits (exact `Fraction` if
   `mpmath` is missing). The flag is resolved iff |r| ≤ 1e-25·(1 + M) for both identities. Every recheck is written to
   `results/MG4c/resolution.json` before any class is read. An unresolved flag stops the driver: NOT SCORED.
7. **Budget.** `results/MG4c/budget_ledger.csv` holds one row per finished run (integration wall time, as in MG4). Before
   each job starts, the driver stops if the ledger has reached 20 core-h, or if the finished plus in-flight runs have
   reached 102. In-flight jobs finish and are recorded. Classification and rechecks in the main process are not in
   the ledger (seconds per run).
8. **Stops and failures.** A stop during Part 2 keeps every completed step in `driver_log.json` (`refined`); the
   evaluator reports the brackets from the recorded steps. A crashed run stops the driver without a record; it is
   rerun after investigation, and the crash is reported.
9. **Evaluator.** It recomputes every class from the stored fields with the registered thresholds and stops if a class
   differs from the driver's log. Statuses are applied as in proposal §5, with *tested b* meaning every b whose five
   runs all exist, plus b = 29. Depth and λ against T (as in MG4's `analyse`: smoothed depth at the particle centre, λ
   from the unsmoothed φ_X profile) are computed only for all-single b at which T is an eligible reference.
10. **Platform.** Each run record and ledger row carries the machine, OS, Python and numpy versions.
11. **Disclosure.** While writing this specification, Claude started a cloud trial of the Part 3 T control (a registered
    run) to check the 1e-6 match ahead of the PC run. It was stopped after about 3 minutes, before it finished, because
    running it before this specification was committed would break rule 1. It produced no output, no record and no field
    file, and its scratch folder was deleted. No other MG4c dynamics have run on a registered graph; the tests use small
    unregistered graphs or synthetic fields.

Result: **W1 (single-particle points found) and O1 (the weighted operator keeps D0 p0 single at b = 29).** Run on the
PI's PC (Windows, Python 3.13.2, 9 workers), 2026-10-31, from main at 494403a. All 102 planned runs completed (Part 1:
40, Part 2: 60, Part 3: 2), 11.33 core-h, below the 20 core-h cutoff. No float64 balance flags, so no recheck was needed.
Records: `results/MG4c/` (`evaluation.json`, `driver_log.json`, `runs/`, `fields/`, `budget_ledger.csv`) and
`results/MG4c_run_log.txt`. Claude recomputed every class, threshold, smoothed minimum, max|dφ/dt|, balance and the T
control independently from the field files (`src/recompute_mg4c.py`); all agree with the evaluator.
- **Window (Case A).** The four D runs are all single at 13 tested b, from 27.2656 to 28.7266. The b = 29 point
  (MG4) is all multiplied.
  - Upper edge bracket: [28.7266, 28.7305]. It is set by one run, D1 p1, which multiplies from 28.7305 up. D0 p0, D0 p1
    and D1 p0 stay single through 28.75.
  - Lower edge bracket: [27.2617, 27.2656]. It is set by one run, D1 p0, classed decayed from 27.2617 down. At 27.2617
    that run had not settled (smoothed min −4.90 against θ −5.12; max|dφ/dt| 4.9e-2 at t = 300). At 27.0, D1 p1 also
    decays. D0 p0 and D0 p1 stay single down to 27.0.
  - Every D class is the same under θ_b and under θ_29.
- **Threshold.** T is an eligible reference from 27.375 to 28.75. Below that, T decays completely (smoothed min about
  1e-15), so the classes at 27.3125 and below use the fixed θ_29. **The lower edge is therefore conditional on the
  fallback detector, and D keeps a single particle at b values where the T lattice does not.**
- **Shape against T** (eligible all-single b only): smoothed depth ratio D/T 0.90–1.14, and λ ratio 0.96–1.01.
- **Part 3.** All construction checks pass. The T control is single and matches MG4's T reference field to 4.5e-13
  (criterion 1e-6). D0 p0 with the cotangent operator is single at b = 29, where it multiplied with the unweighted
  operator.
- **Reading (proposal §5, limited to what was tested).** At b ≈ 27.27–28.73, D graphs keep a single particle under
  MG4's detector and initial condition, so b = 29 lies above D's single-particle range. That range is shifted below
  T's (27.33 < b < 29.20, MG4 Part 0) mainly at the upper edge, by about 0.47. For D0 p0 at b = 29, the unweighted
  operator matters. MG5 design may consider these b values. Both edges each rest on a single run, and the lower edge on
  an unsettled run, so the edge values are estimates within their brackets, not sharp properties of D graphs.
Replication (ChatGPT): **reproduced** (review MG012, 2026-10-31: `review/MG012_chatgpt.md`; files
`review/MG012_replication/`). All 102 runs and 1,122 observations were re-analysed. W1, both brackets and the weighted
single outcome reproduce exactly, and all balances pass. Response: `review/MG012_response.md`. Two items from the
review are recorded here:
- **Deviation (construction check wording).** Proposal §4 says L_w equals L_G on T "entrywise to 1e-12". The implemented
  check used 1e-12 × max|L_G| (= 1.5e-10). The pre-run specification did not state this. The measured max entrywise
  difference is 3.5e-12 (relative 2.4e-14). It passes the implemented check and fails a literal absolute reading.
  - In exact arithmetic both operators are identical on T (every entry 25 or −150; mpmath check in the response). The
    3.5e-12 is float64 rounding in the cotangent arithmetic.
  - The registered T control (final φ_X to 1e-6) passed at 4.5e-13.
  - Whether O1 stands with this deviation recorded is for the PI (sign-off below).
- **Limits of interpretation.**
  - 63 of 72 unweighted single classifications are stationary (max|dφ/dt| < 1e-8). The other 9 range from 2.2e-8 to
    1.4e-2, the largest next to the edges.
  - The weighted D0 p0 single is not stationary (8.3e-6).
  - The edges are finite-time (t = 300) detector boundaries for this initial condition, not established bifurcation
    points.
PI sign-off (design): Richard Benear, 2026-10-31 16:25 UTC, decision card in the project thread ("Approve and
register"), after ChatGPT's three reviews of drafts 1–3 (the round-3 reviewer supports registering draft 3 unchanged).
PI sign-off (result): Richard Benear, 2026-10-31 20:24 UTC, decision card in the project thread ("Accept O1, deviation
recorded"), after ChatGPT review MG012. **W1 + O1 accepted, with the T-operator tolerance deviation recorded above. MG4c
CLOSED.**

---

### MG5 — A particle on an evolving Delaunay graph (Phase 4, no coupling)
Status: **REGISTERED 2026-10-31 (MG014e)**, with PI approval. **No MG5 run of any kind preceded this commit.** The
computation behind the design was a re-reading of MG4c's stored fields, a solver timing and accuracy pilot on an
unregistered graph with no outcome read, and right-hand-side evaluations of the stored MG4c states (proposal §10 and
header, including the disclosure of one pilot script on a registered field that was stopped unread and its output
deleted). The registered design is `spec/MG5_proposal.md` **draft 3**, with its content frozen as committed at
`e2a25e3`. Any change after this point needs a versioned amendment.
- Drafts 1 and 2 are preserved as superseded (`spec/MG5_proposal_v1.md`, `_v2.md`).
- Reviews: `review/MG014_chatgpt.md`, `review/MG014_chatgpt_r2.md`, `review/MG014_replication/` (`review.md`,
  `review_r2.md`, checks). Responses: `review/MG014_response.md`, `review/MG014_response_r2.md`. The round-2 reviewer
  supports implementation and Part 0 once its two corrections are made; both are in draft 3.

Summary of the registered design (the proposal governs where this summary is shorter):
- **Rule (rules amendment v1.2, MG5 only):** field-blind spatial birth–death on the Delaunay triangulation. Nodes die at
  rate ν per node and are born uniformly at rate νρ (ρ = 18 901 / A). Events are applied synchronously in ticks Δ (§3.2
  boundary timeline). Births get barycentric interpolation of the boundary field on the pre-edit mesh. Ephemeral nodes
  are excluded. The mesh is rebuilt once per tick.
- **PI choices:** interpolating transfer; Part B included; unweighted operator with each panel graph's calibrated h.
- **Solver:** ROS2, dt = 0.05, with one Jacobian and LU per tick. Δ = 1 (ν ≤ 0.02, frozen) and Δ = 0.2 (ν = 0.2).
- **Part 0** (unregistered graph):
  - solver check against BDF (≤ 1e-3);
  - step check, dt 0.05 against 0.025, exact replay;
  - tick check, Δ 0.2 against 0.1, exact replay of events;
  - for the two replay checks, strict count agreement, plus centre ≤ 0.1 and depth ≤ 2% while single;
  - fallback and stop rules as in proposal §6.
- **Part A (b = 28, t = 500):** starts from MG4c's four stored b = 28 states. Frozen controls 4, plus ν ∈ {0.002, 0.02,
  0.2} × 4. θ = −4.457582.
- **Part B (b = 29, t = 300):** the same states with φ_Y shifted by −0.9/14. Frozen 4 and ν = 0.2 × 4. θ = −5.117354.
- **Statuses:** COMPLETE / FAILED-NUMERICAL (including any unresolved balance, jump-accounting or interpolation-audit flag)
  / INCOMPLETE, applied before classes. Classes: multiplied (split/remote) > dissolved > single > unresolved. Mobility is
  pinned or mobile at 1.5.
- **Outcomes:**
  - INVALID-CONTROL; NOT SCORED / INCOMPLETE.
  - Gate MG5 (survival): PASS / PARTIAL / FAIL, with ≥ 3/4 single per rate.
  - Mobility: PINNED / MOBILE / MIXED.
  - Part B: A1 / A2 / A3 / A-void, with A1 read as "recovery under rewriting plus transfer".
- **Data (PI: keep the repo small):**
  - Committed: replay specifications with event-list hashes, sparse check records (nodes below −1.0 with row indices),
    per-tick accounting, run records, ledger and evaluation.
  - On the PI's PC, in a committed SHA-256 manifest: event lists, full fields and final states.
- **Budget:** run cap 24; 40 core-h launch cutoff, with in-flight runs finishing. Projected about 20 core-h; restated after
  Part 0 times the complete implementation. Runs on the PI's PC.

Pre-run implementation specification: appended here, dated, before the run it affects.

**Pre-run implementation specification (MG015, 2026-11-01; committed before any Part 0 or registered MG5 run).**
Code: `src/evolve.py` (mesh, ROS2, events, tick edits, audits, detection, classification, flag rechecks),
`src/run_mg5.py` (engine `run_ticks`, registered driver), `src/run_mg5_part0.py` (Part 0), `src/evaluate_mg5.py`,
`src/recompute_mg5.py` (independent). Tests: `tests/test_mg5.py`, 11 tests, all pass on the cloud container (Python 3.13,
numpy 2.5, scipy 1.16).
1. **Order.** `python src/run_mg5.py part0` runs first and writes `results/MG5/part0/decisions.json`. The registered
   driver refuses to start unless that file says `proceed: true`, and takes dt and Δ_fast (Δ at ν = 0.2) from it. Part 0's
   outcome and the restated budget are appended here before Part A starts. Parts A and B then run as one batch of 24,
   longest (ν = 0.2) first.
2. **Streams.** Run seeds: `seed_for(5, b_idx, rate_idx, 1, 2·slot + placement)`, with b_idx 0 (Part A) or 1 (Part B),
   and rate_idx the index of ν in [0, 0.002, 0.02, 0.2]. The same seed drives the event list and the choice of audited
   births. Part 0: graph `seed_for(5, 9, 0, 0, 0)` (18 901 uniform points; h from its homogenized corrector) and events
   `seed_for(5, 9, 3, 1, 0)`. No registered study uses b_idx 9.
3. **Start states.** MG4c's D graph is rebuilt (`run_mg4c.graph_for`). The start mesh must have identical edges (asserted).
   Ids are MG4c's row indices; h is the panel's exact calibrated value. Part B adds S₀(28) − S₀(29) from `steady_state`
   (φ_Y only, −0.9/14).
4. **Tick loop (proposal §3.2).** For k = 1..K: one Jacobian and LU at the state on M_(k−1); Δ/dt ROS2 steps (exact
   integer); failure test; float64 balances (with recheck); check if t_k is a multiple of 10; snapshot if t_k is a
   multiple of 50 or t_end; then tick k's events (times in (t_(k−1), t_k], `searchsorted` right); then the Mesh is rebuilt
   and a new ROS2 system is formed. There are no edits after the final boundary.
5. **Detection.** MG4's `measure_graph.smoother` and `particles` on the current positions. The t = 0 detection on the
   start state gives c(0); a start that is not single is classed unresolved. The sparse record keeps nodes with smoothed φ_X
   < −1.0. `tests/test_mg5.py` checks that sparse re-detection equals full detection, including torus wrap and ties.
6. **Classes.**
   - Multiplied: ≥ 2 particles at two consecutive checks. *Split* vs *remote* is measured at the first check with ≥ 2
     particles: it is the distance, from the last single-particle centre (or c(0)), of the second-nearest particle (≤
     4.0: split).
   - Single requires each step ≤ 2.0. Mobility uses the unwrapped track.
   - D_p is the least-squares fit through the origin of MSD = 4 D_p τ, for lags 10–100.
7. **Flags (all resolved in the run, before the next tick; the original flag and its recheck are kept in the run
   record).**
   - *Balance:* MG4's float64 test. The recheck is mpmath at 40 digits on the operator rebuilt from edges; resolved iff
     |r| ≤ 1e-25 (1 + M).
   - *Jump accounting:* |I_after − (I_before − R + B)| ≤ 1e-10 (1 + |I_before| + |R| + |B|) per species. The recheck uses
     `math.fsum`; resolved iff ≤ 1e-13.
   - *Interpolation audit:* 20 effective births per tick (all if fewer), using `LinearNDInterpolator` on the points within
     1.0. Flag above 1e-9. The recheck is exact-rational over every effective birth of the tick; resolved iff ≤ 1e-12.
   - Any unresolved flag, a non-finite value or |φ| > 1e3 stops the run as FAILED-NUMERICAL.
8. **Replay specification.** The event list is generated in the worker. Its SHA-256 (columns t, type, id, x, y,
   little-endian, in that order) is compared with a regeneration in a fresh Python subprocess. A mismatch writes an
   INCOMPLETE record and the run does not start.
9. **Data.**
   - Committed: `results/MG5/runs/<key>.json` (status, flags, replay specification, driver summary, platform, code
     revision) and `results/MG5/checks/<key>.npz` (sparse records, with offsets, and the per-tick table), plus
     `budget_ledger.csv`, `local_manifest.csv`, `driver_log.json`, `evaluation.json`, `part0/`.
   - PC-only, in `results/MG5/local/` (gitignored): event lists and snapshots (float32 fields and smoothed φ_X every 50
     units; float64 at t_end), each with its SHA-256 in the manifest.
10. **Budget.** One ledger row per finished run (wall time). Before each start, the driver stops if the ledger has
    reached 40 core-h, or if finished plus in-flight runs have reached 24.
11. **Evaluation.** `evaluate_mg5.py` re-detects every check from the committed sparse records with the registered θ. It
    stops if any count sequence or class differs from the driver's record, then applies statuses, validity and outcomes
    (proposal §8; frozen control "pinned" means max d < 0.1). `recompute_mg5.py` repeats the classification with separately
    written detection code.
12. **Disclosure.** While writing this code, dynamics were run only on small unregistered graphs in the tests: a 6 × 5
    torus of 900 points for the ROS2 order test, and a 10 × 10 torus of 3 450 points, ν = 0.2 for 20 time units, for the
    end-to-end test, which asserts only the run status, accounting and record round-trips. The registered start states
    were loaded, and their meshes compared, without integration. No Part 0 run and no registered run has been made.

PI sign-off (design): Richard Benear, 2026-10-31 23:17 UTC, decision card in the project thread ("Approve and
register"), after ChatGPT's two reviews (MG014 r1 and r2) and draft 3's corrections.


**Part 0 outcome and restated budget (MG016, recorded 2026-11-01, after Part 0 and before any Part A or B run).**
Part 0 ran on the PI's PC (RB_LibraryPC, Windows 11 10.0.26200, Python 3.13.2, numpy 2.2.4, scipy 1.16.1) from main
`7ae167a`, 2026-10-31 17:52–20:11 local (8 366 s wall, at most two runs in parallel). Records: `results/MG5/part0/`
(`decisions.json`, four check files) and `results/MG5_part0_log.txt`. Unregistered graph `seed_for(5, 9, 0, 0, 0)`,
N = 18 901, h = 0.2000797; check events `seed_for(5, 9, 3, 1, 0)`, 756 401 events, SHA-256 `2e0490b2…1e313`.
- **Check 1 (solver):** ROS2 dt = 0.05 against BDF over 50 units, max|Δφ| = 1.52e-4 (criterion ≤ 1e-3). **Pass.**
- **Settle:** frozen to t = 300, COMPLETE.
- **Check 2 (step, dt 0.05 vs 0.025, ν = 0.2, Δ = 0.2, 100 units):** counts 1 at all 11 checks in both runs; largest
  centre difference 5.6e-4 (≤ 0.1), largest relative depth difference 0.09% (≤ 2%). **Pass.**
- **Check 3 (tick, Δ 0.2 vs 0.1, dt 0.05):** counts 1 at all 11 checks in both runs; largest centre difference 0.060
  (≤ 0.1), largest relative depth difference 1.34% (≤ 2%). **Pass**, with less margin than check 2.
- Flags in all four check runs: interpolation audit at most 3.6e-15, jump-accounting residual at most 8.8e-14, none
  raised.
- **Decision (registered policy, no PI choice needed):** proceed with dt = 0.05 and Δ = 0.2 at ν = 0.2 (Δ = 1 otherwise).
  No parameter changes.
- Independent recompute (Claude, cloud, from the committed check files with `recompute_mg5.detect`): the same counts,
  centre differences and depth differences as above.

**Restated budget.** Measured cost: 1.41 s per time unit frozen, 27.8 s per unit at ν = 0.2 (5.55 s per tick, Δ = 0.2),
about 6.7 s per unit at ν = 0.002 and 0.02 (estimate). Projection:

| Runs | Number | core-h |
|---|---|---|
| Part A frozen (t = 500) | 4 | 0.78 |
| Part A, ν = 0.002 and 0.02 | 8 | 7.42 |
| Part A, ν = 0.2 | 4 | 15.42 |
| Part B frozen (t = 300) | 4 | 0.47 |
| Part B, ν = 0.2 | 4 | 9.25 |
| **Total** | **24** | **33.3** |

This is about 1.7 times the design figure (20 core-h): the full tick (event handling, hash checks, audits, rebuild and
records) costs about twice what the design pilot timed. It stays under the registered **40 core-h launch cutoff**, and
the run cap stays at 24. Part 0 timed at most two runs at once; with 9 workers, per-run times may rise. The cutoff counts
finished runs only, and the expensive ν = 0.2 runs are launched while most of the cheap runs are finished or in flight,
so the cutoff is not expected to stop any run; if it does, the remaining runs are INCOMPLETE as registered. Expected wall
time with 9 workers: about 4–5 hours (one Part A ν = 0.2 run is about 3.9 h). Committed storage: at most about 11 MB.

**Budget and execution amendment 1 (MG017, recorded 2026-11-01, before the relaunch of Parts A and B; approved by the PI).**
What happened in the first launch (PI's PC, main `594415f`, 9 workers, started 06:23 local; an earlier start at 06:17
died after about 9 s with no run record and no ledger row):
- The driver stopped at about 13:55 local on a console control signal (`^C` in its output, no traceback, no Windows
  error or out-of-memory event). The cause is unknown. No one has read a field or a check from the interrupted runs.
- Two runs had finished, both COMPLETE: `A_nu0.02_D0_p0` and `A_nu0.02_D0_p1` (11 136 s and 11 279 s wall; ledger 6.23
  core-h). Their records are kept (branch `mg5-interrupted-1`, commit `9c678e0`) and they are not rerun.
- Nine runs were in flight and are lost: the eight ν = 0.2 runs (each at t = 250 of 500 or 300) and `A_nu0.02_D1_p0`
  (t = 200). They have no ledger row. About 65 core-hours of machine time went into them, off the ledger. Their
  snapshots are listed with SHA-256 in `results/MG5/interrupted_local_manifest.csv` and are not used. The driver resumes
  at run level only, so these runs restart from t = 0 with the same seeds and event lists; their results are
  unaffected.
- **Slowdown:** runs took about 3.5 times the Part 0 rate (ν = 0.02: about 3.1 h against 0.93 h; ν = 0.2: about 98 s per
  unit against 27.8). Part 0 timed at most two runs at once. With 9 workers, RAM was nearly full (3.2 GB of 31.7 GB
  free; commit charge 56.9 of 61.3 GB), workers kept most of their 1.7–2.1 GB out of RAM, and each used about 0.4 of a
  core while the machine was about 70% idle. Each worker also ran its own unrestricted OpenBLAS thread pool.
- **Error in the MG016 statement (Claude):** MG016 said the expensive ν = 0.2 runs would launch after most cheap runs.
  The driver did the opposite: it sorted by descending ν, so the frozen controls were last in the queue. With the 40
  core-h launch cutoff, which counts finished runs only, the frozen controls would most likely never have started, and
  the gate could not have been scored.

Changes (execution and budget only; no criterion, threshold, rule, seed or run changes):
1. **Order:** frozen controls first (A then B), then ν = 0.2, then ν = 0.02 and 0.002.
2. **Workers:** 5 (`MODELG_WORKERS=5`). The driver sets `OPENBLAS_NUM_THREADS`, `OMP_NUM_THREADS` and
   `MKL_NUM_THREADS` to 1 unless already set. These affect speed only; the registered solver is the same, and the tests
   pass unchanged.
3. **Launch cutoff:** 80 core-h on the ledger (was 40). The run cap stays at 24. Projection for the 22 remaining runs at
   Part 0 rates: about 31 core-h (about 37 with the 6.23 already on the ledger), about 6–7 h of wall time with 5 workers.
   The cutoff leaves room for some remaining contention.
4. **Disclosure** in the final report: ledger core-hours plus the 65 off-ledger core-hours of the lost runs.

PI sign-off (MG017): Richard Benear, 2026-11-01 20:58 UTC, decision card in the project thread ("Amend and relaunch").


**Result (MG018, recorded 2026-11-02; reviewed by ChatGPT (MG018) and signed off by the PI).** Relaunch on the PI's PC (RB_LibraryPC,
Windows 11, Python 3.13.2, numpy 2.2.4, 5 workers), amendment MG017 order. All 24 runs COMPLETE, no flag raised in the
status sense (evaluator ran without a mismatch). Ledger 39.96 core-h (cutoff 80, run cap 24; the cutoff did not bind),
plus about 65 off-ledger core-h from the lost first launch (MG017 item 4). Records: `results/MG5/runs/`, `checks/`,
`evaluation.json`, `driver_log.json`, `local_manifest.csv`, `budget_ledger.csv`.
- **Validity:** Part A frozen controls all single and pinned (max d ≤ 1.3e-13 < 0.1). Valid.
- **Gate MG5 (survival): PASS.** Single in 4/4 runs at ν = 0.002, 0.02 and 0.2.
- **Mobility:** ν = 0.002 PINNED (4 pinned; max d 0.35–1.37); ν = 0.02 PINNED (3 pinned, 1 mobile; max d 1.06–1.52);
  ν = 0.2 MOBILE (3 mobile, 1 pinned; max d 1.14–2.85). Threshold 1.5. D_p (descriptive): 1.2e-4–5.8e-4, 6.8e-4–1.0e-3,
  7.9e-4–2.3e-3 for the three rates. The expectation that mobility rises with ν is met in direction only; ν = 0.02
  has one run just above 1.5 (1.52).
- **Part B: A1.** Frozen controls multiplied in 4/4; ν = 0.2 single in 4/4. Read as "recovery under rewriting plus
  transfer"; MG5 does not separate changing disorder from the smoothing effect of repeated interpolation (§8).
- **Independent check:** `recompute_mg5.py` (separately written detection) agrees with the evaluator on all 24 classes
  and mobility statuses.
- **Not shown (per §8):** survival under token noise, other rules or geometry-changing rules; that the motion is
  Brownian; a general survival probability (4 runs per rate, 2-graph panel).
- **Core-h by group:** Part A frozen 2.54; ν = 0.002 3.64; ν = 0.02 8.19; ν = 0.2 15.08; Part B frozen 1.45; ν = 0.2 9.07.
- **Total compute:** 39.96 core-h on the ledger plus about 65 core-h from the lost first launch, about 105 core-h, before
  Part 0 (2.3 h wall) and setup and replay work outside the run timer.

**Review MG018 (ChatGPT, 2026-11-02, `review/MG018_chatgpt.md`; response `review/MG018_response.md`).** ChatGPT
independently reproduced, from the committed records, every class, track (≤ 1.8e-15), D_p, the gate (PASS) and Part B (A1).
It also reproduced the Part 0 comparisons, the jump accounting over 19 984 edited ticks, and all 16 event lists (25 856 362
events, hashes match). Recommendation: accept PASS and A1 within the registered scope. No blocking error. Three
nonblocking items, all accepted:
1. `recompute_mg5.py` now also checks tracks, D_p and every aggregate outcome, and exits 1 on any difference. Rerun on
   the MG5 records: everything agrees.
2. **Mobility labels are registered descriptive labels.** One ν = 0.02 run is mobile at 1.523 (0.023 above 1.5), and one
   ν = 0.2 run is at 1.598. Part 0's refinement differences do not bound every production displacement. The labels stand,
   with no change to the threshold. They do not establish a sharp transition or a monotonic mobility law.
3. **Survival is not invariance of shape.** In the fast Part A runs, depth relative to the matched frozen control ranges
   from about 0.78 to 1.26 across checks. Descriptive only.
Residual limits (review): no independent re-integration of the evolving-graph ODE, no reconstruction of smoothing from
full fields, and no check of the PC-only field files against the manifest.

**MG5 CLOSED 2026-11-02: Gate PASS; Part B A1.** PI sign-off: Richard Benear, 2026-11-02 13:39 UTC, in the project
thread ("I sign off on MG5"), after ChatGPT's review MG018 and the response.

### MG6 — G-coupled spatial growth around a particle (Phase 5, first campaign)
Status: **REGISTERED 2026-11-02 (MG019c)**, with PI approval. **No MG6 run of any kind preceded this commit.** The
registered design is `spec/MG6_proposal.md` **draft 3**, with its content frozen as committed at `106b3e2` (only its
status lines changed since). Drafts 1 and 2 are preserved as superseded (`spec/MG6_proposal_v1.md`, `_v2.md`). Reviews:
ChatGPT MG019 r1 and r2 (`review/MG019_chatgpt*.md`); responses `review/MG019_response*.md`.

Summary (the proposal governs):
- **Rule (rules amendment v1.3, MG6 only):** dominating candidate birth process at rate νρ₀F, F = e^0.4, with Exp(ν)
  lifetimes for initial nodes and candidates. Acceptance u < exp(κφ_G*)/F, using the boundary field on the pre-edit mesh.
  Rejected candidates are ghosts; ephemeral candidates are skipped. MG5's tick-synchronous transfer is kept. All arms of
  a start share one stream.
- **Arms:** b = 28, ν = 0.2, Δ = 0.2, dt = 0.05, t_end = 200, four MG4c starts. DEC (κ = 0); SH2 and SH4 (shifted
  field, |s| ≥ 6); CP2 and CP4. 20 runs.
- **Bound:** the maximum nodal φ_G is checked against 0.4/κ at every boundary. A breach is INVALID-BOUND and stops for
  the PI.
- **Gate MG6 (κ = 2, relative density response):** window t = 100–200. Run-level D₀(1.5) is compared in paired contrasts
  CP2 − SH2 and CP2 − DEC. **DEMONSTRATED** iff both contrasts have all four negative and a one-sided paired t ≤ −2.353,
  under the stated conditional working model. Otherwise **NOT DEMONSTRATED** (inconclusive), or **NOT SCORED
  (feedback)** if a run in CP2, SH2 or DEC is not single through the window.
- **Feedback (CP2, CP4):** SURVIVES (≥ 3/4 single) or CHANGED.
- **Reported separately:** the absolute sign of CP2's D₀(1.5); the κ = 4 contrasts; the D_c baseline; the quasistatic
  predictions; M₋, M₊ and M; the annular profiles.
- **Falloff:** expected UNRESOLVED in this campaign, stated before any run.
- **Part 0** (unregistered graph `seed_for(6, 9, 0, 0, 0)`):
  - acceptance-rate and determinism;
  - density equilibration;
  - step check with a fixed edit history;
  - fully coupled tick check;
  - bound margin;
  - transient;
  - quadrature;
  - timing with 5 concurrent runs (8 since MG020a).
  Each has registered formulas, fallbacks and stops. Outcomes, and any change of dt, Δ, φ_cap, grid or window, are
  written here before any registered run.
- **Budget:** run cap 24 and a provisional 40 core-h launch cutoff, restated after Part 0. Order: DEC, SH2, SH4, CP2,
  CP4. 5 workers (8 since MG020a) with one BLAS thread, launched with no console attached. Full accounting, including off-ledger work.
- **Deferred:** particle-mass scaling (b-scan), falloff, a curvature reading, the shuffled null, and ν = 0.02.

PI sign-off (design): Richard Benear, 2026-11-02 14:45 UTC, decision card in the project thread ("Register draft 3"),
after ChatGPT's two reviews (MG019 r1 and r2), which support implementation and Part 0.

**Pre-run implementation specification (MG020, 2026-11-02; committed before any Part 0 or registered MG6 run).**
1. **Code.**
   - `src/mg6.py`: the candidate process, shift, coupled state, quadrature grid and per-check measurements.
   - `src/quad6_check.py`: the separately written in-run quadrature cross-check.
   - `src/run_mg6.py`: the driver and tick loop `run_ticks6`.
   - `src/run_mg6_part0.py`: Part 0.
   - `src/evaluate_mg6.py`: the evaluator.
   - `src/recompute_mg6.py`: the independent offline recompute.
   - `tests/test_mg6.py`: 10 tests, all passing in the cloud.
   - Reused from MG5 unchanged: `evolve.Mesh`, `Ros2`, `apply_tick`, the interpolation audit, the jump accounting,
     detection, sparse records and `classify_run`.
2. **Streams and replay.**
   - Candidates: `gen_candidates(N0 = 18 901, ν = 0.2, t_end, Lx, Ly, seed)`. The draw order is fixed: life0 (N0), the
     count n ~ Poisson(ν·N0·F·t_end), tb (sorted uniforms on [0, t_end]), x, y, u, ell. Candidate j has id N0 + j.
   - SHA-256 is taken over N0 (`<i8`) followed by the six arrays (`<f8`) in that order. It is regenerated in a fresh
     subprocess before each run, and a mismatch makes the run INCOMPLETE.
   - Seeds: run candidates `seed_for(6, 0, 0, 1, 2·slot + pl)`; shift `seed_for(6, 0, 1, 1, 2·slot + pl)`.
   - Part 0 seeds: graph `seed_for(6, 9, 0, 0, 0)`; check 1 `(6, 9, 3, 1, 0)`; check 2 `(6, 9, 3, 1, 2)`; checks 3–4
     `(6, 9, 3, 1, 1)`; timing `(6, 9, 3, 1, 10 + i)` for i = 0–4.
   - Candidate lists are kept on the PC only, under the manifest.
3. **Tick edits.**
   - Tick k covers (t0, t1] with t1 = k·Δ and t0 = t1 − Δ.
   - Deaths: alive nodes whose death time lies in (t0, t1].
   - Candidates: those with tb in (t0, t1] (searchsorted, side right). Ephemeral iff tb + ℓ ≤ t1; ephemeral candidates are
     skipped.
   - Acceptance: strict u < p, with p = exp(κ φ_G*)/F, and φ_G* barycentric on M_(k−1) at x (+ s, wrapped). At κ = 0,
     p = 1/F without a field lookup.
   - Accepted candidates get death time tb + ℓ.
   - The edits are passed to MG5's `apply_tick` as an event array: deaths first, then accepted births in candidate order.
   - There are no edits after the final boundary.
   - Imposed acceptances (`forced`) exist only for Part 0 check 3.
4. **Bound.**
   - Before thinning, at every boundary: max nodal φ_G of the boundary field ≤ 0.4/κ (κ > 0). A failure gives
     INVALID-BOUND, and the flag records the time, the node and its position.
   - Candidate assertion: max p ≤ 1 + 1e-15. The per-tick max φ_G is recorded.
5. **Checks.**
   - Checks are taken at t = 0 and every 10 units, on M_(k−1) and the boundary field, before edits.
   - The measurements are centred on the deepest detected particle (MG5 detector, θ_28) when there is at least one.
     They are NaN when there is none.
   - The evaluator uses a run's window only if every window check has exactly one particle.
6. **Quadrature details.**
   - P_c uses the quadrature disc area Σw in place of πR² in its background term, so DEC is exactly zero (§4.3 (iii)).
     Σw differs from πR² by ≤ 0.1% (check (ii)).
   - Annular densities and predictions divide by the exact annulus area. Radial φ_G means divide by the quadrature area.
   - Cross-check (2× finer grid, point mask, its own triangulation): |a − b| ≤ 0.01·max(|a|, |b|) for M₋ and M₊, and
     the same plus 0.5 nodes for P₀(1.5) and P_c(1.5). A miss is a flag; a second miss on re-evaluation makes the run
     FAILED-NUMERICAL.
7. **Statuses.** FAILED-NUMERICAL > INVALID-BOUND > INCOMPLETE > COMPLETE. The flags, rechecks and stops are as in MG5
   (balance, jump accounting, interpolation audit), plus the bound and the quadrature cross-check.
8. **Data.**
   - Committed: `runs/<key>.json` and `checks/<key>.npz`. The npz holds:
     - sparse records, centres and particle counts;
     - every §4.2 measurement per check (`m_*`) and the cross-check values (`x_*`);
     - the per-tick table;
     - the packed accepted bitmap.
   - Also committed: the ledger, `local_manifest.csv`, `driver_log.json`, `evaluation.json` and `part0/`.
   - PC-only: `results/MG6/local/` (candidate lists, and fields every 50 units in float32, float64 at t_end).
9. **Budget and order.** The driver launches DEC, SH2, SH4, CP2, CP4 (asserted by a test). It stops launching at 24 runs,
   or when the ledger reaches the cutoff in `part0/decisions.json` (40 core-h unless restated). It runs 8 workers (MG020a; was 5) with one
   BLAS thread each, launched with no console attached.
10. **Evaluation.**
    - Y = the mean of D₀(1.5) over the window checks (both ends included).
    - Paired contrasts by start; t = mean/(sd/2); the conditional 95% interval uses t₀.₉₇₅,₃ = 3.182.
    - DEMONSTRATED needs both contrasts all negative and t ≤ −2.353. The absolute sign uses the same rule on CP2's Y.
    - observed/predicted = Y_D0/Y_P0 per run.
    - Falloff: the arm's mean coarse radial φ_G profile is interpolated linearly between annulus midpoints to find the
      first zero crossing. Eligible annuli have inner edge ≥ 1.0 and outer edge ≤ that crossing. The exponent is
      UNRESOLVED unless at least 3 eligible annuli have |t| ≥ 2.353; otherwise both fits use `curve_fit` with run-level
      SEs.
    - The evaluator stops if re-detection from the sparse records differs from the driver.
    - `recompute_mg6.py` repeats all of this independently, as stated in proposal §8, and exits 1 on any difference.
11. **Part 0 details.**
    - Checks 1–2 use the static field in fixed-field mode: no ODE, with φ_G re-evaluated at node positions after each
      edit.
    - The settle is MG5's frozen `run_ticks` (Δ = 1, dt = 0.05, t = 300) from the soliton at the torus centre.
    - Checks 3 and 4 share the (Δ = 0.2, dt = 0.05) run; the forced dt/2 replay follows it.
    - Fallbacks are as in proposal §7.
    - The timing runs 8 concurrent coupled runs (MG020a; was 5) for 20 units from the settled state; projected core-h = 20 × t_end ×
      mean s/unit / 3600.
    - MG5's ROS2-against-BDF solver check is not repeated: the solver is unchanged, and MG6 §7 does not list it.
12. **Disclosure.**
    - The tests ran dynamics only on small unregistered graphs (10 × 10 torus, 3 450 points).
    - One smoke test ran on a **registered start state** (D0_p0, b = 28) with an **unregistered** candidate seed (123):
      κ = 4 CP for 20 units, plus the start of a κ = 2 shifted run, which was stopped unread. These values were seen
      (t = 0, 10, 20):
      - one particle at every check;
      - D₀(1, 1.5, 2) = (6.6, 3.1, 8.3), (−51.4, −41.9, 8.3), (−52.4, −60.9, −52.7);
      - P₀(1, 1.5, 2) = (−56.6, −56.2, −19.3), (−58.9, −77.5, −51.5), (−66.8, −99.0, −85.5);
      - max φ_G 0.0475, max p 0.80, 440 s.
      The registered stream differs, nothing registered was changed in response, and the run is not part of any result.
    - In the cloud, Part 0 checks 1, 2 and the quadrature checks (i)–(iii) were exercised on the unregistered Part 0
      graph, to test the code. Their outcomes are recorded only by the PC's Part 0.
13. **Worker count (MG020a, PI decision 2026-11-02 15:36 UTC, before Part 0).** Richard chose 8 workers instead of 5
    (decision card "Use 8 workers for MG6 instead of 5?"). The PC has 16 logical CPUs and 32 GB; 9 workers thrashed memory
    in MG5, so the count stops at 8. This is a speed setting only: random streams come from `seed_for` and do not depend
    on the worker count, and each worker keeps one BLAS thread. Proposal §7 item 8 ("5 concurrent coupled runs") is
    therefore run with 8 concurrent runs, so that the restated budget is measured at the production worker count
    (`N_TIMING = 8` in `run_mg6_part0.py`, seeds `seed_for(6, 9, 3, 1, 10 + i)`, i = 0–7). Contention can raise per-run
    wall time and hence ledger core-h; the restated budget after Part 0 accounts for it. Claude watches free memory during
    the timing; if the PC runs short of memory, or the projected wall time is no better than with 5 workers, Claude
    reports it and the PI decides whether to drop to 6 workers before any registered run.

**Part 0 outcome (MG020b, 2026-11-02, PI's PC, code `410696d`, 8 workers; records in `results/MG6/part0/`, commit
`5d351dc`): STOP, step or tick check failed twice.** Checks 1–3 passed. The tick check failed part (a) both times
(Δ 0.2 vs 0.1: centre up to 0.54 at four checks; Δ 0.1 vs 0.05: centre 0.36 at one check; limit 0.15) and passed part (b)
both times (|ΔD̄₀| = 1.3 and 2.3 nodes against limits of 37). Checks 5–8 were not reached. This result is not re-scored.

**Amendment MG021 (`spec/MG6_amendment_1.md`; PI decision 2026-11-02 21:18 UTC; nothing run under it). Draft 1 text,
superseded by draft 2 below:** Δ = 0.1 and dt = 0.05 for all MG6 runs; part (b) governs the tick check for this campaign, with part (a)'s
failure recorded as a design flaw (single-snapshot centroid noise of about 0.1 against a 0.15 tolerance; diagnostic in
`review/MG020_centre_diagnostic.py`). Checks 5–8 are completed for Δ = 0.1 by `python src/run_mg6.py part0-mg021`, with
the registered criteria, into `results/MG6/part0/decisions_mg021.json`; `decisions.json` is kept. Provisional cutoff
80 core-h (was 40); run cap 24; the budget is restated here from check 8 before any registered run. Gate, arms, seeds,
measurements and evaluation are unchanged.

**Amendment MG021, draft 2 (after ChatGPT review MG021, `review/MG021_chatgpt.md`; response `review/MG021_response.md`;
supersedes the draft 1 text above).** Δ = 0.1, dt = 0.05. Numerical validation is accepted **by exception**: only the
instantaneous-centre clause of tick part (a) for Δ 0.1 vs 0.05 is waived; equal counts, depth ≤ 3% and part (b) remain
required (and are met). This is not a passed tick check; motion and D_p stay descriptive. The transient (§7.6) uses the
Δ 0.1 series with P̄₀ from the Δ 0.05 run. `part0-mg021` requires a hash-verified final state, a COMPLETE settle and
COMPLETE timing runs on the selected grid, the bound margin over timing runs and terminal boundaries, and a projection
≤ 80 core-h; otherwise it writes `proceed = false`. Cutoff 80 core-h, run cap 24.

**Part 0 completion and restated budget (MG021c/d, 2026-11-02/03).** `python src/run_mg6.py part0-mg021` ran on the PI's
PC (code `39f3a88`, 8 workers, 15:28–16:08 local; records `results/MG6/part0/decisions_mg021.json`, timing check records,
commit `0248e78`). Waiver condition met; final state hash-verified. Check 5: max φ_G 0.0534 (records), terminal maxima
0.0438 and 0.0456, timing runs 0.0472–0.0490, all ≤ 0.07. Check 6: B₁ = −147.44, B₂ = −136.44, no drift; window
t = 100–200, t_end = 200. Check 7: ok, grid 398 × 402. Settle COMPLETE (767 s). Check 8: 8 concurrent runs, all COMPLETE
at t = 20, 81.3–82.1 s per time unit; **projected 90.73 core-h**, above the 80 core-h cutoff, so the record says
`proceed = false` (PI decision required).
**PI decision (2026-11-03 00:14 UTC, decision card): cutoff raised to 110 core-h.** Recorded in
`results/MG6/part0/budget_decision_mg021.json`; the driver lifts only the budget stop, and only when the projection is
within that cutoff. Run cap 24, 8 workers, order DEC, SH2, SH4, CP2, CP4; about 13 h of wall time.

**MG6 result (MG022, 2026-11-03; corrected after ChatGPT review MG022, see below).** 20/20 runs COMPLETE on the PI's PC (code `c795770`, 8 workers,
Windows 11, Python 3.13.2, numpy 2.2.4), 83.09 **production ledger** core-h, no INVALID-BOUND or FAILED; records commit `09b72c2`.
`evaluate_mg6.py` → `results/MG6/evaluation.json`; `recompute_mg6.py` agrees on the PC (full check, candidate hashes
included) and in the cloud (`--cross-platform`; see below).
- **Gate (κ = 2): DEMONSTRATED.** Paired Y = D̄₀(1.5) contrasts by start (D0p0, D0p1, D1p0, D1p1):
  CP2 − SH2 = (−67.6, −66.5, −71.5, −62.5), mean −67.1, t = −36.2, 95% CI (−73.0, −61.2);
  CP2 − DEC = (−67.6, −66.5, −71.5, −62.0), mean −66.9, t = −34.1, 95% CI (−73.2, −60.7). All negative, t ≤ −2.353.
  Intervals are conditional on the registered four-start model.
- **Absolute deficit (CP2):** Y = (−66.9, −62.3, −64.1, −54.9), mean −62.05, t = −24.3: negative. Observed/predicted (Y_D0/Y_P0) =
  1.01, 0.94, 1.01, 0.83.
- **Feedback:** CP2 and CP4 SURVIVE (4/4 single each). Mobility (descriptive): CP2_D1_p0 and CP4_D1_p1 mobile, the rest
  pinned; all 12 null-arm runs single and pinned.
- **Descriptive:** CP4 − SH4 = −150.5 (t = −44.1); Dc contrasts agree with D₀ (Dc CP4 − DEC = −149.5, t = −43.7, added
  after review MG022); SH2 − DEC = +0.14 and SH4 − DEC = +0.05
  (the shifted nulls are near-identical to DEC around the particle, as designed).
- **Falloff:** UNRESOLVED at both κ, as pre-declared. CP2: first zero crossing of the mean φ_G profile at r = 1.65; CP4
  (added after review MG022, as §4.4 requires): r = 1.95. In each arm only the annulus [1, 1.5] is eligible (paired
  t = −8.31 for CP2, −23.84 for CP4); three are required.
- **Cost:** production ledger 83.09 core-h (20 runs, cutoff 110). Off-ledger, not on the ledger and only partly measured:
  the original Part 0 (11 662 s wall, at most 3 concurrent runs: ≤ 9.7 core-h); the MG021 completion (2 418 s wall:
  timing about 3.6 core-h, settle 0.2, quadrature small; about 4 core-h); evaluations and replays on the PC and in the
  cloud, and the cloud Part 0 exercises, unmeasured.
- **Reading:** the G-coupled birth rule produces a node deficit around the particle that matches the static prediction
  (P₀) within 0–17%, and the particle survives the feedback. This shows the coupled process does what it was designed to
  do without destroying the particle; it does not test a mechanism beyond that rule.
- **Cross-platform note:** in the cloud (Linux, numpy 2.5.3 and 2.2.4) the candidate lists regenerate with identical
  N0, life0, tb, x, y, u but a different `ell` hash (PC per-field check, CP2_D0_p0), because numpy's exponential sampler
  calls the C library's log on a small fraction of draws. `recompute_mg6.py --cross-platform` reports this and audits the
  rest: rebuilt node counts, detections, classes, window means and the gate all agree for all 20 runs. On the PC the full
  check, hashes included, passes.

**Review MG022 (ChatGPT, `review/MG022_chatgpt.md`; response `review/MG022_response.md`).** Strict replay on Windows
(numpy 2.2.4): all 20 candidate hashes match; node sets, counts and detections rebuilt at 420 checks; primary statistics,
absolute deficit, feedback and mobility labels reproduced. Recommends accepting DEMONSTRATED and SURVIVES (both κ), no
rerun. Three record corrections, all made: the κ = 4 falloff and Dc CP4 − DEC added to the evaluator (both falloffs
UNRESOLVED); `recompute_mg6.py` now compares the stored contrast vectors, t statistics, intervals and the absolute
result; rounding and cost labels corrected above. Gate outcomes unchanged.

**MG6 CLOSED 2026-11-03: gate DEMONSTRATED, feedback SURVIVES (κ = 2 and 4), falloff UNRESOLVED (both κ).** PI sign-off:
Richard Benear, 2026-11-03 13:24 UTC, in the project thread ("I sign off on MG6"), after ChatGPT review MG022 and the
corrections in its response. Scope as reviewed: the deficit is the local node count under the imposed birth rule; MG6
establishes neither intrinsic curvature nor a long-range force.

### MG7 — Mass scaling of the coupled density deficit (b-scan; Phase 5, second campaign)
Status: **REGISTERED 2026-11-03 (MG023b)**, with PI approval. **No MG7 run of any kind preceded this commit.** The
registered design is `spec/MG7_proposal.md` **draft 2**, with its content frozen as committed at `5b44cde` (only its
status lines changed since). Draft 1 is preserved as superseded (`spec/MG7_proposal_v1.md`). Review: ChatGPT MG023
(`review/MG023_chatgpt.md`, static replication `review/MG023_replication/`); response `review/MG023_response.md`. The
multiway/foliation gate formerly listed as MG7 is renumbered MG8.

Summary (the proposal governs):
- **Everything not stated is MG6 as registered and amended (MG021),** code pinned at `c795770`: coupled candidate rule,
  F = e^0.4, full-field bound check, Δ = 0.1, dt = 0.05, ν = 0.2, t_end = 200, checks every 10 units, MG6 §4
  measurements and quadrature.
- **b = 27.5 and 28.5** (new runs); MG6's recorded b = 28 CP2 and SH2 runs (`09b72c2`) are the midpoint, reported only.
- **Start fields:** MG4c's own final array at the requested b, no shift, SHA-256-checked against a committed manifest;
  the ODE's b asserted; b in every run key. Detector threshold fixed at θ = −4.45758 (MG4c b-specific thresholds as a
  descriptive check).
- **Arms:** CP2 and SH2 per endpoint, four starts: 16 runs. One candidate stream `seed_for(7, 0, 0, 1, 2·slot + pl)` and
  one shift `seed_for(7, 0, 1, 1, 2·slot + pl)` per start, shared across both b and both arms.
- **Gate MG7 (directional endpoint result):** c_i(b) = Y_i(CP2, b) − Y_i(SH2, b), Y = window mean (t = 100–200) of
  D₀(1.5); Δ_i = c_i(28.5) − c_i(27.5). **LARGER AT HIGHER b** iff all four Δ_i < 0 and the one-sided paired t ≤ −2.353
  (df = 3). Otherwise **NOT DEMONSTRATED**; **NOT SCORED (feedback)** if any endpoint CP2 or SH2 run is not single at
  every window check; **INCOMPLETE** (precedes all others; PI decides) if any of the 16 runs is not COMPLETE.
- **Reported, not scored:** the three contrasts c(27.5), c(28), c(28.5) with MG6's replication rule; within-start
  trajectories; paired residuals Δ − ΔP (ΔP = change in predicted CP2 minus predicted SH2) with start-field and in-run
  predictions, no agreement label; Fieller ratio; descriptive slope against M₋; feedback per endpoint; mobility; falloff
  (UNRESOLVED expected, declared now); early/late window means (t = 100–140, 160–200).
- **Part 0** (unregistered MG6 Part 0 graph settled at each endpoint; streams `seed_for(7, 9, 3, 1, ib)`): endpoint input
  hashes, single detection and bound (STOP on failure); forced-history step check dt 0.05 vs 0.025 (equal counts, depth
  ≤ 2%, |ΔD̄₀| ≤ 2 nodes); free-acceptance tick check Δ 0.1 vs 0.05 (equal counts, depth ≤ 3%, |ΔD̄₀| ≤ 6 nodes; centre
  recorded, not scored); bound margin ≤ 0.14 at every boundary; every Part 0 run COMPLETE. Outcomes go to
  `results/MG7/part0/decisions.json` and here before any registered run. Validation budget about 25 core-h, cap 30.
- **Budget:** about 74–75 production ledger core-h; cutoff 90, run cap 22, 8 workers, order SH2 then CP2. Spare slots only
  for reruns of runs INCOMPLETE for a non-physical reason; a physical outcome is never replaced. A projection above 90
  needs a recorded PI decision.

PI sign-off (design): Richard Benear, 2026-11-03 13:55 UTC, in the project thread ("register you plan as is"), after
ChatGPT review MG023 and its adoption in draft 2; no second review round.

**Pre-run implementation specification (MG024, 2026-11-03; committed before any Part 0 or registered MG7 run).**
1. **Code.** `src/run_mg7.py` (driver, start-field loader), `src/run_mg7_part0.py` (Part 0), `src/evaluate_mg7.py`,
   `src/recompute_mg7.py`, `tests/test_mg7.py` (8 tests). MG6's engine is reused unchanged: `run_mg6.run_ticks6`,
   `save_checks6`, `summarize6`, `mg6.Coupled`, `mg6.Grid`, `mg6.measure`. A test checks that the reused files
   (`mg6.py`, `run_mg6.py`, `run_mg5.py`, `evolve.py`, `quad6_check.py`, `graph_ode.py`, `lattice.py`, `graphs.py`,
   `runner.py`, `run_mg4c.py`, `run_mg6_part0.py`) are byte-identical to the MG6 production revision `c795770`.
2. **Start fields.** `results/MG7/endpoint_manifest.csv` (committed now) lists the SHA-256 of the eight MG4c files
   `results/MG4c/fields/b{27.5, 28.5}_D{0,1}_p{0,1}.npz`. `start_state7(slot, pl, b)` refuses any b other than 27.5 or
   28.5, checks the file's SHA-256 against the manifest, loads `final` (float64, 3 × 18901) with no shift, checks it equals
   the file, and takes the graph and h from MG4c (`run_mg4c.graph_for`, independent of b). The reaction parameter is
   built by `params_for(b)`, which asserts b. Each run record stores b, the parameter b, the file path, its SHA-256 and
   `shift_applied = false`; the evaluator and the recompute stop if any of these disagree with the spec.
3. **Runs.** Keys `{SH2,CP2}_b{27.50,28.50}_D{slot}_p{pl}`; launch order SH2 (27.5, then 28.5), then CP2. κ = 2, Δ = 0.1,
   dt = 0.05, t_end = 200, checks every 10, MG6 grid 398 × 402, θ = −4.457581675528346. Candidate stream
   `seed_for(7, 0, 0, 1, 2·slot + pl)` and shift `seed_for(7, 0, 1, 1, 2·slot + pl)`, shared across b and arms.
4. **Launch guards.** The driver refuses to launch unless the eight endpoint files match the manifest and
   `results/MG7/part0/decisions.json` has `proceed = true` with Δ 0.1, dt 0.05. If Part 0 stops only on the projected
   cost, a PI decision recorded in `results/MG7/part0/budget_decision.json` (`pi_approved`, `cutoff_core_h`) lifts that
   stop and nothing else. Cutoff 90 production ledger core-h, cap 22 (every started run counts), 8 workers, one BLAS thread.
5. **Spare slots.** The driver never reruns a run that has a record. A non-physical INCOMPLETE (crash, interruption) is
   rerun only after the PI moves its record aside, with the reason recorded; a physical outcome is never replaced.
6. **No silent overwrites (MG021 lesson).** The driver moves any existing output it would replace to
   `results/MG7/superseded/` with its SHA-256 in the name and a line in `moves.csv`. Part 0 refuses to start, or to write,
   if an output already exists.
7. **Part 0.** As proposal §7, in this order:
   - check 1 on the registered start files (static: hashes, array equality, single at θ, max nodal φ_G ≤ 0.14,
     P₀^start(1.5) and M₋ recorded with κ = 2, no shift);
   - settle the MG6 Part 0 graph (`seed_for(6, 9, 0, 0, 0)`) frozen at each b (Δ 1, dt 0.05, t = 300), which must be
     COMPLETE and single at t = 300;
   - for each b, with stream `seed_for(7, 9, 3, 1, ib)` (ib 1 = 27.5, 2 = 28.5), CP2 runs of 100 units: Δ 0.1/dt 0.05
     (A), Δ 0.05/dt 0.05 (C), then Δ 0.1/dt 0.025 with A's accepted edits imposed (B);
   - step check A vs B: equal counts at every check, relative depth ≤ 2%, |D̄₀ difference| ≤ 2 nodes (D̄₀ = mean
     D₀(1.5) over t = 50–100); tick check A vs C: equal counts, depth ≤ 3%, |D̄₀ difference| ≤ 6 nodes; centre
     differences recorded, not scored;
   - bound: max nodal φ_G ≤ 0.14 at every boundary of the six coupled runs;
   - every Part 0 run COMPLETE;
   - Part 0 cap 30 core-h: after A and C, if the spent wall time plus 2 × 2 × mean(A wall) exceeds 30, STOP before B;
   - projection = 16 × 200 × max(MG6 production mean s/unit over CP2 and SH2 = 83.5, mean s/unit of the two A runs) /
     3600. `proceed = true` only if every check passes and the projection ≤ 90; otherwise STOP with the reason
     (a check failure: PI decides; a projection above 90: "PI decision required").
8. **Evaluation.** `evaluate_mg7.py` re-detects every check from the sparse records (stops on any count or class
   difference), reads MG6's b = 28 CP2 and SH2 records as the midpoint, and applies §5: INCOMPLETE first, then NOT
   SCORED (feedback), then the gate on Δ_i = c_i(28.5) − c_i(27.5) (all four negative and one-sided t ≤ −2.353). Reported:
   the three contrasts with MG6's replication rule, per-start trajectories, residuals e = Δ − ΔP with ΔP from the
   start-field (t = 0 check) and in-run (window-mean) predictions of CP2 minus SH2, each with mean and conditional 95%
   interval, a paired Fieller interval for mean Δ / mean ΔP (unbounded if the denominator's leading coefficient is ≤ 0),
   the slope of c against window-mean M₋ over 12 CP2 runs (descriptive), feedback and mobility per endpoint, falloff
   per endpoint (UNRESOLVED expected), single-in-window at the MG4c b-specific thresholds, and early (t = 100–140) and
   late (160–200) window means with |P₀^start| as reference. `recompute_mg7.py` checks all scored quantities
   independently (strict on the records' platform, `--cross-platform` elsewhere).
9. **Check 1 computed in the cloud while writing the code (no dynamics; disclosed):** all eight files match, every start
   is single at θ, max nodal φ_G 0.0564; mean P₀^start(1.5) −35.00 (27.5) and −44.78 (28.5), as in ChatGPT's static
   replication. The PC run of check 1 is the record.
10. **Not committed:** `results/MG7/local/`, `results/MG7/superseded/`, `results/MG7/part0/settled_*.npz` (SHA-256 in
    the records).

**MG7 Part 0 outcome (MG024a, 2026-11-03, PI's PC, code `ca2c513`, 8 workers; records in `results/MG7/part0/`,
commit `9b20219`): PROCEED.** Launched 07:29 local, finished 10:42 (11 549 s wall, 8.63 core-h against the 30 cap).
Re-checked in the cloud from the committed check records (window means recomputed from `m_D0`).

| Check | Result |
|---|---|
| 1. Endpoint inputs | ok: 8 files match the manifest; all starts single; max nodal φ_G 0.0564; mean P₀^start(1.5) −35.00 (27.5), −44.78 (28.5) |
| Settles (t = 300) | COMPLETE and single at both b (SHA-256 of the PC-only settled states in `decisions.json`) |
| 2. Step check, dt 0.05 vs 0.025 (forced) | ok: \|ΔD̄₀\| 0.17 (27.5) and 0.17 (28.5), limit 2; depth ≤ 1.06%; counts equal (1 at every check); centre 0.005 / 0.0003 (recorded) |
| 3. Tick check, Δ 0.1 vs 0.05 | ok: \|ΔD̄₀\| 0.5 (27.5) and 1.0 (28.5), limit 6; depth ≤ 1.46%; counts equal; centre 0.076 / 0.101 (recorded, not scored) |
| 4. Bound margin | ok: max nodal φ_G 0.0621 ≤ 0.14 over the six coupled runs |
| 5. Completion and timing | all eight Part 0 runs COMPLETE; projection 16 × 200 × 83.48 s / 3600 = **74.2 core-h** (MG6 production rate; Part 0's 40.3 s/unit at 4 concurrent runs is lower and not used) |

Production configuration as registered: Δ 0.1, dt 0.05, t_end 200, window t = 100–200, grid 398 × 402, cutoff 90
core-h, cap 22, 8 workers. The Part 0 deficits (D̄₀ ≈ −45 at 27.5 and −73 at 28.5 on the unregistered graph, CP2 alone,
no null) are validation quantities, not MG7 data, and are not interpreted.

**MG7 result (MG025, 2026-11-04; registered runs on the PI's PC 2026-11-03 12:44–21:53 local, code `25f1007`,
8 workers; records `282726f`).** 16/16 runs COMPLETE, 73.01 production ledger core-h (projected 74.2; cutoff 90).
`evaluate_mg7.py` (PC and cloud agree) → `results/MG7/evaluation.json`.

- **Gate MG7: NOT SCORED (feedback).** The registered rule applies because endpoint runs are not single through the
  window: at b = 27.5 the particle dissolved in 7 of 8 runs. No run multiplied; no count returned from 0.

  | b = 27.5 | D0_p0 | D0_p1 | D1_p0 | D1_p1 |
  |---|---|---|---|---|
  | SH2: first check with no particle | t = 160 | 90 | 60 | 110 |
  | CP2: first check with no particle | t = 20 | single throughout (mobile) | 20 | 140 |

- **Feedback:** b = 27.5 **CHANGED** (CP2 1/4 single); b = 28.5 **SURVIVES** (CP2 4/4 single, all pinned). All 8
  b = 28.5 runs are single and pinned.
- **Interpretation (corrected after ChatGPT review MG025):** b = 27.5 failed this finite-panel survival criterion
  through t = 200 in both tested arms; one CP2 run survived there. No survival edge was searched for and monotonic
  survival in b was not tested, so no new edge is estimated. The SH2 losses show that aligning the birth bias with the
  particle is not needed for the loss; SH2 still has field-dependent births elsewhere and MG7 has no DEC arm, so uniform
  rewriting is not isolated as the cause. The result is compatible with a narrower robust-survival region under
  rewriting than MG4c's frozen-graph window (27.27–28.73). MG7 Part 0's unregistered 27.5 run stayed single for 100
  units. The mass-scaling question remains open.
- **Threshold sensitivity (registered, stored by the evaluator per run as `alt_threshold`):** at MG4c's b-specific
  thresholds every run has the same single-in-window verdict. In the seven dissolved runs the final smoothed-X minimum
  is at most 1.34e-8 in magnitude (ChatGPT MG025), so the loss does not depend on the threshold.
- **Unavailable (undefined inputs):** the b = 27.5 contrast, the endpoint change Δ and its t, both residuals, the
  Fieller ratio, the twelve-run slope against M₋, and falloff at b = 27.5. Nothing is substituted for them.
- **Descriptive at b = 28.5** (the registered §5 items whose inputs exist; collected by `src/describe_mg7.py`, code
  added after the data → `results/MG7/descriptive.json`): CP2 − SH2 = (−75.1, −71.0, −74.8, −67.2), mean −72.0,
  t = −38.6 (all four negative). Within-window late − early (t = 160–200 minus 100–140) D₀ means: CP2 +19.4, −6.0, +13.8,
  +16.4; SH2 +19.2, −15.0, −20.6, +5.4; contrast +0.2, +9.0, +34.4, +11.0 nodes. These are temporal changes over the
  fixed window, not a gate or a drift test, and they caution against reading a few-node difference as an equilibrium
  effect. MG6 at b = 28.0: mean −67.1; per-start 28.5 − 28.0 differences (−7.5, −4.5, −3.3, −4.6) use different
  candidate streams and are not the registered test. Window-mean M₋ (CP2): 1.15 at 28.5 against 0.91 at 28.0.
- **Corrections (MG025 response):** the earlier text here said the lower survival edge "lies above 27.5" and that
  threshold sensitivity was "not computed"; both were wrong and are replaced above.
- **Deviation (audit code, after the data):** `recompute_mg7.py` crashed on checks with no particle (no centred
  measurement is stored there, so `m_N` is NaN). Fixed to require n_parts = 0 at such checks and skip the count
  comparison; strict mode now reports a candidate-hash failure without crashing. No scoring rule changed. In the cloud
  (`--cross-platform`; all 16 candidate hashes differ on Linux, as in MG6) everything agrees. Strict rerun on the PI's PC (2026-11-04, main `04ad4e2`): candidate hashes match and
  everything agrees (exit 0); `evaluation.json` unchanged by a re-evaluation.

**MG7 CLOSED 2026-11-04: NOT SCORED (feedback); feedback CHANGED at b = 27.5, SURVIVES at b = 28.5.** PI sign-off:
Richard Benear, 2026-11-04 14:38 UTC, in the project thread ("i sign off on MG7"), after ChatGPT review MG025 and the
interpretation corrections in its response. Scope: b = 27.5 failed the finite-panel survival criterion under rewriting
in both arms; no survival edge is estimated; the mass-scaling question remains open.

### MG9 — Local onset type and front arrest against the particle window (deterministic; Phase 2 follow-up)
Status: **REGISTERED 2026-11-04 (MG026c)**, PI decision card 16:40 UTC ("Register draft 3"). **No MG9 run of any kind
preceded this commit.** The registered design is `spec/MG9_proposal.md` **draft 3**, with its content frozen as
committed at `7810b8b` (only its status lines changed since). Drafts 1 and 2 are preserved as superseded
(`spec/MG9_proposal_v1.md`, `_v2.md`). Reviews: ChatGPT MG026 and MG026 r2 (`review/MG026_chatgpt*.md`, static and
synthetic checks in `review/MG026_replication/`); responses `review/MG026_response*.md`. Disclosed prior computation:
one timing run of the b = 29 MG003 reference; ChatGPT's static scripts re-run; one synthetic envelope-kernel check.

Summary (the proposal governs):
- **P1 (local onset type):** cubic coefficient L₃ of the real cosine amplitude at the D02 mode-34 threshold (frozen
  formula, §5), validated by the V2 fitted-growth runs. Outcomes in order: NOT EVALUATED, DEGENERATE, MIXED,
  SUBCRITICAL / SUPERCRITICAL. Either sign is a valid finding.
- **P2 (primary):** the front-arrest interval on D02 (grid b = 28.50–30.00, step 0.05; ordered classification §6.3;
  refinement, robustness and preparation statuses §§6.1, 6.4) against the rechecked window benchmarks (§6.5). Outcomes
  in order: UNRESOLVED, NO ARREST OBSERVED, OUTSIDE, COINCIDES (tolerance 0.10, both edges resolved), OVERLAPS; flags
  CENSORED, ANOMALOUS, BENCHMARK-SENSITIVE, MESH-SENSITIVE. A positive result reads "front arrest consistent with a
  pinning interpretation", never "snaking".
- **Descriptive:** α, L₃ terms, sweep persistence bracket, E_pat(b), interface velocities, one Landau-style figure.
- **Budget:** run cap 230 integrations (210 planned + 20 reserve), cutoff 10 core-h, cloud container; stop priority as
  §9.
- **Next:** code and tests (`src/landau_mg9.py`, `src/run_mg9.py`, `src/evaluate_mg9.py`, `tests/test_mg9.py`), with the
  synthetic envelope and classification validations passing before any Model G run, and a pre-run implementation
  specification appended here.


**MG9 amendment 1 (MG027, 2026-11-04; ADOPTED: PI decision card 16:59 UTC, "Approve amendment 1"):** `spec/MG9_amendment_1.md`. The draft 3 §6.2
envelope failed its registered synthetic validation (worst edge error 0.27 λ against λ/10); it is replaced by complex
demodulation at k = 2π/λ with a double λ-boxcar low-pass (worst 0.036 λ). Synthetic only; no Model G MG9 run.

**MG9 pre-run implementation specification (MG027, 2026-11-04; written before any MG9 run).** Code: `src/landau_mg9.py`,
`src/run_mg9.py`, `src/evaluate_mg9.py`; tests `tests/test_mg9.py` (16 pass, cloud: Python 3.13.18, numpy 2.5.3,
scipy 1.18.1). The production code revision is the commit that merges this specification. The proposal (as amended)
governs; these items fix what it leaves open.
1. **Platform.** All stages run in the cloud container, one process, stages resumable at run level. The platform and
   package versions are recorded in `results/MG9/evaluation.json`.
2. **Ledger.** One count per integration call. A sweep extension (to 3000) or front extension (to 1800) continues the same
   count. Wall time per integration is summed as core-s. The cap (230) and cutoff (10 core-h) are checked before each new
   integration, and the driver exits with code 2 on a stop.
3. **Completeness.** `mg9_ode` returns `complete` only if `sol.success` and the time reached equals t_end (relative
   1e-9). V2 integrates in 100-unit chunks; any incomplete chunk marks the run incomplete, and so V2 fails.
4. **V2.** A(t) is sampled every 1 unit; the stop check applies at each sample. The fit is linear least squares on
   y = 1/A² of the integrated form (proposal §5), over t ≥ 20 up to the stop. Usable means stop_t ≥ 40. σ₃₄ is the
   largest real part of the eigenvalues of M(b, κ₃₄) at the problem-(i) threshold (about 1e-15 by construction). The
   identifiability ratio uses |σ₃₄|. The tolerance repeat compares the middle amplitude's base and tight fits.
5. **Sweep start.** A₀ = 0.5 times cos(k₃₄x) times the X-normalized leading eigenvector of M(30.50, κ₃₄). No slaved
   terms are added.
6. **Preparation classification.** PREP-UNRESOLVED if incomplete or the residual is above 1e-6 after the extension;
   PREP-LOST if max |φ_X| < 0.1; PERSISTS if the cell count equals round(N h / λ) (34, or 68 on D02L); otherwise
   PREP-CHANGED. Cell count and E_pat use the amendment 1 envelope.
7. **Off-sweep preparation source.** The sweep value above b is the smallest sweep value greater than b. If that sweep
   state is not PERSISTS, its status is propagated to b and nothing is settled. A stage that needs a sweep state not yet
   run raises an error and records nothing.
8. **Tracking.** At t = 0, the tracked patch is the component containing node N/2. Afterwards it is the component of
   largest node overlap with the previous one. After a fill sample the previous set is the whole ring. If the tracked
   patch has no overlap with any component while components exist, the sample is "untracked", and the run is
   UNRESOLVED-INTERNAL. BREAKUP means at least two components each overlap the previous patch over at least λ. REMOTE
   means a component has zero overlap. Interfaces are the linear-interpolation crossings of E = ½ E_pat between each
   end node of the tracked component and its outside neighbour, unwrapped by the nearest-image rule. Outward
   displacement is −(x_L − x_L(t₁)) and x_R − x_R(t₁), with t₁ the first sample at t ≥ 100.
9. **Velocities.** Least-squares slope of the outward displacement over the samples with t ≥ T − 200.
10. **Phase convention.** The cell minimum nearest node N/2 uses the same prominence and separation rule as the cell
    count. The phase-shift repeat is an exact Fourier translation of all three fields by λ/2. The D01 interpolation is
    zero-padding Fourier interpolation with the Nyquist term split.
11. **Refinement order.** The lower edge first, then the upper, each up to two bisections. A midpoint already run is not
    re-run. Robustness repeats then run at the final inner points. On D02L the tiled state is settled to t = 300 (one
    count) and the patch half-width is 50.
12. **ANOMALOUS.** Set if any point below the arrest interval is INVADES or another non-U class other than RETREATS or N,
    or any point above is anything but INVADES or U.
13. **Benchmarks.** b_s "exists" iff φ_X at node N/2 is below −1 at t = 600. b_r "single" iff exactly one interior
    local minimum is below −4 at t = 800. MG003's bisection treated zero minima as single; the recheck follows the
    registered "exactly one", and its verified endpoints make this explicit. Bisection midpoints are rounded to 1e-6.
    The b_r benchmark uses zero initial potentials and the standard seed, as in MG003.
14. **MESH-SENSITIVE.** Each D01 side's edge is the adjacent pair (M or N outside, A inside) among its points. Edge
    positions are bracket midpoints, on both meshes. h = 0.2 uses the rechecked benchmark midpoints; h = 0.1 uses
    MG003's (28.97821, 29.4119), not rechecked.
15. **P2 inputs.** P2 is UNRESOLVED if any grid point lacks a record, or `benchmarks.json` is absent. The order of
    §6.6 is applied in `evaluate_mg9.p2_outcome`, tested on synthetic class sequences.

**MG9 RESULT (MG028, 2026-11-04; cloud, Python 3.13.18, numpy 2.5.3, scipy 1.18.1; code `f7ef9e9` + docs `54b2a08`):
P1 SUBCRITICAL; P2 UNRESOLVED.** 150 integrations, 0.49 core-h (cap 230, cutoff 10). Evaluator
`results/MG9/evaluation.json`; records `results/MG9/runs/`, states `results/MG9/states.npz`.
- **P1 = SUBCRITICAL.** L₃ (real cosine amplitude, original time): (i) D02 mode 34 at b = 30.04577: **+0.013653**
  (terms +0.01563, +0.01311, −0.01508; no near-degeneracy: ordinary partial cancellation, far from the DEGENERATE threshold [wording corrected per ChatGPT MG028]); (ii) +0.014162; (iii) +0.013128; (iv) +0.012602. Code check
  passed. V2 passed: fitted L₃ = 0.013650 (A₀ = 0.1) and 0.013641 (A₀ = 0.2), within 0.1% of (i); A₀ = 0.05 was not
  identifiable (1/A² changed 20%, below 25%) although its fit was 0.013653; tight-tolerance repeat 0.0136500.
  Independent recompute (finite-difference Taylor coefficients of the raw reaction terms): 0.013653.
- **Sweep (descriptive).** The 34-cell pattern PERSISTS from 30.50 down to 29.22. At 29.20 it was still changing at
  t = 3000 (PREP-UNRESOLVED, 30 cells, residual 1.5e-5); from 29.18 to 28.62 a stationary 30-cell pattern (PREP-CHANGED);
  PREP-LOST from 28.60. Primary persistence bracket [29.20, 29.22]; no 34-cell state regained.
- **Benchmarks (rechecked, longer horizons):** b_s ∈ [28.928907, 28.929688], b_r ∈ [29.395313, 29.396094]; shifts
  +0.0006 and −0.0010 from MG003; not BENCHMARK-SENSITIVE.
- **A3 grid (D02).** 28.50–29.15: N (UNAVAILABLE-LOST at 28.50–28.60, UNAVAILABLE-CHANGED at 28.65–29.15); 29.20: U
  (preparation unresolved); 29.25, 29.30, 29.35: U (UNRESOLVED-SLOW after extension to T = 1800: both interfaces moved
  outward 2.22–2.24 over [100, 1800], between λ/2 and λ, decelerating to about λ/8000 per unit time, residual 7.6e-3);
  29.40–30.00: INVADES (all 13). No ARRESTED point, so no refinement or robustness repeats were run, and the D01 comparison
  is NOT EVALUATED.
- **P2 = UNRESOLVED** (§6.6 rule 1/2: no A point and some U points). Independent recompute of every grid label from the
  raw tracks agrees.
- **Descriptive, not scored:** the change from slow creep to invasion lies between b = 29.35 and 29.40, which brackets
  the rechecked b_r (29.3953–29.3961). Below 29.22 the 34-cell pattern does not persist under this preparation, so the
  front test cannot reach the lower window edge (b_s ≈ 28.93) with the registered wavelength.
- **Deviation (after the data, record only):** `evaluate_mg9.py` now writes the platform into `evaluation.json` (pre-run
  spec item 1); no scoring change.

**MG9 CLOSED 2026-11-04: P1 SUBCRITICAL; P2 UNRESOLVED.** PI sign-off: Richard Benear, 2026-11-04 19:21 UTC, on the
project decision card ("Sign off and close MG9"), after ChatGPT result review MG028 and its response
(`review/MG028_response.md`). Scope: P1 is the local onset type at the Turing threshold only; P2 found no finite-time
arrest, and the lower window edge was not testable with the registered 34-cell wavelength. The stochastic barrier S(b)
remains deferred to a separately registered MG9b.

### MG9c — Late-time fate of MG9's slow fronts (deterministic; Phase 2 follow-up to MG9)
Status: **REGISTERED 2026-11-04 (MG029c)**, PI decision card 19:56 UTC ("Register draft 2" with the §9 defaults). **No
MG9c run of any kind preceded this commit.** The registered design is `spec/MG9c_proposal.md` **draft 2**, content frozen
at `372336b` (draft 1 kept as `spec/MG9c_proposal_v1.md`). Reviews: ChatGPT MG029 and MG029 r2
(`review/MG029_chatgpt*.md`; responses `review/MG029_response*.md`). Disclosed prior knowledge: exponential fits of the
committed MG9 tracks (proposal §0). MG9's own classifications are unchanged.
- **Q1 (primary):** at b = 29.25, 29.30, 29.35 (D02, patches prepared from the 34-cell branch, MG9 splice), late-time
  classification on windows from t = 1800 (T = 12,000, extension 30,000): SETTLED / ADVANCES / RECEDES / INVADES /
  RETREATS / UNRESOLVED-*. Robustness of a SETTLED point: D01, D02L (half-width 25), ±5% amplitude perturbations;
  PASS / FAIL / UNRESOLVED. Outcome order: UNRESOLVED, SLOW INVASION, BOUNDED, MIXED. Prediction (informed by §0):
  BOUNDED.
- **Q2 (if 29.35 is SETTLED-ROBUST):** three bisections between 29.35 and 29.40; interval comparison with b_r ∈
  [29.395313, 29.396094], tolerance 0.01: NEAR b_r / BELOW / ABOVE / UNRESOLVED RELATIVE TO THE TOLERANCE.
- **Descriptive:** b = 29.22; width trials (D02L half-width 50; D02 half-width 25 + λ/4).
- **Cap 60 integrations; cutoff 2 core-h; cloud.** Replay of MG9's [0, 1800] tracks to 1e-6 required (REPLAY-MISMATCH
  otherwise). Synthetic tests and a pre-run implementation specification come before any run.


**MG9c pre-run implementation specification (MG029d, 2026-11-04; written before any MG9c run).** Code:
`src/run_mg9c.py`, `src/evaluate_mg9c.py`; tests `tests/test_mg9c.py` (20 pass; no Model G integration, a fake solver
checks segmentation, snapshots and the ledger). MG9's `run_mg9.py` and `landau_mg9.py` are imported unchanged
(identical to production code `f7ef9e9`).
1. **Platform.** Before the first run the versions must equal MG9's recorded ones (Python 3.13.18, numpy 2.5.3,
   scipy 1.18.1); checked in this container. Every run record stores its platform.
2. **Prepared states** come from `results/MG9/states.npz`: the key named in each MG9 grid record's `prep` field
   (`prep_29.2500`, `sweep_29.30`, `prep_29.3500`), and `sweep_29.22`. Splice: `phase_convention`, then
   `splice_patch` (half-width 25). E_pat is the envelope mean of the prepared φ_X. It equals the MG9 record to 1e-12
   (tested).
3. **Segmentation (replay).** Main, D01 and D02L front runs integrate in separate `mg9_ode` calls over [0, 600],
   [600, 1800], [1800, 3600], [3600, 6000], [6000, 12,000] (extension [12,000, 15,000], [15,000, 30,000]), each
   starting from the previous final state, with t_eval every 5 units from t_a + 5. The first two segments are exactly
   MG9's `front_run` calls. Perturbation runs use [0, 1500], [1500, 3000], [3000, 6000]. A run counts once in the
   ledger, extension included.
4. **Replay check.** All 361 samples on [0, 1800] are compared with the MG9 record: identical statuses and
   |Δx_L|, |Δx_R| ≤ 1e-6. The E_pat difference is recorded. A failure is REPLAY-MISMATCH (Q1 UNRESOLVED). Tolerances are
   not changed after a mismatch; the cause is investigated.
5. **Windows.** D_k is the outward displacement between the samples at the window ends. s_k = D_k / |W_k|. v is the
   least-squares slope of outward position over samples with t ≥ T − 1000. The speed-decrease condition applies to each
   interface separately (|s₃| < |s₂| < |s₁|, or |D₃| ≤ λ/400). The cell-count range uses samples with t ≥ W₁ start.
6. **Rules 2–3 as MG9's code.** Fill or loss is checked from t = 0. An "untracked" sample at any time gives
   UNRESOLVED-INTERNAL. BREAKUP/REMOTE and λ/2 reversals are checked on samples with t ≥ 100.
7. **Robustness.** It runs only for main runs labelled SETTLED with a passed replay, plus 29.22 (descriptive; it does
   not enter Q1).
   - (a) D01: `fourier_resample` to N = 1000, then settle to 1000 (extension to 3000 if the residual exceeds 1e-6);
     `prep_classify` must give PERSISTS (34 cells).
   - (b) D02L: tile twice, integrate to 300; PERSISTS = complete, residual ≤ 1e-6, 68 cells. This preparation is shared
     with the descriptive half-width-50 trial (d).
   - (c) The main run's final field (its snapshot at T), multiplied by 1.05 and 0.95; E_pat of the main run; a fresh
     tracker. PASS = SETTLED and both final interfaces within λ/10 of the main run's final interfaces.
8. **Q2.** Midpoints are rounded to 6 decimals. Each is prepared from the nearest sweep state above (`carry_concentration`,
   settle to 1000, extension to 3000, `prep_classify`). The final settled point gets (a)–(c). If all three midpoints
   invade, the bracket is [29.35, last midpoint] [documentation slip corrected per ChatGPT MG030; the code already did this, and Q2 was not run], with 29.35's robustness already established. The interval
   comparison rounds both differences to 1e-9 before applying the inclusive ±0.01 limits.
9. **Records.**
   - Per-run JSON in `results/MG9c/runs/`: the full track, the initial-state SHA-256, E_pat, success, time reached,
     residual, and snapshot names with SHA-256.
   - Snapshots are saved as `.npy` during the run (gitignored), then packed into `results/MG9c/snapshots.npz` with
     `snapshots_manifest.json` (committed).
   - Settled preparation states are saved the same way.
10. **Order and stops.** `q1` → `robust` → `q2` → `desc` → `evaluate_mg9c.py` → `pack`. The cap (60) and cutoff
    (2 core-h) are checked before each integration; on a stop the driver exits with code 2 and the PI decides.

**MG9c RESULT (MG030, 2026-11-04; cloud, Python 3.13.18, numpy 2.5.3, scipy 1.18.1; code `cee5498`): Q1 UNRESOLVED;
Q2 not eligible.** 11 integrations, 0.23 core-h (cap 60, cutoff 2). Evaluator `results/MG9c/evaluation.json`; records
`results/MG9c/runs/`; snapshots `results/MG9c/snapshots.npz` (+ SHA-256 manifest). PI chose "Start the runs now"
(decision card 20:40 UTC).
- **Replay:** all three main runs reproduce MG9's [0, 1800] tracks exactly (max |Δx| = 0, 361/361 samples, E_pat
  identical).
- **Q1 labels:** 29.25, 29.30, 29.35 → **UNRESOLVED-SLOW** (rule 6 after the extension to T = 30,000). Each fails
  only the residual criterion: final residual 1.084e-5, 1.092e-5, 1.085e-5 against the registered 1e-5. Every other
  SETTLED condition holds:
  - |D₃| over [15,000, 30,000] is 0.008–0.009, against λ/40 = 0.074;
  - |v| over the last 1000 units is 4.3–4.5e-7, against λ/40,000 = 7.4e-5;
  - speeds decrease for both interfaces;
  - the cell count is constant: 18, 18, 19.
- **So Q1 = UNRESOLVED** (§7 order 1). Robustness was not run for these points (it requires SETTLED), and Q2 is not
  eligible.
- **Descriptive, not scored (29.22):** SETTLED (residual 9.45e-6). Robustness was ROBUSTNESS-UNRESOLVED:
  - D01: preparation PREP-CHANGED (the 34-cell pattern does not persist on h = 0.1 at 29.22), so the trial was UNRESOLVED;
  - D02L, half-width 25: PASS;
  - ±5% perturbations: PASS (edges return within 0.0023).
  
  Width trials: both SETTLED. The D02 half-width 25 + λ/4 start settled at the same edges as half-width 25, to
  within 1.6e-4 (0.000157 and 0.000147; the width differs by about 1e-5) [corrected per ChatGPT MG030; was "within 1e-4"]. The D02L half-width 50 start settled at ±55.7.
- **Descriptive, not scored (shape of the late motion, from the committed tracks).** The patch width is nearly
  constant, at 58.1–58.3 ≈ 19.8 λ: over [15,000, 30,000] its mean rate of change is 1.2–2.5e-8 per unit time in
  absolute value (interval averages, not bounds at every sample) [corrected per ChatGPT MG030; draft said "after
  t ≈ 6000 … frozen … less than 4e-8", which fails over [6000, 30,000], where the rates are 4.4e-8 to 1.0e-7]. The two interfaces move in the same
  direction (the outward displacements have opposite signs), at about 5e-7 per unit time and slowing (from 1e-6 at
  t ≈ 4000 to 4.7e-7 at 27,500). So the late motion is predominantly translation of the patch, not growth. ChatGPT MG030's post-result projection of the final RHS onto translation leaves a 1.4–1.6% remainder; this is not proof of exact rigid motion. Its size is
  consistent with the final residual of about 1.1e-5. No re-scoring follows from this reading.
- **Independent recompute** (separate script on the raw tracks and snapshots): all labels, D₃, v and cell counts agree.
  Residuals recomputed from the 30,000 snapshots agree with the records.
- **ChatGPT result review MG030** (`review/MG030_chatgpt.md`; replication `review/MG030_replication/`): recommends
  accepting Q1 UNRESOLVED and Q2 NOT ELIGIBLE. Every label, the 50 snapshot hashes, the residuals (to 1e-13) and the
  ledger were reproduced independently. Two descriptive corrections and one documentation slip were accepted (marked
  above). Response `review/MG030_response.md`. Awaiting PI sign-off.

**MG9c CLOSED 2026-11-04: Q1 UNRESOLVED; Q2 NOT ELIGIBLE.** PI sign-off: Richard Benear, 2026-11-04 22:29 UTC, on the
project decision card ("Sign off and close MG9c"), after ChatGPT result review MG030 and its response
(`review/MG030_response.md`). Scope: at b = 29.25, 29.30 and 29.35 the patches stopped growing, but they were not
stationary by the registered residual criterion (1.08–1.09e-5 vs 1e-5). An exploratory post-result diagnostic
(`results/MG9c/diagnostics/`) indicates a coherent slow translation of nearly constant-width patches; its cause is
not identified. Any moving-frame or drift-cause criterion would be a new, separately registered question.


### MG9d — Cause of MG9c's late drift (bounded deterministic diagnostic; follow-up to MG9c)

**REGISTERED 2026-11-05 (PI decision card, 02:00 UTC: "Register; runs after code review"). No MG9d integration
preceded this registration.** The registered design is `spec/MG9d_proposal.md` **draft 3**, with content frozen at
`69ad83e` (drafts 1 and 2 kept as `spec/MG9d_proposal_v1.md` and `_v2.md`, superseded). It includes the PI's pre-run
decisions P0–P2 (01:50 UTC) and the endpoint-average speed estimator (§3, after ChatGPT MG031 r3).

- Reviews: ChatGPT MG031 (draft 1), r2 (draft 2) and r3 (draft 3), with responses `review/MG031_response*.md`.
- Pre-run implementation specification: `spec/MG9d_prerun.md` (MG032), on the same commit.
- Design: seven runs from the Nyquist-filtered b = 29.30, t = 30,000 MG9c state, T = 10,000, no extension. R0
  baseline, R1 tight solver, R2a–c sub-grid shifts, R3 reflection, R4 D01.
- Rules: validity §3; C1–C5 and the ordered overall reading §4; INCONCLUSIVE allowed.
- Budget: cap 10 integrations, cutoff 1 core-h, cloud.
- **Runs start only after ChatGPT's implementation review and a separate PI go.** `results/MG9d/REGISTERED`, which
  `run_mg9d.py runs` requires, is created at that go.
- Nothing here changes MG9c's classifications or thresholds.

**MG9d RESULT (MG033, 2026-11-05; cloud, Python 3.13.18, numpy 2.5.3, scipy 1.18.1; code `9407f59`): MESH EFFECT
SUPPORTED (on D02; the finer-mesh comparison is unavailable because R4 is INVALID).** 7 integrations, 0.084 core-h
(cap 10, cutoff 1). PI run go: thread message 02:12 UTC. Records are in `results/MG9d/runs/`; the evaluation is
`results/MG9d/evaluation.json`; snapshots are `results/MG9d/snapshots.npz` (+ manifest); the figure is
`paper/figures/mg9d_drift.png`. Speeds are endpoint averages of the alignment displacement d (v_f = 5e-9; resolved at 1.5e-8).

| Run | a [2000, 6000] | b [6000, 10,000] | v [2000, 10,000] | Valid |
|---|---:|---:|---:|---|
| R0 baseline | −3.74e-7 | −3.16e-7 | −3.45e-7 | yes |
| R1 tight solver | −3.74e-7 | −3.16e-7 | −3.45e-7 | yes |
| R3 reflected (raw; adjusted = −raw) | +3.74e-7 | +3.16e-7 | +3.45e-7 | yes |
| R2a +h/4 | −1.41e-6 | −1.41e-6 | −1.41e-6 | yes |
| R2b +h/2 | +5.54e-7 | +6.54e-7 | +6.04e-7 | yes |
| R2c +3h/4 | +1.23e-6 | +1.12e-6 | +1.18e-6 | yes |
| R4 h = 0.1 | −2.7e-11 | −4e-13 | −1.4e-11 | **no (item 5)** |

- **C1 NUMERICS-AGREE.** R1 matches R0 to about 1e-11 in every slope.
- **C3 SYMMETRIC.** The reflection-adjusted R3 equals R0.
- **C2 ALIGNMENT-DEPENDENT.** All three shifted arms DIFFER from R0. The speed depends on the patch's sub-grid
  position, in size and in sign: +h/4 moves left four times faster than R0, while +h/2 and +3h/4 move right.
- **C4 DECREASING.** R0 has ρ = b/a = 0.845, and R1 gives the same label.
- **C5 R4 INVALID / LOST-IDENTITY.** R4 failed item 5: S_ref(T) = 5.7e-3 against a limit of 1e-3. On the finer mesh
  the patch kept changing shape after t′ = 2000 (width 58.316 → 58.323, settling by about 6000). Its identity checks
  passed (18 cells throughout, ΔW = 0.12 < λ/2).
- The overall reading is rule 4, MESH EFFECT SUPPORTED, with C4 DECREASING alongside.

**Independent recompute** (from the raw snapshots, by Newton iteration on the Fourier cross-correlation, not via
`run_mg9d.align`): every displacement agrees with the record within 9e-10, and C1–C4 come out the same. All validity
values for R0–R3 are within their limits. The tightest consistency check is R2b's C_phi2: 1.9e-4 off d against a
band of 9.8e-4.

**Descriptive only (not scored):**
- Read together, the four D02 arms are consistent with a lattice pinning landscape. *Hypothesis, not measured
  [corrected per ChatGPT MG033]:* the patch may slide toward a preferred sub-grid position (R0 left; +h/4 left fast;
  +h/2 and +3h/4 right), with R0 slowing. No pinning landscape, preferred position or eventual arrest was measured, and
  lattice-dependent transient dynamics over a finite horizon are not excluded.
- On D01, after its shape adjustment, R4's centre does not move (|v| ≈ 1e-11, far below resolution; its digits are
  not a resolved speed), and its final residual is 2.1e-9. That is well below MG9c's 1e-5 stationarity threshold,
  which the D02 patches missed. R4 is INVALID by item 5 only (shape, S_ref = 5.7e-3); it is not rescored, and the
  validity reference time is not moved after the data.
- R0 ends at residual 7.2e-6, below MG9c's 1e-5 threshold, while its motion stays resolved [added per ChatGPT MG033].
  This is consistent with slowing, but MG9d starts from a filtered state; it does not change MG9c's classification.
- These observations support the mesh-effect reading but are not registered results. R4 is INVALID under the
  registered rule. The early-transient caveat (ChatGPT MG031 r3) applied to R4: its adjustment lasted until about
  6000, not 2000.
- One refinement is a sensitivity statement, not a continuum limit. MG9c's classifications are unchanged.

**MG9d CLOSED 2026-11-05: MESH EFFECT SUPPORTED (on D02; finer-mesh comparison unavailable, R4 INVALID by shape;
C4 DECREASING alongside).** ChatGPT result review MG033 independently replicated all records and labels and found no
result-blocking error. Response: `review/MG033_response.md`, with wording corrections marked "[corrected per ChatGPT
MG033]". PI sign-off: Richard Benear, thread message 2026-11-05 02:32 UTC ("the run was accepted by ChatGPT ... Let's
proceed to MG9b now"). MG9c stays closed as scored. Asymptotic pinning is left to optional, separately registered
work. MG9b proceeds on explicitly lattice-conditional terms, with its own deterministic event-confirmation checks.


### MG9b — Does the noise-activation slope S change with b? (stochastic; lattice-conditional)

**REGISTERED 2026-11-05 (PI decision card, 03:01 UTC: "Register draft 3 after merge"). No MG9b run preceded this
registration.** The registered design is `spec/MG9b_proposal.md` **draft 3**, with content frozen at `185aceb` (drafts
1 and 2 kept as `spec/MG9b_proposal_v1.md` and `_v2.md`, superseded). It keeps the PI's five §8 decisions (02:39 UTC).

- Reviews: ChatGPT MG034 (draft 1) and r2 (draft 2), with responses `review/MG034_response.md` and `_r2.md`.
- Question: under a linear model in b, does S change with b? NEGATIVE / POSITIVE LINEAR TREND / NOT RESOLVED (§5),
  with per-b slopes, the linear-model check, branch-specific calibration labels and the exploratory rule.
- Design: MG3 protocol at τ = 6.667e-4, b = 29.20 and 29.35 plus MG3b's b = 29.00; Part A placement scan; Part B three
  cells from Ω_T and the ordered candidates Ω_T − 15, − 30, + 15, + 30; MG3b batching (≥ 20 events, cap 64).
- Controls V-b1 (reference validity; failure stops that b), V-b2 (late-window confirmation, completion required),
  V-b3 (tight-tolerance confirmation).
- Budget: run caps in §6 (534 stochastic runs plus controls), about 35 core-h projected, cutoff 50 core-h on the PI's
  PC; synthetic calibration in the cloud, cap 4 core-h.
- **Order before any production run:** code and pre-run specification; synthetic calibration (cloud); ChatGPT review
  of the calibration and implementation; then a separate PI go for the PC run.

- **Amendment 1 (ADOPTED 2026-11-05, PI card 03:27 UTC, after ChatGPT MG034 r3; before any production run):** budget truncation of
  rate-data collection at any stage makes a directional result EXPLORATORY (study-level flag), not only a BUDGET-SHORT
  cell; a stop that cuts only controls stays CONTROL-INSUFFICIENT. `spec/MG9b_amendment_1.md`.

- **RESULT (2026-11-05; PC run, records `5a54f7a`; `spec/MG9b_result.md`): NOT RESOLVED.** Branch B3, calibrated; no
  qualifiers. S₁ = −0.023 per unit Ω per unit b, 95% profile CI [−0.094, +0.048]; LR linear-model check p = 0.30.
  Per-b slopes: 0.064 (29.00, MG3b), 0.048 (29.20), 0.057 (29.35). V-b1, V-b2, V-b3 PASS at both b. 280 runs, ledger
  26.28 core-h, no budget stop. Recorded deviation: two tight-tolerance V-b3 confirmations at 29.20 took 5.4 h and 4.4 h
  (median 18 s); both completed and agreed; no time limit was registered. Descriptive: Ω_T at 29.35 is at the top of
  the grid (250). Independently recomputed from the raw runs (`src/check_mg9b.py`).
- **CLOSED 2026-11-05 (PI decision card, 14:03 UTC): NOT RESOLVED.** ChatGPT MG035 reproduced the fit, controls and
  ledger; two wording corrections applied. The event-shape pattern (single-minimum share falls with lower Ω at both b;
  post hoc, pooled p ≈ 0.04) is archived as a hypothesis for a possible separately registered study.

### MG8 — Scheduler sensitivity of the lattice approximation (multiway / foliation; Central Question 5)
Status: **REGISTERED 2026-11-05 (MG036d)**, with PI approval (decision card, 20:02 UTC: adopt all §9 recommendations and
register). **No MG8 run of any kind preceded this commit.** The registered design is `spec/MG8_proposal.md` **draft 4**,
with its content frozen as committed at `3ab4ec8` (only its status line changed since). Drafts 1–3 are preserved as
superseded (`spec/MG8_proposal_v1.md`, `_v2.md`, `_v3.md`). Reviews: ChatGPT MG036 rounds 1–3 (`review/MG036_chatgpt*.md`);
responses `review/MG036_response*.md`. Round 3 found no remaining design issue after its one correction.

Summary (the proposal governs):
- **Setting:** MG2/MG3 1D lattice, P17, b = 29.0, h = 0.2, L = 100, τ₀ = 6.667e-4; reference soliton start.
- **Schedulers:** J (registered synchronous kernel), R (rightward sweep, random seam each step), P (random sequential),
  L (mirror check only). One snapshot of all species per node update; changes applied together.
- **B1 (primary, deterministic expected counts):** J and three seam replicates of R at τ₀/2^k, k = 0–6, T = 300; L checks;
  J and R continued to T = 600 at k = 6. 32 runs.
- **Production step (§4.2):** FATE (SAME/CHANGED), SHAPE (INSUFFICIENT RESOLUTION / UNSETTLED / MATERIAL ≥ 1e-2 / SMALL /
  NEGLIGIBLE < 1e-3), DRIFT (INSUFFICIENT RESOLUTION / UNSETTLED / DRIFTING ≥ 1e-6 / NOT DRIFTING AT THIS HORIZON); summary
  SENSITIVE / PRACTICALLY INSENSITIVE / MINOR / INSUFFICIENT RESOLUTION / UNRESOLVED; SEAM-DEPENDENT if replicates disagree.
- **Refinement (§4.3):** CONVERGES TO J / BELOW MEASUREMENT RESOLUTION (not convergence) / NOT RESOLVED BY THIS LADDER.
- **B2 (descriptive):** J, R, P at Ω = 1000, 8 runs each with independent streams, plus one L mirror; T = 500. Fate,
  drift difference (Welch), stochastic shape (bootstrap), conditional on eligibility.
- **Interpretation:** the primary conclusion is the scheduler sensitivity of this lattice approximation. Central Question
  5 stays open whatever MG8 finds.
- **Part 0** (priced separately, cap 3 core-h): diffusion-only operator benchmark with the frozen refinement and
  conditional-order tests; J regression; R/L mirror; clipping; synthetic validation; timing. STOP on failure.
- **Budget:** 57 runs, cap 64; production cutoff 20 core-h (projected 5–18); PC, 8 workers; BUDGET-SHORT rules as §6.

Result (MG038, 2026-11-06; `results/MG8/MG8_results.md`): all implementation checks PASS. **Production step: SENSITIVE**
(all three R replicates lose the particle at t ≈ 8; J keeps it). **Refinement: NOT RESOLVED BY THIS LADDER** (eligible
k = 3–6 only; slope 1.000, residuals and horizon pass, but E at the smallest eligible τ is 0.125 of the largest, above the
0.1 limit). B2 (descriptive): losses J 7/8, R 8/8, P 8/8 (R and P within 20 units); drift and shape NOT ESTIMATED.
55 records, 3.09 core-h, cloud. Independent recompute agrees on all 106 quantities.
Replication (ChatGPT MG038, 2026-11-06): results accepted; recompute reproduced (106 quantities) plus 187 detection counts,
both mirrors, all stochastic loss classifications, the intervals and the ledger. Four reporting corrections made
(`review/MG038_response.md`): J losses all after t = 100 (MG2 horizon) and Ω = 1000 insufficient for T = 500; particle
"unchanged" withdrawn; causal reading provisional (diagnostic archived); seam spread small, not zero; full B2 counts
discarded (retention deviation).
PI sign-off: registration 2026-11-05 20:02 UTC; deviation D1 and production go 2026-11-05 23:37 UTC; **closure 2026-11-06 02:11 UTC** (registered results with the MG038 corrections). Central Question 5 stays open.


### MG10 — Explicit A/B reservoirs (DRAFT; MG039, 2026-11-06)

**Current status (MG043): first deterministic batch REGISTERED MG042 and completed; result UNRESOLVED; PI closure pending. Historical drafting and approval records follow.**

Design: [spec/MG10_proposal.md](spec/MG10_proposal.md), Draft 1. NOT REGISTERED; no scientific or development trajectory run.
Proposed sequence: A formal state/reaction/accounting model; B implementation qualification; C deterministic recovery of the fixed-reservoir baseline; D finite-capacity estimation; E causal supply interventions; F recomputation and PI closure.
Only A/B/C are the first proposed scope. D/E require separately frozen arms and budgets. Numerical margins, sample sizes, run caps and compute ceiling remain to be finalized before PI registration.
Continuous-time dynamics remain the proposed model; synchronous stepping is a numerical approximation. Finite A/B, inert Z/W and explicit chemostat inputs are modeling hypotheses, not a physical derivation of the vacuum.
Geometry-derived supply, classical correlations and an actual quantum-entanglement extension are later candidate studies, not authorized gates.

MG039b authorization, 2026-11-06: Richard approved MG10 Draft 1 decisions 1-4; the working proposal is now Draft 2, preserving Draft 1. Numerical limits still to be finalized before registration are not retroactively invented. Prefer local PC execution when runs could pressure assistant usage, with at most eight workers and one numerical-library thread per worker. Richard determines external review per run/batch. No MG10 run has occurred.

MG040 (2026-11-06): numerical candidate and compute ceilings specified in spec/MG10_run_plan.md. Static depletion arithmetic motivated the PI-approved addition of C=100000 and 1000000. First batch is deterministic: 64 solver trajectories, cap 4 core-hours plus 0.25 qualification/timing core-hour, local <=8 workers. No trajectory run, no measured forecast, and no production registration yet.

MG041 (2026-11-06): PI explicitly waived external review and authorized proceeding. Isolated availability qualification FAILED under MG040 check 4 (four of 21 Euler paths), independently recomputed using Decimal arithmetic. Production stopped as required; no recovery outcome or registration. Proposed qualification correction: spec/MG10_qualification_correction.md. Evidence: results/MG10/qualification; internal record: review/MG041_internal_check.md. All other implementation/qualification and timing remain pending.

MG042 authorization (2026-11-06): PI adopted spec/MG10_qualification_correction.md (“yes, go ahead”) and authorized continuing qualification, benchmark and conditional production. External review remains waived. The original MG041 failure is preserved. Registration will be prospective after qualification and measured projection.

MG042 REGISTRATION (2026-11-06, before production): MG10-A/B/C deterministic batch REGISTERED under PI authorization; external review waived; adopted MG041 qualification correction. MG10_run_plan.md version 2 and MG10_prerun.md govern. Qualification PASS; projected 6488.182472 worker-seconds (1.802273 core-hours) <=10800 allowed with reserve, 64 fixed jobs, <=8 workers, ceiling 14400 seconds, per job 600. Durable raw retention verified. Full source/spec/reference hashes in results/MG10/registration.json. Earlier Draft 3/run-plan Draft 1 preserved. No scientific production preceded this registration.

MG043 RESULT (2026-11-06): UNRESOLVED. All 64 registered trajectories completed; all accounting/oracle checks pass; separate recomputation agrees on all arm labels and batch verdict. Top capacities C=100000 and 1000000 recover for both starts, but C=1 homogeneous E2/O error 0.4996807 exceeds 0.001. Production 771.0576318 worker-seconds (0.21418268 core-hours), no timeout/budget stop. A posthoc raw-data check also identifies roundoff detections on exactly uniform depleted backgrounds; no detector or verdict is revised. Full raw arrays, ledger, evaluation and recomputation retained in results/MG10/production. Report: results/MG10/MG10_results.md; same-author internal check: review/MG043_internal_check.md. PI closure pending.


### MG10-C1 — Depleted homogeneous stability and resolved minima (MG044)

REGISTERED prospectively, 2026-11-06, under PI instruction to proceed with the next study. spec/MG10C1_proposal.md governs. Original MG10 UNRESOLVED result is unchanged. Ten jobs: two reduced references, two complete 250-mode linear-response scans, six nonlinear paths (zero and two prescribed cosine amplitudes, E2/O+), with prospectively frozen mode selection and linear-window rules. Detector qualification uses 23 fixed fixtures and a numerical prominence floor, not a particle discovery claim. Six development qualification tests and separate raw checks PASS. Six T5 timing cases project 462.838314 worker-seconds; production ceiling1800, qualification360, local<=8 workers with library threads1. No scientific production preceded registration. Full hashes/denominator/retention in results/MG10C1/registration.json. External review not requested with execution instruction. Closure remains PI decision.

MG045 RESULT (2026-11-06): ten of ten MG10-C1 jobs complete, no timeout/budget stop/replacement, production116.5758434 seconds (0.03238218 core-hours). HOMOGENEOUS REFERENCE RESOLVED (reduced disagreement8.1783e-10; uniform E2 error1.5568e-7, exact spatial uniformity). CONTROLLED AMPLIFICATION SUPPORTED (prespecified scan selects mode32; two amplitudes, both solvers, worst eligible relative error0.03295% versus5%). RESOLVED-MINIMUM QUALIFICATION PASS (23 fixtures; all301 uniform E2 states have zero new detections). Independent raw checker agrees; full data retained. Formal large tangent gains are not physical amplitudes or particle populations. Original MG10 UNRESOLVED result unchanged. Report: results/MG10C1/MG10C1_results.md. Same-author check: review/MG045_internal_check.md. PI closure pending.

MG046 PI CLOSURE (2026-11-06): Richard: “I sign off on MG10-c1. proceed to next.” MG10-C1 CLOSED with its three MG045 outcomes: HOMOGENEOUS REFERENCE RESOLVED; CONTROLLED AMPLIFICATION SUPPORTED; RESOLVED-MINIMUM QUALIFICATION PASS. Original MG10 remains UNRESOLVED and is not rescored. PI authorizes the next study; MG10-S1 finite-reservoir stochastic sampler qualification precedes a full finite-capacity particle survey.


### MG10-S1 — Finite-reactant stochastic sampler (MG047)

REGISTERED prospectively2026-11-06 under PI instruction to proceed after MG10-C1 sign-off. spec/MG10S1_proposal.md governs. Sixteen small-system groups (SSA and three leap maxima for four fixtures), R32768/group in64 fixed chunks512. Primary96 marginal CDF comparisons at finest step with simultaneous99% DKW/union-bound intervals and .05 equivalence margin. Coarse distances descriptive. Bounded joint B/X exits and A conversion; no clipping. Eight qualification tests and separate raw checks PASS; 1024 unique chunk seeds and full source/spec hashes frozen in results/MG10S1/registration.json. Measured forecast1032.297616 worker-seconds fits new3600 ceiling with25% reserve; qualification77.2601245 seconds under360; local<=8 workers, threads1. No scientific group preceded registration. Original MG10 unchanged; MG10-C1 PI CLOSED. Full lattice/particle claims outside scope.

MG048 RESULT (2026-11-06): QUALIFIED ON TESTED SMALL SYSTEMS. All16 groups x32768 paths complete; all96 finest-step marginal CDF simultaneous99% upper bounds below .05 (maximum .0361552912; empirical maximum .0107727051). Exact integer ledgers/nonnegativity/probability safety pass; separate checker agrees on all primary bounds and label, hashes, inputs, seeds and cost. Production176.2166985 worker-seconds (0.04894908 core-hours), no timeout/budget stop/replacement; all raw checkpoints/ledgers retained. Qualification limited to these cases/marginals, not full joint law or actual h.2/high-Omega lattice particle dynamics. Report: results/MG10S1/MG10S1_results.md; same-author check review/MG048_internal_check.md. MG10-S1 PI closure pending; original MG10 unchanged; MG10-C1 PI CLOSED.

MG049 PI CLOSURE (2026-11-06): Richard: “i sign off on S1. Proceed to the next.” MG10-S1 CLOSED as QUALIFIED ON TESTED SMALL SYSTEMS. Original MG10 remains UNRESOLVED; MG10-C1 remains CLOSED. Next authorized study is MG10-S2: actual h=.2/Omega100,1000 spatial-noise and lattice sampler qualification, prospectively specified before production.


### MG10-S2 — Actual-lattice spatial-noise qualification (MG050)

REGISTERED prospectively2026-11-06 under PI proceed instruction after S1 closure. spec/MG10S2_proposal.md governs; Draft1/calibration preserved. N500 H/S Omega100/1000 one-step laws versus independent snapshot oracle, and short N16 exact SSA comparisons at actualh.2. Revised R16384/group follows a static precision correction before scientific data; alpha/margins/physics unchanged. Fixed88 jobs/32 groups/524288 paths,1280 unique chunk seeds. Primary192+96 selected spectral-power/cross-product CDF comparisons, separate alpha.005 families and .08 margin; both must qualify. Eight qualification tests and separate raw check PASS. Rescaled measured forecast3062.691866 worker-seconds fits new7200 ceiling with reserve; qual102.7158449 seconds under900, maxjob600, local<=8 workers/threads1. Full reconstructible event/population shards retained. No scientific group preceded registration. Original MG10 unchanged; S1 CLOSED.

MG051 RESULT (2026-11-06): QUALIFIED WITHIN REGISTERED S2 SCOPE. All88 jobs/32 groups x16384 samples complete. Snapshot192 CDF bounds and short96 bounds meet.08 (maxima .0573459/.0507281). Exact integer reconstruction, actual population hashes, initial means, noise matrices, seeds and independent merging-CDF calculations agree. All1280 reconstructible raw shards retained and uploaded in data batches. Production1304.9744223 seconds (0.3624929 core-hours), no timeout/budget stop/replacement, maximum poll gap.251052 seconds. Coarser short steps descriptive; no long-time particle/noise-physics claim. Report results/MG10S2/MG10S2_results.md and same-author review/MG051_internal_check.md. S2 closure pending; S1 CLOSED MG049; original MG10 unchanged.


MG052 PI CLOSURE (2026-11-06): Richard: “I sign off on S2 results. Proceed with the next step.” MG10-S2 CLOSED as QUALIFIED WITHIN REGISTERED S2 SCOPE. S1 and MG10-C1 remain CLOSED; original MG10 remains UNRESOLVED. Next step is MG10-D1 finite-supply estimation design, implementation qualification and disclosed development timing. A new long-horizon production matrix and budget require prospective registration; no scientific production yet.

MG052 DEVELOPMENT DISPOSITION: MG10-D1 Draft2 prepared;20 implementation checks and separate36-archive reconstruction/CDF audit pass. Three synthetic measurement checks pass after a disclosed measurement-only roundoff correction; sampler source unchanged. Timing20 classes projects41.5715core-hours including25% allowance, about5.20hours on8workers. Original48-hour reserve check fails and is preserved; proposed60-hour ceiling passes planning reserve, pending PI approval. Eighty candidate T60 jobs; no production registered or run. Full source/stream audit and development raw retained.


MG053 AUTHORIZATION / REGISTRATION: Richard: “I approve of the run.” MG10-D1 REGISTERED prospectively under spec/MG10D1_prerun.md and approved Draft2. Eighty fixed T60 paths, N500/h.2/Omega1000, five supply arms, H/S, two steps, R4 each. Twenty-one driver checks and synthetic complete/incomplete end-to-end pipeline PASS; earlier qualification retained. Forecast41.682623core-hours <=45 allowed reserve; approved60-core-hour production cap, <=8workers/threads1,3600seconds/job. All80 seeds unique/no collisions against2376 collected earlier/development values; actual stocks, source/spec/reference hashes and retention frozen in registration.json. No scientific production before registration commit/push. External review not requested with execution approval; same-author checks only.


MG054 D1 RESULT / FOLLOW-UP: COMPLETE DESCRIPTIVE PILOT as registered:80/80 complete,33.333427core-hours, no timeout/budget stop, separate integer/raw/480-summary audit PASS. All80 raw archives retained. Posthoc HIGH fine-step A-source deficit0.300-0.304% against autonomous analytic expectation; coarse agrees. HIGH recovery interpretation UNRESOLVED pending separately registered numerical diagnostic/correction; original verdict/data not rescored. PI: “looks ok, proceed with the next steps.” Record as acknowledgment/authorization to investigate and archive, not evidence that source bias is resolved.


MG055 DIAGNOSTIC / CORRECTED REGISTRATION: Source-law plan registered e79fa73 before samples. Revised2 primary+4 supplemental scalar laws qualify; independent80-digit Decimal/histogram and raw checks PASS. Original/NumPy primary CDF bounds also pass the deliberately limited margin; direct legacy fine P(0) error and mean deficits documented, no retrospective CDF FAIL. Twelve low-stock regressions and4 corrected-stream matches PASS; synthetic32/32 and31/32 reporting/check PASS. New32-path CHEMO/HIGH T60 study REGISTERED before production; forecast15.475784core-hours <=18 planning reserve, new24-hour cap within original60 (combinedworst57.583427), local8workers/threads1,3600seconds/job. Source/spec/reference hashes, unique seeds and retention frozen. D1 remains unchanged and fully archived with caveat; no new scientific trajectory yet.


MG056 D2 RESULT: COMPLETE DESCRIPTIVE PILOT32/32, production13.373943core-hours under24; elapsed1.850166hours; no failure/timeout/cap. Independent raw/seed/source/integer reconstruction and192-summary check PASS. Old0.30% HIGH fine-source deficit absent; new offsets+0.0051%/+0.0062% explicitly reported (approximate intervals do not all cover expectation), not equivalence. Fine HIGH/CHEMO mean-X differences under0.06% descriptive; seeded contrast/particle identity and numerical equivalence unresolved withn4. Full raw retention/upload follows; no old result rescored. PI closure pending.


MG057 PI CLOSURE (2026-11-07): Richard: “I approve closure for D2.” MG10-D2 CLOSED as COMPLETE DESCRIPTIVE PILOT under its registered scope. All32 paths and independent raw/source/seed/integer/192-summary checks completed; full archive uploaded and verified. Residual source offsets, n4 uncertainty, numerical shape sensitivity and unqualified particle identity/persistence remain explicit. Closure is not stochastic equivalence or a particle-survival verdict. Original D1 and earlier outcomes remain unchanged; no new production registered or launched by closure.


MG058 REGISTRATION: MG10-P1 single-reference X-motif measurement qualification, authorized next step after D2 closure. Frozen68 frame fixtures,15 temporal fixtures,7 invalid-API assertions; independently recompute by dense Gaussian convolution/direct template dots and separate conservative tracking. Prospective criteria in spec/MG10P1_proposal.md; all raw synthetic inputs/expectations/seeds retained. One worker/threads1;300-second total owned budget,240 fixture/60 oracle timeouts. No stochastic production or D1/D2 rescoring; qualified sampled field continuity would not establish particle lineage/survival. Registration source/reference/matrix hashes pushed before execution.


MG059 P1 RESULT: QUALIFIED ON REGISTERED SYNTHETIC FIXTURES. All90 assertions (68 frame/15 temporal/7 API) pass; independent dense-convolution/direct-score and temporal oracle agrees on511 raw frame assessments. Total owned charge2.016197000seconds under300, no failure/timeout/stop, no criteria revision/retry. Scope is sampled single reference-like X continuity, not physical lineage/particle survival; actual model-noise and coupled-field/depleted-shape generalization remain open. No D1/D2 rescoring or new stochastic trajectory. Full raw decision archive retained/upload follows; PI closure pending, no next production registered.


MG060 PI CLOSURE (2026-11-07): Richard: “I approve,” in response to the MG10-P1 closure request. MG10-P1 CLOSED as QUALIFIED ON REGISTERED SYNTHETIC FIXTURES. All90 assertions and independent511-frame dense/direct recomputation completed; all16 raw archives uploaded and retention verified. X-only, fixed-reference, sampled-continuity and noise/lineage limitations remain explicit; closure is not a physical particle/survival or real-noise false-positive verdict. Earlier outcomes preserved. No next production registered or launched.


MG061 P2 PREPARATION: New model-noise observability/calibration draft prepared, scientific production NOT REGISTERED.1024 candidate paths/8 cells atT.25,128 percell, two prepared supply states, H/S and two steps. Proposed5% H reference-match cap via4 simultaneous95% exact upper bounds; S retention/coupled-field observability descriptive. State/projection/CP and independent8-path timing-ledger checks PASS after preserved metadata-only preparation correction. Forecast2.168291core-hours/~16.26minutes on8, candidate4-hour cap with reserve. New population/criterion/budget require PI decision and final production freeze; no scientific run/old rescoring.


MG062 P2 AUTHORIZATION/REGISTRATION: Richard: “I approve the plan. Proceed.” Approved1024-path/8-cell P2 protocol and4-hour cap registered before production under spec/MG10P2_prerun.md. Original preparation/P1/D2 source unchanged; full synthetic complete1024/incomplete1023 pipeline and separate4096/4092 frame/metric/exact-bound audits PASS. All1024 path streams and128 warmups distinct/disjoint, no collision against2464 recorded prior/development values. Forecast2.168291core-hours <=3 reserve threshold;14400worker-second cap,8workers/threads1,180seconds/job. Primary4 H endpoint boundsalpha.0125,5% cap; no family score on incomplete batch. Complete/partial prefix raw retention and automatic audit prepared. Registration commit/push before scientific samples.


MG063 P2 RESULT: SHORT-HORIZON UNSEEDED MATCH CAP SUPPORTED ON TESTED CONDITIONS.1024/1024 paths,128/128jobs; four H cells0/128 each, simultaneous>=95% one-sided upper.0336552 below.05; S128/128 each descriptive. Independent raw/source/seed/state/integer/4096-instrument-frame and40 metric-summary-array audit PASS. No failure/timeout/cap. Elapsed14.223561minutes, production1.814736core-hours under4. G global contrast-noise proxy<1 does not prevent known-phase ensemble projection; all observability limits retained. No physical particle/survival/stationarity or CTMC equivalence claim. Full raw archive upload follows; PI closure pending.


MG064 PI CLOSURE (2026-11-07): Richard: “I sign off of P2.” MG10-P2 CLOSED as SHORT-HORIZON UNSEEDED MATCH CAP SUPPORTED ON TESTED CONDITIONS. All1024 paths, independent4096-frame/raw/metric/exact-bound audit and verified128-file raw archive retained. Four primary H groups0/128 with simultaneous upper3.3655%; seeded retention descriptive. All short-time, prescribed-state, noise/observability and physical-lineage/survival limitations preserved. No prior outcome rescored; no next production registered or launched.


MG065 P3 PREPARATION: Draft64-path FULL baseline temporal estimation, H/S,two steps,R16,T25/quarter grid;20-unit initial-S visibility and any-H episode frequencies descriptive, no inherited5% cap or physical survival label. Fixed/all-path vs conditional X-aligned0..20 averages distinguish movement/selection and temporal correlation. All18 preparation assertions PASS; independent8 episode/4 averaging cases and12 raw archives agree; no failed qualification/old-result rescoring. Forecast14.147329core-hours/~1.768hours on8; proposed NEW24-hour cap awaiting PI approval. No scientific production registered or started.


MG066 P3 APPROVAL/REGISTRATION: Richard: “i approve.” Approved64-path FULL baseline temporal-estimation/T25,20-unit sampled events and24-core-hour ceiling registered before science. R16 H/S/two steps,101quarter checkpoints.18 prep assertions,4 T1 stream matches, full/partial algebraic pipeline and separate6464/6363-frame/direct-average/exact-stat audit PASS; allowance377.296947seconds <900.128 scientific/warmup values unique/disjoint, no collision against3620 recorded prior/development. Forecast14.147329core-hours <=18 reserve; local8threads1,1800seconds/path, owned cap86400. No fraction success threshold/physical survival verdict; initial S vs any H episodes and conditional averaging preserved. Source/spec/reference/stream/retention freeze committed/pushed before launch.


MG067 P3 RESULT: COMPLETE TEMPORAL ESTIMATION,64/64 jobs/paths, no failure/timeout/cap; independent raw/source/seed/state/integer6464-frame/direct/statistical audit PASS,96 temporal averages. S16/16 initial20 event each step, H0/16 any20 episode each, per-cell descriptive intervals retained (no5%/pooled claim). Fine one initial signal censors first22.75 after reaching20, not physical death. All32 S aligned-eligible. Mean residual ratios G4.4-4.5/X~2/Y~2.6, reference projections~.98-.99 descriptive; G per-path residual still material, no new classifier/equivalence. Elapsed90.296155minutes, production11.460091core-hours under24. Full raw upload/verification follows; PI closure pending, no next production registered.


MG068 PI CLOSURE / REPORTING CONVENTION (2026-11-07): Richard requested a one-paragraph project-level summary in every stage result, supplied the desired P3 paragraph, and approved P3 conditional on adding it. Exact paragraph added near the top of P3 results; condition satisfied. MG10-P3 CLOSED as COMPLETE TEMPORAL ESTIMATION. Registered data/verdict/uncertainties/censoring/physical-claim limits and verified archive unchanged. New prospective summary convention recorded in RESEARCH_PROTOCOL.md and AGENTS.md; no prior result rescored or automatically rewritten. No next production registered.


MG069 E1 PREPARATION: DRAFT spatial A-source intervention pilot qualified, not scientific registration.32 T5 parents/128 T20 forks, local/remote/uniform/sham, H/S/two steps/8 blocks. Primary fixed-origin X seed-by-location contrast subtracts unseeded forced background.15 assertions/3 source ensembles and separate35 raw checks PASS; no failure or tuning. Forecast22.642613core-hours/~2.83hours on8, proposed NEW36-hour cap awaiting PI. Preparation185.646991worker-seconds under900;320 candidate and55 development seeds distinct/disjoint from2726 recorded prior values. All-one kernel identity supports collapsed untouched/sham. Final full production pipeline/controller/audit qualification and registration after approval. No physical death/reservoir-origin claim or old rescoring.


MG070 E1 APPROVAL/REGISTRATION: Richard: “I approve the E1 study. Proceed running it.” Approved MG069 design/NEW36-core-hour cap unchanged.32 fixed-T5 parents then128 T20 gated forks, two steps/8 blocks/H/S/four arms; primary fixed-origin X late15..20 seed-by-location interaction, nominal97.5% t intervals per step, no mortality/equivalence claim.20 held-development stream replays match exactly; complete128/missing127/available-prefix pipelines and separate dense-adjoint projection/profile/contrast/statistical audits PASS; missing step counts7/8. Total preparation allowance335.932328seconds<900. Exact320 streams, source/spec/reference hashes, masks and budgets frozen in registration before science; forecast22.642613core-hours<27 reserve.8workers/threads1, parent600/branch1800second timeouts, aggregate phase charge129600 ceiling and full raw retention. Commit/push before launch. No external review requested or old rescoring.


MG071 E1 RESULT: COMPLETE INTERVENTION ESTIMATION,32/32 parents and128/128 branches, no failure/timeout/budget stop;160 raw/fork/source/projection/profile/statistical audits PASS. Primary fixed-origin X seed-by-location contrast coarse-0.833949 nominal97.5%CI[-1.102872,-0.565027], fine-0.856842 CI[-1.001798,-0.711885], n8/step; registered direction-consistent model-pilot evidence, all block contrasts negative. Local seed-excess~.02 vs remote~.86-.87; uniform projection increases but no particle-health inference. All intervention S arms fail endpoint reference match including remote/uniform; primary field response is not mortality.19.808474core-hours under36, ~153.19minutes through audit, no replacement/rescore. Full archive publication/verification follows; PI closure pending, no next production registered.


MG072 PI CLOSURE (2026-11-08): Richard: “I sign off”, following verified publication of E1 results/full archive. MG10-E1 CLOSED as COMPLETE INTERVENTION ESTIMATION.32/32 parents and128/128 branches complete, separate160-file raw/fork/projection/statistical audit PASS; all160 published archives match audited hashes. Registered direction-consistent model-pilot resource-location interaction retained; small-n approximate intervals, fixed-origin measurement, finite horizon and physical-mortality/reservoir-origin/entanglement limits unchanged. No criteria/data rescore, no new scientific production authorized.


MG073 E2 DEVELOPMENT: Richard authorized developing E2. DRAFT reversible local A shutoff/restoration design,32 fresh fixed-T5 parents/128 T25 forks, pulse ends post5, four arms/H-S/two steps/n8 balanced blocks. Primary local-restore versus continuous-off seed interaction; secondary movement/shape/history diagnostics, no physical revival/quantum memory claim. Frozen preparation048d705 preceded8 short qualification/16 timing paths.11 assertions/7 field fixtures and24 raw/dense/statistical audits PASS, no failure/retry/tuning. Controlled157.545613seconds, total allowance162.907576<900. Forecast29.558218core-hours/~3h42 on8, proposed NEW40core-hour cap below30 reserve; PI final design/budget/public production retention pending. No science or old rescoring; full production implementation/raw-pipeline qualification and registration remains after approval.


MG074 E2 REGISTRATION: Approved32-parent/128-branch/T25 recovery protocol and40-core-hour cap/public result+full raw retention registered prospectively.20 held writer streams match every snapshot/statistic; complete128/missing127/prefix127 and empty-parent-stop raw-file pipelines pass separate integer/dense movement/profile/paired/Student audits. Missing primary counts7/8, empty0/0; no fabrication. Total preparation336.239390seconds<900.320 candidate/warm values unchanged, unique/disjoint with40 development and previous registrations. Forecast29.558218core-hours<=30 reserve.144000aggregate worker-second cap,8workers/threads1,600parent/2400branch timeouts, all parents before branches, prior phase charge carried over, whole-batch stop/no retry. Commit/push exact registration before science; no old result rescored.


MG075 DIAGNOSIS/AMENDMENT DRAFT: Original incomplete E2 preserved and all40 raw published hashes verified.9 held old/new timing pairs bit-identical, E2/old ratio.932..1.017; current CPU-bound probes do not establish past sustained-slowdown hardware cause. New admission-only forecast stop drains active jobs; hard cap/timeout/failure still terminate, synthetic tests PASS. Exact-prefix guard replays held pulse and rejects altered input.137.006255diagnosticseconds<300,9.494611controlseconds<30. Original19628.563985worker-seconds retained; proposed same-seed7replays+120unstarted, reuse32parents/1completed. Newforecast95.763055aggregate core-hours/~11h17remaining8workers; proposed144aggregate cap/5400seconds job DRAFT pending PI. No scientific restart, no rescore/new seed/denominator change. Final continuation orchestration/mixed-provenance audit must qualify before registration/science.


MG076 E2-C1 REGISTRATION: PI approved144total core-hour continuation cap inclusive of original19628.563984799664worker-seconds,5400/job,forecast admission-drain.32reused parents/1reused complete branch,7exact replay+120unstarted original seeds; unchanged128matrix/endpoints/arms/windows/n8 inference. Final mixed-provenance full/missing/prefix160-file/direct-statistical/source/fork/replay/aggregate-charge pipeline PASS,77.777412seconds<300; original fake tamper rejected. Original real E2 source/data/report preserved. Exact127job registration/source/reference/original-file manifest committed/pushed before science.518400cap,8workers/threads1,.25poll,hard cap/timeout/failure still terminate,no retry. Public E2 results/full-raw publication already approved. Forecast95.763055aggregate/~11h17remaining, uncertain fine throughput. No original rescoring/new seed/outcome selection.


MG077 E2-C1 RESULT: COMPLETE RESTORATION ESTIMATION,32parents/128branches;127continuation jobs complete with no admission/hard stop/timeout/failure.7exact prefix replay proofs,160raw/fork/source/direct-movement/statistical audits PASS. Primary coarse+.102771 nominal97.5%CI[.025840,.179701], fine+.140728 CI[.063955,.217502], n8/step: registered direction-consistent restoration response. Restored fixed-origin seed excess.11/.15 vs sham.95/.93; raw restoration largely unseeded background, not full recovery/physical identity/quantum memory. Aggregate28.363830core-hours incloriginal5.452379, continuation174.610150minutes, finishedNov8~12:02Mountain. Conservative11h17forecast exceeded actual runtime; earlier throughput variation cause remains unproven. Original stopped E2 not rescored/overwritten. Public full archive publication/verification follows, already PI-approved; scientific closure pending.


MG078 PI CLOSURE (2026-11-08): Richard: “i sign off”, following completed E2-C1 result and verified public full archive. MG10-E2-C1 CLOSED as COMPLETE RESTORATION ESTIMATION.32parents/128branches,127continuation jobs,7exact original-prefix replay proofs and160-file raw/field/statistical audit PASS; full published raw hashes verified. Registered modest direction-consistent restoration interaction retained; full-recovery/equivalence/physical identity/quantum-memory limitations unchanged. Aggregate28.363830core-hours including original stopped attempt preserved. Original E2 remains INCOMPLETE/immutable; no frozen criteria/data rescored, no next scientific production authorized.


MG079 ENVIRONMENTAL NETWORK DEVELOPMENT: PI authorized building test. MG11-N1 DRAFT fixed symmetric/equal-volume A/B graph and finite closed contacts built, no active-lattice change/replenishment. Exact small-system SSA and continuum/linear limits;8grouped assertions/15physical raw plus deterministic direct checks PASS, source2cf1358 frozen before development, no failure/retry/tuning. Planned32jobs/131072paths,3analytic endpoint laws+independent full-tiny SSA population laws, frozen DKW margins/families.64prospective+15dev streams distinct/disjoint. Development19.686674worker-seconds, total allowance24.906960<900. Forecast.090118core-hours/~1minute8workers, proposed NEW1hour cap and public results/full raw retention pending PI;195035136byte array bound. No scientific production/old rescore or particle/full-lattice stochastic/native-vacuum/entanglement claim; final production/pooling/incomplete audit registration remains after approval.


MG080 N1 APPROVAL/REGISTRATION: Richard approved all model/criteria/1hour/public retention items.32canonical jobs/131072paths plus32warmups,64candidate streams unchanged; 6analytic bounds D+sqrt(log24000/(2*32768))<=.03 and28two-sample bounds D+2sqrt(log224000/(2*16384))<=.08, combined>=99.9%nominal DKW family coverage. Formal evidence/source unchanged. Five held streaming/writer replays exact; full131072 passing, deliberate valid-accounting distribution failure, missing/prefix127104 and empty algebraic pipelines/independent audits PASS,18.769082seconds<300. No production seed drawn during qualification.3600worker-second cap,180/job,8threads1/.25poll,forecast324.424602seconds,2GBfree. Frozen registration/code/criteria/seeds/formal/reference hashes published before samples; no old rescore or broad particle/vacuum/quantum claim.


MG081 N1 RESULT: QUALIFIED WITHIN REGISTERED N1 SCOPE,32/32jobs,131072scientific paths+32excluded recorded warmups.6/6analytic bounds max.014669808<.03;28/28tiny two-implementation population bounds max.048666983<.08. Formal8grouped checks/source preserved. Separate64raw accounting/graph/seed/joint-generator/sortedCDF/statistical/owned-charge audit PASS. No failure/timeout/cap/admission pause/replacement/criteria amendment.105.096124worker-seconds=.029193core-hours<1cap,18.769seconds elapsed through audit. Full raw/config/counters/probability/audit archives prepared for approved public publication/hash verification. Small endpoint/formal/linear scope only, no full nonlinear particle/long-lattice/vacuum/energy/quantum claim. Old MG10 records unchanged; PI closure pending, no next production automatic.


MG082 PI CLOSURE (2026-11-08): Richard: “i sign off on MG11-N1”. MG11-N1 CLOSED as QUALIFIED WITHIN REGISTERED N1 SCOPE.32/32jobs,131072scientific paths+32excluded recorded warmups,6/6analytic and28/28two-sample criteria supported,64raw separate audit and full published hash-verified archive accepted. Formal/linear/small-system endpoint limits remain unchanged; no nonlinear full-lattice particle, vacuum/energy/entanglement claim. Source/data/registered criteria/earlier outcomes unchanged, no next scientific production authorized.


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

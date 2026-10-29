# Gates

Each gate's prediction, pass criterion and null control are filled in and
committed **before** the run that tests it. Status: PENDING / PASS / FAIL / REVISED.

| Gate | Phase | Question | Status |
|---|---|---|---|
| MG1 | 1 | Continuum PDE reproduces a published stationary Model G soliton | PASS by PI exception (registered FAIL on 1D stationarity preserved) |
| MG2 | 2 | Stochastic token rules converge to the same-lattice mean-field solution as Ω → ∞ | **FAIL** as registered (C3, C4); convergence criteria C1, C6 PASS — see MG2 result |
| MG3 | 2 | Particles nucleate spontaneously from intrinsic fluctuations at finite Ω | **FAIL** as registered (C2); C1, C3, C4, V5 PASS — see MG3 result |
| MG4 | 3 | Particles survive irregular connectivity; graph Laplacian ≈ ∇² at soliton scale | PENDING |
| MG5 | 4 | Behavior of solitons on an evolving hypergraph (persist / drift / dissolve) | PENDING |
| MG6 | 5 | G-coupled growth produces a volume-growth deficit; falloff and mass scaling | PENDING |
| MG7 | 6 | Multiway branch divergence and foliation dependence | PENDING |

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
tests `tests/test_mg3.py` (15 pass), `tests/test_track.py`.
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
| C3 (activated scaling) | PASS | k̂ = 4.05e−6, 4.63e−5, 1.58e−4 at Ω = 150, 120, 100 (5, 10, 11 events). Weighted fit S = 0.0726 per unit Ω, 95% CI [0.0516, 0.0936]; R² = 0.993. With 3 cells the t-based sensitivity interval (1 df) includes 0 |
| C4 (distance from the edge, Ω = 120) | PASS | k̂ = 4.63e−5, 1.52e−4, 2.73e−3 at b = 29, 29.16, 29.3 (monotone); k(29.3)/k(29) = 58.9, 95% CI [10.4, 229]. Most b = 29.3 runs (29/32) nucleated during burn-in |
| V5 (τ/2, Ω = 120) | PASS | ratio 0.883, 95% CI [0.365, 2.221] (14 vs 10 events) |
| N1 (b = 28.5, Ω = 100) | no confirmed event | 16 runs to T = 500, 0 persistent detections |
| N2 | PASS | mean map stays homogeneous (exact start 0; 10⁻³ perturbation decays to 1e−11) |

Independent recomputation from the raw run files (rates, the C3 fit, the C4 interval, the C2(a) counts) agrees.

Interpretation (Claude; for review). Spontaneous nucleation from intrinsic noise is **observed and confirmed** (C1).
Its rate is consistent with an activated law ln k ≈ c − S·Ω over Ω = 100–150 (C3). It rises toward the window interior
(C4). It is insensitive to halving τ (V5). The FAIL is C2(a): about 30% of persistent detections (tracks lasting ≥ 20
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

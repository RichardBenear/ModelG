# Gates

Each gate's prediction, pass criterion and null control are filled in and
committed **before** the run that tests it. Status: PENDING / PASS / FAIL / REVISED.

| Gate | Phase | Question | Status |
|---|---|---|---|
| MG1 | 1 | Continuum PDE reproduces a published stationary Model G soliton | PASS by PI exception (registered FAIL on 1D stationarity preserved) |
| MG2 | 2 | Stochastic token rules converge to the same-lattice mean-field solution as Ω → ∞ | PROPOSED (awaiting PI sign-off) |
| MG3 | 2 | Particles nucleate spontaneously from intrinsic fluctuations at finite Ω | PROPOSED (awaiting PI sign-off) |
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
Protocol: `spec/phase2_protocol_v1.0.md` (lattice, algorithm, measurement definitions, budget).
1D periodic ring, h = 0.2, L = 100. Rules v1.1. Tau-leaping with τ = 1.333e-3 (τ/2 if V2 fails).
Validations V1 (tau-leap vs exact SSA) and V2 (τ-halving) must pass before any gate run.
Operating points: **b = 29** (published) and **b_mid = 29.16** (centre of the h = 0.2 single-soliton
window [28.929, 29.397]). Deterministic references: `results/MG003/`.

### MG2 — Ω → ∞ convergence of the token rules to the lattice mean field
Status: PROPOSED 2026-10-27 (MG003).

Setup: start each run from the lattice-ODE soliton at that b (counts = round(Ω·(S₀ + φ*))), T = 100,
discard t < 10, record every 0.5. Ω ∈ {10, 30, 100, 300, 1000, 3000, 10000}, 6 runs per Ω per b.

Prediction: the aligned mean profile converges to the same-lattice ODE soliton with bias O(1/Ω).
Far-field fluctuations match the linear-noise approximation (LNA).

Pass criteria:
- C1: E(Ω = 10⁴) < 1% at both b.
- C2: bias scaling. A fit E_b ∝ Ω^(−α) over Ω ≥ 100 gives α ∈ [0.7, 1.3] at both b. If E_b ≤ σ_E at
  ≥ 2 of the Ω values ≥ 100, leaving < 3 usable points, C2 is INCONCLUSIVE (reported, not a fail).
- C3: far-field Ω·Var(φ_S) for S = G, X, Y within 10% of the LNA at every Ω ≥ 100.
  LNA targets: b = 29: 15.59 / 112.7 / 2.73; b = 29.16: 15.59 / 121.9 / 2.84.
- C4: persistence at b_mid. For Ω ≥ 100, every run keeps exactly one particle for all t ∈ [0, 100]
  (no decay, no replication). At b = 29, persistence is reported as a measurement (noise-induced decay
  near the window edge), not a criterion.
- C5: seeded nucleation (b = 29, homogeneous start, R_seed rule, detection by t = 60, 20 runs per cell).
  At seed scale 1.5 (1.58× the lattice threshold 0.949), P(nucleate) ≥ 0.9 at each of Ω = 300, 1000,
  3000. At seed scale 0.5, P ≤ 0.1 at Ω = 3000. Exploratory: P(Ω) at scale 1.0.
- C6: lattice independence (secondary). On the h = 0.1 lattice, b = 29, Ω = 10⁴, 4 runs, E < 1% against
  the h = 0.1 lattice ODE.

Gate passes if C1, C3, C4, C5 and C6 pass and C2 is not FAIL.

Null controls:
- N1: homogeneous start, b = 29, Ω = 10⁴, 4 runs, T = 100 → no particle detected.
- N2 (mis-specified reference; tests the metric's discriminating power): E of the b = 29 runs computed
  against the b_mid reference must plateau near 9.46% (the deterministic difference) and fail C1.

Exploratory (reported): particle diffusion coefficient D_p(Ω) (expected ∝ 1/Ω); noise-induced decay
events at b = 29.

### MG3 — Spontaneous nucleation from intrinsic noise
Status: PROPOSED 2026-10-27 (MG003).

Setup: b = 29, homogeneous start (independent Poisson counts at the steady state), no seed, T = 500.
Particle detection as defined in the protocol §5 (smoothed φ_X < −4.51, persisting ≥ 20).
Rate = first-nucleation events per unit length per unit time. This is a censored-exponential maximum-
likelihood estimate over runs; only the first event per run is used, since later events are seeded by
the first particle's tails.

Pre-registered adaptive design:
- Stage A (scan): Ω descending through {300, 100, 50, 30, 20, 15, 10, 7, 5, 3, 2, 1}, 4 runs each.
  Stop after the first Ω at which all 4 runs nucleate before t = 50.
- Stage B (rates): every Ω from the grid at which Stage A had 1–3 of 4 runs nucleate. If fewer than 3
  such values, add geometric midpoints between the bracketing grid values until there are 3.
  16 runs each, T = 500.

Prediction: nucleation is an activated (large-deviation) process with ln k ≈ c − S·Ω.

Pass criteria:
- C1: spontaneous nucleation exists: ≥ 10 first-nucleation events at some Ω (Stage A + B).
- C2: what nucleates is a Model G soliton. The aligned mean profile 30 time units after detection
  (runs where it persists) has E ≤ max(10%, 2·E_MG2(Ω)) against the lattice-ODE soliton, and λ is
  within 5% of the lattice value.
- C3: activated scaling. A weighted fit of ln k vs Ω over the Stage B values (each with ≥ 5 events)
  gives slope −S with the 95% CI excluding 0 and R² ≥ 0.9.
- C4: distance from the window edge. At the median Stage-B Ω, rates at b ∈ {29.0, 29.16, 29.3}
  (16 runs each) increase monotonically, and k(29.3)/k(29.0) > 1 at 95% confidence.

Gate passes if C1–C4 pass.

Null control:
- N1: b = 28.5 (below the lattice b_s, where no stationary soliton exists), smallest Stage-B Ω, 16 runs,
  T = 500 → zero persistent particles. This checks the detector at comparable noise.

Interpretation of a FAIL (pre-stated):
- No events at any Ω ≥ 1 → the barrier is too high for intrinsic noise at the lattice scale; Ω_len ≤ 5
  needed.
- Events that aren't solitons (C2 fail) → noise drives a different structure (e.g. noisy Turing
  patches).
- Non-activated scaling → nucleation isn't a single-barrier process here.

Exploratory (reported): annihilation rate (soliton → homogeneous) at Stage B Ω; structure factor S(k)
before first nucleation vs the LNA; multi-particle growth after the first event.

Configs (both gates): `configs/MG2_*.json`, `configs/MG3_*.json` (written at sign-off).
Result:
Replication (ChatGPT):
PI sign-off (criteria):
PI sign-off (result):

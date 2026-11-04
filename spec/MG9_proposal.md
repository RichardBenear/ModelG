# MG9 — Local onset type of the Model G Turing instability, and front arrest against the particle window. Design proposal, DRAFT 2

*Stage MG026a. Claude, 2026-11-04. Status: **DRAFT 2, not registered.** Draft 1 is kept as `spec/MG9_proposal_v1.md`
(superseded). This draft answers ChatGPT's review MG026 (`review/MG026_chatgpt.md`; response `review/MG026_response.md`).
All seven items are accepted.*

*PI decisions (Richard Benear, 2026-11-04 16:23 UTC, project thread): "I agree with the PI decision suggestions it
made", i.e. ChatGPT's MG026 recommendations on scope, P2 tolerance, detector, budget and the stochastic follow-up. They
are written in below and listed in §10.*

*Computation disclosed (item 7). Before draft 1, Claude timed one existing reference integration: the b = 29 MG003
single-soliton reference (`lattice_ode_1d(P17, N, h, 400)`, default seed, N = 500 and N = 1000), about 8 s and 6 s in the
cloud container. Since draft 1, Claude re-ran ChatGPT's static script (`review/MG026_replication/design_checks.py`), with
identical output. No MG9 nonlinear integration, cubic-coefficient calculation, continuation or front run has been
made.*

*Naming: MG8 stays reserved for the multiway/foliation gate. A stochastic barrier study, if any, is a separate MG9b.*

## 1. Question

> **(P1) Is the primary Turing onset of Model G locally subcritical or supercritical as b crosses b_c?
> (P2) Over what range of b do fronts between a periodic pattern and the uniform state arrest, and how does that range
> compare with the measured single-particle window?**

**Why it matters.** The single particle exists only in a narrow band of b just below the Turing threshold b_c. In the
standard picture of a subcritical instability, a patterned state coexists with the stable uniform state below
threshold, and fronts between the two arrest over a finite range of the control parameter. Inside that range, localized
states can be stable. If Model G shows both, b plays the role of the quadratic coefficient r(b) in a Landau normal form,
and the particle window is consistent with a pinning interpretation.

**What each part can and cannot show (item 1).**

- P1 decides only the **local onset type** of the selected primary branch, from the cubic coefficient under the
  normalization of §5. It does not by itself establish a stable saturated pattern, a fold, a pinning interval or
  snaking: with L₃ > 0 the cubic equation has no saturating term, and the other state needs higher-order terms or the
  full system. Either sign is a valid finding.
- A validated **supercritical** result (L₃ < 0) does not contradict the patterns already known below b_c. They could lie
  on a disconnected branch, or on a branch that turns back farther from onset.
- P2 is a comparison of **measured finite-time front arrest** with the **operational** single-particle window. Agreement
  supports consistency with a pinning interpretation. It does not demonstrate homoclinic snaking, which needs the
  organization of localized branches and their successive folds; that is not computed here. Burke & Knobloch (2006)
  separate the small-amplitude normal form from the large-amplitude localized states in the same way. A negative P2
  leaves the local reduction (P1) intact.
- No global variational structure has been established for Model G. A potential for a reduced real amplitude does not
  supply physical energy units, and does not prove a global potential for the three coupled fields (Burke & Dawes 2011
  on non-variational extensions). MG9 tests no Higgs identification and converts nothing to vacuum energy.

## 2. What is already known

These are disclosed so the items are not mistaken for blind tests. All on the 1D lattice, reference parameters (eq. 17)
with b varied.

| Quantity | Value | Provenance |
|---|---|---|
| Turing threshold, continuous k (continuum and h = 0.2, 0.1 lattices) | 30.04407748 | ChatGPT MG026 static check; the lattice relation reaches the same optimum Laplacian eigenvalue |
| Critical k, continuous | 2.09928 (continuum), 2.11501 (h = 0.2), 2.10315 (h = 0.1) | same |
| Threshold, L = 100 box | 30.04577 (h = 0.2, mode 34, k = 2.13628); 30.04757 (h = 0.1, mode 33) | same; 30.0458 in MG003 is the h = 0.2 box value |
| Lower window edge b_s | h = 0.2: 28.9286–28.9287; h = 0.1: 28.9782 | `existence_window.py`: seed-free continuation from the b = 29 soliton, "exists" if φ_X(0) < −1 at t = 300 |
| Upper window edge b_r | h = 0.2: 29.3965–29.3969; h = 0.1: 29.4117–29.4121 | `replication_edge.py`: standard seed, "single" if exactly one minimum with φ_X < −4 at t = 400 |
| Deterministic seed threshold | 0.949 (b = 29), 0.887 (29.16), 0.763 (29.5) | MG003 references |
| Above b_r the seeded particle self-replicates until the domain fills | | MG003 |

Consequences:

1. A space-filling pattern is already known for b_r < b < b_c. P1's outcome is therefore not blind in practice, but its
   local sign is a separate fact (§1) and either sign is reported as found.
2. The window edges are finite-time, detector-specific benchmarks. Their narrow bisection brackets do not make them
   equally precise stationary bifurcation points (item 6). §6.5 rechecks their time sensitivity before P2 uses them.
3. The seed-threshold trend is already known over 29–29.5. It is a protocol-specific proxy for access to the particle's
   basin, not a stochastic action. Its extension is deferred (PI decision, §10).

## 3. Model, numerics and frozen domains

- Equations: dimensionless Model G (eqs. 7a–7c), reference set eq. 17 (a = 14, g = 0.1, p = q = 1, s = u = w = 0,
  d_x = 1, d_y = 12), b varied. Fields are potentials φ = (G − G₀, X − X₀, Y − Y₀) about the uniform state at the
  current b. G₀ = X₀ = 140/9 do not depend on b; Y₀ = b/X₀ does.
- **Carrying a field to a new b (item 4).** Concentrations are preserved: φ_Y(new) = φ_Y(old) + Y₀(b_old) − Y₀(b_new),
  with φ_G and φ_X unchanged. This is recorded with each run.
- Integrator: `lattice.lattice_ode_1d`, BDF, rtol 1e-8, atol 1e-10, called **with seed_scale = 0 explicitly in every
  MG9 integration** (item 3). The default adds a localized forcing pulse even when y0 is given. A test checks this.
- Solver success: `sol.success` must be true. A failed integration, after the routine's own first-step retries, is
  classified UNRESOLVED, never as a dynamical outcome.
- Stationarity: max |dφ/dt| (the routine's returned residual) ≤ 1e-6 unless stated.

**Frozen domains (item 2).** Every MG9 integration uses one of these:

| ID | h | N | L | Critical mode | k | Box threshold b_c^box | Use |
|---|---|---|---|---|---|---|---|
| D02 | 0.2 | 500 | 100 | 34 | 2.136283 | 30.04577 | V2, sweep, A3 primary |
| D02L | 0.2 | 1000 | 200 | 68 | 2.136283 | 30.04577 | A3 long-domain repeats |
| D01 | 0.1 | 1000 | 100 | 34 | 2.136283 | (computed; reported) | A3 mesh comparison |

D01 holds the **physical pattern wavelength fixed** at λ = 100/34 = 2.941 (mode 34) rather than its own critical mode
(33). The mesh comparison therefore changes only h, not the pattern (item 2). λ/h = 14.7 at h = 0.2 and 29.4 at h = 0.1.

## 4. Part A1 — linear coefficient (descriptive)

σ_max(b) for b ∈ [28.5, 30.5], step 0.05, on: the continuum with continuous k; the h = 0.2 and h = 0.1 lattices with
continuous k (Laplacian eigenvalue −(4/h²) sin²(kh/2)); and the D02 and D01 boxes (discrete modes).

- **V1:** reproduce the §2 thresholds: 30.04408 for continuous k, 30.04577 for D02, 30.04757 for h = 0.1 at its own
  mode 33, each to 1e-4. A continuous-k calculation is never compared with a box value.
- Report α = dσ/db at threshold in original time t, on each problem, and σ of mode 34 on D01 as a function of b.
  ChatGPT's static values (0.1981 for D02 mode 34, 0.1906 for h = 0.1 mode 33) serve as a cross-check to 1e-3.

## 5. Part A2 — cubic coefficient and local onset type (P1)

**Normalization (frozen).** Let M(κ) = J − κ D, with J the reaction Jacobian at the uniform state, D = diag(1, 1, 12),
and κ the positive Laplacian eigenvalue (k² in the continuum, (4/h²) sin²(kh/2) on a lattice). At b = b_c on the problem
in question, φ is the right null vector of M(κ_c), normalized so its X component is 1. ψ is the left null vector,
normalized so ψ·φ = 1. The pattern is a **real cosine amplitude** in **original time**:

φ(x, t) = A(t) cos(kx) φ + A² [w₀ + w₂ cos(2kx)] + …,  dA/dt = σ(b) A + L₃ A³ + O(A⁵).

A is the amplitude of the cos(kx) part of φ_X. With a complex envelope convention (A e^{ikx} + c.c.) the coefficient would
be 4L₃. Only the real cosine form is used here.

**Formula (frozen).** The nonlinearity of Model G is the single term X²Y, entering the X equation with + and the Y
equation with −, so with e = (0, 1, −1):

Q(u, v) = e [Y₀ u_X v_X + X₀ (u_X v_Y + u_Y v_X)],  C(u) = e u_X² u_Y.

The slaved second-order fields are w₀ = −½ M(0)⁻¹ Q(φ, φ) and w₂ = −½ M(κ(2k))⁻¹ Q(φ, φ), and

**L₃ = ψ · [2 Q(φ, w₀) + Q(φ, w₂) + ¾ C(φ)].**

L₃ > 0 means locally **SUBCRITICAL**, L₃ < 0 locally **SUPERCRITICAL**.

**Problems.** (i) D02 at mode 34, b = 30.04577: this is primary, because V2 can check it directly. (ii) The h = 0.2
lattice with continuous k_c. (iii) The continuum with continuous k_c. (iv) D01 at mode 34. All are reported. The primary
sign is (i). If the signs disagree across (i)–(iv), P1 is reported as MIXED with all values.

**Code check (kept from draft 1).** The same routine, applied to the Swift–Hohenberg equation
u_t = r u − (1 + ∂²)² u + c₂ u² − u³, must reproduce the classical onset-type boundary c₂² = 27/38 (subcritical above)
to 1e-6. This check is normalization-independent. V2 checks the normalization.

**V2: fitted cubic coefficient on D02 (item 3).**

- Initial field: A₀ cos(kx) φ + A₀² [w₀ + w₂ cos(2kx)] at b = 30.04577 (the full critical eigenmode with its slaved
  mean and second harmonic), seed_scale = 0. A₀ ∈ {0.05, 0.1, 0.2}.
- Measurement: A(t) = ψ · û₃₄(t), where û₃₄ is the cosine coefficient of mode 34 of the potential fields (2/N Σ φ_j
  cos(k x_j), per component), sampled every 1 time unit. t = 0–20 is discarded as transient.
- Fit: dA/dt = σ₃₄ A + L₃ A³, with σ₃₄ fixed at the computed linear growth rate of mode 34 on D02 at that b (not
  fitted). It is done as a least-squares fit of the integrated form over t ∈ [20, t_end]. t_end is the first time A
  reaches 2A₀ or A₀/2, or t = 3000.
- Identifiability conditions, all required: |L₃| A₀² ≥ 100 |σ₃₄| at the chosen b; the relative change in 1/A² over the
  fit window is at least 25%; and a repeat of the A₀ = 0.1 run at rtol 1e-10, atol 1e-12 changes the fitted L₃ by at
  most 2%. If the 25% condition fails for every A₀, the set is doubled once to {0.1, 0.2, 0.4} (pre-declared).
- **V2 passes** if the fitted L₃ from every A₀ in the final set has the same sign as the derived (i), each lies within
  10% of it, and the fits converge with A₀ (the spread across A₀ is at most 10% of the mean).
- If V2 fails or cannot be identified, P1 is **NOT EVALUATED**, and the derived values are reported as unvalidated.

**P1 outcome.** SUBCRITICAL or SUPERCRITICAL (validated, from (i)); MIXED; or NOT EVALUATED. No outcome is called a
contradiction or an error (item 1).

**Pattern persistence sweep (item 4; descriptive).** This produces the periodic states that A3 needs, and a finite-time
persistence bracket. It does not locate a fold.

- Start: on D02 at b = 30.5, the mode-34 field with A₀ = 0.5, integrated to t = 1000. It must be stationary with
  34 minima (prominence rule of §6.2). Otherwise the sweep start is UNRESOLVED and reported, and A3 cannot run.
- Step b down by 0.02 to 28.50 (101 steps), carrying the field by the concentration rule (§3), and integrating each step
  up to t_settle = 1000. If the residual is still > 1e-6 at 1000, the step is extended once to 3000.
- Each step is classified: **PERSISTS** (stationary, 34 cells); **CHANGED** (stationary, a different cell count, e.g. a
  phase slip or wavelength change; recorded and the sweep continues from that state, but A3 uses only 34-cell states);
  **LOST** (max |φ_X| < 0.1; the sweep stops); **UNRESOLVED** (not stationary at 3000, or solver failure; the sweep stops).
- Reported: the bracket [last PERSISTS, first LOST or UNRESOLVED] as a **finite-time persistence bracket** under this
  protocol, compared with b_s descriptively. A true fold would need steady-state continuation with a phase condition and
  a stability check. That is deferred.
- Each stored 34-cell state is the "periodic state at b" for A3. A3 grid points at 28.50, 28.60, …, 30.00 coincide with
  sweep steps (16 points). The other 15 (28.55, 28.65, …) use the sweep state at b + 0.01, carried to b and settled to
  t = 1000 under the same rules.

## 6. Part A3 — front arrest against the particle window (P2, primary)

### 6.1 Set-up
- On D02, the periodic state at b (34 cells; §5 sweep), translated by a whole number of nodes so that a cell minimum sits
  at the node nearest x = 0 (the **phase convention**), and multiplied by a window
  w(x) = ½ [tanh((x + 25)/ℓ) − tanh((x − 25)/ℓ)] with ℓ = λ/2. The uniform state (φ = 0) is used elsewhere. The patch is
  50 units = 17 cells.
- If no 34-cell periodic state exists at b (the sweep was LOST or UNRESOLVED above b, or CHANGED there), the point is
  classified **NO-BRANCH**. It counts as not arrested; an interval edge set by NO-BRANCH is labelled BRANCH-LIMITED.
- Integrate to T = 600, with an extension rule in §6.3.

### 6.2 Observables (item 5)
- **Primary: interface positions.** Envelope E(x, t) = √2 × the local RMS of φ_X about its local mean, over a centred
  15-node window (periodic). E_pat(b) is the spatial mean of E for the periodic state at b. The patterned region is where
  E ≥ ½ E_pat. The left and right interfaces x_L(t), x_R(t) are the outermost crossings of ½ E_pat, located by linear
  interpolation and unwrapped on the periodic domain. They are sampled every 5 time units and tracked separately; no
  mirror symmetry is assumed. Outward displacement is positive for both.
- **Secondary: cell count n(t)**, the number of local minima of φ_X with prominence ≥ ½ E_pat and separation ≥ λ/2.
- The pattern **fills** the box if no point has E < ½ E_pat. It is **lost** if no point has E ≥ ½ E_pat.

### 6.3 Classification (frozen before any run)
Let Δ_L, Δ_R be the interface displacements over [100, T]; R_L, R_R their ranges (max − min) over [100, T]; v_L, v_R
their least-squares velocities over the last 200 units; r_T the final residual.

- **INVADES:** Δ_L ≥ λ and Δ_R ≥ λ, or the pattern fills the box.
- **RETREATS:** Δ_L ≤ −λ and Δ_R ≤ −λ, or the pattern is lost.
- **ARRESTED** (finite-time arrest): |Δ_L|, |Δ_R| ≤ λ/2; R_L, R_R ≤ λ; |v_L|, |v_R| ≤ λ/2000 per unit time;
  r_T ≤ 1e-5; and n changes by at most 1 over [100, T].
- **Extension:** if after T = 600 the point is not INVADES or RETREATS and not ARRESTED with margin (|Δ| ≤ λ/4 and
  |v| ≤ λ/4000 for both interfaces), it is extended once to T = 1800 and reclassified with the same rules.
- **UNRESOLVED:** everything else after extension. That includes asymmetric motion (one interface moving and one not),
  motion followed by reversal (range above λ), oscillation, a cell-count change of 2 or more with interfaces fixed
  (flagged INTERNAL), solver failure, and a missing interface (fronts meeting) before t = 100.
- **Robustness repeats.** At the two outermost ARRESTED grid points (lowest and highest b), run (a) the same set-up on
  D02L with a 100-unit patch, and (b) D02 with the patch centred on a cell maximum instead of a minimum (half-period
  phase shift). An outermost ARRESTED point stands only if both repeats are also ARRESTED. Otherwise it becomes
  UNRESOLVED.
- Interface velocities v(b) at non-arrested points are reported descriptively. The cubic amplitude equation alone does
  not supply the invaded state, so no front-speed prediction is made.

### 6.4 Grid and refinement
- **Grid (h = 0.2):** b = 28.50 to 30.00, step 0.05 (31 points).
- **Transition brackets:** between adjacent grid points where the classification changes from ARRESTED to INVADES or
  RETREATS (or NO-BRANCH). Each of the two brackets is bisected twice (to 0.0125), each midpoint classified by the same
  rules (4 runs). A bracket whose refinement point is UNRESOLVED is kept at its unrefined width.
- **Mesh comparison (D01, h = 0.1):** for each refined h = 0.2 bracket [b₁, b₂], the points b₁ − 0.05, b₁, b₂, b₂ + 0.05
  (8 points). The D01 periodic state at each b is the D02 state interpolated by Fourier interpolation (mode-34 phase
  preserved) and settled on D01 to t = 1000 (stationary, 34 cells, or NO-BRANCH). If the 8 points do not bracket a
  transition at h = 0.1, up to 2 further points per edge are added, stepping outward by 0.05. h = 0.1 is compared with
  its own benchmarks (b_s = 28.978, b_r = 29.412).
- **MESH-SENSITIVE flag:** set if, between h = 0.2 and h = 0.1, the arrest-interval width changes by more than half, or
  either edge offset from its own window benchmark changes by more than 0.10. It is a flag, not an attribution:
  wavelength resolution, phase, front-speed resolution and truncation error all change with h, as well as lattice
  locking (item 6).

### 6.5 Benchmark recheck (item 6)
Before P2 is scored, the h = 0.2 window edges are recomputed by the same MG003 rules with longer horizons: b_s with
t = 600 (instead of 300) and b_r with t = 800 (instead of 400), 8 bisection steps each, with seed_scale = 0 for the b_s
continuation as in MG003. The b_r seeded runs keep MG003's standard seed, which is a defined part of that benchmark. P2
uses the **rechecked** values, and both are reported. If either rechecked edge moves by more than 0.05, P2 is labelled
**BENCHMARK-SENSITIVE**, in addition to its outcome.

### 6.6 P2 outcomes (h = 0.2; tolerance 0.10, PI decision)
Let [l₁, l₂] and [u₁, u₂] be the refined lower and upper transition brackets, and b_s*, b_r* the rechecked benchmarks.

| Label | Criterion |
|---|---|
| **COINCIDES** | both brackets are resolved, and each lies entirely within ±0.10 of its benchmark: [l₁, l₂] ⊂ [b_s* − 0.10, b_s* + 0.10] and [u₁, u₂] ⊂ [b_r* − 0.10, b_r* + 0.10] |
| **OVERLAPS** | both brackets are resolved, the arrest interval intersects [b_s*, b_r*], and COINCIDES fails |
| **NO ARREST OBSERVED** | every grid point is classified (none UNRESOLVED), and none is ARRESTED. This is "not observed on this grid at this resolution", not an empty pinning interval |
| **OUTSIDE** | both brackets are resolved, and the arrest interval does not intersect [b_s*, b_r*] |
| **UNRESOLVED** | an edge lies outside the scan, a grid point adjacent to an edge is UNRESOLVED, or there are several disconnected arrest intervals |

**ANOMALOUS** is added to any label if the ordering is not RETREATS/NO-BRANCH below, ARRESTED, INVADES above.
**BRANCH-LIMITED**, **BENCHMARK-SENSITIVE** and **MESH-SENSITIVE** are added where they apply. The wording of any
positive result is "front arrest consistent with a pinning interpretation of the window", never "snaking demonstrated".

## 7. Not registered here: stochastic barrier S(b)

This remains a possible separate MG9b, with its own design. The deterministic seed threshold is only a
protocol-specific basin-access proxy. A stochastic action depends on the noise law and the escape path as well as the
deterministic dynamics (PI decision, §10).

## 8. Descriptive outputs (no criteria)

- **D1:** α (A1, all problems); L₃ on problems (i)–(iv); the sweep's persistence bracket; E_pat(b); and the interface
  velocities v(b).
- **D2:** one Landau-style figure: σ(b), the persistence bracket, the arrest interval and the particle window on one
  b-axis.
- Deferred from draft 1 (PI decision): the seed-threshold extension to new b values.

## 9. Runs, cost and records (item 7)

| Item | Integrations |
|---|---|
| V2 (3 amplitudes + 1 tolerance repeat) | 4 |
| Sweep start + 101 steps (extensions to 3000 are part of the same integration count) | 102 |
| A3 grid (h = 0.2) | 31 |
| Settling for A3 points between sweep steps | 15 |
| Extensions to T = 1800 | part of the same integrations |
| Robustness repeats | 4 |
| Bracket refinement | 4 |
| D01: settle + front, 8 points | 16 |
| Benchmark recheck (2 × 8 bisection + 1 reference) | 17 |
| **Planned total** | **193** |
| Reserved: V2 doubling (3), D01 extra points (8 points, 16 integrations), refinement settles (4), re-runs after solver failure (4) | 27 |
| **Run cap** | **220** |

Refinement midpoints that need a settle use the reserve. If the cap would be exceeded, the run stops and the PI decides.

- Cost: about 8 s per 400 time units for a single-particle field (the disclosed timing, header). Pattern fields with many active cells are
  slower, so 10× is assumed: the sweep (102 × up to 1000 units) about 1.5 core-h; A3 (about 55 runs × 600–1800 units)
  about 1.5 core-h; the rest under 0.5 core-h. **About 3–4 core-h. Cutoff 8 core-h.** It runs in the cloud container in
  resumable stages, each under the 2-hour background limit. No PC time.
- Records, committed: per-run JSON (b, domain ID, initial-state provenance, carry rule, classification, x_L(t), x_R(t),
  n(t), E_pat, residuals, solver status); the A1/A2 numbers; the V2 fits; and the stored periodic states (34-cell fields
  on N = 500, about 120 states × 12 kB) as one `.npz` under 5 MB, under `results/MG9/`.
- Code: `src/landau_mg9.py` (A1, A2, Swift–Hohenberg check, V2 fit), `src/run_mg9.py` (sweep, A3, refinement, D01,
  benchmark recheck; resumable), `src/evaluate_mg9.py`, and `tests/test_mg9.py`. The tests cover: seed_scale = 0
  enforced; the concentration carry rule; the Swift–Hohenberg 27/38 boundary; the V1 thresholds; the interface tracker
  on synthetic fronts (moving, fixed, asymmetric, reversing); and the classification table on synthetic tracks.
- **Pre-run implementation specification** (CLAUDE.md rule 1): anything left open here is fixed in GATES.md before the
  first run.

## 10. Decisions (adopted from ChatGPT MG026 by the PI, 2026-11-04 16:23 UTC)

| Decision | Adopted |
|---|---|
| Scope | Deterministic Part A only. P1 is the local onset type; P2 is the relation between measured front arrest and the operational single-particle window. Any snaking claim is deferred to branch continuation. |
| P2 tolerance | 0.10, provisionally, with explicit edge brackets, targeted refinement and UNRESOLVED (§6.4–6.6). |
| Detector | A direct envelope/interface measurement for motion (primary); unsmoothed minima for cell counts (secondary). Definitions, thresholds and phase convention frozen in §6. |
| Budget | P1 and P2 have priority. Every integration is counted, with validation runs reserved (§9). D2 of draft 1 (seed-threshold extension) is deferred. |
| Stochastic follow-up | MG9b stays separate (§7). |

Open for the PI: whether to send draft 2 to ChatGPT for a round 2 before registration (recommended: the frozen
classification rules in §6.3 and the V2 identifiability conditions are new and worth a second reading), or to register
it as is.

## 11. What each outcome would mean

- **SUBCRITICAL and COINCIDES:** the onset is locally subcritical, and fronts arrest over the same b-range in which a
  single particle survives. That is consistent with b acting as the control parameter of a subcritical Landau-type
  instability, with the window as a pinning region. This is the strongest statement MG9 can support. It concerns the
  normal form and finite-time front dynamics, not a global potential, and not a Higgs identification.
- **SUBCRITICAL and OVERLAPS:** the picture is right in kind, but the single particle's window is not simply the
  front-arrest range. Its edges would need their own explanation, for example the single-dip branch's own folds.
- **SUPERCRITICAL (validated), with any P2:** the primary onset is locally supercritical. The patterns known below b_c
  would lie on a branch not connected to onset in the simple way, and any front arrest is reported as found.
- **NO ARREST OBSERVED or OUTSIDE:** the window is not explained by front arrest at this resolution. The Landau reading
  of b reduces to "b controls the stability of the uniform state", and P1 stands on its own.
- **UNRESOLVED, or any flag:** reported as such, with what would resolve it (longer runs, continuation, finer meshes).

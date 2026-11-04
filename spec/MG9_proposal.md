# MG9 — Local onset type of the Model G Turing instability, and front arrest against the particle window. Design proposal, DRAFT 3 (REGISTERED)

*Stage MG026b. Claude, 2026-11-04. Status: **REGISTERED 2026-11-04 (MG026c), draft 3, content frozen at `7810b8b`; PI decision card 16:40 UTC.** No MG9 run of any kind preceded registration. Drafts 1 and 2 are kept as
`spec/MG9_proposal_v1.md` and `spec/MG9_proposal_v2.md` (superseded). This draft answers ChatGPT's round 2
(`review/MG026_chatgpt_r2.md`; response `review/MG026_response_r2.md`). All three main items and all seven
implementation points are accepted and written into the text below, not left to the pre-run specification.*

*PI decisions (Richard Benear, 2026-11-04 16:23 UTC, project thread): adopt ChatGPT MG026's recommendations on scope,
P2 tolerance, detector, budget and the stochastic follow-up (§10).*

*Computation disclosed:*
- *Before draft 1, Claude timed one existing reference integration: the b = 29 MG003 single-soliton reference
  (`lattice_ode_1d(P17, N, h, 400)`, default seed, N = 500 and 1000), about 8 s and 6 s in the cloud container.*
- *Claude re-ran ChatGPT's static scripts (`review/MG026_replication/design_checks.py`, `draft2_checks.py`), with
  identical output. These are synthetic or linear checks only.*
- *No MG9 nonlinear integration, cubic-coefficient calculation, continuation or front run has been made.*

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

**What each part can and cannot show.**

- P1 decides only the **local onset type** of the selected primary branch, from the cubic coefficient under the
  normalization of §5. It does not by itself establish a stable saturated pattern, a fold, a pinning interval or
  snaking: with L₃ > 0 the cubic equation has no saturating term, and the other state needs higher-order terms or the
  full system. Either sign is a valid finding.
- A validated **supercritical** result (L₃ < 0) does not contradict the patterns already known below b_c. They could lie
  on a disconnected branch, or on a branch that turns back farther from onset.
- P2 is a comparison of **measured finite-time front arrest** with the **operational** single-particle window. Agreement
  supports consistency with a pinning interpretation. It does not demonstrate homoclinic snaking, which needs the
  organization of localized branches and their successive folds; that is not computed here (cf. Burke & Knobloch 2006).
  A negative P2 leaves the local reduction (P1) intact.
- No global variational structure has been established for Model G. A potential for a reduced real amplitude does not
  supply physical energy units, and does not prove a global potential for the three coupled fields (cf. Burke & Dawes
  2011). MG9 tests no Higgs identification and converts nothing to vacuum energy.

## 2. What is already known

These are disclosed so the items are not mistaken for blind tests. All on the 1D lattice, reference parameters (eq. 17)
with b varied.

| Quantity | Value | Provenance |
|---|---|---|
| Turing threshold, continuous k (continuum and h = 0.2, 0.1 lattices) | 30.04407748 | ChatGPT MG026 static check; the lattice relation reaches the same optimum Laplacian eigenvalue |
| Critical k, continuous | 2.09928 (continuum), 2.11501 (h = 0.2), 2.10315 (h = 0.1) | same |
| Box onset, L = 100 | 30.04577 (h = 0.2, mode 34); 30.04757 (h = 0.1, mode 33) | same; 30.0458 in MG003 is the h = 0.2 box value |
| Box onset, L = 200, h = 0.2 | 30.044468 (mode 67); selected mode 68 has threshold 30.045770 | ChatGPT MG026 r2 |
| Lower window edge b_s | h = 0.2: 28.9286–28.9287; h = 0.1: 28.9782 | `existence_window.py`: seed-free continuation from the b = 29 soliton (fields carried as potentials), "exists" if φ_X(0) < −1 at t = 300 |
| Upper window edge b_r | h = 0.2: 29.3965–29.3969; h = 0.1: 29.4117–29.4121 | `replication_edge.py`: standard seed, "single" if exactly one minimum with φ_X < −4 at t = 400 |
| Deterministic seed threshold | 0.949 (b = 29), 0.887 (29.16), 0.763 (29.5) | MG003 references |
| Above b_r the seeded particle self-replicates until the domain fills | | MG003 |

Consequences:

1. A space-filling pattern is already known for b_r < b < b_c. P1's outcome is therefore not blind in practice, but its
   local sign is a separate fact (§1) and either sign is reported as found.
2. The window edges are finite-time, detector-specific benchmarks. Their narrow bisection brackets do not make them
   equally precise stationary bifurcation points. §6.5 rechecks their time sensitivity before P2 uses them.
3. The seed-threshold trend is already known over 29–29.5. It is a protocol-specific proxy for access to the particle's
   basin, not a stochastic action. Its extension is deferred (PI decision, §10).

## 3. Model, numerics and frozen domains

- Equations: dimensionless Model G (eqs. 7a–7c), reference set eq. 17 (a = 14, g = 0.1, p = q = 1, s = u = w = 0,
  d_x = 1, d_y = 12), b varied. Fields are potentials φ = (G − G₀, X − X₀, Y − Y₀) about the uniform state at the
  current b. G₀ = X₀ = 140/9 do not depend on b; Y₀ = b/X₀ does.
- **Carrying a field to a new b.** Concentrations are preserved: φ_Y(new) = φ_Y(old) + Y₀(b_old) − Y₀(b_new), with
  φ_G and φ_X unchanged. This is recorded with each run. **Exception:** the b_s benchmark recheck (§6.5) reproduces
  MG003's convention exactly (potentials carried unchanged). That keeps it a pure time-horizon recheck.
- **Seed forcing.** Every MG9 integration sets seed_scale = 0 explicitly. The routine's default adds a localized
  forcing pulse even when y0 is given. **The only exception** is the b_r benchmark recheck (§6.5), which keeps MG003's
  standard seed (seed_scale = 1) because the seed defines that benchmark. A test checks both settings.
- **Solver wrapper.** MG9 uses its own wrapper, `mg9_ode`, with the same right-hand side, Jacobian, BDF settings
  (rtol 1e-8, atol 1e-10, max_step 0.5) and first-step retries as `lattice.lattice_ode_1d`. It additionally returns
  `success`, the time reached and the retry used. `lattice_ode_1d` itself is unchanged, because earlier studies use it.
  A test shows that the wrapper reproduces `lattice_ode_1d` to 1e-12 on the b = 29 reference. **A run is complete only
  if success is true and the time reached equals the requested end time.** Anything else is a numerical failure,
  never a dynamical outcome, even though the routine returns the last available state.
- Stationarity: max |dφ/dt| (the routine's residual) ≤ 1e-6 unless stated.
- Thresholds are computed at full double precision (brentq, xtol 1e-12). The rounded values in tables are for reading
  only.

**Frozen domains.**

| ID | h | N | L | Selected pattern mode | k | Selected-mode threshold | Box onset (first unstable mode) | Use |
|---|---|---|---|---|---|---|---|---|
| D02 | 0.2 | 500 | 100 | 34 | 2.136283 | 30.045770 | same (mode 34) | V2, sweep, A3 primary |
| D02L | 0.2 | 1000 | 200 | 68 | 2.136283 | 30.045770 | 30.044468 (mode 67) | A3 long-domain repeats |
| D01 | 0.1 | 1000 | 100 | 34 | 2.136283 | computed in A1; reported | 30.04757 (mode 33) | A3 mesh comparison |

The physical pattern wavelength is fixed at λ = 100/34 = 2.941 on every domain. λ/h = 14.7 at h = 0.2 and 29.4 at
h = 0.1.

## 4. Part A1 — linear coefficient (descriptive)

σ_max(b) for b ∈ [28.5, 30.5], step 0.05, on: the continuum with continuous k; the h = 0.2 and h = 0.1 lattices with
continuous k (Laplacian eigenvalue −(4/h²) sin²(kh/2)); and the three domains (discrete modes).

- **V1:** reproduce, to 1e-4, the continuous-k threshold 30.04408 and each domain's selected-mode threshold and box onset
  from §3. A continuous-k calculation is never compared with a box value.
- Report α = dσ/db at threshold in original time t on each problem. ChatGPT's static values (0.1981 for D02 mode 34,
  0.1906 for h = 0.1 mode 33) serve as a cross-check to 1e-3.

## 5. Part A2 — cubic coefficient and local onset type (P1)

**Normalization (frozen).** Let M(κ) = J − κ D, with J the reaction Jacobian at the uniform state, D = diag(1, 1, 12),
and κ the positive Laplacian eigenvalue (k² in the continuum, (4/h²) sin²(kh/2) on a lattice). At b = b_c on the problem
in question, φ is the right null vector of M(κ_c), normalized so its X component is 1. ψ is the left null vector,
normalized so ψ·φ = 1. The pattern is a **real cosine amplitude** in **original time**:

φ(x, t) = A(t) cos(kx) φ + A² [w₀ + w₂ cos(2kx)] + …,  dA/dt = σ(b) A + L₃ A³ + O(A⁵).

A is the amplitude of the cos(kx) part of φ_X. With a complex-envelope convention (A e^{ikx} + c.c.) the coefficient would
be 4L₃. Only the real cosine form is used here.

**Formula (frozen; checked algebraically by ChatGPT MG026 r2).** The nonlinearity of Model G is the single term X²Y,
entering the X equation with + and the Y equation with −, so with e = (0, 1, −1):

Q(u, v) = e [Y₀ u_X v_X + X₀ (u_X v_Y + u_Y v_X)],  C(u) = e u_X² u_Y.

The slaved second-order fields are w₀ = −½ M(0)⁻¹ Q(φ, φ) and w₂ = −½ M(κ(2k))⁻¹ Q(φ, φ), and

**L₃ = ψ · [2 Q(φ, w₀) + Q(φ, w₂) + ¾ C(φ)].**

L₃ > 0 means locally **SUBCRITICAL**, L₃ < 0 locally **SUPERCRITICAL**.

**Problems.** (i) D02 at mode 34, at its selected-mode threshold: this is primary, because V2 can check it directly.
(ii) The h = 0.2 lattice with continuous k_c. (iii) The continuum with continuous k_c. (iv) D01 at mode 34. All are
reported.

**Code check.** The same routine, applied to the Swift–Hohenberg equation u_t = r u − (1 + ∂²)² u + c₂ u² − u³, must give
the cosine-amplitude coefficient 19c₂²/18 − 3/4 (ChatGPT MG026 r2) at c₂ = 0, 0.5, 1 to 1e-8, and so the boundary
c₂² = 27/38. V2 checks the Model G normalization.

**V2: fitted cubic coefficient on D02.**

- Initial field: A₀ cos(kx) φ + A₀² [w₀ + w₂ cos(2kx)] at b = the D02 mode-34 threshold (the full critical eigenmode with
  its slaved mean and second harmonic), seed_scale = 0. Amplitude set S = {0.05, 0.1, 0.2}.
- Measurement: A(t) = ψ · û₃₄(t), where û₃₄ is the cosine coefficient of mode 34 of the potential fields (2/N Σ φ_j
  cos(k x_j), per component), sampled every 1 time unit. t = 0–20 is discarded as transient.
- Integration stops at the first time A reaches 2A₀ or A₀/2, or at t = 3000. A run that stops before t = 40 (fewer than
  20 usable samples) is **unusable**.
- Fit: dA/dt = σ₃₄ A + L₃ A³, with σ₃₄ fixed at the computed linear growth rate of mode 34 on D02 at that b (not fitted).
  It is done as a least-squares fit of the integrated form over [20, stop].
- **Identifiability**, per usable run: |L₃_fit| A₀² ≥ 100 |σ₃₄|, and the relative change in 1/A² over the fit window
  is at least 25%.
- **Doubling rule (pre-declared, at most once):** if no run in S is both usable and identifiable, S becomes
  {0.1, 0.2, 0.4}.
- **Tolerance repeat:** the middle amplitude of the final S (0.1, or 0.2 after doubling) is repeated at rtol 1e-10,
  atol 1e-12. The fitted L₃ must change by at most 2%.
- **V2 passes** if at least two amplitudes of the final S are usable and identifiable; their fitted L₃ values have the
  sign of the derived (i) and lie within 10% of it; their spread is at most 10% of their mean; and the tolerance repeat
  passes.

**P1 outcome (evaluated in this order; the first that applies is the outcome).**
1. **NOT EVALUATED** if the code check or V2 fails.
2. **DEGENERATE** if |L₃(i)| < 1e-3 × (|ψ·2Q(φ,w₀)| + |ψ·Q(φ,w₂)| + |ψ·¾C(φ)|), i.e. the terms cancel to better than
   0.1%, so the sign is not robust.
3. **MIXED** if the signs of (i)–(iv) disagree. All values are reported.
4. **SUBCRITICAL** or **SUPERCRITICAL**, from the validated sign of (i).

No outcome is called a contradiction or an error.

**Pattern persistence sweep (descriptive; also prepares A3's states).** This yields a finite-time persistence bracket
under this protocol. It does not locate a fold. A true fold would need steady-state continuation with a phase condition
and a stability check, which is deferred.

- Start: on D02 at b = 30.50, the mode-34 field with A₀ = 0.5, integrated to t = 1000. It must be complete, stationary
  and have 34 cells (cell rule of §6.2). Otherwise the sweep start is PREP-UNRESOLVED, and P2 is UNRESOLVED.
- 30.50, 30.48, …, 28.50 is 101 values: the start plus 100 downward steps. At each step, the field is carried by the
  concentration rule (§3) and integrated up to t_settle = 1000, extended once (as the same integration) to 3000 if the
  residual is still above 1e-6.
- Each step is classified, in this order:
  - **PREP-UNRESOLVED:** not complete, or not stationary at 3000.
  - **PREP-LOST:** complete and max |φ_X| < 0.1.
  - **PERSISTS:** stationary with 34 cells.
  - **PREP-CHANGED:** stationary with a different cell count.
- **The sweep never stops early.**
  - After PREP-CHANGED it continues from the changed state.
  - After PREP-LOST it continues from the uniform state. Every later step is then PREP-LOST unless the uniform state
    destabilizes, which cannot happen below b_c.
  - After PREP-UNRESOLVED it continues from the last complete stationary state, carried to the next b.
- **Reported:** the full status sequence; the **primary persistence bracket**, [lowest b of the unbroken PERSISTS run
  that starts at 30.50, the next b below it]; and any later b values that regain 34 cells, listed separately as
  "regained".
- A stored 34-cell state (PERSISTS, initial run or regained) is the **prepared periodic state at b** for A3. A regained
  state is flagged as such in A3's records.

## 6. Part A3 — front arrest against the particle window (P2, primary)

### 6.1 Preparation statuses
- **Prepared state at b.** At sweep values (b = 28.50, 28.52, …) the stored state is used. At any other b (grid points
  28.55, 28.65, …; refinement midpoints; D01 comparison points), the stored state at the nearest sweep value above b is
  carried to b and settled on D02 to t = 1000. This settle is classified like a sweep step (PREP-UNRESOLVED,
  PREP-LOST, PERSISTS or PREP-CHANGED). Every settle is counted in §9.
- **Front test availability at b:**
  - **AVAILABLE:** a PERSISTS state exists at b.
  - **UNAVAILABLE-LOST** or **UNAVAILABLE-CHANGED:** the preparation gave a resolved PREP-LOST or PREP-CHANGED. The
    front test is not run. This describes this preparation path only. It does not show that no 34-cell state exists
    at b.
  - **UNRESOLVED:** the preparation was PREP-UNRESOLVED. It counts as an unresolved front classification at b.
- Set-up when available: on D02, the prepared state translated by a whole number of nodes so that a cell minimum sits
  at the node nearest x = 0 (the **phase convention**), and multiplied by a window
  w(x) = ½ [tanh((x + 25)/ℓ) − tanh((x − 25)/ℓ)] with ℓ = λ/2. The uniform state (φ = 0) is used elsewhere. The patch is
  50 units = 17 cells.

### 6.2 Observables
- **Envelope (physical-width kernel).** K_λ is a boxcar average over a centred window of physical width exactly λ,
  discretized at each h with full weight for nodes wholly inside the window and fractional weights for the two end
  nodes, so the weights sum to λ/h (periodic wrap). Local mean m = K_λ[φ_X]; envelope E = √(2 K_λ[(φ_X − m)²]). For a
  pure cosine of wavelength λ this gives E = 1 at every point and phase, apart from discretization.
- **Envelope validation (synthetic, both D02 and D01, before any Model G run):**
  - For a unit cosine at 10 phase offsets spanning one period, E lies within [0.98, 1.02] everywhere.
  - For tanh-windowed cosine patches with known edges, at 10 phase offsets, the tracked interfaces lie within λ/10 of
    the true edges.
  - Each failure blocks the study until the kernel is fixed and re-registered.
- E_pat(b) is the spatial mean of E for the prepared state at b.
- **Patterned set and the tracked patch.** P(t) = {x : E ≥ ½ E_pat}, split into connected components on the periodic
  domain.
  - Components separated by a gap shorter than λ are merged (small holes are ignored).
  - Components shorter than λ are discarded.
  - The **tracked patch** at t = 0 is the component containing x = 0. At each later sample, it is the component with
    the largest overlap with the previous tracked patch.
  - **BREAKUP:** the tracked patch splits into two components that each have overlap at least λ.
  - **REMOTE:** a component not overlapping the tracked patch appears.
- **Interfaces.** x_L(t) and x_R(t) are the two ends of the tracked patch, located by linear interpolation of E = ½ E_pat,
  sampled every 5 time units, and unwrapped. Outward displacement is positive for both. They are tracked separately.
  - **Fill:** the tracked patch has no gap of length λ or more.
  - **Loss:** no component remains.
- **Cell count n(t)** (secondary): local minima of φ_X with prominence ≥ ½ E_pat and separation ≥ λ/2.

### 6.3 Classification (one ordered procedure; the first rule that applies is the classification)
Notation over the window [100, T]: Δ_L, Δ_R are the interface displacements; v_L, v_R are the least-squares velocities
over the last 200 units; r_T is the final residual. A **reversal** of an interface means there are times t₁ < t₂ < t₃ in
[100, T] with (x(t₂) − x(t₁)) and (x(t₃) − x(t₂)) of opposite sign, each of magnitude at least λ/2.

1. **UNRESOLVED-NUMERICAL:** the run is not complete (§3), or the tracked patch cannot be identified at t = 0.
2. **Terminal events (from t = 0, any time up to T):**
   - fill first occurring at t_f and holding to T → **INVADES** (t_f recorded, even if t_f < 100);
   - loss first occurring at t_l and holding to T → **RETREATS** (t_l recorded);
   - fill or loss that later reverses (recovery) → **UNRESOLVED-TEMPORARY**.
3. **Disqualifiers on [100, T]:** BREAKUP or REMOTE → **UNRESOLVED-INTERNAL**; a reversal of either interface →
   **UNRESOLVED-REVERSAL**.
4. **Sustained motion:**
   - Δ_L ≥ λ and Δ_R ≥ λ → **INVADES**;
   - Δ_L ≤ −λ and Δ_R ≤ −λ → **RETREATS**;
   - any other case with max(|Δ_L|, |Δ_R|) ≥ λ (one front moving and one not, or drift of the whole patch) →
     **UNRESOLVED-ASYMMETRIC**.
5. **Arrest:** |Δ_L|, |Δ_R| ≤ λ/2; |v_L|, |v_R| ≤ λ/2000 per unit time; r_T ≤ 1e-5; and the cell-count range
   max n − min n over [100, T] is at most 1 → **ARRESTED**.
6. **Otherwise** → extension (below) or **UNRESOLVED-SLOW**.

**Extension.** Every run first goes to T = 600.
- A classification from rules 1–4 at T = 600 is final.
- ARRESTED at T = 600 with margin (|Δ_L|, |Δ_R| ≤ λ/4 and |v_L|, |v_R| ≤ λ/4000) is final.
- Anything else (ARRESTED without margin, or rule 6) is continued as the same integration to T = 1800 and reclassified
  by rules 1–6 over [100, 1800]. Rule 6 at T = 1800 gives UNRESOLVED-SLOW.

Synthetic tests of this procedure (before any Model G run) include tracks that satisfy more than one draft-2 clause
(ChatGPT's 2λ-out, λ-back example must give UNRESOLVED-REVERSAL), early fill and early loss, temporary loss and recovery,
both directions of asymmetric motion, breakup and remote components, and slow drift.

Interface velocities v(b) at non-arrested points are reported descriptively. No front-speed prediction is made.

### 6.4 Grid, refinement, robustness and mesh comparison
- **Grid (D02):** b = 28.50 to 30.00, step 0.05 (31 points).
- **Front classes.** Each b gets one class: **A** (ARRESTED), **M** (INVADES or RETREATS; resolved motion), **N**
  (UNAVAILABLE-LOST or UNAVAILABLE-CHANGED; no front test), or **U** (any UNRESOLVED).
- **Arrest interval.** The maximal run of A points. Several disconnected runs make P2 UNRESOLVED (§6.6).
- **Edges.** An edge of the arrest interval is:
  - **resolved** if the neighbouring grid point outside it is M;
  - **censored** if that neighbour is N, giving a branch-limited bound (the arrest interval reaches at least to the
    last A point);
  - **unresolved** if that neighbour is U, or the interval reaches the end of the scan.
- **Refinement.** Each resolved edge is bisected twice (step 0.0125) between its A point and its M neighbour. Each
  midpoint is prepared (§6.1), then classified. An A midpoint moves the inner end of the bracket; an M midpoint moves
  the outer end. A midpoint that is N or U ends refinement on that side, and the edge becomes censored or unresolved
  accordingly.
- **Robustness (follows the final edges).** After refinement, the outermost A point on each side of the final interval
  gets two repeats:
  - (a) **D02L:** the prepared D02 state at that b tiled twice (exact, since L = 200 = 2 × 100 with the same h), settled
    on D02L to t = 300 (complete, stationary, 68 cells), then a 100-unit patch with the same splice;
  - (b) **D02 with a half-period phase shift:** the prepared state translated by exactly λ/2 by Fourier interpolation
    (all three fields, mode content preserved). It is not re-settled, because re-settling could relock the phase. The
    patch is centred on a cell maximum. The shifted field is not exactly stationary on the lattice, which is disclosed.

  The point stays A only if both repeats are ARRESTED. If either is not, the point becomes U, and that edge becomes
  **unresolved**. The interval is not re-searched.
- **Mesh comparison (D01, h = 0.1).** For each final resolved h = 0.2 bracket [b₁, b₂], the points b₁ − 0.05, b₁, b₂ and
  b₂ + 0.05 (up to 8 points). At each point, the D02 prepared state (§6.1) is interpolated to D01 by Fourier
  interpolation, preserving the mode-34 phase, then settled on D01 to t = 1000 (classified like a sweep step), then
  front-tested on D01 with the same rules.
  - If a side's points do not bracket a transition, one further point is added on that side, stepping outward by 0.05
    (at most 2 extra points in total).
  - D01 is compared with its own benchmarks (b_s = 28.978, b_r = 29.412).
  - **MESH-SENSITIVE flag:** set if the arrest-interval width at h = 0.1 differs by more than half from h = 0.2's, or
    either edge offset from its own benchmark changes by more than 0.10. It is a flag, not an attribution to lattice
    locking: wavelength resolution, phase, front-speed resolution and truncation error all change with h. It reads
    NOT EVALUATED if D01 brackets could not be formed.

### 6.5 Benchmark recheck
Before P2 is scored, the h = 0.2 window edges are recomputed by MG003's rules, with longer horizons only:

- **b_s:** MG003's convention exactly: y0 = the seeded b = 29 reference at t = 400; potentials carried unchanged;
  seed_scale = 0; "exists" if φ_X(0) < −1. The horizon is t = 600 instead of 300. The bisection bracket is frozen at
  [28.80, 29.00].
- **b_r:** standard seed (seed_scale = 1, the §3 exception); "single" if exactly one minimum with φ_X < −4. The horizon
  is t = 800 instead of 400. The bisection bracket is frozen at [29.30, 29.50].
- **Procedure:** both bracket endpoints are verified first. b_s must not exist at 28.80 and must exist at 29.00. At 29.30
  the run must be single, and at 29.50 it must not be. If an endpoint fails, that benchmark is BENCHMARK-UNRESOLVED and
  P2 is UNRESOLVED. Then 8 bisection steps follow (final width 0.00078). The final brackets [b_s⁻, b_s⁺] and
  [b_r⁻, b_r⁺] are carried into §6.6.
- If either rechecked bracket's midpoint moves by more than 0.05 from MG003's, P2 is flagged **BENCHMARK-SENSITIVE**.

### 6.6 P2 outcomes (h = 0.2; tolerance 0.10, PI decision; the first that applies)
Let [l₁, l₂] and [u₁, u₂] be the final lower and upper edge brackets (l₂ and u₁ are arrested points; l₁ and u₂ are the
adjacent M points).

| Order | Label | Criterion |
|---|---|---|
| 1 | **UNRESOLVED** | the sweep start failed; a benchmark is BENCHMARK-UNRESOLVED; the arrest intervals are disconnected; or an edge is unresolved |
| 2 | **NO ARREST OBSERVED** | no grid point is A, and none is U. This means "not observed on this grid, among available points, at this resolution"; it is not an empty pinning interval. N points are listed |
| 3 | **OUTSIDE** | the arrest interval [l₂, u₁] does not intersect [b_s⁺, b_r⁻] |
| 4 | **COINCIDES** | both edges resolved; [l₁, l₂] ⊂ [b_s⁺ − 0.10, b_s⁻ + 0.10] and [u₁, u₂] ⊂ [b_r⁺ − 0.10, b_r⁻ + 0.10] |
| 5 | **OVERLAPS** | otherwise (the interval intersects the window, and COINCIDES fails or is not available) |

- **CENSORED** is added when an edge is censored. A censored edge can never give COINCIDES. Its bound is reported as
  "arrest extends at least to b = …, branch-limited".
- **ANOMALOUS** is added if the ordering is not M or N below, A, M above.
- **BENCHMARK-SENSITIVE** and **MESH-SENSITIVE** are added where they apply.
- The wording of any positive result is "front arrest consistent with a pinning interpretation of the window", never
  "snaking demonstrated".

## 7. Not registered here: stochastic barrier S(b)

This remains a possible separate MG9b, with its own design. The deterministic seed threshold is only a
protocol-specific basin-access proxy. A stochastic action depends on the noise law and the escape path as well as the
deterministic dynamics (PI decision, §10).

## 8. Descriptive outputs (no criteria)

- **D1:** α (A1, all problems); L₃ on problems (i)–(iv) with its three terms; the sweep status sequence and primary
  persistence bracket; E_pat(b); and the interface velocities v(b).
- **D2:** one Landau-style figure: σ(b), the persistence bracket, the arrest interval and the particle window on one
  b-axis.
- Deferred from draft 1 (PI decision): the seed-threshold extension to new b values.

## 9. Runs, cost, stop rule and records

| Item | Integrations |
|---|---|
| V2: 3 amplitudes + 1 tolerance repeat | 4 |
| Sweep: start + 100 steps (extensions to 3000 belong to the same integration) | 101 |
| A3 grid fronts (D02) | 31 |
| Settles for the 15 off-sweep grid points | 15 |
| Refinement: 2 edges × 2 midpoints × (settle + front) | 8 |
| Robustness: 2 points × (D02L settle + D02L front + phase-shifted front) | 6 |
| D01: 8 points × (D02 settle + D01 settle + D01 front) | 24 |
| Benchmark recheck: 2 × (2 endpoints + 8 steps) + 1 seeded b = 29 reference | 21 |
| **Planned total** | **210** |
| Reserve: V2 doubling (3 amplitudes + 1 repeat) | 4 |
| Reserve: D01 extra points (2 × 3) | 6 |
| Reserve: re-runs after numerical failure | 10 |
| **Run cap** | **230** |

- Extensions (sweep to 3000, fronts to 1800) continue the same integration and count once.
- **Cost.** About 8 s per 400 time units for a single-particle field (the disclosed timing). Pattern fields are assumed
  10× slower:
  - the sweep (101 × up to 3000 units) about 1.5–2 core-h;
  - fronts (about 50 × 600–1800 units) about 1.5 core-h;
  - the rest under 1 core-h.

  **Expected 4–5 core-h. Cutoff 10 core-h.** It runs in the cloud container in resumable stages, each under the
  2-hour background limit. No PC time.
- **Stop rule.** If the run cap or the cutoff would be exceeded, the run stops and the PI decides. Work proceeds in this
  priority order:
  1. A1, the code check, V2 and P1;
  2. the sweep;
  3. the benchmark recheck;
  4. the A3 grid;
  5. refinement and robustness;
  6. the D01 comparison.

  Anything not reached is reported NOT RUN. P2 is UNRESOLVED if stage 3, 4 or 5 is incomplete. MESH-SENSITIVE is NOT
  EVALUATED if stage 6 is incomplete.
- **Records, committed:**
  - per-run JSON: b; domain ID; preparation provenance; carry rule; seed setting; success and time reached;
    classification and the rule that gave it; x_L(t), x_R(t), n(t); components; E_pat; residuals;
  - the A1/A2 numbers and the V2 fits;
  - the stored 34-cell states (N = 500, about 120 states × 12 kB), as one `.npz` under 5 MB in `results/MG9/`.
- **Code:** `src/landau_mg9.py` (A1, A2, the Swift–Hohenberg check, V2 fit); `src/run_mg9.py` (`mg9_ode` wrapper,
  sweep, preparation, A3, refinement, robustness, D01, benchmark recheck; resumable); `src/evaluate_mg9.py`; and
  `tests/test_mg9.py`. The tests cover:
  - seed settings and the wrapper's agreement with `lattice_ode_1d`;
  - the carry rules;
  - the Swift–Hohenberg coefficient;
  - the V1 thresholds;
  - the envelope validation of §6.2;
  - the classification procedure of §6.3 on synthetic tracks;
  - the edge and outcome logic of §§6.4–6.6 on synthetic class sequences.
- **Pre-run implementation specification** (CLAUDE.md rule 1): anything still open is fixed in GATES.md before the
  first run.

## 10. Decisions (adopted from ChatGPT MG026 by the PI, 2026-11-04 16:23 UTC)

| Decision | Adopted |
|---|---|
| Scope | Deterministic Part A only. P1 is the local onset type; P2 is the relation between measured front arrest and the operational single-particle window. Any snaking claim is deferred to branch continuation. |
| P2 tolerance | 0.10, provisionally, with explicit edge brackets, targeted refinement and UNRESOLVED (§6.4–6.6). |
| Detector | A direct envelope/interface measurement for motion (primary); unsmoothed minima for cell counts (secondary). Definitions, thresholds and phase convention frozen in §6. |
| Budget | P1 and P2 have priority; every integration is counted, with reserves (§9). Draft 1's seed-threshold extension is deferred. |
| Stochastic follow-up | MG9b stays separate (§7). |

Open for the PI: register draft 3, or send it to ChatGPT for a short round 3. ChatGPT r2 said that once its three items
and the implementation details are fixed, it has no remaining objection to the scientific scope. All are fixed in this
text. Claude recommends registering, with ChatGPT's next look being the pre-run specification and code.

## 11. What each outcome would mean

- **SUBCRITICAL and COINCIDES:** the onset is locally subcritical, and fronts arrest over the same b-range in which a
  single particle survives. That is consistent with b acting as the control parameter of a subcritical Landau-type
  instability, with the window as a pinning region. This is the strongest statement MG9 can support. It concerns the
  normal form and finite-time front dynamics, not a global potential, and not a Higgs identification.
- **SUBCRITICAL and OVERLAPS (or CENSORED):** the picture is right in kind, but the single particle's window is not
  simply the front-arrest range, or one edge is limited by the preparation path. The edges would need their own
  explanation, for example the single-dip branch's own folds.
- **SUPERCRITICAL (validated), with any P2:** the primary onset is locally supercritical. The patterns known below b_c
  would lie on a branch not connected to onset in the simple way, and any front arrest is reported as found.
- **NO ARREST OBSERVED or OUTSIDE:** the window is not explained by front arrest at this resolution. The Landau reading
  of b reduces to "b controls the stability of the uniform state", and P1 stands on its own.
- **DEGENERATE, MIXED, NOT EVALUATED, UNRESOLVED, or any flag:** reported as such, with what would resolve it (longer
  runs, continuation, finer meshes).

# MG9 — Is b the control parameter of a subcritical Landau-type instability? Design proposal, DRAFT 1

*Stage MG026. Claude, 2026-11-04. Status: **DRAFT, not registered.** No MG9 computation of any kind has been run. The
only numbers quoted below are already committed (MG001, MG003). This draft is for PI decisions (§10) and ChatGPT
design review.*

*Origin: project thread, 2026-11-04. The PI asked how b maps to physics, then to vacuum energy, and relayed ChatGPT's
suggestion of a Higgs-like effective potential V_eff(A; b) = ½ r(b) A² + ¼ u A⁴ + …. Claude and ChatGPT both proposed
testing whether such a potential actually emerges from Model G, rather than assuming it. The PI asked for this draft
("yes, draft that", 16:03 UTC).*

*Naming: MG8 stays reserved for the multiway/foliation gate. This is MG9. Part B (the stochastic barrier scan, §7) is
not part of this registration and would be registered separately as MG9b.*

## 1. Question

> **Does b act as the control parameter of a subcritical, Landau-type pattern instability, with the single-particle
> window as that instability's pinning region?**

In plain terms: the single particle exists only in a narrow band of b just below the Turing threshold b_c. A standard
mathematical picture predicts exactly that kind of band. In a **subcritical** instability (the cubic Landau coefficient
has the "wrong", destabilising sign), a patterned state coexists with the stable uniform state below threshold, and
fronts between the two do not move over a finite range of the control parameter (the **pinning region**). Inside the
pinning region, localized states such as single dips are stable. This is "homoclinic snaking", studied in detail for
the Swift–Hohenberg equation (Burke & Knobloch 2006, 2007) and in reaction–diffusion systems.

If Model G fits this picture, then b plays the role the quadratic coefficient r(b) plays in a Landau (and, by analogy,
Higgs) potential, and the window edges b_s and b_r have a precise meaning. If it doesn't fit, the Higgs-like reading of b
is not supported at this level, and that is worth knowing before anyone builds on it.

**What MG9 does not test.** It does not test any Higgs identification. That would need symmetry breaking and
mass-generating couplings, which Model G does not supply (ChatGPT's point, adopted). It does not assume Model G has an
energy functional; it does not (the system is driven and dissipative). The "potential" here is the normal form of the
amplitude equation near onset, which exists for non-variational systems too. It does not convert anything to physical
units or to vacuum energy.

## 2. What is already known (disclosed so the predictions are not mistaken for blind tests)

All on the 1D lattice, h = 0.2, L = 100, reference parameters (eq. 17) with b varied:

| Quantity | Value | Source |
|---|---|---|
| Turing threshold b_c | 30.0458 (h = 0.2); 30.0476 (h = 0.1); 30.04408 continuum | `results/MG003/existence_window.txt` |
| Lower window edge b_s | 28.9286–28.9287 (h = 0.2); 28.9782 (h = 0.1) | same |
| Upper window edge b_r | 29.3965–29.3969 (h = 0.2); 29.4117–29.4121 (h = 0.1) | `results/MG003/replication_edge.txt` |
| Deterministic seed threshold | 0.949 (b = 29), 0.887 (b = 29.16), 0.763 (b = 29.5) | `results/MG003/lattice_reference*.txt` |
| Above b_r, the particle self-replicates until the domain fills | MG003 | `replication_edge.py` |

Two consequences for honesty:

1. A space-filling patterned state is already known to exist at b between b_r and b_c (b < b_c). So "the instability is
   subcritical" is **already strongly indicated**. MG9's subcriticality item (P1) is therefore a check that the
   weakly nonlinear theory agrees with what the simulations show, not a discovery.
2. The seed threshold already falls with b across three values. A barrier proxy that decreases with b is therefore not a
   new prediction over 29–29.5. MG9 reports it (D2) but does not score it.

The genuinely open item is **P2**: whether the range of b where pattern fronts are pinned coincides with the measured
single-particle window. Nobody has computed that for Model G as far as this project knows, and the answer could go
either way (§6).

## 3. Model and numerics

- Equations: dimensionless Model G, eqs. 7a–7c, reference set eq. 17 (a = 14, g = 0.1, p = q = 1, s = u = w = 0,
  d_x = 1, d_y = 12), b varied. Same as every Phase 2 study.
- Geometry: periodic 1D lattice, h = 0.2 (primary), with h = 0.1 repeats where stated. 1D is primary because the window
  edges are measured precisely there and the 1D amplitude equation has no quadratic (hexagon) term, so its cubic
  coefficient decides super- versus subcriticality on its own.
- Integrator: the existing deterministic lattice ODE (`lattice.lattice_ode_1d`, BDF, rtol 1e-8, atol 1e-10), with a
  continuum cross-check from `mg_continuum` where stated. No noise anywhere in MG9.
- Everything is deterministic, so every number is reproducible to solver tolerance on any machine.

## 4. Part A1 — linear coefficient r(b)

Compute the largest growth rate σ_max(b) of the uniform state (lattice dispersion relation, Laplacian eigenvalue
−(4/h²) sin²(kh/2) in place of −k²) over b ∈ [28.5, 30.5] in steps of 0.05, on both the continuous-k relation and the
discrete modes of L = 100.

Report: b_c (validation V1: reproduces 30.0458 to 1e-3), the critical wavenumber k_c, and the slope
**α = dσ_max/db at b_c**. In the Landau form, r(b) = −σ_max(b) ≈ −α (b − b_c). α is the number that fixes how fast the
"potential" tilts with b. It is descriptive, with no pass criterion.

## 5. Part A2 — cubic Landau coefficient (P1)

Derive the 1D amplitude equation at onset by the standard multiple-scale expansion:

dA/dT = α (b − b_c) A + ξ² ∂²A/∂X² + L₃ |A|² A,

computed numerically from the Jacobian, the second- and third-order reaction terms (the X²Y nonlinearity) and the
second-order responses at wavenumbers 0 and 2k_c. Sign convention: **L₃ > 0 means subcritical** (the cubic term
destabilises).

- Computed on the continuum relation and on the h = 0.2 lattice.
- Validation V2: integrate a single small-amplitude mode (A₀ = 1e-3 to 1e-2 relative to the uniform state) at b = b_c on a
  domain holding a whole number of critical wavelengths, and fit the early growth dA/dt ≈ L₃ A³. The fitted L₃ must
  agree with the derived L₃ to within 10%, with the same sign. If V2 fails, P1 is reported as NOT EVALUATED.

**P1 (prediction): L₃ > 0 (subcritical), on both the continuum and the lattice.**

Also reported for P1: the fold b_f of the periodic pattern branch, found by deterministic continuation. Start from a
saturated pattern at b = 30.5 on a domain holding an integer number of critical wavelengths, step b down by 0.02,
integrate each step to stationarity from the previous state, and record the last b at which the pattern survives.
Expected: b_f < b_c (already indicated by §2). Reported against b_s: theory places the pinning region inside
[b_f, b_c], so b_f ≤ b_s is expected. The ordering b_f ≤ b_s is reported, not scored.

## 6. Part A3 — front pinning against the particle window (P2, primary)

**Set-up.** On the L = 100 periodic lattice, place the saturated periodic pattern at the same b (taken from the A2
continuation branch) in half the domain and the uniform state in the other half, joined by a smooth tanh window of
width one wavelength. The two fronts are mirror images. Integrate deterministically to t = 600.

**Measurement.** Count the pattern cells (minima of φ_X below −0.5·|A*_s(b)|, the production detector rule applied to
the unsmoothed field; or the smoothed-detector rule if the PI prefers, §10) every 10 time units. Fit the cell count
n(t) over t ∈ [100, 600].

- **INVADES:** n increases by at least 2 over the fit window.
- **RETREATS:** n decreases by at least 2, or the pattern disappears.
- **PINNED:** |n(600) − n(100)| ≤ 1, and the front position (outermost cell) moves less than one wavelength.

**Grid.** b = 28.70 to 29.70 in steps of 0.05 (21 values), h = 0.2. Then the same grid at h = 0.1, but only the six
values nearest each measured window edge (12 runs), to show whether lattice discreteness is doing the pinning.

**P2 (prediction):** the pinned interval [b_p⁻, b_p⁺] (bounded by the grid points where the classification changes)
agrees with the single-particle window.

| Label | Criterion (h = 0.2) |
|---|---|
| **COINCIDES** | each pinned-interval edge is within 0.10 of its window edge: \|b_p⁻ − 28.929\| ≤ 0.10 and \|b_p⁺ − 29.397\| ≤ 0.10 |
| **OVERLAPS** | the pinned interval is non-empty and intersects [28.929, 29.397], but COINCIDES fails |
| **DISJOINT** | the pinned interval is empty, or does not intersect the window |

INVADES above the pinned interval and RETREATS below it are part of the expected picture. Any other ordering (for
example, invading below and pinned above) is recorded as **ANOMALOUS** regardless of the label.

**Why this could fail.** Several things could break the correspondence:

- The snaking theory's pinning region belongs to fronts of an extended pattern. A single dip has its own existence
  range, which can differ.
- Front pinning on a lattice also has a purely numerical cause: fronts lock to the grid, as in Peierls–Nabarro
  pinning. The h = 0.1 repeat separates this from intrinsic pinning. If the pinned interval shrinks by more than half
  at h = 0.1, the result is flagged **LATTICE-PINNED**.
- Model G is three-component and non-variational, so fronts can also drift or oscillate.

The tolerance 0.10 is about 20% of the window width (0.47). It is the PI decision in §10(b).

## 7. Part B (not registered here) — the stochastic barrier S(b)

MG3b measured the nucleation barrier S = 0.064 per unit Ω (≈ 0.047–0.082) at b = 29. A Landau picture predicts that S
falls as b rises toward b_c, because the critical nucleus shrinks to zero amplitude at a subcritical onset. Measuring S at
a second b needs a new MG3b-style rate scan (at least three Ω values at the new b). That is expensive, and near either
window edge noise drives loss or replication (MG2 C4), which complicates detection.

Proposal: register Part B as MG9b **only if** P1 holds and P2 is COINCIDES or OVERLAPS, with its own pre-run
specification and budget. MG9 produces the deterministic inputs Part B would need: the critical-nucleus proxy (D2) and the
fronts.

## 8. Descriptive outputs (no criteria)

- **D1:** α, k_c, L₃, b_f, on both h values and the continuum.
- **D2:** the deterministic seed threshold θ(b) at b ∈ {28.95, 29.0, 29.1, 29.2, 29.3}, by bisection with the existing
  MG003 seed (judged at t = 150, as in `lattice_reference.py`). Validation: θ(29.0) reproduces 0.9492–0.9493. This is
  the deterministic barrier proxy. Its trend over 29–29.5 is already known (§2).
- **D3:** front velocities v(b) outside the pinned interval (cells per unit time × λ/2), for comparison with the
  amplitude-equation estimate near b_c.
- **D4:** a plain-language Landau diagram: σ_max(b), the pattern branch, the pinned interval and the particle window on
  one b-axis.

## 9. Cost, runs and records

- A1, A2 and V2: seconds to minutes.
- A2 continuation: about 80 b-steps × up to 300 time units on a small domain, under 0.3 core-h.
- A3: 33 runs × 600 time units on N = 500 (h = 0.2) or N = 1000 (h = 0.1). A single-particle reference run takes about
  8 s per 400 time units at either size (timed in the cloud container, 2026-11-04, b = 29 reference only). A half-domain
  pattern has many more active cells, so allow 10× that: about 2 core-min per run, roughly 1 core-h in total.
- D2: 5 bisections × about 14 runs, well under 0.5 core-h.
- **Total: about 1–2 core-h. Run cap 120 deterministic runs. Cutoff 6 core-h.** It runs in the cloud container. No
  stochastic runs, no PC time needed.
- Records: a per-run JSON (b, h, classification, n(t), front position, max |dφ/dt|), committed, plus the A1/A2 numbers
  in `results/MG9/`. Final fields are small (N ≤ 1000 × 3) and can be committed as one `.npz` under 5 MB.
- Code: `src/landau_mg9.py` (A1, A2, V2), `src/run_mg9.py` (continuation, A3, D2), `src/evaluate_mg9.py`, tests in
  `tests/test_mg9.py`. These include a test that the derived L₃ reproduces the known cubic coefficient of a textbook
  case (the Brusselator or Swift–Hohenberg) before it is applied to Model G.

## 10. Decisions for the PI

(a) **Scope.** Register Part A only (recommended), with Part B as a separate conditional MG9b. The alternative is to
register both now, with Part B's b-values fixed in advance.

(b) **P2 tolerance.** 0.10 in b (recommended, about 20% of the window). A tighter 0.05 would make COINCIDES a sharper
claim but risks failing on grid resolution alone, since the grid step is 0.05.

(c) **Cell detector for A3.** Unsmoothed φ_X with the −0.5·|A*_s(b)| threshold (recommended, since there is no noise to
smooth), or the production smoothed detector (σ = 0.4) for continuity with MG3.

(d) **Review.** Send this draft to ChatGPT for a design review before registration (recommended, because P2's
set-up and the L₃ derivation both have choices a second reader should check).

## 11. What each outcome would mean

- **P1 holds and P2 COINCIDES:** b behaves as the control parameter of a subcritical Landau instability, and the
  particle window is its pinning region. In that picture, b − b_c plays the role of the Landau or Higgs quadratic
  coefficient, and particles are the localized states of a vacuum held just below its pattern threshold. This is the
  strongest statement MG9 can support. It is a statement about the dynamics' normal form, not a Higgs identification.
- **P1 holds, P2 OVERLAPS:** the Landau picture is right in kind, but the single particle's window is not simply the
  front-pinning range. The edges need their own explanation, for example the single-dip branch's own folds.
- **P2 DISJOINT or ANOMALOUS, or LATTICE-PINNED:** the particle window is not explained by snaking on this lattice. The
  Higgs-like reading of b gets no support beyond "b controls stability".
- **P1 fails (L₃ < 0, validated):** a contradiction with the known coexistence below b_c. It would point to an error in
  the derivation or to higher-order terms dominating, and is reported as such before anything else.

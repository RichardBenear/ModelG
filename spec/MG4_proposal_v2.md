# MG4 — Particles on frozen irregular graphs (Phase 3): design proposal, DRAFT 2

**SUPERSEDED** by draft 3 (`spec/MG4_proposal.md`). Kept unchanged below as the record.

Status: **DRAFT 2** (Claude, 2026-10-30, stage MG009). Draft 2 applies the PI's scope decision (2026-10-30 13:41 UTC,
project thread: "2D graphs, noise-free only"): draft 1's stochastic Part B is removed and deferred, and the primary h
returns to 0.2. Draft 1 is in the git history of this file (commit on PR #12). Not registered. Nothing in it binds until the PI approves it after
ChatGPT's design review and it is entered in `GATES.md`. **No MG4 dynamics of any kind has been run.** The only
computation behind this draft is the graph statistics in "Design-stage numbers" (degree counts and edge-length moments of
random point sets, with no reaction or diffusion), used to size the time step and the cost.

Roadmap question (README, `InitialPlan.md`): *do particles survive irregular connectivity, and does the graph Laplacian
behave like ∇² at the soliton's scale?*

## 1. Conclusions of the design, first
1. **MG4 runs in 2D, not 1D or 3D.** In 1D the frozen rules make an irregular graph pointless: hopping uses only
   connectivity (rules v1.1, D_S: one rate per incident edge), and every connected nearest-neighbour graph on a ring *is*
   the regular ring. 3D is about 30–100× dearer per stochastic run at the resolution needed. In 2D, b = 29 lies well
   inside the continuum window (b_s = 27.34; MG1's 2D soliton is stationary at b = 29), so MG4 avoids the 1D near-edge
   fragility that dominated MG2 and MG3.
2. **The dynamics see only the graph.** Node positions are used for two things only: building the graph, and measuring
   (detection, profiles, distances). One global spacing h per graph turns hop counts into lengths. It is calibrated from
   the graph's long-wavelength diffusion, not from the soliton (§4).
3. **Comparator, not continuum, is the reference.** Each irregular graph is compared with a regular triangular lattice of
   the same effective spacing, run with identical code. This separates "irregular" from "2D" and from "discretized".
   The MG1 continuum 2D soliton is reported alongside.
4. **Noise-free only (PI decision).** Part 0 (graph construction and calibration, automatic rules, frozen output) and
   Part A (deterministic graph ODE, the Ω → ∞ limit that MG2 showed the token rules converge to; criteria C1–C3). Survival
   under token noise is **not** tested by MG4. It is deferred to a later, separately registered study (§6). MG4's verdict
   therefore applies to the mean-field limit of the rules on these graphs.
5. **Wolfram hypergraph snapshots are deferred to a separate study (MG4b).** They have no embedding, need a choice of
   rule from the PI, and need graph-intrinsic detection and dimension estimates. MG4's intrinsic spectral test (C3) is
   written so that it carries over.

## 2. Graph families (all periodic, on a square torus of side L_box = 20 ≈ 6.5 λ₀)
| Code | Graph | Role | Parts |
|---|---|---|---|
| **T** | Triangular lattice (degree 6) | Regular comparator and positive control | 0, A |
| **D** | Delaunay triangulation of a periodic Poisson point set (planar, mean degree 6, degree 3–13) | **Primary irregular family** | 0, A |
| **R** | Random geometric graph of a periodic Poisson point set, connection radius set for mean degree 12, largest component only | Second irregular family (denser, non-planar) | 0, A |
| **N3** | Random 6-regular graph with D's node count (no geometry) | Negative control for C3's discriminating power | 0 (A descriptive) |
| **W** | D's graph with finite-volume (Voronoi) edge weights | Optional diagnostic: separates "irregular points" from "unweighted rule". **Needs a rules amendment** (per-edge rates), so it is exploratory and deterministic only, and runs only if the PI includes it | A (exploratory) |

Construction is fully seeded (stream in §8). D and R: 4 independent realizations each. T: one graph. N3: 2.
The point density of D and R is set so that their calibrated h (§4) equals the target, which gives different node counts
from T. (This matters only for a later noisy study, which would match noise per unit area, Ω·ρ.)

## 3. Operating point and resolution
- b = 29 (published), eq. 17 parameters otherwise, rules v1.1 unchanged.
- **Primary h = 0.2** for every criterion, the Phase 2 spacing (triangular spacing 0.163, about 15 nodes per λ₀). Draft 1
  used 0.3 only to make the stochastic part affordable; with that part removed, the finer spacing costs little.
- **h = 0.3 resolution check:** Part A repeated on T and on one D realization (8 positions), descriptive only.
- **Positive control (Part 0), automatic:** on T at h = 0.2, measure the deterministic single-particle window [b_s, b_r]
  by bisection (as `src/existence_window.py`, to ±0.005). If b = 29 is not inside it with margin ≥ 0.1 on both sides,
  or the T particle's core φ_X depth is more than 20% off the MG1 continuum 2D value (−13.6), MG4 stops for a PI decision
  (no criterion is evaluated).
- Box check (Part 0): T particle at L_box = 20 against L_box = 30; radial φ_X profile difference over r < 6 must be
  < 0.5%, otherwise L_box = 30 for everything and the budget is re-projected before Part A.

## 4. Calibrating h from the graph (Part 0, automatic)
With one hop rate r per edge, the graph Laplacian is r(A − D_deg) (A adjacency, D_deg degree). It is symmetric and
conserves tokens on any graph. Its long-wavelength behaviour is set by the **homogenized diffusivity**, which for an
irregular network is below the naive second-moment estimate (the arithmetic mean of each node's Σ_j |Δx_ij|²/4, the
"Voigt" bound). For D the per-node second moment has coefficient of variation 0.72 (Design-stage numbers), so the two
can differ materially.
- **Rule:** for each graph, solve the periodic corrector problem (two sparse linear solves, one per axis) for the
  effective diffusivity tensor D_eff of the unit-rate Laplacian, and set h² = trace(D_eff)/2. The anisotropy
  |D_eff,xx − D_eff,yy| / trace is reported.
- The naive second-moment h is reported next to it.
- T's h is exact (h² = 1.5 a²). Applying the corrector to T must return this to 1e-10 (implementation test).
- N3 has no embedding, so it gets no calibration. It uses D's h, so its rates match the primary family.

## 5. Part A — deterministic graph ODE (criteria C1–C3)
Equations: the rules' Ω → ∞ mean field on the graph, with ∇² replaced by (A − D_deg)/h². Stiff integration as
`lattice_ode_1d` (BDF, sparse Jacobian), t ∈ [0, 300].

**Initial condition (survival test):** the MG1 continuum 2D stationary soliton, its radial profiles interpolated onto node
positions around a centre c. This does not depend on the graph. Centres: 8 per D or R realization (32 per family), drawn
uniformly on the torus from the registered stream. T: 4 centres (a node, an edge midpoint, a triangle centre, one random).
**Secondary (seeded formation, descriptive):** the published seed (eq. 18), scale 1.0 (1.75× the 2D continuum threshold),
at the same centres, homogeneous start.

**Measurement (embedding used here only):** φ_X smoothed with a Gaussian of σ = 0.4 over node positions (normalized
weights; the MG2 detector width). Detection threshold −0.5·|A*_s(T)|, where A*_s(T) is the smoothed depth of the T particle.
Particle centre: the smoothed minimum, refined by a weighted centroid of the nodes below threshold. Radial profiles:
mean φ_S in rings of width 0.25 around the centre, out to r = 8.

**Outcome classes** (the §5.4 logic, without noise): at the 11 check times t = 200, 210, …, 300,
- *single particle*: exactly one detection at every check, identity-linked (within 1.0 between checks);
- *decayed*: no detection at t = 300;
- *multiplied*: two or more detections at any check;
- *unresolved*: anything else.

**C1 — survival.** Per family (D and R separately), the fraction of the 32 runs ending as *single particle*.
Clopper–Pearson 95% interval.
- PASS if the lower bound ≥ 0.75 (needs ≥ 28 of 32);
- FAIL if the upper bound < 0.75 (≤ 18 of 32);
- INCONCLUSIVE otherwise.
Reported with it: the class counts, the drift |c_final − c| (pinning), and max |∂φ/∂t| at t = 300.

**C2 — it is the same particle.** Over the family's *single particle* runs, against T's particle (mean of its 4 runs):
- (a) relative L2 difference of the radial φ_X profile over r < 6, area-weighted (r dr);
- (b) core depth φ_X at the centre (smoothed and unsmoothed);
- (c) λ by MG3's core estimator (2 × mean zero spacing of the radial φ_X, 1 < r < 7).
PASS if the **median** run has (a) < 10%, (b) within 10% and (c) within 5%, and at least 80% of runs meet (a). FAIL if the
median run misses any of the three. INCONCLUSIVE otherwise, or if fewer than 8 runs are *single particle*.
Why 10%: on the 1D ring the h = 0.2 lattice particle is 8% off the continuum. The tolerance lets disorder cost about as
much as discretization already does, and no more.

**C3 — the Laplacian is ∇² at the soliton's scale.** Intrinsic and embedding-free, so it carries to MG4b.
- Compute the graph Laplacian's eigenvalues μ_m up to 1.5 k*² (k* = 2π/λ₀ ≈ 2.04; about 200 eigenvalues for L_box = 20,
  by shift-invert Lanczos; about 24 000 nodes for D at h = 0.2).
- Weyl test: in 2D the continuum count is N(μ) = A·μ/(4π), with A = L_box². Define k_eff(k) = √(4π N_graph(k²)/A).
- **Criterion:** |k_eff(k*)/k* − 1| ≤ 0.025 (≈ ±5% in k²) for every realization of the family.
- PASS if every realization meets it. FAIL if the family's mean misses it. INCONCLUSIVE otherwise.
- Descriptive: the plane-wave spectral function at |k| = k* (16 directions): its median, compared with k*², and its width.
  The width measures how far modes at the soliton scale are from plane waves (scattering). T's width is near 0.
- **Discriminating power (N3):** C3 must FAIL on the random regular graph. Its spectrum has no small eigenvalues apart from
  the constant mode. If N3 passed, C3 would not discriminate, and C3 would be INCONCLUSIVE for every family.

**Checks that apply to every Part A run (implementation, not criteria):**
- The homogeneous state is an exact fixed point on any graph (rows of A − D_deg sum to 0). Unseeded runs from it stay at
  round-off level to t = 300.
- **Linear stability needs no new computation:** σ depends on the graph only through each Laplacian eigenvalue μ ≤ 0,
  and MG001 showed max Re σ < 0 (−0.195) for every μ = −k² ≤ 0 at b = 29. So the homogeneous state is linearly stable on
  every graph in MG4 (Othmer–Scriven). A 10⁻³ random perturbation must decay (one run per family).
- **Integral identity on graphs:** for a symmetric Laplacian, a stationary particle has Σ_i φ_G = Σ_i φ_X = 0 (node sums;
  `model_g_published.md` §9.2). Check: |Σ φ_S| / Σ |φ_S| < 10⁻⁶ at t = 300 in every *single particle* run.

## 6. Deferred: survival under token noise (draft 1's Part B)
Removed by the PI's scope decision. Draft 1 proposed stochastic runs on T and D at matched noise per unit area
(Ω_T = 1000), 8 runs each to T = 100, with graph-kernel validations (bit-identity with `tau_leap_euler_1d` on a cycle,
the mean map, exact SSA on a small graph, a clip pilot), at about 40 core-hours. A later study can take it up as written,
using Part A's stationary particles as its starting states and reference. Part A saves every final field so that this
remains possible without re-running.

## 7. Gate outcome
- **PASS:** the positive controls hold (T exists in the window, box check, C3 on T), N3 fails C3 as predicted, and C1, C2,
  C3 PASS for D and R.
- **FAIL:** any of C1–C3 FAILs for D or R. Each family's results are reported separately.
- **INCONCLUSIVE:** a positive control fails at both h, N3 passes C3, any criterion is INCONCLUSIVE, or a budget stop.
- W (if included), the h = 0.3 check and the seeded-formation runs are reported and enter no criterion.

**Prediction (stated before any run; not a criterion).** D and R pass C1–C3. Disorder pins the particle (drift under one
λ₀ from the starting centre). N3 fails C3. Least certain: C2 on D. Its per-node diffusivity varies by 72%, and
its disorder correlation length (about h) is only about 6× smaller than λ₀, so homogenization at the core scale is not
guaranteed.

Pre-stated interpretation of failures:
- C1 fails, with *decayed* or *multiplied* outcomes: disorder shifts the existence window enough to push b = 29 out of it
  locally. The window on D would then be measured as a follow-up.
- C1 passes and C2 fails: particles exist but are a different object (graph-specific), not Model G's.
- C3 fails: the graph is not ∇²-like at the soliton scale, and survival (if any) does not rest on continuum behaviour.

## 8. Implementation, streams and records
- New code (tests first): `src/graphs.py` (builders, periodic Delaunay by 3×3 tiling, corrector, spectra),
  `src/graph_ode.py`, `src/measure_graph.py`, the drivers
  `src/run_mg4.py` and `src/evaluate_mg4.py`, and `tests/test_mg4.py`.
- Streams: `runner.seed_for` with gate id 5. Stage index: Part 0 graphs 0, Part A centres 1, timing 6. b_idx 0; om_idx = family (T 0, D 1, R 2, N3 3, W 4); replicate = realization × 100 + centre index.
- Records: every Part A final field (float32) and its outcome class, and the fields at the 11 check times.
  Part 0 writes `results/MG4/part0/decisions.json` (h per graph, L_box, node counts, degree maxima, T window).
  The Part A driver refuses to start without it.
- Everything a choice leaves open is appended to the registration, dated, before the run it affects.

## 9. Budget (run caps bind; core-hour ceilings are stop triggers checked before every batch)
| Item | Runs (cap) | Projected core-h | Ceiling |
|---|---|---|---|
| Part 0 (graphs, corrector, spectra, T window bisection, box check) | about 30 ODE runs | 3 | 8 |
| Part A at h = 0.2 (T 4 + D 32 + R 32, both initial conditions; N3 4; W 32 if included) | 140 (172) | 10 | 25 |
| Part A h = 0.3 check (T 4 + D 8) | 12 | 1 | 3 |
Total projected about 14 core-hours. These are rough: the 2D graph ODE's cost (sparse BDF, about 50 000–70 000 unknowns
for D and R at h = 0.2) has not been measured, so the first Part 0 runs time it and the projection is re-stated, in the
registration, before Part A starts. All of it is deterministic and can run in the cloud container or on the PI's machine.

## 10. Design-stage numbers (graph statistics only; scratch, not committed results)
Periodic Poisson sets on a 20 × 20 torus, 3 realizations each:
- Delaunay: mean degree 6.00, minimum 3, maximum 12–13. Per-node second moment Σ_j|Δx|²: coefficient of variation 0.72,
  maximum/mean 8.3. At the naive h, ρ·h² = 2.39 (T: 1.73), so D needs about 1.38× T's nodes. At h = 0.3 that is about
  10 600 nodes (T: 7 700).
- RGG at mean degree 12: minimum degree 0–2, maximum 24–27. Hence "largest component only". At a fixed h it needs about
  3× T's nodes.
- At h = 0.2 on the 20 × 20 torus: T about 17 300 nodes, D about 23 900.

## 11. Not in MG4 (later, separately registered)
- **MG4b:** frozen Wolfram hypergraph snapshots (rule chosen by the PI), graph-intrinsic detection, dimension estimation.
- Survival under token noise (§6); 3D graphs; nucleation (MG3-type rates) on graphs; the existence window on D; gradient drift and depinning (v = −0.0599 m
  needs a mobile particle, so pinning matters for MG6).
- D1–D3 stay deferred. The h = 0.3 / 0.2 comparison here says something about D2's question in 2D, without replacing it.

## 12. Questions for the PI and the reviewer
1. Scope: decided by the PI (2D, noise-free only; Wolfram snapshots as MG4b).
2. Include W (finite-volume weights)? It needs a rules v1.2 amendment allowing per-edge rates, used only in this
   exploratory arm.
3. The thresholds (C1 0.75, C2 10% / 10% / 5% with 80%, C3 2.5% in k) are new and set without data. MG3's
   80% threshold is the cautionary case. Are they the right questions at the right strictness?
4. Is homogenized calibration of h (§4) the right choice over the naive second moment? It is the physically correct
   long-wavelength diffusivity, but it makes C3's agreement at small k hold by construction.
5. Is the continuum-soliton initial condition (survival) the right primary test, with seeded formation descriptive?

# MG4 — Particles on frozen irregular graphs (Phase 3): design proposal, DRAFT 1

Status: **DRAFT 1** (Claude, 2026-10-30, stage MG009). Not registered. Nothing in it binds until the PI approves it after
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
4. **Three parts.** Part 0 (graph construction and calibration, automatic rules, frozen output). Part A (deterministic
   graph ODE, the Ω → ∞ limit that MG2 showed the token rules converge to; cheap; carries criteria C1–C3). Part B
   (stochastic tokens at one finite Ω; carries C4; about 40 core-hours).
5. **Wolfram hypergraph snapshots are deferred to a separate study (MG4b).** They have no embedding, need a choice of
   rule from the PI, and need graph-intrinsic detection and dimension estimates. MG4's intrinsic spectral test (C3) is
   written so that it carries over.

## 2. Graph families (all periodic, on a square torus of side L_box = 20 ≈ 6.5 λ₀)
| Code | Graph | Role | Parts |
|---|---|---|---|
| **T** | Triangular lattice (degree 6) | Regular comparator and positive control | 0, A, B |
| **D** | Delaunay triangulation of a periodic Poisson point set (planar, mean degree 6, degree 3–13) | **Primary irregular family** | 0, A, B |
| **R** | Random geometric graph of a periodic Poisson point set, connection radius set for mean degree 12, largest component only | Second irregular family (denser, non-planar) | 0, A |
| **N3** | Random 6-regular graph with D's node count (no geometry) | Negative control for C3's discriminating power | 0 (A descriptive) |
| **W** | D's graph with finite-volume (Voronoi) edge weights | Optional diagnostic: separates "irregular points" from "unweighted rule". **Needs a rules amendment** (per-edge rates), so it is exploratory and deterministic only, and runs only if the PI includes it | A (exploratory) |

Construction is fully seeded (stream in §8). D and R: 4 independent realizations each. T: one graph. N3: 2.
The point density of D and R is set so that their calibrated h (§4) equals the target, which gives different node counts
from T. Part B therefore matches **noise per unit area** (Ω·ρ, with ρ the node density), not Ω per node.

## 3. Operating point and resolution
- b = 29 (published), eq. 17 parameters otherwise, rules v1.1 unchanged.
- **Primary h = 0.3** for every criterion (triangular spacing 0.245, about 10 nodes per λ₀). This is chosen for Part B's
  cost: per stochastic time unit the cost scales roughly as h⁻⁴, and h = 0.2 would multiply Part B by about 3.4.
- **h = 0.2 resolution check:** Part A repeated on T and on one D realization (8 positions), descriptive only.
- **Positive control (Part 0), automatic:** on T at h = 0.3, measure the deterministic single-particle window [b_s, b_r]
  by bisection (as `src/existence_window.py`, to ±0.005). If b = 29 is not inside it with margin ≥ 0.1 on both sides,
  or the T particle's core φ_X depth is more than 20% off the MG1 continuum 2D value (−13.6), the primary h becomes 0.2.
  Part A then runs at h = 0.2, and Part B stops for a PI budget decision before it starts.
- Box check (Part 0): T particle at L_box = 20 against L_box = 30; radial φ_X profile difference over r < 6 must be
  < 0.5%, otherwise L_box = 30 for everything and the budget is re-projected before Part B.

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
  by shift-invert Lanczos).
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

## 6. Part B — stochastic tokens (criterion C4)
Runs on T and on D (2 of the 4 realizations, 4 centres each). Initial state counts = round(Ω·(S₀ + φ*)), where φ* is that
graph's own Part A stationary particle at that centre (MG2's start). Runs whose Part A outcome was not *single particle*
are replaced by the next registered centre.
- **Noise level:** Ω_T = 1000 per node on T; on D, Ω_D = 1000·ρ_T/ρ_D, so Ω·ρ (tokens per unit area at unit
  concentration) is equal. Ω_D is rounded to 3 significant figures.
- 8 runs per family, T = 100, snapshots every 0.5.
- **Kernel:** the unsplit multinomial tau-leap of MG3 generalized to graphs (each token picks one of: its reactions, a hop
  along any incident edge, or staying). τ = 2.5 × 10⁻⁴ for both families. Total exit probability at the highest-degree D
  node (degree 13) stays below 0.6. The clip rule and V3 threshold (10⁻⁴) apply; Part 0 records the measured degree maxima.
- Tracking: MG2's rules (loss, birth, replication within 2λ*, merge) with embedding distances. Detection as Part A,
  except that the threshold uses that graph's own Part A depth.

**C4 — survival under noise.** The number of the 8 D runs that keep their original particle as the only particle for all
t ∈ [0, 100].
- PASS if ≥ 7 of 8, provided T has ≥ 7 of 8;
- FAIL if ≤ 4 of 8 while T has ≥ 7 of 8;
- INCONCLUSIVE otherwise, including when T itself has < 7 of 8 (then the problem is 2D noise, not irregularity).
Ω = 1000 is chosen where MG2 saw no events on the ring at Ω ≥ 1000. Under matched noise, a failure on D can then be
attributed to the graph.

**Exploratory (reported, no criterion):** the particle's mean-squared displacement and Ω·D_p on T and D. Pinning by
disorder would show as a lower D_p or a plateau on D. Far-field per-node Ω·Var against the graph LNA (exact on T;
computed from the Part 0 spectrum on D, for modes up to the computed cutoff). The positions at which particles pin.

**Kernel validations (all before the first Part B run; failure stops Part B for a PI decision):**
- VG1: on the cycle graph C₅₀₀ with neighbour order (left, right), the graph kernel is **bit-identical** to
  `tau_leap_euler_1d` for a fixed seed.
- VG2: the Ω → ∞ mean map equals forward Euler of the graph ODE (expected increments, exactly; and numerically at Ω = 10⁸ on
  a 200-node D graph, as V2).
- VG3: exact SSA against tau-leap on a 50-node D graph at Ω = 3 and 30 (V1's procedure).
- VG4: clip pilot on D at Ω_D, T = 20 (V3's procedure).
- Timing: two 20-unit runs per family on the production machine. Their measured cost rescales the Part B projection.

## 7. Gate outcome
- **PASS:** the positive controls hold (T exists in the window, box check, C3 on T), N3 fails C3 as predicted, and C1, C2,
  C3 PASS for D and R, and C4 PASS for D.
- **FAIL:** any of C1–C4 FAILs for D or R. Each family's results are reported separately.
- **INCONCLUSIVE:** a positive control fails at both h, N3 passes C3, any criterion is INCONCLUSIVE, or a budget stop.
- W (if included), the h = 0.2 check and the seeded-formation runs are reported and enter no criterion.

**Prediction (stated before any run; not a criterion).** D and R pass C1–C3. Disorder pins the particle (drift under one
λ₀, lower D_p on D than on T). C4 passes. N3 fails C3. Least certain: C2 on D. Its per-node diffusivity varies by 72%, and
its disorder correlation length (about h) is only about 6× smaller than λ₀, so homogenization at the core scale is not
guaranteed.

Pre-stated interpretation of failures:
- C1 fails, with *decayed* or *multiplied* outcomes: disorder shifts the existence window enough to push b = 29 out of it
  locally. The window on D would then be measured as a follow-up.
- C1 passes and C2 fails: particles exist but are a different object (graph-specific), not Model G's.
- C3 fails: the graph is not ∇²-like at the soliton scale, and survival (if any) does not rest on continuum behaviour.
- C4 fails on D only: noise and disorder together destabilize a particle that each alone does not.

## 8. Implementation, streams and records
- New code (tests first): `src/graphs.py` (builders, periodic Delaunay by 3×3 tiling, corrector, spectra),
  `src/graph_ode.py`, `src/graph_lattice.py` (numba kernel on CSR adjacency), `src/measure_graph.py`, the drivers
  `src/run_mg4.py` and `src/evaluate_mg4.py`, and `tests/test_mg4.py`.
- Streams: `runner.seed_for` with gate id 5. Stage index: Part 0 graphs 0, Part A centres 1, Part B 2, validations 3,
  timing 6. b_idx 0; om_idx = family (T 0, D 1, R 2, N3 3, W 4); replicate = realization × 100 + centre index.
- Records: every Part A final field (float32) and its outcome class; every Part B run record as in MG2, plus the field at
  t = 100. Part 0 writes `results/MG4/part0/decisions.json` (h per graph, primary h, L_box, node counts, degree maxima,
  τ check). The Part A and Part B drivers refuse to start without it.
- Everything a choice leaves open is appended to the registration, dated, before the run it affects.

## 9. Budget (run caps bind; core-hour ceilings are stop triggers checked before every batch)
| Item | Runs (cap) | Projected core-h | Ceiling |
|---|---|---|---|
| Part 0 (graphs, corrector, spectra, T window bisection, box check) | about 30 ODE runs | 2 | 6 |
| Part A at h = 0.3 (T 4 + D 32 + R 32, both initial conditions; N3 4; W 32 if included) | 140 (172) | 6 | 15 |
| Part A h = 0.2 check (T 4 + D 8) | 12 | 2 | 5 |
| VG1–VG4 and timing | — | 1 | 3 |
| Part B (T 8 + D 8, T = 100) | 16 | 40 | 60 |
Total projected about 51 core-hours, mostly Part B (about 5–6 h wall with 9 workers). Part 0 and Part A are deterministic
and can run in the cloud container. Part B runs on the PI's machine. If the primary h falls back to 0.2, Part B needs a
new PI budget decision (about 3.4× the cost).

## 10. Design-stage numbers (graph statistics only; scratch, not committed results)
Periodic Poisson sets on a 20 × 20 torus, 3 realizations each:
- Delaunay: mean degree 6.00, minimum 3, maximum 12–13. Per-node second moment Σ_j|Δx|²: coefficient of variation 0.72,
  maximum/mean 8.3. At the naive h, ρ·h² = 2.39 (T: 1.73), so D needs about 1.38× T's nodes. At h = 0.3 that is about
  10 600 nodes (T: 7 700).
- RGG at mean degree 12: minimum degree 0–2, maximum 24–27. Hence "largest component only", and R is deterministic only.
  At a fixed h it needs about 3× T's nodes, and its time step would be about 3× smaller.
- Part B cost: scaling MG3b's measured 0.89 core-s per time unit (500 nodes, τ = 6.667e-4), with about 2× per node-step
  for the extra hop channels, gives about 2.0 (T) and 2.8 (D) core-hours per 100-unit run at τ = 2.5e-4.

## 11. Not in MG4 (later, separately registered)
- **MG4b:** frozen Wolfram hypergraph snapshots (rule chosen by the PI), graph-intrinsic detection, dimension estimation.
- 3D graphs; nucleation (MG3-type rates) on graphs; the existence window on D; gradient drift and depinning (v = −0.0599 m
  needs a mobile particle, so pinning matters for MG6).
- D1–D3 stay deferred. The h = 0.3 / 0.2 comparison here says something about D2's question in 2D, without replacing it.

## 12. Questions for the PI and the reviewer
1. Scope: 2D embedded graphs now, with Wolfram snapshots as MG4b (recommended)? Keep or drop Part B (about 40 core-h)?
2. Include W (finite-volume weights)? It needs a rules v1.2 amendment allowing per-edge rates, used only in this
   exploratory arm.
3. The thresholds (C1 0.75, C2 10% / 10% / 5% with 80%, C3 2.5% in k, C4 7 of 8) are new and set without data. MG3's
   80% threshold is the cautionary case. Are they the right questions at the right strictness?
4. Is homogenized calibration of h (§4) the right choice over the naive second moment? It is the physically correct
   long-wavelength diffusivity, but it makes C3's agreement at small k hold by construction.
5. Is the continuum-soliton initial condition (survival) the right primary test, with seeded formation descriptive?

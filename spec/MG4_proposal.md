# MG4 — Particles on frozen irregular graphs (Phase 3): design proposal, DRAFT 3

Status: **DRAFT 3** (Claude, 2026-10-30, stage MG009a). Draft 3 answers ChatGPT's design review
(`review/MG009_chatgpt.md`; response `review/MG009_response.md`). Draft 2 applied the PI's scope decision (2026-10-30 13:41 UTC,
project thread: "2D graphs, noise-free only"): draft 1's stochastic Part B is removed and deferred, and the primary h
returns to 0.2. Drafts 1 and 2 are preserved as `spec/MG4_proposal_v1.md` and `_v2.md`. Not registered. Nothing in it binds until the PI approves it after
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
   rule from the PI, and need graph-intrinsic detection and dimension estimates. MG4's spectral-count benchmark (C3(a)) is
   intrinsic, so it carries over.

## 2. Graph families and the torus
**Torus (draft 3, review item 3).** An equilateral triangular lattice has no exactly square periodic cell. Every family
therefore uses the same **rectangular torus commensurate with T**: n_x columns at spacing a and n_y rows (n_y even) at
spacing a√3/2, where a = h_t/√1.5 is the triangular spacing for the target h_t. At h_t = 0.2: n_x = 122, n_y = 142, giving
L_x = 19.92 and L_y = 20.08 (aspect ratio 1.008). The box check uses n_x = 184, n_y = 212 (30.05 × 29.98). The h = 0.3
check uses n_x = 82, n_y = 94 (20.09 × 19.94). The random point sets are drawn on the same L_x × L_y torus. Nothing is
strained.

| Code | Graph | Role |
|---|---|---|
| **T** | Triangular lattice on that torus (degree 6) | Regular comparator and positive control |
| **D** | Delaunay triangulation of a periodic Poisson point set (planar, mean degree 6, degree 3–13) | **Primary irregular family** |
| **R** | Random geometric graph of a periodic Poisson point set, radius set for mean degree 12, largest component only | Second irregular family (denser, non-planar) |
| **N3** | Random 6-regular graph with D's mean node count (no geometry) | Negative control for C3(a) |

**Benchmark panel (review item 1).** D and R each have a **fixed panel of 4 graphs × 8 placements**, drawn once from the
registered stream. Results are reported per graph and summed over the panel. The criteria are explicit counts on this
panel. **No confidence statement is made about a population of random graphs**: the 8 placements on one graph are
clustered, not 32 independent trials. T has one graph with 4 placements; N3 has 2 graphs.

W (finite-volume weights) is **deferred** (§11). Doing it properly needs node volumes and weighted conservation, not just
edge weights.

## 3. Operating point and resolution
- b = 29 (published), eq. 17 parameters otherwise, rules v1.1 unchanged.
- **Primary h_t = 0.2** for every criterion, the Phase 2 spacing (triangular spacing 0.163, about 15 nodes per λ₀).
- **h_t = 0.3 resolution check:** Part A repeated on T and on the first D panel graph (8 placements). Descriptive only.
- **Positive control (Part 0), automatic:** on T at h_t = 0.2, measure the deterministic single-particle window
  [b_s, b_r] by bisection (as `src/existence_window.py`, to ±0.005). MG4 stops for a PI decision, with no criterion
  evaluated, if b = 29 is not inside the window with margin ≥ 0.1 on both sides, or if the T particle's core φ_X depth is
  more than 20% off the MG1 continuum 2D value (−13.6).
- **Box check (Part 0):** the T particle at 20 × 20 against 30 × 30. The radial φ_X profile difference over r < 6 must be
  < 0.5%. Otherwise the 30 × 30 torus is used for everything, and the budget is re-projected before Part A.

## 4. Calibrating h from the graph (Part 0, automatic)
With one hop rate r per edge, the graph Laplacian is r(A − D_deg), where A is the adjacency matrix and D_deg the degree
matrix. It is symmetric and conserves node sums on any graph. Its long-wavelength behaviour is set by the **homogenized
diffusivity tensor**, which for an irregular network lies below the naive second-moment estimate (the "Voigt" bound).
For D, the per-node second moment has coefficient of variation 0.72.
- **h for each graph:** solve the periodic corrector problem (two sparse linear solves) for the effective diffusivity
  tensor D_eff of the unit-rate Laplacian, with eigenvalues e₁ ≥ e₂. Set h² = (e₁ + e₂)/2, and use that graph's own h in
  its hop rate.
- **Reported (review item 3):** the full tensor, its eigenvalues and principal angle, the anisotropy (e₁ − e₂)/(e₁ + e₂),
  and the naive second-moment h.
- **Reaching the target h:** the node density sets the resolution. In Part 0, four *pilot* realizations per family (not
  in the panel) estimate the constant c_F = ρ·h². Each panel graph has N = round(c_F·L_x·L_y / h_t²) points. A panel graph
  is accepted if |h/h_t − 1| ≤ 0.03. Otherwise it is replaced by the next draw from the stream, and the replacement is
  recorded. The R threshold applies to the largest component, and its node fraction is reported.
- **Tests:** on T the corrector must return h² = 1.5a² and e₁ = e₂ to 1e-10.
- N3 has no embedding. It uses the mean h of the D panel.

## 5. Part A — deterministic graph ODE (criteria C1–C3)
Equations: the rules' Ω → ∞ mean field on the graph, with ∇² replaced by (A − D_deg)/h². Stiff integration as in
`lattice_ode_1d` (BDF, sparse Jacobian), over t ∈ [0, 300].

**Initial condition (survival test):** the MG1 continuum 2D stationary soliton, its radial profiles interpolated (cubic)
onto node positions around a centre c. This does not depend on the graph. Placements: 8 per panel graph, drawn uniformly
on the torus. T: 4 placements (a node, an edge midpoint, a triangle centre, one random).
**Secondary (seeded formation, descriptive):** the published seed (eq. 18), scale 1.0 (1.75× the 2D continuum threshold),
at the same centres, from the homogeneous state.

### 5.1 Measurement specification (review item 4; the embedding is used here only)
Distances are periodic Euclidean distances on the torus.
- **Smoothing:** s_i = Σ_j w_ij φ_X,j / Σ_j w_ij, with w_ij = exp(−d_ij²/(2σ²)), σ = 0.4, cut off at d = 4σ.
- **Threshold:** θ = −0.5·|A*_s(T)|, where A*_s(T) is the smoothed depth of the T particle (mean over its 4 placements).
- **Candidates:** nodes with s_i < θ and s_i ≤ s_j for every node j within distance 1.0. A tie goes to the lower node
  index.
- **Particles:** candidates within 1.5 (about λ₀/2) of each other are merged, by single linkage. Each cluster is one
  particle, located at its deepest node. Its centre is refined to the weighted centroid of the nodes within 0.6 of that
  node, with weights max(θ − s_i, 0).
- **Linking between checks:** a particle is linked to the nearest particle at the previous check within 1.0. Otherwise it
  is a new identity.
- **Radial profiles:** mean φ_S of the nodes in each ring around the centre. Bins are [0, 0.5), then width 0.25 out to
  r = 8. An empty bin is left out of every norm, and the number of empty bins with r < 6 is reported for each run.
- **λ:** zero crossings come from linear interpolation between consecutive non-empty bin means of φ_X; λ is 2 × the mean
  spacing of the crossings in 1 < r < 7. Fewer than 2 crossings there make λ *missing* for that run.

### 5.2 Outcome classes (checks at t = 200, 210, …, 300), applied in this order of precedence
1. *multiplied*: two or more particles at any check;
2. otherwise *decayed*: no particle at t = 300;
3. otherwise *single particle*: exactly one particle at every check, with one identity throughout;
4. otherwise *unresolved* (for example, no particle at an intermediate check, or a broken link).

**Survival is not stationarity (review item 5).** Each run also reports max |∂φ/∂t| at t = 300 (from the ODE right-hand
side). A run is labelled *stationary* if that is < 10⁻⁸ (the 2D MG1 value was 1.4 × 10⁻¹⁰). The class does not depend on
the stationary label. The panel counts of *single particle and stationary* are reported alongside C1.

### 5.3 Criteria
**C1 — survival on the benchmark panel.** For each family (D and R separately), count the *single particle* runs.
- PASS if the panel total is ≥ 28 of 32 **and** every graph has ≥ 6 of 8;
- FAIL if the panel total is ≤ 20 of 32 **or** some graph has ≤ 3 of 8;
- INCONCLUSIVE otherwise.
These are benchmark counts. In draft 2, "≥ 28 of 32" was wrongly presented as a 95% Clopper–Pearson lower bound of 0.75;
the bound at 28 of 32 is 0.710 (review item 1). The per-graph counts, the class counts, the drift |c_final − c| (pinning)
and the stationary counts are reported.

**C2 — the particle's shape matches T's.** Over the family's *single particle* runs, against T's particle (mean of its 4
placements, at the same bin centres):
- (a) relative L2 difference of the radial φ_X profile over r < 6, area-weighted (r dr), non-empty bins only;
- (b) core depth φ_X at the centre bin (smoothed and unsmoothed);
- (c) λ as in §5.1. A *missing* λ counts as not meeting (c).
- PASS if the **median** run meets (a) < 10%, (b) within 10% and (c) within 5%, **and** at least 80% of runs meet (a);
- FAIL if the median run misses any of the three;
- INCONCLUSIVE otherwise, or if fewer than 8 runs are *single particle*.
The per-graph medians are reported. Why 10%: on the 1D ring the h = 0.2 lattice particle is 8% off the continuum, so the
tolerance lets disorder cost about as much as discretization already does.

**C3 — the Laplacian at the soliton's scale (review item 2).** Draft 2's single Weyl count is narrowed and joined by an
operator test. C3 has two parts, and both must PASS.

- **C3(a) — spectral-count benchmark** (intrinsic, so it carries over to MG4b).
  - Compute the eigenvalues μ_m of the graph Laplacian up to 1.5 k*², where k* = 2π/λ₀ ≈ 2.04. That is about 200
    eigenvalues, found by shift-invert Lanczos.
  - The reference is the **exact finite-torus count** N_c(k) = #{q = (2πm/L_x, 2πn/L_y) : |q| ≤ k}, not Weyl's
    asymptotic formula. On the square 20 × 20 torus N_c(k*) = 137 while Weyl gives 132.5; on the 19.92 × 20.08 torus N_c(k*) is also 137.
  - Statistic: the mean of N_graph(k²)/N_c(k) over 21 equally spaced k in [0.9 k*, 1.1 k*]. Averaging over the band
    smooths out the staircase.
  - PASS if the statistic is within 1 ± 0.05 (about ±2.5% in k) for every panel graph. FAIL if the panel mean is outside
    that range. INCONCLUSIVE otherwise.
  - The same statistic on T is reported. T must pass as a positive control.
  - It **measures only the count of modes**, not isotropy or localization.
- **C3(b) — action on the soliton, coarse-grained.**
  - For each placement, take f = the MG1 continuum profile of each field S ∈ {G, X, Y} interpolated onto the nodes (as
    for the initial condition). Compare the graph's (A − D_deg)f/h² with the analytic continuum ∇²f = f″ + f′/r, both
    evaluated at the nodes.
  - Smooth both with the §5.1 kernel and compute e_S = ‖S(L_G f) − S(∇²f)‖ / ‖S(∇²f)‖ over nodes with r < 6.
  - PASS if the panel median of e_S is ≤ 10% for each of G, X and Y. FAIL if any field's median exceeds 20%.
    INCONCLUSIVE otherwise.
  - Reported with it: the values on T (its discretization error), per-graph medians, and the unsmoothed e_S. The
    unsmoothed values are expected to be large on D.
- **Descriptive (no criterion): plane-wave spectra.**
  - Uses **allowed periodic wavevectors only**: every q with |q| ∈ [0.95 k*, 1.05 k*]. Arbitrary directions would add
    seam discontinuities.
  - For each q: the distribution of e_q's weight over the computed eigenvectors (its median against |q|² and its width),
    plus the weight above the 1.5 k*² cutoff.
  - The spread of the medians across directions is the isotropy measure. T's width and spread should be near 0.
- **Discriminating power (N3):** C3(a) must FAIL on both N3 graphs, whose spectra have no small eigenvalues apart from
  the constant mode. If N3 passed C3(a), C3 would be INCONCLUSIVE for every family. C3(b) is not defined on N3, which has
  no embedding.

### 5.4 Implementation checks (every run; not criteria)
- The homogeneous state is an exact fixed point on any graph, since the rows of A − D_deg sum to 0. Unseeded runs stay at
  round-off level up to t = 300.
- **Linear stability needs no new computation.** σ depends on the graph only through each Laplacian eigenvalue μ ≤ 0, and
  MG001 found max Re σ < 0 (−0.195) for every μ = −k² ≤ 0 at b = 29. So the homogeneous state is linearly stable on every
  graph in MG4 (Othmer–Scriven). A 10⁻³ random perturbation must decay (one run per family).
- **Node-sum balances, with their time-derivative terms (review item 5).** Let I_S = Σ_i φ_S,i. For a symmetric
  Laplacian with q = p = 1 and g = 0.1, these hold at all times:
  - dI_G/dt = −I_G + 0.1 I_X;
  - d(I_G + I_X + I_Y)/dt = −0.9 I_X.
  - At t = 300 in every run, both residuals, with the dI/dt terms summed from the ODE right-hand side, must be < 10⁻⁹
    relative to the sum of the magnitudes of their terms. This tests conservation in the implementation.
  - The zero sums I_G = I_X = 0 are expected **only for stationary runs**. For those, |I_S| / Σ|φ_S| is reported, and it
    must be < 10⁻⁶. For non-stationary runs, I_S is reported together with the derivative terms, as in MG1.

## 6. Deferred: survival under token noise (draft 1's Part B)
Removed by the PI's scope decision. Draft 1 proposed stochastic runs on T and D at matched noise per unit area
(Ω_T = 1000), 8 runs each up to T = 100. It came with graph-kernel validations: bit-identity with `tau_leap_euler_1d` on a
cycle, the mean map, exact SSA on a small graph, and a clip pilot. It cost about 40 core-hours. A later study can take it
up as written, using Part A's stationary particles as its starting states and reference. Part A saves every final field,
so this stays possible without re-running.

## 7. Gate outcome
- **PASS:** the positive controls hold (T inside the window, the box check, and C3(a) and C3(b) on T), N3 fails C3(a) as
  predicted, and C1, C2 and C3 PASS for both D and R.
- **FAIL:** any of C1–C3 FAILs for D or R. Each family's results are reported separately.
- **INCONCLUSIVE:** a positive control fails, N3 passes C3(a), any criterion is INCONCLUSIVE, or a budget stop.
- The h = 0.3 check, the seeded-formation runs and the plane-wave spectra are reported and enter no criterion.

**Prediction (stated before any run; not a criterion).** D and R pass C1–C3. Disorder pins the particle: it drifts less
than one λ₀ from its starting centre. N3 fails C3(a). Least certain are C2 and C3(b) on D: its per-node diffusivity varies
by 72%, and its disorder correlation length (about h) is only about 15× smaller than λ₀, so homogenization at the core
scale is not guaranteed.

Pre-stated reading of failures (softened per the review; each is a possible reading, not a diagnosis):
- C1 fails: irregular connectivity changes whether a single particle persists from this start. A local shift of the
  existence window is one possible cause among others (pinning-induced splitting, different basins). Separating them
  would need a follow-up study, for example the window measured on D.
- C1 passes and C2 fails: stationary localized solutions of graph Model G exist, but they differ measurably in shape from
  the comparator's. They are still solutions of Model G on that graph. What fails is agreement with the regular-lattice
  particle.
- C3 fails: the graph operator departs from ∇² at the soliton scale by more than the registered tolerances, in mode count
  (a) or in coarse-grained action (b). Whether any survival rests on continuum-like behaviour is then left open.

## 8. Implementation, streams and records
- New code (tests first): `src/graphs.py` (builders, periodic Delaunay by 3×3 tiling, corrector, spectra),
  `src/graph_ode.py`, `src/measure_graph.py`, the drivers
  `src/run_mg4.py` and `src/evaluate_mg4.py`, and `tests/test_mg4.py`.
- Streams: `runner.seed_for` with gate id 5. Stage index: Part 0 graphs 0, Part A centres 1, timing 6. b_idx 0; om_idx = family (T 0, D 1, R 2, N3 3); replicate = graph × 100 + placement index (pilot graphs 90–93,
  replacement draws 50 onward).
- Records: every Part A final field (float32) and its outcome class, and the fields at the 11 check times.
  Part 0 writes `results/MG4/part0/decisions.json` (h and D_eff tensor per graph, torus, node counts, degree maxima,
  replacements, T window, box check).
  The Part A driver refuses to start without it.
- Everything a choice leaves open is appended to the registration, dated, before the run it affects.

## 9. Budget (run caps bind; core-hour ceilings are stop triggers checked before every batch)
| Item | Runs (cap) | Projected core-h | Ceiling |
|---|---|---|---|
| Part 0 (graphs, corrector, spectra, T window bisection, box check) | about 30 ODE runs | 3 | 8 |
| Part A at h = 0.2 (T 4 + D 32 + R 32, both initial conditions; N3 4) | 140 | 10 | 25 |
| Part A h = 0.3 check (T 4 + D 8) | 12 | 1 | 3 |
Total projected about 14 core-hours. These are rough: the 2D graph ODE's cost (sparse BDF, about 50 000–70 000 unknowns
for D and R at h = 0.2) has not been measured, so the first Part 0 runs time it and the projection is re-stated, in the
registration, before Part A starts. All of it is deterministic and can run in the cloud container or on the PI's machine.

## 10. Design-stage numbers (graph statistics only; scratch, not committed results)
Periodic Poisson sets on a 20 × 20 torus, 3 realizations each (computed for draft 1; the torus is now 19.92 × 20.08):
- Delaunay: mean degree 6.00, minimum 3, maximum 12–13. Per-node second moment Σ_j|Δx|²: coefficient of variation 0.72,
  maximum/mean 8.3. At the naive h, ρ·h² = 2.39 (T: 1.73), so D needs about 1.38× T's nodes. At h = 0.3 that is about
  10 600 nodes (T: 7 700).
- RGG at mean degree 12: minimum degree 0–2, maximum 24–27. Hence "largest component only". At a fixed h it needs about
  3× T's nodes.
- At h = 0.2 on the 20 × 20 torus: T about 17 300 nodes, D about 23 900.

## 11. Not in MG4 (later, separately registered)
- **MG4b:** frozen Wolfram hypergraph snapshots (rule chosen by the PI), graph-intrinsic detection, dimension estimation.
- **W:** D's graph with finite-volume (Voronoi) weights, to separate "irregular points" from "unweighted rule". It
  needs node volumes, weighted conservation and a rules amendment (review).
- Survival under token noise (§6); 3D graphs; nucleation (MG3-type rates) on graphs; the existence window on D; gradient drift and depinning (v = −0.0599 m
  needs a mobile particle, so pinning matters for MG6).
- D1–D3 stay deferred. The h = 0.3 / 0.2 comparison here says something about D2's question in 2D, without replacing it.

## 12. Questions for the PI and the reviewer
Draft 2's questions are settled: the scope is decided by the PI; W is deferred; the reviewer supports homogenized
calibration, the continuum-profile survival test (with seeded formation descriptive) and C2's tolerances.
Open for the review of draft 3:
1. Are C1's panel counts (≥ 28 of 32 with ≥ 6 of 8 per graph; FAIL at ≤ 20 or a graph ≤ 3) and C3(b)'s 10% / 20% bounds
   at the right strictness?
2. Is the band-averaged count ratio over [0.9, 1.1] k* with ±5% a sound version of C3(a)?

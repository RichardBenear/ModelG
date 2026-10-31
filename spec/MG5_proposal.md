# MG5 — A particle on an evolving Delaunay graph (Phase 4, no coupling): design proposal, DRAFT 1

*Stage MG014. Claude, 2026-10-31. Status: **DRAFT, not registered.** No MG5 run of any kind has been made. The only
computation behind this draft is (i) a re-reading of MG4c's stored fields (particle positions, §2) and (ii) a solver
timing and accuracy pilot on an **unregistered** graph, with no particle outcome read (§9, which also discloses one
stopped script). Items marked **[PI]** are decisions for
the PI; items marked **[R]** are points for the reviewer.*

## 1. Question
The roadmap's Gate MG5 (`InitialPlan.md`, Phase 4): *tokens ride on a graph that is being rewritten by a spatial rule;
does the particle persist, drift, or dissolve, and why?* This is the last step before Phase 5 (MG6), where the local G
concentration will modulate the spatial rule's firing rate.

MG5 asks, on 2D Delaunay graphs at a b where the **frozen** graph keeps a single particle (MG4c):

1. **Survival:** does a settled single particle stay a single particle when the graph under it is continually
   rewritten, at rewriting rates from slow to fast compared with the particle's own relaxation?
2. **Mobility:** does it stay where it is (pinned, as on the frozen graph), or does it move? If it moves, how does its
   effective diffusion coefficient depend on the rewriting rate?
3. **Annealing (Part B):** at b = 29, where every frozen D run multiplied (MG4) but the regular lattice T stayed single,
   does fast rewriting restore a single particle? That is, does a graph whose disorder keeps changing behave like the
   averaged (regular) medium rather than like one frozen disordered graph?

Question 2 matters for MG6. LaViolette's gravity picture needs a particle that can move down a G gradient. On frozen D
the particle is pinned (§2), so a mobile particle has to come from somewhere, and an evolving graph is the natural
candidate. Question 3 is the most direct test so far of the idea that particles live in an *emergent* space whose
microstructure never stops changing.

MG5 is deterministic in the fields (Ω → ∞ mean field of rules v1.1, as MG4 and MG4c). The only randomness is the
rewriting itself. Token noise stays deferred.

## 2. What MG5 builds on
- **MG4c (W1):** on unweighted Delaunay graphs, the four D runs are all single for b ≈ 27.27–28.73. **b = 28.0** sits near
  the middle of that range (0.73 from either edge) and inside T's window (27.33–29.20). At b = 28.0, T is an eligible
  threshold reference: θ_28 = ½ × (−8.915163) = **−4.457582** (`results/MG4c/evaluation.json`). All four D runs at
  b = 28.0 are single and stationary (max|dφ/dt| ≤ 8e-13), with smoothed depth ratio to T 0.91–1.02.
- **MG4c stored the float64 final state** (G, X, Y) of every run. MG5 starts from those settled particles, so no time is
  spent on formation and the frozen-graph outcome (ν = 0) is already known exactly.
- **Pinning on frozen D (design-stage reading of MG4c's fields; descriptive).** At b = 27.5, 28.0 and 28.5, each D run's
  particle sits 0.23–0.76 from its placement at t = 200 and at exactly the same place at t = 300 (to 3 decimals). The
  initial condition is a continuum soliton centred on the placement, so the particle slides a short way to a nearby
  favourable spot of the graph and stops. At b = 28.0 the four shifts are 0.45, 0.68, 0.35 and 0.24.
- **MG4c O1:** the cotangent operator keeps D0 p0 single at b = 29. MG5 keeps the unweighted operator (rules v1.1 D_S
  unchanged) for its primary runs; see §10 for the weighted option.

## 3. The spatial rule (rules amendment v1.2, MG5 only)
### 3.1 The rule
Space is the Delaunay triangulation of a set of points on the 19.92 × 20.08 torus. Its triangles are ternary hyperedges;
its edges are the binary relations {u, v} along which D_S hops (rules v1.1 §1, unchanged). Two events fire, independently
of the fields:

| ID | Event | Rate |
|---|---|---|
| S1 (death) | a node is deleted; the triangulation is re-formed over the remaining points | ν per node |
| S2 (birth) | a node is created at a uniformly random position; the triangulation is re-formed | ν ρ per unit area |

ρ = N₀ / A with N₀ = 18 901, the D panel's node count (A = 19.92 × 20.08). "Re-formed" means the Delaunay triangulation of
the new point set. In hypergraph terms each event is a local rewrite: an insertion is a 1→3 Pachner move followed by 2↔2
flips, a deletion is the reverse. The flips are chosen by the Delaunay (empty-circle) criterion, which reads the node
positions.

**Honest note on "Wolfram-style".** Wolfram's models rewrite purely combinatorially. This rule uses node positions to
choose its moves. The purely combinatorial version (uniformly random flips) is known to drive 2D triangulations to a
non-flat geometry of Hausdorff dimension 4 (dynamical triangulations / the Brownian map), which would destroy the 2D
space the particle lives in, and our detector needs an embedding. That route belongs after MG4b (graph-intrinsic
detection). MG5 asks the simpler question first: a 2D space that is locally rewritten but stays statistically flat.

### 3.2 Why this rule: the snapshot ensemble does not change
S1 + S2 is a spatial birth–death process whose stationary law is the Poisson point process at density ρ. MG4's D graphs
are Delaunay graphs of N₀ uniform points, so they are (conditioned on N) samples of that same law. Starting from a D
graph, every later snapshot is again a D-type graph: N(t) is Binomial(N₀, e^(−νt)) survivors plus Poisson(N₀(1 − e^(−νt)))
newcomers, with mean N₀ and standard deviation at most √N₀ ≈ 137 (0.7%).

So **any difference from the frozen outcome is caused by the rewriting itself, not by different graph statistics.** The
rule is also the natural uncoupled baseline for MG6: there, the birth rate would become ν ρ f(φ_G), and a particle's
G-well would change the local node density (the "volume") around it.

### 3.3 What happens to the fields (mean-field transfer rule) **[PI]**
- **Death:** the node's values are removed. No other node's values change.
- **Birth:** the new node's (φ_G, φ_X, φ_Y) are the linear (barycentric) interpolation from the three corners of the
  Delaunay triangle that contains the new point, *before* the insertion.

This is the **interpolating** rule: local concentrations are preserved and the node's tokens come from, or return to,
the reservoirs, as a new cell of vacuum appears or disappears at the local concentrations. Consequences, stated before any
run:
- In the homogeneous state (φ = 0 everywhere) the rule does nothing at all: it creates no background perturbation, so it
  cannot nucleate particles. It perturbs the fields only where they are curved (the error of linear interpolation is
  O(h² ∇²φ)).
- The node sums (total tokens) jump at every event by the inserted or deleted node's content. The §5.4 node-sum
  balances therefore hold between events only, and are checked there (§6).

The alternative is a **token-conserving** rule (a deleted node's tokens are shared among its neighbours, a new node takes
an equal share from its neighbours). In the homogeneous state it changes every neighbour by about 1/6 ≈ 17% of the
background concentration at each event, so it acts as a strong noise source everywhere. That mixes Phase 4's question with
MG3's nucleation question. **Recommendation:** interpolating rule for MG5; the conserving rule is a candidate for MG6's
coupling design, where token bookkeeping is the point.

### 3.4 Fixed from MG4 / MG4c
Rules v1.1 reactions and the eq. 17 parameters except b. Unweighted operator L = (A − D_deg)/h² with **h = 0.2 fixed**
(the D panel's calibrated values are 0.20007 and 0.20003; since ρ is held fixed, h is held at the torus target rather
than re-calibrated on every snapshot). MG4's smoother (σ = 0.4), particle detector and merge radii, built on the current
graph at each check.

## 4. Integration
- **Events** are generated as an exact continuous-time process (event times, types, birth positions and the deleted
  node), from a stream `seed_for(GATE=5, ...)`, before and independently of the field integration.
- **Solver (new for MG5).** MG4's BDF cannot be restarted cheaply: in the timing pilot (§9) each restart cost about 75 s
  for one time unit, against about 2 s per time unit on a frozen graph. MG5 therefore uses a fixed-step, second-order,
  linearly implicit Rosenbrock method (ROS2, Verwer et al. 1999, γ = 1 + 1/√2), with step dt = 0.05 and the Jacobian
  evaluated and LU-factorised (nested-dissection order) once per epoch. ROS2 keeps second order with an approximate
  Jacobian, so reusing it within an epoch does not lower the order.
- **Epochs.** Events are applied at epoch boundaries t_k = kΔ, all events in (t_(k−1), t_k] in their time order (each birth
  interpolates on the point set as it stands after the earlier events of the epoch). Between boundaries the graph ODE is
  integrated on the current graph with ROS2.
- **Epoch length:** Δ = min(1.0, 0.05/ν), so that on average at most 5% of the nodes are replaced per epoch, and Δ divides
  10 (the check interval). That gives Δ = 1, 1 and 0.25 at the three rates.
- **Solver check (Part 0, unregistered graph, frozen):** ROS2 against MG4's BDF from the continuum initial condition at
  b = 28 for 50 time units. Criterion: max|Δφ| ≤ 1e-3 over all three fields (the design pilot gave 1.1e-4 after 20 units,
  with clean second-order convergence: 4.6e-4, 1.1e-4 and 1.4e-5 at dt = 0.1, 0.05 and 0.02). If it fails, dt is halved
  and the check repeated once; a second failure stops for the PI.
- **Splitting check (Part 0, unregistered graph):** at the fastest rate, the same event stream is run with Δ and Δ/2 for
  100 time units. Criterion: particle centres agree to ≤ 0.1 and smoothed depths to ≤ 2% at every check. If it fails, Δ
  is halved for all rates and the check is repeated once; a second failure stops for the PI.

## 5. Runs
**Initial condition:** MG4c's stored float64 final state of each of the four D runs at b = 28.0 (D0 p0, D0 p1, D1 p0,
D1 p1), on its own graph. Rewriting starts at t = 0. Each of the four runs at a given ν uses its own event stream.

**Part A (b = 28.0): survival and mobility.**
- Rates ν ∈ {0.002, 0.02, 0.2} per node per unit time. In units of the core: about 250 nodes lie within 1.5 of the
  particle centre, so the core's graph is renewed on time scales 1/ν = 500, 50 and 5. These bracket the frozen pinning
  relaxation (the slide to the pinning spot finishes within 200 time units in MG4c): slow (quasi-static), comparable, and
  fast.
- Horizon t_end = 500 (ν t_end = 1, 10, 100 renewals per node).
- 3 rates × 4 runs = **12 runs**.
- The frozen control (ν = 0) is MG4c itself: single, stationary, pinned. It is not re-run.

**Part B (b = 29.0): annealing.**
- Same four initial states (settled at b = 28.0), with b raised to 29.0 at t = 0, horizon t_end = 300 (MG4's horizon).
- ν = 0 (frozen control, 4 runs) and ν = 0.2 (4 runs): **8 runs**.
- The frozen control is needed because MG4 multiplied from the continuum initial condition, not from a settled b = 28
  particle.

## 6. Measurement and classes
- **Checks** every 10 time units from t = 10 to t_end, on the graph current at the check (just before that boundary's
  events are applied). Smoothed φ_X; particles with the fixed threshold (θ_28 = −4.457582 in Part A, θ_29 = −5.117354 in
  Part B, both registered in MG4c/MG4).
- **Per-run class** (precedence as MG4: multiplied > dissolved > single > unresolved):
  - *multiplied*: ≥ 2 particles at two consecutive checks. Recorded as **split** if the second particle is first seen
    within 4.0 of the original, otherwise **remote**.
  - *dissolved*: 0 particles at two consecutive checks.
  - *single*: exactly one particle at every check, with each step between checks ≤ 2.0 (a link radius suited to a moving
    particle; MG4's 1.0 was for a fixed one).
  - *unresolved*: anything else (for example isolated checks with 0 or 2 particles).
- **Mobility (single runs only).** Displacement d(t) = |c(t) − c(0)| on the torus, unwrapped along the track.
  - *pinned*: max d < 1.5. The largest frozen-D pinning shift seen in MG4c is 0.76; 1.5 is about half a Turing wavelength.
  - *mobile*: max d ≥ 1.5.
- **Particle diffusion coefficient (descriptive, per rate):** D_p from the time-averaged mean squared displacement of the
  single runs, fitted as MSD = 4 D_p τ over lags 10–100, with the four runs' spread reported. Not used in any status.
- **Shape (descriptive):** smoothed depth at each check against the run's own frozen depth (MG4c), median and range.
- **Balances:** the §5.4 node-sum balances of MG4, checked at the end of each epoch's integration (before events); flags
  handled as in MG4c (mpmath recheck before any class is read).

## 7. Outcomes (pre-registered)
**Gate MG5 (survival, Part A).** Per rate, *survives* means ≥ 3 of the 4 runs are single.
- **PASS:** survives at all three rates.
- **PARTIAL:** survives at some rates and not others. Report the largest surviving rate.
- **FAIL:** survives at no rate.
- Each failing rate is described by its classes (dissolved / split / remote / unresolved), which is the "why".

**Mobility (Part A, registered statuses, no pass/fail).** Per surviving rate: **PINNED** if ≥ 3 single runs are pinned,
**MOBILE** if ≥ 3 are mobile, **MIXED** otherwise. Pre-registered expectation (not a criterion): mobility increases with ν.

**Annealing (Part B).**
- **A1 (annealing restores single):** frozen control multiplies in ≥ 3/4 runs and ν = 0.2 is single in ≥ 3/4.
- **A2 (no restoration):** frozen control multiplies in ≥ 3/4 and ν = 0.2 multiplies in ≥ 3/4.
- **A3 (mixed):** frozen control multiplies in ≥ 3/4, ν = 0.2 neither.
- **A-void:** the frozen control does not multiply in ≥ 3/4 runs. Part B then cannot test annealing; the runs are reported
  only.

**How they read (limited to what was tested):**
- PASS: a particle survives local rewriting of the graph that carries it, for this rule, transfer rule, b and horizon.
  It does not show survival under token noise, under other rules, or under a rule that changes the geometry.
- MOBILE: rewriting unpins the particle; MG6 then has a particle that can respond to a gradient.
- A1: fast-changing disorder acts like an averaged medium at this b. A2: the disorder's effect does not average out at
  ν = 0.2.

## 8. Implementation (fixed in GATES.md before any run)
- `src/evolve.py`: the ROS2 stepper; the event process (pure function of the seed), the epoch splitter, barycentric interpolation with
  periodic wrap, and the Delaunay rebuild (reusing `graphs.delaunay`). Tests first: stationarity of N(t) on a small
  torus; interpolation exact for linear fields; homogeneous state left exactly unchanged by events; determinism (same seed,
  same events, independent of worker count); unwrapped displacement across the torus boundary.
- `src/run_mg5.py`: resumable driver with run cap and core-hour launch cutoff, process pool (`MODELG_WORKERS`), UTF-8 log,
  platform record per run. Stores per check: particle list, the smoothed φ_X minimum, node count; and the field (float32)
  every 50 time units plus float64 at the end.
- `src/evaluate_mg5.py` and an independent `src/recompute_mg5.py` (separate detection code, as for MG4).
- Runs on the PI's PC.

## 9. Budget
**Design-stage timing pilot (disclosed).** Run in the cloud container (2 cores) on an **unregistered** Delaunay graph
(18 901 uniform points from a stream not used by any registered study, seed 999), from the continuum initial condition at
b = 28. No particle detection was run and no outcome was read; only wall times and field differences between solvers.
(An earlier pilot script that started from MG4c's registered D0 p0 field was stopped after a few minutes, when Claude noticed it used a registered graph. Its output file was never opened
and has been deleted.)
- Delaunay rebuild about 1.0 s, nested-dissection order 0.05 s, one LU of the 56 703-unknown system 1.5 s.
- BDF restarted every time unit: about 75 s per time unit (about 450 function evaluations per restart). Not viable.
- ROS2 at dt = 0.05: about 0.8 s per time unit of stepping plus about 2.6 s per epoch (LU, Delaunay, interpolation).

**Projected cost (cloud figures; the PI's PC was about 20% faster in MG3):**
| Item | Runs | Per run | Total |
|---|---|---|---|
| Part 0: solver check, splitting check, unit tests | — | — | about 1 core-h |
| Part A, ν = 0.002 and 0.02 (Δ = 1, t = 500) | 8 | 0.47 h | 3.8 core-h |
| Part A, ν = 0.2 (Δ = 0.25, t = 500) | 4 | 1.6 h | 6.2 core-h |
| Part B, frozen (t = 300) | 4 | 0.1 h | 0.4 core-h |
| Part B, ν = 0.2 (t = 300) | 4 | 0.93 h | 3.7 core-h |
| **Total** | **20** | | **about 15 core-h** |

On the PI's PC with 9 workers this is about 2–3 hours of wall time. Registered as a **run cap of 20** plus a **30 core-h
launch cutoff** (in-flight runs finish), with core-hours reported. The projection is restated after Part 0 measures the
real per-epoch cost on the PI's PC, before Part A starts.

## 10. Questions for the PI and the reviewer
1. **[PI] Transfer rule:** interpolating (recommended) or token-conserving (§3.3).
2. **[PI] Part B (annealing at b = 29):** include (recommended; 8 runs) or drop.
3. **[PI] Weighted operator:** MG5 keeps the unweighted operator. Adding a cotangent-operator copy of Part B (4 more
   runs) would show whether annealing and the weighted operator do the same job. Not recommended now; it can follow if
   Part B gives A1.
4. **[R]** Is the 1.5 pinned/mobile radius sensible, given frozen shifts ≤ 0.76 and λ ≈ 3.1?
5. **[R]** Is the epoch splitting (Δ ≤ 0.05/ν, with a Δ/2 check) adequate, and is the switch from BDF to fixed-step ROS2
   (with the Part 0 solver check) acceptable?
6. **[R]** Are 4 runs per rate enough for the 3/4 rules, given that the four runs share two graphs and start settled?

## 11. Not in MG5
Token noise; combinatorial (non-geometric) rewriting and Wolfram hypergraph snapshots (MG4b); growing or shrinking space
(net birth ≠ death); coupling of the rule to the fields (MG6); 3D.

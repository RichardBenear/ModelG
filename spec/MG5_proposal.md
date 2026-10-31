# MG5 — A particle on an evolving Delaunay graph (Phase 4, no coupling): design proposal, DRAFT 2

*Stage MG014b. Claude, 2026-10-31. Status: **DRAFT, not registered.** No MG5 run of any kind has been made. Draft 1 is
preserved as `spec/MG5_proposal_v1.md` (superseded). Review: `review/MG014_chatgpt.md` (with the reviewer's checks in
`review/MG014_replication/`, to be committed by the PI). Response: `review/MG014_response.md`.*

*Computation behind this draft, all of it disclosed: (i) a re-reading of MG4c's stored fields (particle positions, §2);
(ii) a solver timing and accuracy pilot on an **unregistered** graph, with no particle outcome read (§10, which also
discloses one stopped script); (iii) for draft 2, the right-hand side of the four stored MG4c b = 28 states evaluated at
two values of h, and the steady states at b = 28 and 29 (no integration; §5).*

**PI decisions already taken (2026-10-31 21:14 UTC):** interpolating transfer rule; Part B included; unweighted operator
only. Draft 2 keeps all three.

## 1. Question
The roadmap's Gate MG5 (`InitialPlan.md`, Phase 4): *tokens ride on a graph that is being rewritten by a spatial rule;
does the particle persist, drift, or dissolve, and why?* This is the last step before Phase 5 (MG6), where the local G
concentration will modulate the spatial rule's firing rate.

MG5 asks, on 2D Delaunay graphs at a b where the **frozen** graph keeps a single particle (MG4c):

1. **Survival:** does a settled single particle stay single when the graph under it is continually rewritten, at
   rewriting rates from slow to fast compared with the particle's own relaxation?
2. **Mobility:** does it stay where it is (pinned, as on the frozen graph), or does it move? If it moves, how does its
   effective diffusion coefficient depend on the rewriting rate?
3. **Recovery at b = 29 (Part B):** at b = 29, where every frozen D run multiplied (MG4) but the regular lattice T stayed
   single, does fast rewriting give a single particle?

Question 2 matters for MG6: LaViolette's gravity picture needs a particle that can move down a G gradient, and on frozen D
the particle is pinned (§2). Question 3 is the closest test so far of whether a space whose microstructure keeps changing
behaves differently from one frozen disordered graph.

MG5 is deterministic in the fields (Ω → ∞ mean field of rules v1.1, as MG4 and MG4c). The only randomness is the
rewriting. Token noise stays deferred.

## 2. What MG5 builds on
- **MG4c (W1):** on unweighted Delaunay graphs, the four D runs are all single for b ≈ 27.27–28.73. **b = 28.0** sits near
  the middle (0.73 from either edge) and inside T's window (27.33–29.20). At b = 28.0, T is an eligible threshold
  reference: θ_28 = ½ × (−8.915163) = **−4.457582** (`results/MG4c/evaluation.json`). All four D runs at b = 28.0 are single
  and stationary (max|dφ/dt| ≤ 8.6e-13 under their own operators).
- **MG4c stored the float64 final state** (G, X, Y) of every run, which MG5 uses as its starting states.
- **Pinning on frozen D (design-stage reading of MG4c's fields; descriptive).** At b = 27.5, 28.0 and 28.5, each D run's
  particle sits 0.23–0.76 from its placement at t = 200 and at the same place at t = 300 (to 3 decimals). At b = 28.0 the
  shifts are 0.45, 0.68, 0.35 and 0.24.

## 3. The spatial rule (rules amendment v1.2, MG5 only)
### 3.1 Events
Space is the Delaunay triangulation of a point set on the 19.92 × 20.08 torus. Its triangles are ternary hyperedges; its
edges are the binary relations {u, v} along which D_S hops (rules v1.1 §1, unchanged). Two kinds of event occur,
independently of the fields:

| ID | Event | Rate |
|---|---|---|
| S1 (death) | a node is deleted | ν per node |
| S2 (birth) | a node is created at a uniformly random position | ν ρ per unit area |

ρ = N₀ / A with N₀ = 18 901, the D panel's node count. Each node carries a persistent integer id (MG4c's node index for
the starting nodes, then N₀, N₀ + 1, … in birth order).

**Event list.** For each run, the event list (times, types, birth positions, and the id of each deleted node) is generated
in advance as an exact continuous-time process from a stream `seed_for(5, …)`, independently of the fields. A death picks
uniformly among the nodes alive at its time in the continuous-time process. The event list is saved with the run.

In hypergraph terms, an insertion is a 1→3 Pachner move followed by 2↔2 flips, a deletion the reverse, with the flips
chosen by the Delaunay (empty-circle) criterion, which reads node positions.

**Honest note on "Wolfram-style".** Wolfram's models rewrite purely combinatorially. This rule uses node positions to
choose its moves. Uniformly random combinatorial flips are known to drive 2D triangulations to a non-flat geometry of
Hausdorff dimension 4 (dynamical triangulations / the Brownian map), which would destroy the 2D space the particle lives
in. That route belongs after MG4b (graph-intrinsic detection). MG5 asks the simpler question first: a 2D space that is
locally rewritten but stays statistically flat.

### 3.2 Ticks: how events are applied (the registered update rule)
Events are applied **synchronously in ticks** of length Δ (§4). At the tick boundary t_k = kΔ, all events with times in
(t_(k−1), t_k] are applied at once:
1. Every **birth** of the tick gets the linear (barycentric) interpolation of (φ_G, φ_X, φ_Y) from the three corners of
   the triangle containing it **in the tick-start mesh and field** (the triangulation and field at t_k, before any event of
   this tick).
2. A node born and deleted within the same tick is never added.
3. The tick's deleted nodes are removed.
4. The Delaunay triangulation of the resulting point set is built once.

This is a deliberate definition. Interpolating each birth on the mesh as it stands after the earlier events of the same
tick would need a local Delaunay update per event: at ν = 0.2 that is about 3.8 million events per run, and a global
rebuild per event (about 1 s) would cost about 1 000 hours per run (review item 1). The registered rule is therefore the
tick-synchronous rule above, with Δ part of the rule. **The tick check (§6) measures how much the results depend on Δ.**
The Δ → 0 limit is the continuous-time rule.

### 3.3 Why this rule: the snapshot ensemble does not change
S1 + S2 is a spatial birth–death process whose stationary law is the Poisson point process at density ρ. MG4's D graphs
are Delaunay graphs of N₀ uniform points, i.e. (conditioned on N) samples of that law, so every later snapshot is again a
D-type graph. N(t) is Binomial(N₀, e^(−νt)) survivors plus Poisson(N₀(1 − e^(−νt))) newcomers: mean N₀, standard
deviation at most √N₀ ≈ 137 (0.7%). Ticking does not change which points exist at the tick boundaries.

Any difference from the frozen outcome is therefore caused by rewriting plus field transfer, not by different graph
statistics. The rule is also the natural uncoupled baseline for MG6, where the birth rate would become ν ρ f(φ_G).

### 3.4 Field transfer (PI decision: interpolating)
- **Death:** the node's values are removed. No other node's values change.
- **Birth:** the new node's values are the barycentric interpolation in the tick-start mesh (§3.2).

Consequences, stated before any run:
- In the homogeneous state (φ = 0 everywhere) the rule changes nothing, so it cannot nucleate particles by itself.
- It preserves constant and linear fields but systematically changes curved profiles (linear interpolation of a convex
  profile overshoots, error O(h² ∇²φ)). Repeated remapping can therefore affect a particle independently of the changing
  connectivity (review item 4). §8 limits the readings accordingly.
- The node sums (total tokens) jump at every tick by what was removed and added. The jump accounting of §7 records this
  exactly.

### 3.5 Fixed from MG4 / MG4c
Rules v1.1 reactions and the eq. 17 parameters except b. Unweighted operator L = (A − D_deg)/h². MG4's smoother (σ = 0.4,
cutoff 1.6, distance-based), detector and merge radii, applied on the current node positions at each check.

**Spacing h (changed in draft 2, review item 2).** Each run keeps **its panel graph's calibrated h for its whole run**:
h_D0 = 0.20007382779836555 and h_D1 = 0.20002762164692553 (`results/MG4/part0/decisions.json`). Draft 1's h = 0.2 is
dropped. Reason: the stored MG4c states are equilibria of their own operators (max|dφ/dt| 6.3e-13 to 8.5e-13), but under
h = 0.2 the same states have max|dφ/dt| = 0.049, 0.049, 0.018 and 0.018 (checked for draft 2, confirming the reviewer's
figures). Since ρ is fixed, the panel h is as good an estimate of the rewritten graphs' homogenized spacing as 0.2, and it
keeps t = 0 an exact equilibrium of the frozen system. The matched frozen controls (§5) test this directly.

## 4. Integration
- **Solver.** MG4's BDF cannot be restarted cheaply: in the timing pilot (§10) each restart cost about 75 s per time unit,
  against about 2 s per time unit on a frozen graph. MG5 uses the fixed-step, second-order, L-stable linearly implicit
  Rosenbrock method ROS2 (Verwer et al. 1999; γ = 1 + 1/√2), which keeps second order with an approximate Jacobian.
  - Step **dt = 0.05**.
  - The Jacobian is evaluated at the tick-start state and LU-factorised (nested-dissection order) **once per tick**, and
    reused for every step in the tick.
  - On a frozen run (ν = 0) the tick structure is kept (Δ = 1), so the Jacobian is refreshed every time unit.
- **Tick length:** Δ = 1.0 for ν ≤ 0.04, and Δ = 0.2 at ν = 0.2. So νΔ ≤ 0.04: on average at most 4% of nodes are replaced
  per tick. **Alignment:** every Δ used (1, 0.2, and the checks' 0.1) is an integer multiple of every dt used (0.05, 0.025),
  and divides the check interval 10. No partial steps occur; the number of LU factorisations equals the number of ticks.
- **Numerical failure:** any non-finite value, or max|φ| > 1e3, stops the run with status FAILED-NUMERICAL (§7).

## 5. Runs
**Initial condition (Part A):** MG4c's stored float64 final state at b = 28.0 of each of the four D runs (D0 p0, D0 p1,
D1 p0, D1 p1), on its own graph, with its own h. Rewriting starts at t = 0. Each run at a given ν has its own event list.

**Part A (b = 28.0): survival and mobility.** Horizon t_end = 500.
- **Frozen controls (new in draft 2):** ν = 0, 4 runs, with MG5's solver, h and horizon.
- **Rewriting:** ν ∈ {0.002, 0.02, 0.2} per node per unit time, 4 runs each. About 250 nodes lie within 1.5 of a particle
  centre, so the core's graph is renewed on time scales 1/ν = 500, 50 and 5 (slow, comparable with, and fast compared with
  the pinning relaxation seen in MG4c). ν t_end = 1, 10 and 100 renewals per node.
- 4 + 12 = **16 runs**.

**Part B (b = 29.0): recovery.** Horizon t_end = 300 (MG4's).
- **Parameter change at t = 0 (defined in draft 2, review item 2).** Absolute concentrations are kept: the stored
  potentials are shifted by the change of steady state, φ_S → φ_S + S₀(28) − S₀(29). With the eq. 17 parameters G₀ and X₀
  do not depend on b (both 15.5556) and Y₀ rises from 1.8 to 1.8642857, so only **φ_Y is shifted, by −0.0642857** (exactly
  −0.9/14). The particle starts as the settled b = 28 particle in a medium whose b has just changed.
- ν = 0 (frozen control, 4 runs) and ν = 0.2 (Δ = 0.2, 4 runs): **8 runs**.

**Total: 24 runs.**

## 6. Part 0 checks (unregistered graph; before any registered run)
All on a Delaunay graph of 18 901 uniform points from a stream used by no registered study, with its own calibrated h, at
b = 28, started from the continuum soliton and first settled frozen (ROS2) to t = 300. The checks compare fields and
detections. They do not count toward any gate.

1. **Solver check (frozen):** ROS2 (dt = 0.05) against MG4's BDF from the continuum initial condition for 50 time units.
   Criterion: max|Δφ| ≤ 1e-3 over all three fields. (The design pilot gave 1.1e-4 after 20 units, with clean second-order
   convergence: 4.6e-4, 1.1e-4 and 1.4e-5 at dt = 0.1, 0.05 and 0.02.)
2. **Step check (exact replay):** from the settled state, ν = 0.2, Δ = 0.2, for 100 time units, with the **same event list
   and the same ticks**, run at dt = 0.05 and dt = 0.025.
3. **Tick check (exact replay of events, different synchrony):** the same event list at Δ = 0.2 and Δ = 0.1 (both at dt =
   0.05).

**Comparison rule for checks 2 and 3**, at every check time (every 10 units):
- the particle counts must agree at every check;
- while both runs are single: centre distance ≤ 0.1 and smoothed-depth difference ≤ 2%;
- if either run reaches *multiplied* or *dissolved* (§7), the other must reach the same class, with the first
  qualifying check no more than one check (10 units) apart.

**Policy:**
- Failure of check 1 or 2: dt is halved for all runs and the failed check repeated once. A second failure stops for the PI.
- Failure of check 3: Δ is halved at ν = 0.2 (Δ = 0.1 in Parts A and B, so about twice the cost of those runs) and the
  check repeated as Δ = 0.1 against 0.05. A second failure stops for the PI.

Each check's outcome, and any change of dt or Δ, is written to `results/MG5/part0/decisions.json` and to GATES.md before
Part A starts.

## 7. Measurement, statuses and classes
**Checks:** every 10 time units from t = 10 to t_end, on the state at the tick boundary, before that boundary's events.
Smoothed φ_X on the current positions; particles with the fixed thresholds θ_28 = −4.457582 (Part A) and θ_29 = −5.117354
(Part B), both registered in MG4c/MG4.

**Run status (applied before any scientific class):**
- **COMPLETE:** reached t_end, no numerical failure, every balance and jump-accounting check passed (or a flag was resolved
  by the recorded precision recheck, as in MG4c).
- **FAILED-NUMERICAL:** non-finite values, max|φ| > 1e3, or an unresolved balance or jump-accounting flag.
- **INCOMPLETE:** not run or not finished (run cap, launch cutoff, crash). A crashed run is rerun after investigation and
  the crash reported, as in MG4c.
Only COMPLETE runs are classified.

**Per-run class** (precedence multiplied > dissolved > single > unresolved):
- *multiplied*: ≥ 2 particles at two consecutive checks. Recorded as **split** if the second particle is first seen
  within 4.0 of the original, otherwise **remote**.
- *dissolved*: 0 particles at two consecutive checks.
- *single*: exactly one particle at every check, with each step between consecutive checks ≤ 2.0.
- *unresolved*: anything else.

**Mobility (single runs).** Displacement d(t) = |c(t) − c(0)| from the settled starting centre (detected at t = 0 with
the same detector), unwrapped along the track.
- *pinned*: max d < 1.5. *mobile*: max d ≥ 1.5. (1.5 is a declared descriptive threshold: frozen-D pinning shifts are ≤
  0.76 and λ ≈ 3.1.)
- **D_p (descriptive):** time-averaged MSD = 4 D_p τ over lags 10–100, per run, with the four runs' spread.

**Shape (descriptive):** smoothed depth at each check against the run's own frozen control at the same time.

**Balances and jump accounting.**
- Between events: MG4's §5.4 node-sum balances at the end of each tick's integration (float64 check, mpmath recheck of
  any flag before any class is read, as in MG4c).
- **At each tick (new in draft 2, review item 5):** for each species, the node sums I_before (just before the tick's
  events), the removed sum R (deleted nodes' values) and the added sum B (births' interpolated values), and I_after. The
  identity I_after = I_before − R + B is checked to 1e-10 (1 + |I_before| + |R| + |B|).
- **Interpolation audit:** at every tick the driver recomputes 20 of the tick's births (or all, if fewer) with an
  independent interpolator (`scipy.interpolate.LinearNDInterpolator` on a periodically padded point set) and records the
  maximum difference; above 1e-9 is a flag.
- Counts of births, deaths and same-tick birth-deaths are recorded per tick.

**Saved data (review item 5), per run:**
- The event list.
- At **every check**: the positions (float64) and persistent ids of all nodes with smoothed φ_X < −1.0, their smoothed
  φ_X (float64), the node count, and the global minimum of the smoothed φ_X with its position. Detection at any threshold
  θ ≤ −1.0 uses only nodes below θ (a candidate is below θ, and so is every node that enters its centroid weights), so these
  records are sufficient to recompute every class, track and mobility status exactly, without the run.
- Every 50 time units and at t_end: the full state, positions (float64), ids and smoothed φ_X (float32), so that a
  check's records can be reconciled with the full field and the mesh rebuilt.
- At t_end: full float64 state.
- Per tick: the jump-accounting row and the interpolation-audit maximum.
- Estimated size about 15–25 MB per run, about 0.5 GB in all, compressed. **[PI]** Commit everything (as for MG4c's 111
  MB), or commit the check records and accounting and keep the full-field snapshots on the PC with checksums.

## 8. Outcomes (pre-registered)
**Validity first.**
- Part A is **INVALID-CONTROL** if any Part A frozen control is not *single* and *pinned* with max d < 0.1. Then Part A
  is reported but not scored, and it stops for the PI.
- A rate is **NOT SCORED** unless all 4 of its runs are COMPLETE.
- The gate is **INCOMPLETE** unless all three rates are scored.

**Gate MG5 (survival, Part A).** At a scored rate, *survives* means ≥ 3 of the 4 runs are single.
- **PASS:** survives at all three rates.
- **PARTIAL:** survives at some rates and not others. Report the largest surviving rate.
- **FAIL:** survives at no rate.
- Each failing rate is described by its classes (dissolved / split / remote / unresolved).
- Four runs per rate over a fixed panel of two graphs: the result describes this benchmark. It does not estimate a
  general survival probability.

**Mobility (Part A, registered statuses, no pass/fail).** Per surviving rate: **PINNED** if ≥ 3 single runs are pinned,
**MOBILE** if ≥ 3 are mobile, **MIXED** otherwise. Expectation (not a criterion): mobility increases with ν.

**Recovery (Part B).** Both groups of 4 must be COMPLETE, otherwise Part B is NOT SCORED.
- **A1 (recovery under rewriting):** frozen control multiplies in ≥ 3/4 runs and ν = 0.2 is single in ≥ 3/4.
- **A2 (no recovery):** frozen control multiplies in ≥ 3/4, and ν = 0.2 multiplies in ≥ 3/4.
- **A3 (mixed):** frozen control multiplies in ≥ 3/4, ν = 0.2 neither.
- **A-void:** the frozen control does not multiply in ≥ 3/4 runs. Part B then cannot test recovery; the runs are reported
  only.

**How they read (limited to what was tested):**
- **PASS:** a particle survives local rewriting of the graph that carries it, for this rule, transfer rule, b, tick and
  horizon. It does not show survival under token noise, under other rules, or under geometry-changing rules.
- **MOBILE:** rewriting plus transfer unpins the particle, so MG6 has a particle that can respond to a gradient.
- **A1:** at b = 29, rewriting at ν = 0.2 *with interpolating transfer* gives a single particle where the frozen graph
  multiplies. Two explanations are possible and MG5 does not separate them: changing disorder acting like an averaged,
  more regular medium; or the smoothing effect of repeated interpolation on the particle's profile (§3.4). The result is
  reported as "recovery under rewriting plus transfer", not as annealing.
- **A2:** the frozen graph's replication is not removed by rewriting at ν = 0.2 with this rule.

## 9. Implementation (fixed in GATES.md before any run)
- `src/evolve.py`: ROS2 stepper; event-list generator (pure function of the seed); tick application (§3.2); barycentric
  interpolation with periodic wrap; Delaunay rebuild (reusing `graphs.delaunay`); jump accounting; interpolation audit.
  Tests first:
  - stationarity of N(t) on a small torus;
  - interpolation exact for linear fields, and agreement with the independent interpolator;
  - homogeneous state left exactly unchanged by events;
  - jump-accounting identity on random events;
  - determinism (same seed, same event list, independent of worker count);
  - unwrapped displacement across the torus boundary;
  - ROS2 order on a small linear test problem;
  - a frozen settled state preserved by ROS2 (stationary fixed point).
- `src/run_mg5.py`: resumable driver with run cap and core-hour launch cutoff, process pool (`MODELG_WORKERS`), UTF-8 log,
  platform record per run.
- `src/evaluate_mg5.py`, and an independent `src/recompute_mg5.py` that classifies from the saved check records with
  separate detection code (as for MG4).
- Runs on the PI's PC. The pre-run implementation specification is appended to GATES.md before the first registered run.

## 10. Budget
**Design-stage timing pilot (disclosed).** Run in the cloud container (2 cores) on an **unregistered** Delaunay graph
(18 901 uniform points from a stream not used by any registered study, seed 999), from the continuum initial condition at
b = 28. No particle detection was run and no outcome was read; only wall times and field differences between solvers.
(An earlier pilot script that started from MG4c's registered D0 p0 field was stopped after a few minutes, when Claude
noticed it used a registered graph. Its output file was never opened and has been deleted.)
- Delaunay rebuild about 1.0 s, nested-dissection order 0.05 s, one LU of the 56 703-unknown system 1.5 s.
- BDF restarted every time unit: about 75 s per time unit (about 450 function evaluations per restart). Not viable.
- ROS2 at dt = 0.05: about 0.8 s per time unit of stepping plus about 2.6 s per tick (LU, Delaunay, interpolation).

**Projected cost** (cloud figures; the PI's PC was about 20% faster in MG3). Per-tick overhead 2.6 s, one rebuild per tick
(§3.2); frozen runs skip the rebuild (1.5 s per tick).
| Item | Runs | Per run | Total |
|---|---|---|---|
| Part 0 checks (§6) | — | — | about 2 core-h |
| Part A frozen controls (Δ = 1, t = 500) | 4 | 0.3 h | 1.2 core-h |
| Part A, ν = 0.002 and 0.02 (Δ = 1, t = 500) | 8 | 0.47 h | 3.8 core-h |
| Part A, ν = 0.2 (Δ = 0.2, t = 500) | 4 | 1.9 h | 7.7 core-h |
| Part B frozen (Δ = 1, t = 300) | 4 | 0.2 h | 0.8 core-h |
| Part B, ν = 0.2 (Δ = 0.2, t = 300) | 4 | 1.15 h | 4.6 core-h |
| **Total** | **24** | | **about 20 core-h** |

If the tick check forces Δ = 0.1 at ν = 0.2, those 8 runs cost about twice as much (total about 33 core-h). Registered
as a **run cap of 24** plus a **40 core-h launch cutoff** (in-flight runs finish), with core-hours reported. On the PI's PC
with 9 workers this is about 3–4 hours of wall time. The projection is restated after Part 0 measures the per-tick cost on
the PI's PC.

## 11. Reviewer questions from draft 1 (answered in review MG014)
- The 1.5 mobility radius is acceptable as a declared descriptive threshold from the settled starting centre.
- Four runs per rate are adequate for this fixed-panel benchmark, without claiming a general survival probability.
- ROS2 and the tick approach are acceptable conditionally on the implementation, controls and validation now in §§3–7.

**Remaining PI question:** data storage (§7, last item).

## 12. Not in MG5
Token noise; combinatorial (non-geometric) rewriting and Wolfram hypergraph snapshots (MG4b); growing or shrinking space
(net birth ≠ death); coupling of the rule to the fields (MG6); a weighted-operator arm; separating "annealing" from
"interpolation smoothing" (a possible follow-up only if Part B gives A1); 3D.

# MG6 — G-coupled spatial growth (Phase 5): design proposal, DRAFT 2

*Stage MG019a. Claude, 2026-11-02. Status: **DRAFT 2, not registered.** It revises draft 1 (preserved as
`spec/MG6_proposal_v1.md`) after ChatGPT's review MG019 (`review/MG019_chatgpt.md`; my response
`review/MG019_response.md`). No MG6 run of any kind has been made.*

*PI decisions (Richard Benear, 2026-11-02 14:13 UTC, project thread): adopt ChatGPT's recommendations for all eight §9
decisions of draft 1. They are written into this draft and listed in §10.*

*Computation behind this draft: the static calculations on MG4c's four stored b = 28 start states only (`spec/MG6_design_calc.py`
and ChatGPT's `review/MG019_replication/static_checks.py`). No integration and no MG6 dynamics.*

## 1. Question
MG5 (CLOSED, PASS) showed that a settled particle survives continual field-blind rewriting of the Delaunay graph under it,
and is mobile at ν = 0.2. MG6 couples the rule to the particle: **the local G value modulates the birth rate**. The
first campaign asks one focused question:

> **Does the G-dependent birth rule produce the predicted embedded node-density response around the particle, and does
> the particle survive, or change, under that feedback?**

Not asked here (§9): curvature as such, a gravitational falloff law, and scaling with particle mass.

**What a positive result means (stated now).** With a rule that sets the local target density to ρ₀ f(φ_G), a density
deficit in the G-well is largely built in. A demonstrated response therefore **validates that the implemented coupled
process does what it is meant to do, with the particle present and responding**. It is not a physical discovery by itself.
The informative parts are the feedback on the particle, the size of the response against its quasistatic prediction, and
the contrast with the nulls.

## 2. What MG6 keeps from MG5
Rule v1.2's spatial setting, with tick-synchronous application (Δ = 0.2 at ν = 0.2) and interpolating transfer of φ to
newborn nodes on the pre-edit mesh M_(k−1). Also: ROS2 with dt = 0.05, panel-calibrated h, MG4c's b = 28 start states, MG5's
detection, θ_28, classes, smoothing, balance and jump accounting and interpolation audit, and the MG017 execution settings
(5 workers, one BLAS thread). Everything not stated below is as in MG5.

## 3. The coupled birth–death process (rules amendment v1.3, MG6 only)
### 3.1 Construction (review item 1)
A **dominating candidate process** with independent lifetimes. Per run, with persistent integer IDs:
- **Initial nodes:** each of the N₀ start nodes gets a lifetime ℓ ~ Exp(ν). It dies at time ℓ if ℓ < t_end.
- **Candidates:** a Poisson process on the torus × [0, t_end] at rate ν ρ₀ F per unit area per unit time, with
  ρ₀ = N₀/A and F = f_max (§3.2). Each candidate has a position x (uniform), a time t_b, a mark u ~ U(0, 1), a lifetime
  ℓ ~ Exp(ν) and a fresh ID.
- **Acceptance:** at the tick boundary that applies a candidate's birth, it is **accepted iff u < f(φ_G*(x)) / F**. Here
  φ_G* is the boundary field barycentrically interpolated on the pre-edit mesh M_(k−1), exactly as for the transferred
  value, and evaluated at x (or at x + s for the shifted null, §3.3). An accepted candidate becomes a node, with the
  transferred field. It dies at t_b + ℓ if that is before t_end.
- **Rejected candidates are ghosts.** They never enter the graph, and their deaths are ignored. Their marks and lifetimes
  are still drawn, so the stream does not depend on acceptance.
- **Ephemeral candidates:** a candidate born and dying within the same tick interval never exists at a boundary. As in MG5,
  it is skipped, and no acceptance is evaluated.
- **Tick ordering** is as in MG5: at boundary k, deaths in the interval are applied, then births. Births are accepted on
  M_(k−1), and one Delaunay rebuild follows. Survivors keep their order; accepted births are appended in candidate order.

With f ≡ 1 (decoupled arm) the accepted births form a Poisson process at rate ν ρ₀ with uniform positions. Deaths are at
rate ν per node. This is MG5's rule in law. It is **not** MG5's event list: hash-identical reproduction of MG5's generator
is not claimed, and it is not needed (§3.4).

**Common random numbers.** The initial lifetimes and the candidate stream depend only on the start state (slot, placement),
not on the arm. The five arms of one start therefore share one stream, and differ only in which candidates are accepted.
This pairs the arms for the contrasts in §5.

### 3.2 Coupling function and bound (decisions 1 and 2)
- **F1: f(φ_G) = exp(κ φ_G)**, with κ ∈ {0, 2, 4}. κ = 2 is the **primary response arm**; κ = 4 is a **stronger
  feedback (stress) arm**. Neither the panel nor the ratio of the two responses is a small-κ linearity test. The expected
  core birth factors are about 0.58–0.60 (κ = 2) and 0.33–0.36 (κ = 4). The predicted core-deficit ratio κ4/κ2 is about
  1.70–1.73, not 2.
- **Dominating bound.** φ_cap = 0.10 (the start states have φ_G ≤ 0.049). F = exp(4 · 0.10) = e^0.4 ≈ 1.492 for **all** arms
  (κ = 0 included), so all arms share one stream. The κ = 0 acceptance probability is 1/F ≈ 0.67.
- **Breach policy.** If any acceptance evaluation has f(φ_G*)/F > 1, that is φ_G* > φ_cap·4/κ, the run stops with status
  **INVALID-BOUND**. The record keeps the tick, the position and the value. The run is not rerun with a changed bound, and
  any breach stops the campaign for the PI.

### 3.3 Arms and the randomized null (decisions 3 and 5)
Five arms × four starts = **20 runs**, all at b = 28, ν = 0.2, Δ = 0.2, dt = 0.05, t_end = 200:

| Arm | Coupling |
|---|---|
| DEC | κ = 0 (decoupled; rerun under the MG6 pipeline) |
| SH2, SH4 | κ = 2, 4 evaluated on the **shifted** field φ_G*(x + s) |
| CP2, CP4 | κ = 2, 4, coupled at x |

**Shift s:** one fixed periodic offset per start, shared by SH2 and SH4. It is drawn from its own stream, uniform on the
torus and conditioned on a minimum-image length of at least 6 (rejection sampling, first accepted draw). The induced
density well then sits at least 6 from the particle at the start. The shuffled null is deferred (§9).

**MG5 as an external check (not an arm).** DEC's class and descriptive D_p are compared with MG5's ν = 0.2, b = 28 runs.
Same starts, different streams: this is a compatibility check, not independent evidence.

### 3.4 Seeds
- **Initial lifetimes and candidates:** `seed_for(6, 0, 0, 1, 2·slot + pl)`.
- **Shift:** `seed_for(6, 0, 1, 1, 2·slot + pl)`.
- **Part 0 graph:** `seed_for(6, 9, 0, 0, 0)`. **Part 0 streams:** `seed_for(6, 9, 3, 1, k)`.

The candidate list is a pure function of its seed. It is regenerated and hash-checked in a fresh process before each run,
as in MG5.

## 4. Measurements (review items 2, 3 and 5)
### 4.1 What is measured, and what it is called
The primary quantity is the **embedded node-density deficit** around the particle, adopted as a **volume proxy**. Discs
and annuli are in the embedding coordinates. It is **not** shown to be intrinsic curvature: that would need a separate
metric/volume argument or a graph-distance measurement, which this campaign does not make (§9).

### 4.2 Per-check records (every 10 time units, on the boundary state, before that boundary's edits)
Around the detected particle centre c(t) (MG5 detector, θ_28), using periodic minimum-image distances:
- **Cumulative counts:** N(<R) for R ∈ {1.0, **1.5**, 2.0, 3.0, 5.0, 8.0}. Two baselines:
  - **D₀(R) = N(<R) − ρ₀πR²** (absolute change, reference density ρ₀ = N₀/A);
  - **D_c(R) = N(<R) − ρ̄(t)πR²** (redistribution relative to the current mean, ρ̄(t) = N(t)/A).
- **Annular counts and densities:** annuli with edges {0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0, 8.0}. Densities
  n_j = count_j / area_j, with δn_j relative to both ρ₀ and ρ̄(t).
- **Radial φ_G:** the annulus means of φ_G on the quadrature grid (§4.3).
- **G-well integrals** (decision 4), over the whole torus by quadrature:
  - **M₋ = ∫ max(−φ_G, 0) dA** (well);
  - **M₊ = ∫ max(φ_G, 0) dA** (halo);
  - signed **M = ∫ φ_G dA**.
  M₋ is an operational **G-well integral**, not an established gravitational mass. In the coupled arms ∫φ_G dA need not
  vanish, because Σφ_G = 0 holds on the node measure and the mesh is deliberately non-uniform. All three are reported.
- **Quasistatic predictions** (§4.3): P₀(R) and P_c(R), and the annular predictions.
- The MG5 records: sparse φ_X, centres, classes, depth, node count, and the per-tick accounting.

### 4.3 Quadrature and the prediction
**Coordinate-area quadrature:** a fixed regular grid of spacing 0.05 on the torus (cell area dA = 0.0025). φ_G is
barycentrically interpolated on the current mesh at each grid point. Node-equal-area weights (A/N per node) are **not** used:
the nodes are deliberately sampled non-uniformly, and that would bias the integral.

**Quasistatic prediction** of the coupled density, given the measured φ_G at the check:
- P₀(R) = ∫_disc ρ₀ (f(φ_G^(eval)(x)) − 1) dA;
- P_c(R) = ∫_disc ρ₀ f dA − πR² · (1/A) ∫_torus ρ₀ f dA;
- annular versions likewise.

Here φ_G^(eval) is φ_G at x for CP arms, at x + s for SH arms, and f ≡ 1 for DEC. This is the stationary target for a
**prescribed static field**. For a field that changes, the mean density obeys ∂ₜn = ν[ρ₀ f − n] and depends on field
history, so the prediction is **quasistatic and approximate**. It is also built from the particle's own reshaped field, so
agreement with it does **not** rule out feedback. The ratio observed/predicted is reported descriptively, not scored.

### 4.4 Falloff (descriptive; review item 5)
ΔN(R) is cumulative, and can show a plateau, a zero crossing or cancellation. A fitted power of a cumulative count is not
a density exponent (in 2D, dΔN/dR = 2πR·δn). The falloff analysis therefore uses **annular densities**:
- **Quantity:** for each annulus, the run-level paired contrast δn_j(CP2) − δn_j(SH2), with the ρ₀ baseline. κ = 4 is done
  the same way, as a second descriptive set.
- **Fit range:** annuli from r = 1.0 outward, up to the last annulus inside the first zero crossing of the arm's mean
  radial φ_G profile (the well), keeping the sign as measured, with no absolute values.
- **Minimum signal:** at least 3 annuli in the range must each have a paired-contrast |t| ≥ 2.353 (n = 4, df = 3).
  Otherwise the exponent is reported as **UNRESOLVED** and no fit is chosen.
- **Fits:** if resolved, both an exponential A·e^(−r/λ) and a power law A·r^(−p), by least squares on run-level means with
  run-level SEs. Both are reported, with residuals. Neither is presumed. The halo (outside the zero crossing) is reported
  per annulus, not fitted.
- The predicted annular profile (§4.3) is shown alongside.

## 5. Outcomes and statistics (review item 4)
**Time window:** checks at t = 100, 110, …, 200 (11 checks), provisional, subject to Part 0 (§7, item 6).

**Run level first.** For each run, Y = the window mean of D₀(1.5) around its own detected centre. All inference uses the four
run-level values per arm. Snapshots within a run are not treated as independent evidence.

**Paired contrasts.** The arms of one start share a stream (§3.1), so the contrasts are paired by start i:
- c_i^SH = Y_i(CP2) − Y_i(SH2);
- c_i^DEC = Y_i(CP2) − Y_i(DEC).

**Statistical assumption (registered):** the four run-level paired contrasts are independent and approximately normal. With
n = 4, a one-sided paired t-test at α = 0.05 needs t ≤ −2.353 (df = 3).

**Gate MG6 (primary: density response at κ = 2):**
- **NOT SCORED (feedback):** some run of CP2, SH2 or DEC is not single at every window check, so no centre is defined. The
  feedback outcome is then the result.
- **DEMONSTRATED:** both contrasts (vs SH2 and vs DEC) have all four c_i < 0 and the paired t ≤ −2.353.
- **NOT DEMONSTRATED (inconclusive):** otherwise. This does not establish that there is no deficit. No equivalence
  criterion is registered.

**Also reported, not scored:**
- the same contrasts at κ = 4 (CP4 vs SH4, vs DEC) and with the D_c baseline;
- each contrast's mean, its paired SD, and a 95% t interval, which is explicitly conditional on the assumption above;
- observed/predicted at R = 1.5 per run;
- DEC versus SH (SH2 − DEC, SH4 − DEC). A non-zero particle-centred profile in DEC or SH is **not** treated as invalidity.
  The detected centre depends on fields that respond to the graph, so particle-centred sampling need not be spatially
  unbiased. The paired treatment-minus-null contrast is the test.

**Feedback outcome (per arm, CP2 and CP4):**
- **SURVIVES:** at least 3 of 4 runs single (MG5's classes and rule).
- **CHANGED:** otherwise; the classes are reported (multiplied / dissolved / unresolved).

Also reported descriptively: depth relative to the paired DEC run at matched checks, mobility and D_p, and M₋ over time.
Multiplication or dissolution under coupling is a **physical feedback result**, not "no deficit". Such runs are never
discarded, and averages are never taken over survivors only.

**Invalidity** (reserved for failures of construction or implementation):
- INVALID-BOUND (§3.2);
- an unresolved flag (balance, jump accounting, interpolation audit, as MG5; FAILED-NUMERICAL);
- a candidate-list hash mismatch (INCOMPLETE);
- a Part 0 stop.

**Status precedence per run:** FAILED-NUMERICAL > INVALID-BOUND > INCOMPLETE > COMPLETE. The gate is scored only if all 12
runs of CP2, SH2 and DEC are COMPLETE. Otherwise it is INCOMPLETE, and the PI decides.

## 6. Runs, cost, budget and order (decisions 6 and 8)
- **Rate:** ν = 0.2 only. ν = 0.02 is deferred (§9).
- **Order:** DEC, then SH2 and SH4, then CP2, then CP4. Nulls first, so a cutoff can never remove them. The ordering is
  tested in code (a unit test asserts the queue order) before launch.
- **Cost estimate:**
  - MG5 measured 3.77 ledger core-h for a ν = 0.2 run at t = 500, about 27 s per unit with 5 workers.
  - At t = 200 that gives about 1.5 core-h per run.
  - The candidate stream is about 1.49 times MG5's birth events, and quadrature adds about 20 checks × 160 000 point locations
    per run, small against ROS2.
  - **20 runs ≈ 30 ledger core-h.**
- **Budget:** run cap 24 (20 runs + 4 slots of headroom), and a **provisional 40 core-h launch cutoff**. Both are restated
  after Part 0's full-pipeline timing, which must time **5 concurrent runs** (MG5 lesson). The launch cutoff is not a hard
  spending ceiling. The report gives ledger core-h plus all off-ledger work: Part 0, event generation and replay, and any
  interrupted runs.
- **Wall time:** about 6–7 h with 5 workers on the PI's PC, launched with no console attached (as MG5's relaunch).

## 7. Part 0 (before any registered run, on an unregistered graph; review item 6)
On a Delaunay graph from `seed_for(6, 9, 0, 0, 0)` with its own calibrated h, b = 28, settled frozen to t = 300 (as MG5 Part 0).
Each check has a quantitative criterion and a fallback rule.

1. **Acceptance-rate and determinism check.** Static synthetic field φ_G = −0.25·exp(−r²/2) (no ODE), κ = 4, 20 time units.
   - The accepted fraction must be within 4 binomial SD of the expected mean f/F over the candidates.
   - Two runs with the same seed must give identical accepted bitmaps.
   - Failure means an implementation error: stop for the PI.
2. **Density equilibration.** Same static field, κ = 4, 50 time units, graph process only. Over t ∈ [25, 50], for each of the
   annuli inside r = 3, the time-averaged annular density must agree with ρ₀ f (exact for a static field) within 4
   Poisson-plus-autocorrelation SDs. The SD uses the measured lag-1 autocorrelation of the annular count. Failure: stop
   for the PI.
3. **Step check (fixed edit history).** Coupled κ = 4, ν = 0.2, Δ = 0.2, 100 units: run at dt = 0.05 and record its accepted
   edits. Replay at dt = 0.025 with **the same accepted edits imposed** (no re-evaluation of acceptance), which isolates the
   ODE step error.
   - Criterion as MG5: equal counts at every check, centre ≤ 0.1, depth ≤ 2%.
   - Failure: dt halved for all runs and the check repeated once; a second failure stops for the PI.
4. **Tick check (fully coupled).** κ = 4, Δ = 0.2 against Δ = 0.1, both fully coupled. Acceptances may differ, so the
   comparison has two parts:
   - (a) equal counts at every check, centre ≤ 0.15, depth ≤ 3%;
   - (b) window-mean D₀(1.5) (t = 50–100) differing by at most 0.25·|P₀(1.5)|.
   - Failure: Δ = 0.1 for all MG6 runs (about twice the cost) and the check repeated as 0.1 against 0.05; a second
     failure stops for the PI.
5. **Bound check.** In all Part 0 coupled runs, max φ_G* at candidate positions must stay ≤ 0.07 (70% of φ_cap). If it
   does not, the run stops for the PI, who can raise φ_cap before registration and restate F.
6. **Transient and dependence.** From the tick-check run at Δ = 0.2: the D₀(1.5) series over t = 0–100 and its lag-1
   autocorrelation at 10-unit spacing are reported. If the 20-unit block means over t = 60–100 still drift by more than
   0.25·|P₀(1.5)|, the window start moves to t = 150 and t_end to 250 (cost × 1.25). That is the only window change
   allowed, decided before any registered run.
7. **Timing.** 5 concurrent coupled runs for 20 units each; the budget is restated in GATES.md before the runs.

All Part 0 outcomes, and any change of dt, Δ, φ_cap or window, are written to `results/MG6/part0/decisions.json` and to
GATES.md before any registered run.

## 8. Data and recomputation (review item 6)
**Committed** (small), per run:
- the record (status, flags, replay spec with candidate-list SHA-256, shift s, platform, code revision);
- `checks/<key>.npz` with:
  - MG5's sparse records;
  - centres and classes;
  - N(<R), annular counts and areas;
  - ρ₀, N(t) and ρ̄(t);
  - the radial φ_G profile;
  - M₋, M₊ and M, with the quadrature definition (grid spacing and origin);
  - P₀, P_c and the annular predictions;
  - the per-tick accounting;
- `accepted/<key>.npz`: a packed bitmap over the candidate list, in candidate order. About 1.1 M candidates per run, so
  about 140 kB.

Together with the regenerated candidate list, the bitmap reconstructs the full node history. Also committed: the budget
ledger, driver log, local manifest, `evaluation.json` and `part0/`.

**PC-only,** with SHA-256 manifest: candidate lists, full-field snapshots and final states.

**Recomputation:** an independent script (separate detection, counting, quadrature, contrasts and outcome rules, as
`recompute_mg5.py` now does) recomputes everything in §5 from the committed records. It exits nonzero on any difference.

## 9. Not in MG6 (deferred, recorded)
- **Particle-mass scaling (roadmap Phase 5 objective; decision 4).** A b-scan (for example b = 27.5, 28, 28.5) with the
  M₋ definition is **deferred to a separately registered study**. Varying κ tests coupling strength, not scaling with the
  particle. This defers a roadmap objective, and the record says so.
- **Curvature:** an intrinsic interpretation (metric/volume argument or graph-distance balls).
- **Shuffled null** (decision 3), **ν = 0.02** (decision 6).
- Token noise, combinatorial rewriting, 3D, a weighted operator, coupling of deaths, density-dependent h, motion as a
  gate, and multiway/frame tests (MG7).

## 10. Record of the PI decisions (draft 1 §9; ChatGPT MG019 recommendations adopted)
| # | Decision | Adopted |
|---|---|---|
| 1 | Coupling | F1 = exp(κ φ_G), with a registered dominating bound and a breach policy (§3.2) |
| 2 | κ panel | {2, 4}: κ = 2 primary, κ = 4 feedback/stress; not a linearity test |
| 3 | Randomized null | Fixed spatial shift, matched at both κ; shuffled deferred |
| 4 | Mass / b-scan | M₋ = ∫max(−φ_G, 0) dA as an operational G-well integral, plus M₊ and the signed integral; b-scan deferred to a separate study |
| 5 | Decoupled arm | Rerun under the MG6 pipeline with MG6 streams; MG5 as an external compatibility check |
| 6 | Rates | ν = 0.2 only |
| 7 | Thresholds | R = 1.5 primary, other radii descriptive; t_end = 200, window 100–200 provisional (Part 0 item 6); run-level paired contrasts replace the 3σ rule; NOT DEMONSTRATED replaces NO DEFICIT |
| 8 | Budget | 20 runs, cap 24, provisional 40 core-h launch cutoff restated after Part 0; full accounting |

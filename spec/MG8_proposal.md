# MG8 — Does the order in which the rewrite rules are applied change the outcome? Design proposal, DRAFT 1

*Stage MG036. Claude, 2026-11-05. Status: **DRAFT 1, not registered. No MG8 run of any kind has been made.** Requested by
the PI (project thread, 14:09 UTC: "Let's plan MG8 now ... whether different admissible orders of applying the rewrite
rules lead to measurably different outcomes"). The only computation behind this draft is static arithmetic on the
registered parameters (§2.4); no dynamics were run.*

*Naming: MG8 is the multiway/foliation gate of the roadmap (Phase 6). `InitialPlan.md` called it MG7; it was renumbered
MG8 when the mass-scaling study took MG7 (GATES.md, MG7 registration; STAGE_LOG MG023b). It has been held in reserve
since, and every later proposal (MG9, MG9b, MG9c, MG9d) says so.*

## 0. What MG8 is, confirmed against the charter
- `README.md`, Central Question 5: **"Is SQK's preferred frame a non-confluence effect? Quantify multiway branch
  divergence and foliation dependence."**
- `InitialPlan.md`, Phase 6: "quantify how much branches diverge across update orders. Test whether SQK's preferred frame
  shows up as a measurable foliation dependence. This turns the 'non-confluence = absolute ether' observation into a
  result."
- `spec/rules_v1.1.md` §4: "Continuous-time stochastic simulation selects one weighted path through the multiway system
  ... This is a choice of foliation and is examined in Phase 6."
- GATES.md row MG8: "Multiway branch divergence and foliation dependence", PENDING.

So the PI's reading is right: MG8 asks whether different admissible orders of applying the rules give measurably different
outcomes, and, behind that, whether Model G's preferred (ether) frame could come from that order dependence.

## 1. Question
> **Q (primary): On the 1D lattice of MG2/MG3, does a different admissible update order change what the particle does,
> and if it does, does the difference survive as the time step goes to zero?**

The second half matters for Central Question 5. An order effect that disappears as τ → 0 belongs to the discretisation:
in the continuous-time limit the rules pick no frame of their own. An effect that survives refinement would be a frame
defined by the update order itself, which is what "the preferred frame is a non-confluence effect" would require.

## 2. What is known before any run (Part A: structural, no simulation)
These are arguments, not results. They are written down so that the gate tests something they do not already settle,
and ChatGPT can check them.

**A1. Model G has a preferred frame already in the continuum, before any rewriting.** The equations are
∂ₜφ = D∇²φ + f(φ) with a medium at rest. Under a Galilean change of frame x′ = x − vt they acquire an advection term
v·∇φ, and they are not Lorentz invariant either (first order in time, second order in space). A soliton moving at speed v
relative to the medium is therefore not a solution of the same equations as one at rest; the rest frame of the medium is
singled out by the mean-field equations themselves. This is the ether frame of SQK, and it needs no discreteness and no
non-confluence. MG8 can therefore only ask whether update order adds a *second* frame, or changes the first.

**A2. The token rules are confluent except at resource conflicts.** Rules R1–R5 and hopping act on token counts per
node. Two events that are both enabled and remain enabled after the other fires give the same final counts in either
order, because their count changes add. The only non-commuting ("critical") pairs are those where one event uses up a
token the other needs: for example R4 (Y + 2X → 3X) needs n_X ≥ 2, and an X-consuming event at the same node can disable
it when n_X = 2. At the homogeneous state G₀ = X₀ = 15.6, Y₀ = 1.86 (b = 29), so at Ω = 1000 a node holds about 15,600
G, 15,600 X and 1,860 Y tokens, and conflicts need counts of 0–2. So **branches of the multiway system differ almost
entirely by which events fire (the stochastic choice), not by the order in which the same events fire.** Two caveats:
(i) this is state confluence, which is weaker than Wolfram's causal invariance (isomorphic causal graphs): with
indistinguishable tokens the causal graph is not uniquely defined without a token-identity convention; (ii) it says
nothing about what a scheduler does when its *steps* group events, which is A4 and the gate.

**A3. The spatial rewriting of MG5–MG7 is order-independent by construction.** The Delaunay graph is a function of the
point set alone, so the order of births and deaths within a tick does not change the graph. The registered tick
(`src/evolve.py` `apply_tick`) is synchronous: every birth is interpolated on the pre-edit mesh, so the field values do not
depend on order either. An order effect could only appear if births were interpolated sequentially, and then only for
edits in adjacent cavities in the same tick. The frequency of such pairs can be counted from the stored MG5 event lists
without running any dynamics (Part A3 count, §5). MG6's tick check did find that tick *size* moves the particle centre by
0.2–0.5 (MG021); that is a granularity effect of the synchronous tick, not an ordering effect, and it was already
recorded as a weakness of the tick-check design.

**A4. Directional orders let tokens ride the sweep (static arithmetic, the design calculation).** The registered kernel
`tau_leap_euler_1d` is *synchronous*: every node's exits are drawn from the start-of-step counts and applied together
(Jacobi order). In a *sequential* order, nodes are processed one at a time and a token that arrives at a node not yet
processed is processed again in the same step. If the order sweeps left to right, a token that hops right can hop again,
while a token that hops left cannot. Per step the mean extra displacement is h·r²/(1 − r), with r = Dτ/h² the per-step hop
probability in one direction, so each species is advected at

  v_s(τ) = h·r²/((1 − r)·τ) ≈ D²τ/h³.

| species | D | r at τ₀ = 6.667e-4 | v_s at τ₀ | v_s at τ₀/16 | v_s at τ₀/64 |
|---|---|---|---|---|---|
| G | 1 | 0.0167 | 0.085 | 0.0052 | 0.0013 |
| X | 1 | 0.0167 | 0.085 | 0.0052 | 0.0013 |
| Y | 12 | 0.200 | 15.0 | 0.76 | 0.19 |

(h = 0.2, L = 100.) Y is advected about 175 times faster than G and X, so a directional order is a *differential flow*,
not a uniform shift, and at the production step it is not a small perturbation: a Péclet number for Y over a particle
width of about 2 is 2.5 at τ₀ and falls to about 0.03 at τ₀/64. Every v_s is proportional to τ, so the naive expectation is
that the order effect vanishes linearly with τ. The gate tests that expectation; it does not assume it, because
differential flow can destabilise a Turing pattern and pinning by the mesh (MG9d) can hide a small drift.

**A5. What A1–A4 imply for the outcome.** The likely result is "order matters at the production step, and the effect
vanishes as τ → 0", which would answer Central Question 5 in the negative for the token layer: the preferred frame is the
medium's rest frame (A1), not a non-confluence effect. This expectation is disclosed now so it is not mistaken for a
finding later. The value of the study is to measure the size of the effect at the step every earlier result used, and to
check the scaling rather than assume it.

## 3. Design (Part B: the gate)
Everything not stated here is MG2/MG3 as registered: P17 parameters with b = 29.0, h = 0.2, L = 100 (N = 500, periodic),
D_G = 1, D_X = 1, D_Y = 12, the per-token exit probabilities of `tau_leap_euler_1d`, the deterministic reference soliton
from `results/MG003/ref_profiles.npz`, smoothing, detection, tracking (`measure.track`, x0 track only) and
`runner.seed_for` streams.

### 3.1 Schedulers (the admissible orders)
All four use the same per-token exit probabilities and the same sources; they differ only in the order in which nodes are
processed within a step and so in whether tokens that arrive during a step are processed again.
- **J (synchronous, the registered kernel).** All exits drawn from start-of-step counts; arrivals wait for the next step.
  This is the control and the order used by every earlier gate.
- **R (rightward sweep).** Nodes processed in the order s, s+1, …, s−1 (mod N), with the seam s drawn uniformly each
  step so that no fixed node is privileged. A node's exits are drawn from its counts at the moment it is processed,
  arrivals included; R4's per-Y probability uses that moment's n_X.
- **L (leftward sweep).** The mirror image of R. Run with R's stream and the mirrored seam and node order, L must
  reproduce R's trajectory mirrored node for node. **L is an implementation check, not independent data.**
- **P (random sequential).** As R, but the node order is a fresh uniform random permutation each step. P is
  parity-symmetric in distribution, so it can have no mean drift; it is the isotropic alternative order.

Deterministic versions (Ω → ∞) apply the same orders to expected counts, with R4's mean-field form X²Y. P has no single
deterministic map and is run stochastically only.

### 3.2 Arms
- **B1, deterministic ladder:** J and R at τ_k = τ₀/2^k, k = 0, 1, …, 6 (τ₀ to τ₀/64), plus L at k = 0 and k = 6 as mirror
  checks. 16 runs. Start: the reference soliton at its registered position. T = 300.
- **B2, stochastic at the production step:** J, R and P at τ₀, Ω = 1000, 8 runs each, plus one L mirror run on R's first
  stream. 25 runs. Start: the reference soliton (MG2 "soliton" start, `np.rint(Ω·c)`), T = 500. Streams
  `seed_for(8, 0, arm, 0, rep)`, shared across arms for the same rep (arms differ in order, not in seed).
- Ω = 1000 because MG2 recorded no loss or replication events at Ω ≥ 1000, so the J control should stay single and any
  change of fate under R is attributable to the order (decision 3).

### 3.3 Measurements
Per run, every 1 time unit after t = 50:
- **Fate:** single / lost / multiplied, by MG2/MG3 detection; first time it is not single.
- **Centre and velocity:** sub-grid centre of the x0 track (deterministic: the parabolic minimum of φ_X on the unwrapped
  axis); velocity v = least-squares slope of the centre over t ∈ [150, T].
- **Shape difference (deterministic):** E(τ) = min over shift s of ‖φ_R(x + s) − φ_J(x)‖₂ / ‖φ_J‖₂ over all three
  species, using the final fields at the same τ (shift by Fourier interpolation). E measures the distortion the order
  causes even when the mesh pins the particle and v = 0.
- **Settledness (deterministic):** E and v at t = T − 50 and T; NOT SETTLED if E changes by more than 1% (relative) or v
  by more than 10% plus 1e-7.
- **Conflict density (A2, stochastic):** fraction of node-steps in which any species count is ≤ 2 (where a resource
  conflict is possible). Descriptive; predicted ≈ 0 at Ω = 1000.
- **Descriptive (B2):** particle diffusion D_p (MG2 method) for J, P and R-detrended; far-field Ω·Var(X) for J and P.

## 4. Outcomes
**Implementation checks (any failure: STOP, PI decides).** (i) L equals R mirrored to 1e-12 relative (deterministic) and
exactly (stochastic counts); (ii) deterministic J at τ₀ stays at the reference soliton: E against the stored reference
≤ 1e-6 and |v| ≤ 1e-7; (iii) J at τ₀ with the new code reproduces `tau_leap_euler_1d` exactly on a shared stream (one
short run).

**Gate MG8-F (does order matter at the production step?), from B1 at k = 0 and B2:**
- **ORDER-DEPENDENT** if deterministic R at τ₀ is not single at T, or E(τ₀) ≥ 1e-3, or |v_R(τ₀)| ≥ 1e-4; **and** B2
  agrees: at least 3 of 8 R runs are not single while at least 7 of 8 J runs are single, or the 8 R velocities (single
  runs, at least 6) have the deterministic sign in all runs and a two-sided t with |t| ≥ 3.5.
- **ORDER-INDEPENDENT (at this resolution)** if E(τ₀) < 1e-5 and |v_R(τ₀)| < 1e-6, and B2 R matches J in fate (fate
  counts differ by at most one) with |mean v_R| below 2 standard errors.
- **UNRESOLVED** otherwise (including deterministic and stochastic disagreeing).

**Gate MG8-C (does the effect survive τ → 0?), from B1, scored only if MG8-F is ORDER-DEPENDENT:**
Use the levels k at which R is single and settled. Fit log E against log τ over the four smallest such levels.
- **VANISHES** if the fitted order p ≥ 0.8 and E at the smallest level is ≤ 0.1·E at the largest single level.
- **PERSISTS** if p ≤ 0.2 (E roughly flat as τ falls by at least a factor 8).
- **UNRESOLVED** otherwise, or if fewer than four levels are single and settled.
The same rule is reported for |v| as a descriptive check (it may be zero on pinned levels).

**Reported, not scored:** P vs J at τ₀ (fate, D_p, Ω·Var(X), drift, which should be zero), conflict density, the
deterministic R velocity and shape at every level, any level where R loses or multiplies the particle, and the A3 count.

**Mapping to Central Question 5 (decision 6; agreed before any run).**
- ORDER-DEPENDENT + VANISHES: update order defines a frame only through the finite step. In the continuous-time limit the
  token rules pick out no frame of their own, and Model G's preferred frame is the medium's rest frame (A1). Q5 answered
  "no" for the token layer on the lattice.
- ORDER-DEPENDENT + PERSISTS: the update order defines a frame that survives refinement. Q5's hypothesis is supported for
  the token layer; a follow-up on irregular and evolving graphs would be the next study.
- ORDER-INDEPENDENT: no measurable order dependence at the production step; Q5 answered "no" at this resolution.
- Anything UNRESOLVED: Q5 stays open; the record says which part failed.

## 5. Part 0 (before any registered run)
1. **Code checks:** the three implementation checks of §4 on unregistered streams (`seed_for(8, 9, …)`), plus unit tests:
   one-step mean of each scheduler against its expected-count map on a small ring (Ω = 10⁶, 200 replicates); the sweep
   advection v_s of a single species with reactions off, against §2 A4, within 5%.
2. **Timing:** one short run per arm type (deterministic at τ₀/64; stochastic R at τ₀), projected cost restated in
   GATES.md. Projection above the cutoff: `proceed = false` until the PI decides.
3. **A3 count:** from the committed MG5 event lists, the fraction of births per tick whose insertion cavity touches
   another edit of the same tick (static; no dynamics).
Outcomes go to `results/MG8/part0/decisions.json` before any registered run.

## 6. Runs, cost and budget
- **Runs:** 16 deterministic + 25 stochastic = 41, cap 48 (the spare 7 only to rerun a run interrupted for a non-physical
  reason; a physical outcome is never replaced).
- **Cost (estimate, to be restated by Part 0):** B2 at about 1 core-s per time unit per run (MG3 measured 0.89–1.07):
  25 × 500 ≈ 3.5 core-h. B1: the ladder costs 127 τ₀-equivalents per scheduler; the deterministic map has no random draws,
  so at an assumed 0.2–0.4 core-s per τ₀ time unit, 2 × 127 × 300 ≈ 4–9 core-h, plus about 0.5 for the L checks.
  **Total about 8–13 core-h; cutoff 20 core-h.**
- **Where:** the PI's PC, 8 workers, about 1.5–2 hours of wall time. The cloud container (2 cores, jobs killed after 2 h)
  would need about 5–7 hours in several pieces.
- **Order:** Part 0, then B1 (deterministic, decides most of the gate), then B2.

## 7. Data and recomputation
Committed: per-run records (centre track, fate, E and v at checkpoints, settledness, conflict density), final
deterministic fields for B1 (small: 3 × 500 floats each), ledger, `decisions.json`. Stochastic full fields stay on the PC
with a SHA-256 manifest. `evaluate_mg8.py` and an independent `recompute_mg8.py` (E, v, fits and labels from the raw
records). Platform recorded.

## 8. Not in MG8
- The spatial layer as a dynamic test (sequential births on the evolving graph): A3 shows the registered tick is
  order-free, so this would test a scheme no earlier study used; only the static A3 count is made (decision 1).
- Exact Gillespie (the τ → 0 limit itself): at Ω = 1000 the event rate is about 10⁹ per time unit on this ring, so it is
  not affordable; the τ ladder approaches it instead.
- Causal-graph (Wolfram) invariance with token identities: needs a token-identity convention and is a separate question.
- Other orders (red-black, tilted foliations with several slopes): R is the extreme directional case; if R's effect
  vanishes, milder orders are bounded by it.

## 9. Decisions for the PI (Claude's recommendations)
| # | Decision | Recommendation |
|---|---|---|
| 1 | Which rule layer | Token rules on the 1D lattice only; spatial layer by static count (A3). Cheapest, cleanest, and the place where an order effect can actually occur |
| 2 | Orders compared | J (control), R (directional), P (random sequential), L as mirror check only |
| 3 | Ω for the stochastic arm | 1000 (J should stay single, so any fate change is the order's). Ω = 300 adds noise-driven losses (17 of 24 MG2 runs) that would muddy fate |
| 4 | τ ladder | τ₀ to τ₀/64 (7 levels), deterministic only; stochastic only at τ₀ |
| 5 | Budget | Cutoff 20 core-h, cap 48 runs, PC with 8 workers |
| 6 | Q5 mapping | Agree the outcome-to-answer mapping of §4 before any run |
| 7 | Review | ChatGPT design review of this draft before registration, as for MG7 and MG9b |

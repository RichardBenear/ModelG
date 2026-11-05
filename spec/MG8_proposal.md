# MG8 — Does the order in which the rewrite rules are applied change the outcome? Design proposal, DRAFT 2

*Stage MG036a. Claude, 2026-11-05. Status: **DRAFT 2, not registered. No MG8 run of any kind has been made.** Draft 1 is
kept as `spec/MG8_proposal_v1.md` (superseded). This draft answers ChatGPT's review MG036 (`review/MG036_chatgpt.md`;
response `review/MG036_response.md`). The only computation behind it is static arithmetic on the registered parameters
and a numerical check of the one-token step distribution in §2 A4; no dynamics were run.*

*Naming: MG8 is the multiway/foliation gate of the roadmap (Phase 6). `InitialPlan.md` called it MG7; it was renumbered
MG8 when the mass-scaling study took MG7 (GATES.md, MG7 registration; STAGE_LOG MG023b).*

## 0. What MG8 is, confirmed against the charter
- `README.md`, Central Question 5: "Is SQK's preferred frame a non-confluence effect? Quantify multiway branch
  divergence and foliation dependence."
- `InitialPlan.md`, Phase 6: "quantify how much branches diverge across update orders. Test whether SQK's preferred frame
  shows up as a measurable foliation dependence."
- `spec/rules_v1.1.md` §4: the simulation "selects one weighted path through the multiway system ... This is a choice of
  foliation and is examined in Phase 6."
- GATES.md row MG8: "Multiway branch divergence and foliation dependence", PENDING.

## 1. Question and what MG8 can and cannot conclude
> **Q: On the 1D lattice of MG2/MG3, how sensitive is the particle to the order in which nodes are updated within a time
> step, at the production step τ₀ = 6.667e-4, and does the deterministic (expected-count) difference shrink to zero as the
> step is refined?**

**MG8 is a bounded study of scheduler sensitivity of this lattice approximation.** That is its primary conclusion
(review item 2). It does not test coordinate-frame transformations or causal-graph invariance, and it does not establish
stochastic equivalence of schedulers as τ → 0 (the stochastic arm is run only at τ₀). Central Question 5 stays open
whatever MG8 finds; §4.4 states exactly what each outcome does and does not say about it.

## 2. What is known before any run (Part A: structural, no simulation)
These are arguments, not results, written down so that the study tests something they do not already settle.

**A1. Model G already has a preferred frame in the continuum.** ∂ₜφ = D∇²φ + f(φ) describes a medium at rest. Under
x′ = x − vt the equations acquire an advection term v·∂ₓφ, and they are not Lorentz invariant (first order in time, second
in space). The medium's rest frame is singled out by the mean-field equations themselves, with no discreteness and no
update order. (ChatGPT: "substantially sound".) So scheduler effects, whatever their size, are not needed to produce a
preferred frame in Model G; MG8 does not test whether they could produce one.

**A2. Commuting pairs, narrowed (review item 6).** For a *particular pair* of token events that are both enabled and
remain enabled after either fires, the two orders give the same counts, because count changes add (a commuting diamond).
This does not establish global confluence: firing one event changes the propensities of later events (R4's
n_X(n_X−1)n_Y/Ω², every first-order rate), so order changes the probabilities of what follows even far from low counts.
A pair fails to commute outright only when one event disables the other (counts of 0–2). Nothing in MG8 measures
stochastic order sensitivity through this route; draft 1's count-≤2 diagnostic is dropped.

**A3. The registered spatial rewriting is order-free by construction.** The Delaunay graph depends only on the point set,
and the registered tick (`src/evolve.py` `apply_tick`) interpolates every birth on the pre-edit mesh. Order within a tick
therefore changes neither graph nor field in MG5–MG7. MG6's tick-size sensitivity (MG021) is a granularity effect, not an
ordering effect. MG8 makes no spatial-layer measurement (decision 1).

**A4. What a sequential sweep does to one token, diffusion only (corrected; review item 1).** Draft 1 claimed the sweep
advects each species. That was wrong: it counted repeated forward hops and omitted the backward hops a token can make when
it is processed again. Corrected calculation, for scheduler R (left-to-right sweep, §3.1), reactions and sources off,
per-step one-direction hop probability r = Dτ/h², ignoring the seam (§3.1). A token is processed at its node; it stays
(1 − 2r), hops left (r), or hops right (r). A token that hops left lands on a node already processed and waits. A token
that hops right lands on a node not yet processed and is processed again, with the same three options. The one-step
displacement k (in units of h) is therefore

  P(−1) = r, P(0) = 1 − 2r + r², P(m) = r^m(1 − 2r) + r^(m+2) for m ≥ 1,

with moments (checked numerically by direct summation for r = 0.2):

- mean **0** (ChatGPT's argument: each processing sends equal expected numbers left and right);
- variance h²[2r + 2r²/(1 − r)] against 2r h² for J, so the **effective diffusion is D/(1 − r)**, an excess of r/(1 − r);
- third moment 6r²h³/(1 − r)² against 0 for J: an **odd, dispersive** term −κ₃∂ₓ³u in the mean map, with
  κ₃ = h³r²/((1 − r)²τ) ≈ D²τ/h. It is the only parity-breaking term of the diffusion part, and it is the only route by
  which R could make a symmetric particle move.

| species | D | r at τ₀ | excess D at τ₀ | κ₃ at τ₀ | excess D at τ₀/16 | κ₃ at τ₀/16 | excess D at τ₀/64 | κ₃ at τ₀/64 |
|---|---|---|---|---|---|---|---|---|
| G, X | 1 | 0.0167 | 1.7% | 0.0034 | 0.10% | 2.1e-4 | 0.03% | 5.2e-5 |
| Y | 12 | 0.200 | 25% | 0.75 | 1.3% | 0.031 | 0.31% | 0.0075 |

(h = 0.2.) **Reaction re-exposure:** a token processed twice in one step also gets a second chance to react. In J each
token has at most one exit per step; in R about a fraction r of each species gets two. For Y at τ₀ this raises the per-step
R4 conversion by roughly r ≈ 20%. MG9 put the single-particle window at b ∈ [28.93, 29.40], so corrections of this size
could plausibly change the particle at τ₀. All corrections are proportional to τ at fixed h.
**P** (random sequential) re-processes a hopped token with probability about 1/2 in either direction. By parity symmetry
its odd moments vanish in expectation, but its variance and re-exposure are also enhanced; Part 0 computes them exactly.
**Seam:** chains stop at the seam, which is drawn uniformly each step; the exact seam-averaged step operator is computed in
Part 0, and the formulas above are its N → ∞ form.

**A5. Disclosed expectation.** Material scheduler sensitivity at τ₀ is plausible (Y diffusion 25% higher under R), with
deterministic differences shrinking about linearly in τ. A particle drift under R is possible through κ₃ but its size and
sign are not predicted. The study measures these rather than assuming them.

## 3. Design (Part B)
Everything not stated here is MG2/MG3 as registered: P17 with b = 29.0, h = 0.2, L = 100 (N = 500, periodic),
D_G = D_X = 1, D_Y = 12, the per-token exit probabilities and probability-sum rescaling of `tau_leap_euler_1d`, the
reference soliton from `results/MG003/ref_profiles.npz`, smoothing and the MG2 detector (threshold −0.5|A*_s|).

### 3.1 Scheduler semantics (executable specification; review item 5)
**Node update (identical for all schedulers).** Take a snapshot (g, x, y) of the node's counts at the moment of update.
Draw every exit of every species from that one snapshot, in a fixed sequence: G exits (R2, R1r, hop L, hop R), X exits (R2r,
R3, R5, hop L, hop R), Y exits (R4 with probability x(x − 1)τ/Ω², hop L, hop R), then the R1 source. Then apply all of the
node's changes at once: conversions stay at the node, and hops are added to the neighbours' current counts. Species are
never processed one after another against changing counts.
**Schedulers.** They differ only in when the snapshot is taken:
- **J (synchronous; the registered kernel).** All snapshots are start-of-step counts; all changes are applied after the
  last node. For the same stream it must reproduce `tau_leap_euler_1d` exactly.
- **R (rightward sweep).** Draw a seam s uniformly from {0, …, N − 1}; update nodes s, s + 1, …, s − 1 (mod N) in turn, each
  snapshot taken when the node is reached (arrivals included).
- **L (leftward sweep; mirror check only).** Seam N − 1 − s, update order s′, s′ − 1, …; at each node the hop-L and hop-R
  draws swap roles. Started from the mirror image of R's start with R's seed, L must reproduce R's trajectory mirrored
  node for node, exactly (stochastic) or to 1e-12 relative (deterministic).
- **P (random sequential).** Each step a uniform random permutation of the nodes is drawn and nodes are updated in that
  order.
**Random numbers.** One numba generator per run, seeded `seed_for(8, 0, 0, 0, rep)` (the same seed for every arm with the
same rep). Per step the draws are consumed in this order: the order draw (seam for R/L, permutation for P, none for J),
then the node updates in processing order with the fixed per-node sequence above. Because the schedulers consume draws in
different orders, **a shared seed does not couple trajectories across arms; comparisons between arms are distributional.**
The shared seed is bookkeeping only. Part 0 uses `seed_for(8, 9, …)`.
**Deterministic maps (Ω → ∞).** The same orders applied to expected counts: each exit equals count × probability (with the
same rescaling rule), and R4 uses x²y in the mean-field form. Deterministic J must equal `lattice.euler_mean_map` to 1e-12.
Deterministic P is not defined and not run.
**Clipping.** Every scheduler counts node-steps whose exit probabilities summed above 1 and were rescaled, per species.

### 3.2 Arms
- **B1, deterministic ladder (primary):** J and R at τ_k = τ₀/2^k, k = 0, …, 6; L at k = 0 and k = 6 (mirror checks);
  plus a longer-horizon continuation of R at k = 6 from T = 300 to T = 600. Start: the reference soliton, T = 300.
  17 runs.
- **B2, stochastic at τ₀ (descriptive):** J, R and P at Ω = 1000, 8 runs each, plus one L mirror run on R's first seed.
  Start: the MG2 soliton start (`np.rint(Ω·c)`), T = 500. 25 runs.

### 3.3 Measurements
**Deterministic (B1), every unit time:**
- **Fate:** number of detections of the smoothed φ_X; SINGLE if exactly one at every check from t = 50 to T.
- **Centre and velocity:** parabolic sub-grid minimum of smoothed φ_X, unwrapped; v = least-squares slope over
  [T − 150, T]; v′ = the same over [T − 200, T − 50] (settledness).
- **Shape:** for each species, E_s(τ) = min over shift σ of ‖φ_s^R(x + σ) − φ_s^J(x)‖₂ / ‖φ_s^J‖₂ on the final fields at the
  same τ (Fourier shift, σ optimised jointly for all three species); **E = max over species**. E′ is the same at T − 50.
- **Floors:** E_floor(τ) = E between J's own fields at T − 50 and T (J against itself) + 1e-9; v_floor(τ) = |v_J(τ)| + 1e-8.
- **Settled:** |E − E′| ≤ 0.01E + 1e-7 and |v − v′| ≤ 0.1|v| + 1e-7.
**Stochastic (B2):**
- **Fate:** MG2's persistence-qualified classification (`measure.analyze_run` events: loss, replication, birth, merge),
  over the whole run. **TRACK-FAIL** is a separate label: the x0 track ends while some detection lies within 3 units of its
  last position at the next snapshot (a linking failure, not a physical loss).
- **Drift:** per run, slope of the unwrapped x0 track over [150, 500] when the run has no loss/replication/merge/TRACK-FAIL.
- **Shape:** per run, the time-mean aligned profile over [150, 500] (`analyze_run` block profiles); E_stoch of each arm's
  mean profile against J's, with a bootstrap 95% interval over runs (2000 resamples).
- **Also:** D_p (MG2 method), far-field Ω·Var(X), clipping counts.

## 4. Outcomes
### 4.1 Implementation checks (any failure: STOP, PI decides)
J regression (stochastic exact; deterministic 1e-12), R/L mirror (exact; 1e-12), Part 0 benchmark (§5) passed.

### 4.2 Scheduler sensitivity at τ₀ (primary; deterministic B1 at k = 0; separate labels, review item 3)
- **FATE:** CHANGED if R is not SINGLE; SAME otherwise.
- **SHAPE** (only if FATE SAME): MATERIAL if E ≥ 1e-2 (the 1% mean-field tolerance of MG2 C1); NEGLIGIBLE if E < 1e-3 or
  E < 10·E_floor; SMALL otherwise. E must be settled; if not, UNSETTLED.
- **DRIFT** (only if FATE SAME): DRIFTING if settled and |v| ≥ max(10·v_floor, 1e-6); NOT DRIFTING AT THIS HORIZON if settled
  and below that (a pinned particle can hide a parity-breaking force, MG9d); UNSETTLED otherwise.
- **Summary:** SENSITIVE if FATE CHANGED, SHAPE MATERIAL or DRIFTING; PRACTICALLY INSENSITIVE if FATE SAME, SHAPE NEGLIGIBLE
  and NOT DRIFTING; MINOR if the only effect is SHAPE SMALL; UNRESOLVED if anything is UNSETTLED.
The same labels are reported for every ladder level.

### 4.3 Deterministic refinement (MG8-C; scored whatever 4.2 gives; review item 4)
Eligible levels: k = 2, …, 6 (refinement span 16) at which R is SINGLE and settled and E ≥ 10·E_floor.
- **CONVERGES TO J** if (a) at least 4 levels are eligible, spanning at least a factor 8 in τ; the least-squares slope p of
  ln E on ln τ has a 95% t interval (df = n − 2) with lower bound ≥ 0.5; every residual is within ±0.15 in ln E; and
  E(k_max) ≤ 0.1·E(k_min) over the eligible levels; and the longer-horizon continuation changes E at k = 6 by at most
  0.01E + 1e-7; **or** (b) fewer than 4 levels are eligible only because E falls below 10·E_floor, E decreases
  monotonically over k = 0, …, 6, and R is SINGLE and settled at every level.
- **NOT RESOLVED BY THIS LADDER** otherwise, with the reason recorded (R not single at some k ≥ 2, unsettled, slope interval
  too wide, plateau, failed residual or horizon check). A plateau is reported as such; MG8 has no label asserting a
  non-zero limiting discrepancy.
- Velocity at each level is descriptive. Shape convergence does not imply velocity convergence; velocities near the floor
  are reported as "not drifting at this horizon", never as pinned.

### 4.4 Stochastic results (B2; descriptive)
Per arm: fate counts with exact binomial 95% intervals; mean drift with its 95% t interval; E_stoch with its bootstrap
interval; D_p; Ω·Var(X); clipping. An arm is called **equivalent to J within the margin** only when the whole interval lies
inside the margin (drift |v| < 1e-3; E_stoch < 1e-2); an interval that straddles the margin is reported as **imprecise**,
not as equivalence. Minimum detectable drift (planning): with D_p ≈ 1e-4 (MG2: Ω·D_p ≈ 0.1), one run's drift over 350
units has SE ≈ 7.6e-4, so 8 runs give a 95% half-width of about 6e-4.

### 4.5 What the outcomes say about Central Question 5 (review item 2)
- The primary conclusion is always stated as **scheduler sensitivity of this lattice approximation**, with the 4.2 labels.
- CONVERGES TO J says the deterministic expected-count dynamics of R approach J's as τ → 0 at fixed h. It does not show
  stochastic equivalence, causal-graph invariance or frame invariance.
- NOT RESOLVED BY THIS LADDER says nothing about a limiting frame.
- Q5 remains open after MG8. A1 shows that Model G's preferred frame exists without any order dependence; any further
  link between scheduler sensitivity and a physical frame would need a separate mathematical argument and study.

## 5. Part 0 (before any registered run; priced separately)
1. **Exact small-system diffusion benchmark (replaces draft 1's 5% check).** One species, reactions and sources off.
   (a) Exact seam-averaged one-step operators T_R on a ring of N = 12, and the permutation-averaged T_P on N = 6 (all 720
   orders), for r ∈ {0.01, 0.05, 0.2}, built from the per-token chain rules of §3.1. (b) The simulator with 10⁶ tokens per
   node, 2000 replicate steps per case: conservation exact; mean unwrapped displacement zero within 4 standard errors;
   variance and third moment of displacement match T within 4 standard errors; R and P mean one-step maps match T to
   1e-12 in the deterministic version. (c) The deterministic per-unit-time operator difference
   ‖(T_R^n − T_J^n)‖ at fixed elapsed time falls about linearly in τ (ratio 0.4–0.6 per halving, r from 0.2 to 0.0125).
   (d) The infinite-ring formulas of A4 are compared with T_R (descriptive). Any failure in (a)–(c): STOP.
2. **J regression and R/L mirror** on unregistered seeds (as §4.1), plus deterministic J against `euler_mean_map`.
3. **Clipping:** per scheduler, the fraction of rescaled node-steps in a 50-unit run at τ₀ (stochastic, Ω = 1000) and in
   the deterministic k = 0 run; reported.
4. **Synthetic validation of the measurements:** (a) shape: profiles shifted by known sub-grid σ and given a known odd
   deformation of norm 1e-3 and 1e-2 — E recovered within 2% and σ within 1e-3; (b) velocity: profiles moved at 1e-3 and
   1e-5 per unit, with a wrap — v within 2%; (c) tracking and fate: synthetic detection sequences that move, wrap,
   disappear, split, and jump beyond the link distance — classified as moving/single, loss, replication and TRACK-FAIL
   respectively; (d) convergence labels: synthetic E ladders with p = 1 (CONVERGES), p = 0.3, a plateau, and values near
   the floor — the §4.3 rule gives the intended label. Any failure: fix and repeat before production.
5. **Timing** of the sequential kernels (seam and permutation included) at τ₀ (stochastic) and τ₀/64 (deterministic);
   projected production cost restated in GATES.md. Above the cutoff: `proceed = false` until the PI decides.
Part 0 outcomes go to `results/MG8/part0/decisions.json` before any registered run. Part 0 is priced at about 1 core-h,
cap 3 core-h, **not counted** in the production cutoff.

## 6. Runs, cost and budget
- **Runs:** 17 deterministic + 25 stochastic = 42, cap 48 (spares only to rerun a run interrupted for a non-physical
  reason; a physical outcome is never replaced).
- **Cost (estimate; Part 0 restates it):** stochastic at 1–2 core-s per τ₀ time unit (MG3: 0.89–1.07 for J; sequential
  updates and order draws add overhead): 25 × 500 → 3.5–7 core-h. Deterministic: the ladder is 127 τ₀-equivalents per
  scheduler, plus 64 × 300 for the continuation, about 95,000 τ₀ time units in all; with no random draws, at 0.03–0.3 core-s
  each, about 1–8 core-h. **Total about 5–15 core-h; production cutoff 20 core-h.**
- **Order:** Part 0, then B1, then B2.
- **Budget stop:** if the cutoff is reached, the run in progress finishes, nothing new starts, and results are labelled
  BUDGET-SHORT. A stop in B1 leaves 4.2 scored only if k = 0 completed and 4.3 NOT RESOLVED BY THIS LADDER (BUDGET-SHORT);
  a stop in B2 leaves B2 descriptive on the completed runs.
- **Where:** the PI's PC, 8 workers, about 1–2 hours of wall time.

## 7. Data and recomputation (review "validation" list)
Committed: all B1 records and fields (final fields and smoothed φ_X every 10 units, float64; a few MB), per-unit detection
lists, tracks and events for every B2 run, the B2 run-mean aligned profiles, and B2 smoothed φ_X snapshots every 5 units
in float32 (about 5 MB compressed), the ledger and `decisions.json`. Full B2 fields stay on the PC with a SHA-256
manifest (decision 8). `evaluate_mg8.py` and an independent `recompute_mg8.py` recompute every label and interval from the
committed records. Platform recorded.

## 8. Not in MG8
- The spatial layer (A3).
- Exact Gillespie (about 10⁹ events per time unit on this ring at Ω = 1000).
- Stochastic refinement (B2 at smaller τ).
- Causal-graph invariance with token identities, and coordinate-frame (boost) tests.
- Other orders (red-black, tilted foliations). Draft 1's claim that R bounds milder orders is withdrawn (no proof).

## 9. Decisions for the PI (Claude's recommendations)
| # | Decision | Recommendation |
|---|---|---|
| 1 | Which rule layer | Token rules on the 1D lattice only; spatial layer argued in A3, not measured |
| 2 | Orders compared | J (control), R (directional), P (random sequential, stochastic only), L as mirror check only |
| 3 | Ω for the stochastic arm | 1000 (J should stay single, so a fate change under R is the order's) |
| 4 | τ ladder | τ₀ to τ₀/64 (7 levels), deterministic only, plus the T = 600 continuation at τ₀/64 |
| 5 | Budget | Production cutoff 20 core-h, cap 48 runs; Part 0 separate (cap 3); PC with 8 workers |
| 6 | Conclusion framing | Primary conclusion "scheduler sensitivity of this lattice approximation"; Q5 stays open (§4.5) |
| 7 | Stochastic arm status | Descriptive with intervals and stated margins (not a gate) |
| 8 | Committed B2 data | Per-unit detection records plus φ_X snapshots every 5 units (about 5 MB); full fields on the PC |
| 9 | Review | ChatGPT round 2 on this draft before registration |

## 10. Changes from draft 1 (review MG036)
| Review item | Change |
|---|---|
| 1 (blocking) A4 advection wrong | A4 replaced by the exact one-token step distribution: zero mean, D/(1 − r), odd κ₃ term, re-exposure; Part 0 exact benchmark (§5.1) |
| 2 (blocking) Q5 mapping | Primary conclusion is scheduler sensitivity; PERSISTS label removed; Q5 stays open (§1, §4.5); R-bounds claim withdrawn (§8) |
| 3 Separate outcomes, equivalence | FATE / SHAPE / DRIFT labels with margins (§4.2); stochastic arm descriptive with intervals and margins (§4.4) |
| 4 Convergence safeguards | Floors, absolute + relative settledness, span ≥ 8, slope interval, residuals, horizon continuation; "not resolved by this ladder" (§4.3) |
| 5 Semantics and streams | One snapshot per node update; order, seed and draw sequence specified; shared seed declared non-coupling; exact mirror construction (§3.1) |
| 6 A2, fate, tracking | A2 narrowed and the count-≤2 diagnostic dropped; fate = MG2 persistence-qualified classification; TRACK-FAIL separate; synthetic validation (§5.4) |
| Validation list | Clipping diagnostics (§3.1, §5.3); committed snapshots (§7); fresh timing, budget-stop labels, Part 0 outside the cutoff (§5–6) |

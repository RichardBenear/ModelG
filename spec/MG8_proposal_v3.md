# MG8 — Does the order in which the rewrite rules are applied change the outcome? Design proposal, DRAFT 3

*Stage MG036b. Claude, 2026-11-05. Status: **SUPERSEDED by draft 4 (`spec/MG8_proposal.md`) after ChatGPT review MG036 round 3 (`review/MG036_chatgpt_r3.md`). DRAFT 3, not registered. No MG8 run of any kind has been made.** Drafts 1
and 2 are kept as `spec/MG8_proposal_v1.md` and `_v2.md` (superseded). Draft 2 answered ChatGPT's review MG036
(`review/MG036_chatgpt.md`; response `review/MG036_response.md`); this draft answers round 2
(`review/MG036_chatgpt_r2.md`; response `review/MG036_response_r2.md`). Changes from draft 2 are listed in §11. The only computation behind it is static arithmetic on the registered parameters
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
  κ₃ = h³r²/((1 − r)²τ) ≈ D²τ/h. This is the **leading** parity-breaking contribution of the diffusion part; the exact
  map also has higher odd moments, and mixed reaction–diffusion effects of the sequential update (below) can also break
  parity. Any of these could make a symmetric particle move; MG8 measures motion, it does not attribute it.

| species | D | r at τ₀ | excess D at τ₀ | κ₃ at τ₀ | excess D at τ₀/16 | κ₃ at τ₀/16 | excess D at τ₀/64 | κ₃ at τ₀/64 |
|---|---|---|---|---|---|---|---|---|
| G, X | 1 | 0.0167 | 1.7% | 0.0034 | 0.10% | 2.1e-4 | 0.03% | 5.2e-5 |
| Y | 12 | 0.200 | 25% | 0.75 | 1.3% | 0.031 | 0.31% | 0.0075 |

(h = 0.2.) **Reaction re-exposure:** a token processed twice in one step also gets a second chance to react. In J each
token has at most one exit per step; in R about a fraction r of each species gets two. **Heuristic only:** for Y at τ₀
this could raise the per-step R4 conversion by something of order r ≈ 20%, but hopping competes with reaction exits and
the local X-dependent R4 rate changes during the sweep, so the figure is an order of magnitude, not a prediction. MG9 put the single-particle window at b ∈ [28.93, 29.40], so corrections of this size
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

### 3.1 Scheduler semantics (executable specification)
**Node update (identical for all schedulers).** Take a snapshot (g, x, y) of the node's counts at the moment of update.
Draw every exit of every species from that one snapshot, in a fixed sequence: G exits (R2, R1r, hop L, hop R), X exits (R2r,
R3, R5, hop L, hop R), Y exits (R4 with probability x(x − 1)τ/Ω², hop L, hop R), then the R1 source. Then apply all of the
node's changes at once: conversions stay at the node, and hops are added to the neighbours' current counts. Species are
never processed one after another against changing counts.
**Orders.** An *order* o fixes the sequence in which the N nodes are updated within one step:
- **J (synchronous; the registered kernel).** All snapshots are start-of-step counts; all changes are applied after the
  last node. For the same stream it must reproduce `tau_leap_euler_1d` exactly.
- **R_s (rightward sweep with seam s).** Nodes s, s + 1, …, s − 1 (mod N), each snapshot taken when the node is reached
  (arrivals included).
- **L_s (leftward sweep; mirror check only).** Nodes s′, s′ − 1, … with s′ = N − 1 − s, hop-L and hop-R draws swapping
  roles at each node. Started from the mirror image of R's start with R's seed and mirrored seams, L must reproduce R's
  trajectory mirrored node for node, exactly (stochastic) or to 1e-12 relative (deterministic).
- **P_π (sequential in permutation π).**
**Schedulers** draw a fresh order every step: R draws s uniformly from {0, …, N − 1}; P draws π uniformly from all
permutations; J has no draw.
**Three deterministic objects (review r2 item 1).** For each order o, the **conditional map** M_o applies that order to
expected counts (each exit = count × probability, the same rescaling rule, R4 in mean-field form x²y). Then:
(i) M_s for one specified seam; (ii) the **order-averaged one-step map** M̄_R = (1/N)Σ_s M_s (and M̄_P likewise);
(iii) a **random-seam trajectory**: M_{s_n} ∘ … ∘ M_{s_1} with s_1, s_2, … drawn from a seeded stream.
**B1 tests object (iii)**, replicated over three independent seam streams per level (§3.2). It is an expected-count
trajectory under randomly drawn orders, not the order average; object (ii) is computed only in the diffusion-only
benchmark of Part 0 (§5.1), which validates the update rule, not the production experiment. Seam streams at different
levels are independent (the levels take different numbers of steps, so no sequence can be shared); seam realisation
effects are measured by the replicate spread. Deterministic J must equal `lattice.euler_mean_map` to 1e-12 whenever no
rescaling occurs. Deterministic P_π is defined (used in Part 0) but not run in B1.
**Random numbers.** One numba generator per run. Per step the draws are consumed in this order: the order draw (seam or
permutation; none for J), then the node updates in processing order with the fixed per-node sequence above.
- B1 (deterministic; only seams are random): seam stream `seed_for(8, 1, k, 0, j)` for level k and replicate j = 0, 1, 2.
- B2 (stochastic): **independent streams per arm (review r2 item 4)**: `seed_for(8, 0, a, 0, rep)` with a = 0 (J), 1 (R),
  2 (P). The L mirror run uses R's seed. Arms are independent samples; nothing is paired across arms.
- Part 0: `seed_for(8, 9, …)`.
**Clipping.** Every scheduler counts node-steps whose exit probabilities summed above 1 and were rescaled, per species.

### 3.2 Arms
- **B1, deterministic ladder (primary):** at τ_k = τ₀/2^k, k = 0, …, 6: J (one run; no randomness) and R with three seam
  replicates j = 0, 1, 2; T = 300. Mirror checks: L at k = 0 (T = 300) and k = 6 (T = 50). **Horizon continuation:** J and
  R (j = 0) at k = 6 continued from T = 300 to T = 600. Start: the reference soliton. 32 runs.
- **B2, stochastic at τ₀ (descriptive):** J, R and P at Ω = 1000, 8 runs each, plus one L mirror run on R's first seed.
  Start: the MG2 soliton start (`np.rint(Ω·c)`), T = 500. 25 runs.

### 3.3 Measurements
**Deterministic (B1), every unit time:**
- **Fate:** number of detections of the smoothed φ_X; SINGLE if exactly one at every check from t = 50 to T.
- **Centre and velocity:** parabolic sub-grid minimum of smoothed φ_X, unwrapped; v = least-squares slope over
  [T − 150, T]; v′ = the same over [T − 200, T − 50].
- **Shape:** for each species, E_s = min over shift σ of ‖φ_s^R(x + σ) − φ_s^J(x)‖₂ / ‖φ_s^J‖₂, R and J taken **at the same
  time and the same τ** (Fourier shift; σ optimised jointly for all three species by an integer-node search followed by
  golden-section search on [−h/2, h/2] to 1e-10); **E = max over species**. E at T, E′ at T − 50.
- **Per-level summary:** E_k = median of the three replicates' E; seam spread S_k = (max − min)/E_k.
- **Floors:** E_floor(τ_k) = E between J's own fields at T − 50 and T, plus 1e-9; v_floor(τ_k) = |v_J(τ_k)| + 1e-8.
- **Settled (per run):** |E − E′| ≤ 0.01E + 1e-7 and |v − v′| ≤ 0.1|v| + 1e-7.
- **Horizon check (k = 6):** E₆₀₀ = E(R at t = 600, J at t = 600); passes if R is SINGLE through 600, settled at 600 (the
  rule above with T = 600), and |E₆₀₀ − E₃₀₀| ≤ 0.05·E₃₀₀ + 1e-7.
**Stochastic (B2):**
- **Fate:** MG2's persistence-qualified classification (`measure.analyze_run` events: loss, replication, birth, merge),
  over the whole run. **TRACK-FAIL** is a separate label: the x0 track ends while some detection lies within 3 units of its
  last position at the next snapshot (a linking failure, not a physical loss).
- **Eligible run:** no loss, replication, merge or TRACK-FAIL in [0, 500]. Drift and shape use eligible runs only, and are
  reported **conditional on eligibility**, with the eligible count per arm.
- **Drift:** per eligible run, slope of the unwrapped x0 track over [150, 500].
- **Shape:** `analyze_run` with t_discard = 150, T = 500, K = 7 blocks; per-run profile = the snapshot-count-weighted mean of
  its block profiles (weights n_per_block); arm profile = unweighted mean over its eligible runs; E_stoch(arm) = the E
  definition above applied to arm profile against J's arm profile. A J split-half baseline (E between the mean profiles of
  J runs 1–4 and 5–8) is reported as the noise reference.
- **Also:** D_p (MG2 method), far-field Ω·Var(X), clipping counts.

## 4. Outcomes
### 4.1 Implementation checks (any failure: STOP, PI decides)
J regression (stochastic exact; deterministic 1e-12 against `euler_mean_map` without rescaling), R/L mirror (exact;
1e-12), Part 0 benchmark (§5.1) passed.

### 4.2 Scheduler sensitivity at τ₀ (primary; deterministic B1 at k = 0)
Labels are assigned per replicate. Each label set is mutually exclusive and evaluated in the order written.
- **FATE:** CHANGED if R is not SINGLE; SAME otherwise.
- **SHAPE** (if FATE SAME): INSUFFICIENT RESOLUTION if 10·E_floor ≥ 1e-3; else UNSETTLED if not settled; else MATERIAL if
  E ≥ 1e-2 (the 1% mean-field tolerance of MG2 C1), SMALL if 1e-3 ≤ E < 1e-2, NEGLIGIBLE if E < 1e-3.
- **DRIFT** (if FATE SAME): INSUFFICIENT RESOLUTION if 10·v_floor ≥ 1e-6; else UNSETTLED if not settled; else DRIFTING if
  |v| ≥ 1e-6, NOT DRIFTING AT THIS HORIZON otherwise (a pinned particle can hide a parity-breaking force, MG9d).
- **Summary per replicate:** SENSITIVE if FATE CHANGED, SHAPE MATERIAL or DRIFTING; otherwise PRACTICALLY INSENSITIVE if
  SHAPE NEGLIGIBLE and NOT DRIFTING; otherwise MINOR if SHAPE SMALL and NOT DRIFTING; otherwise INSUFFICIENT RESOLUTION if
  either label is; otherwise UNRESOLVED.
- **Production-step result:** the common label if all three replicates agree; otherwise **SEAM-DEPENDENT**, with the three
  labels reported.
The same labels are reported for every ladder level.

### 4.3 Deterministic refinement (MG8-C; scored whatever 4.2 gives)
**Eligible level:** k ∈ {2, …, 6} with all three R replicates SINGLE and settled, S_k ≤ 0.25, and E_k ≥ 10·E_floor(τ_k).
Outcomes, mutually exclusive, evaluated in this order:
1. **CONVERGES TO J** if: at least 4 eligible levels spanning at least a factor 8 in τ; the least-squares slope p of ln E_k
   on ln τ_k is ≥ 0.7, and the lower end of its regression diagnostic interval (t, df = n − 2) is ≥ 0.5; every residual is
   within ±0.15 in ln E; E at the smallest eligible τ ≤ 0.1·E at the largest; and the horizon check passes.
2. **BELOW MEASUREMENT RESOLUTION** if fewer than 4 levels are eligible only because E_k < 10·E_floor(τ_k) at some
   k ≥ 2, and: 10·E_floor(τ_k) < 1e-4 at every k ≥ 2; every replicate at every k = 0, …, 6 is SINGLE and settled; E_k is
   non-increasing in k; E₆ < 1e-4; and the horizon check passes (with E₆₀₀ < 1e-4 in place of the relative clause when E₃₀₀
   is below the floor). This says the discrepancy fell below what MG8 can resolve. **It is not demonstrated convergence.**
3. **NOT RESOLVED BY THIS LADDER** otherwise, with the reason recorded (R not single at some k ≥ 2, unsettled,
   seam-sensitive level, slope or residual failure, plateau, floor too high, failed horizon check). MG8 has no label
   asserting a non-zero limiting discrepancy.
**The regression interval is a diagnostic, not calibrated uncertainty about the convergence order (review r2 item 6):** the
residuals are deterministic discretisation and settling errors, not independent noise, so no probability model is
claimed. Velocity at each level is descriptive; shape convergence does not imply velocity convergence, and velocities
near the floor are reported as "not drifting at this horizon", never as pinned.

### 4.4 Stochastic results (B2; descriptive)
Per arm: fate counts with exact binomial 95% intervals; eligible counts. Against J (arms are independent):
- **Drift:** Welch 95% interval for mean(v_arm) − mean(v_J) over eligible runs.
- **Shape:** E_stoch(arm) with a percentile bootstrap 95% interval from independent resampling of eligible runs within each
  arm (J resampled alongside), 2000 resamples, seed `seed_for(8, 2, 0, 0, 0)`; reported beside the J split-half baseline.
- Requires at least 5 eligible runs in the arm and in J; otherwise NOT ESTIMATED.
- **Equivalent to J within the margin** only when the whole interval lies inside the margin (|Δv| < 1e-3; E_stoch < 1e-2);
  an interval that crosses the margin is **imprecise**, not equivalence.
- Planning: with D_p ≈ 1e-4 (MG2: Ω·D_p ≈ 0.1), one run's drift over 350 units has SE ≈ 7.6e-4, so with 8 runs per arm
  the difference interval has a half-width of about 8e-4, close to the 1e-3 margin: drift equivalence may well be
  imprecise. This is accepted for a descriptive arm.

### 4.5 What the outcomes say about Central Question 5
- The primary conclusion is always stated as **scheduler sensitivity of this lattice approximation**, with the 4.2 labels.
- CONVERGES TO J says the deterministic expected-count trajectories of R under random seams approach J's as τ → 0 at fixed
  h. It does not show stochastic equivalence, causal-graph invariance or frame invariance.
- BELOW MEASUREMENT RESOLUTION says only that the discrepancy fell below what MG8 resolves.
- NOT RESOLVED BY THIS LADDER says nothing about a limiting frame.
- Q5 remains open after MG8. A1 shows that Model G's preferred frame exists without any order dependence; any further
  link between scheduler sensitivity and a physical frame would need a separate mathematical argument and study.

## 5. Part 0 (before any registered run; priced separately)
1. **Diffusion-only operator benchmark (review r2 item 3).** One species, reactions and sources off, h = 1, D = 1, so
   r = τ.
   (a) **Conditional maps as matrices.** For R on a ring of N = 12 (each of the 12 seams), P on N = 6 (all 720
   permutations) and J, build T_o by applying the production deterministic code's M_o to each basis vector e_j (column j).
   Build the same matrices independently from the per-token chain rules (separate code); agreement 1e-12. Average them
   explicitly: T̄_R = (1/12)Σ_s T_s, T̄_P = (1/720)Σ_π T_π.
   (b) **Unwrapped displacement (augmented benchmark).** Propagate labelled cohorts (start node, unwrapped displacement d,
   |d| ≤ N) through each M_o to obtain A_o[j, d], the expected fraction of tokens starting at j with displacement d; average
   over orders as in (a). Then the stochastic simulator in benchmark mode, with cohort labels, runs one step from a
   **localised** population of 10⁶ tokens at node 0 with a random order each replicate, 2000 replicates per case,
   r ∈ {0.01, 0.05, 0.2}. Pass: conservation exact; mean displacement within 4 SE of 0; variance and third moment within
   4 SE of Ā's; chi-square goodness of fit to Ā with p ≥ 1e-3 (cells with expected count < 5 pooled).
   (c) **Refinement (frozen construction).** N = 12, elapsed time t = 1, τ ∈ {0.2, 0.1, 0.05, 0.025, 0.0125} with n = 1/τ
   steps; Δ(τ) = ‖T̄^n − T_J^n‖_F (Frobenius norm of the 12 × 12 averaged operators; N = 6 for P). Pass: Δ decreases at every
   halving, and the last two ratios Δ(τ/2)/Δ(τ) lie in [0.4, 0.6]. Descriptive: the same for one random-seam product
   T_{s_n}⋯T_{s_1} (seed `seed_for(8, 9, 5, 0, 0)`).
   (d) The infinite-ring formulas of A4 against Ā (descriptive).
   Any failure in (a)–(c): STOP.
2. **J regression and R/L mirror** on unregistered seeds (as §4.1).
3. **Clipping:** per scheduler, the fraction of rescaled node-steps in a 50-unit run at τ₀ (stochastic, Ω = 1000) and in a
   deterministic k = 0 run; reported.
4. **Synthetic validation of the measurements:** (a) shape: profiles shifted by known sub-grid σ and given a known odd
   deformation of norm 1e-3 and 1e-2 — E recovered within 2% and σ within 1e-3; (b) velocity: profiles moved at 1e-3 and
   1e-5 per unit, with a wrap — v within 2%; (c) tracking and fate: synthetic detection sequences that move, wrap,
   disappear, split, and jump beyond the link distance — classified as single, loss, replication and TRACK-FAIL
   respectively; (d) labels: synthetic E ladders and floors that should give each §4.2 and §4.3 label, including
   CONVERGES, a plateau, BELOW MEASUREMENT RESOLUTION, a large floor (INSUFFICIENT RESOLUTION / not resolved) and a
   seam-sensitive level — every label as intended; (e) the stochastic eligibility, weighting and bootstrap code on
   synthetic run records with known answers. Any failure: fix and repeat before production.
5. **Timing** of the sequential kernels (seam and permutation included) at τ₀ (stochastic) and τ₀/64 (deterministic);
   projected production cost restated in GATES.md. Above the cutoff: `proceed = false` until the PI decides.
Part 0 outcomes go to `results/MG8/part0/decisions.json` before any registered run. Part 0 is priced at about 1 core-h,
cap 3 core-h, **not counted** in the production cutoff.

## 6. Runs, cost and budget
- **Runs:** 32 deterministic + 25 stochastic = 57, cap 64 (spares only to rerun a run interrupted for a non-physical
  reason; a physical outcome is never replaced).
- **Cost (estimate; Part 0 restates it):** stochastic at 1–2 core-s per τ₀ time unit (MG3: 0.89–1.07 for J; sequential
  updates and order draws add overhead): 25 × 500 → 3.5–7 core-h. Deterministic: J and three R replicates over the ladder
  are 4 × 127 × 300 τ₀ time units, plus the two continuations (2 × 64 × 300) and the L checks, about 195,000 τ₀ time units;
  with no token draws, at 0.02–0.2 core-s each, about 1–11 core-h. **Total about 5–18 core-h; production cutoff 20 core-h.**
- **Order:** Part 0, then B1 (k = 0 first, then the ladder, then the continuations), then B2.
- **Budget stop:** if the cutoff is reached, runs in progress finish, nothing new starts, and results are labelled
  BUDGET-SHORT. A stop in B1 leaves 4.2 scored only if all k = 0 runs completed, and 4.3 NOT RESOLVED BY THIS LADDER
  (BUDGET-SHORT); a stop in B2 leaves B2 descriptive on the completed runs.
- **Where:** the PI's PC, 8 workers, about 1–3 hours of wall time.

## 7. Data and recomputation (review r2 item 5)
Committed (small):
- **B1:** unit-time centres and detection counts for every run; full three-species fields (float64) at every checkpoint
  any decision uses — t = 50, 100, …, T for every run (so T − 50 and T are included) and t = 550, 600 for the
  continuations — for R, L and the J references alike (about 2–3 MB in all); the seam-stream seeds, the generator
  specification and the SHA-256 of each realised seam sequence (sequences are regenerated, not stored); alignment
  settings (§3.3).
- **B2:** per-unit detection lists, tracks and events; per-run block profiles with their snapshot counts; per-run drift;
  seeds; the bootstrap seed and resample count; smoothed φ_X every 5 units in float32 (about 5 MB compressed).
- The ledger and `decisions.json`.
PC only: full B2 fields, with a committed SHA-256 manifest (decision 8). `evaluate_mg8.py` and an independent
`recompute_mg8.py` recompute every label, interval and fit from the committed records. Platform recorded.

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
| 4 | τ ladder | τ₀ to τ₀/64 (7 levels), deterministic only, three seam replicates of R per level, plus the T = 600 continuations of J and R at τ₀/64 |
| 5 | Budget | Production cutoff 20 core-h, cap 64 runs; Part 0 separate (cap 3); PC with 8 workers |
| 6 | Conclusion framing | Primary conclusion "scheduler sensitivity of this lattice approximation"; Q5 stays open (§4.5) |
| 7 | Stochastic arm status | Descriptive with intervals and stated margins (not a gate); independent streams per arm |
| 8 | Committed B2 data | Per-unit detection records plus φ_X snapshots every 5 units (about 5 MB); full fields on the PC |
| 9 | Review | A short ChatGPT round 3 check of the items 1–3 fixes before registration (ChatGPT asked for items 1–3 to be resolved before registration and 4–6 frozen before the implementation review) |

## 10. Changes from draft 1 to draft 2 (review MG036)
| Review item | Change |
|---|---|
| 1 (blocking) A4 advection wrong | A4 replaced by the exact one-token step distribution: zero mean, D/(1 − r), odd κ₃ term, re-exposure; Part 0 exact benchmark (§5.1) |
| 2 (blocking) Q5 mapping | Primary conclusion is scheduler sensitivity; PERSISTS label removed; Q5 stays open (§1, §4.5); R-bounds claim withdrawn (§8) |
| 3 Separate outcomes, equivalence | FATE / SHAPE / DRIFT labels with margins (§4.2); stochastic arm descriptive with intervals and margins (§4.4) |
| 4 Convergence safeguards | Floors, absolute + relative settledness, span ≥ 8, slope interval, residuals, horizon continuation; "not resolved by this ladder" (§4.3) |
| 5 Semantics and streams | One snapshot per node update; order, seed and draw sequence specified; shared seed declared non-coupling; exact mirror construction (§3.1) |
| 6 A2, fate, tracking | A2 narrowed and the count-≤2 diagnostic dropped; fate = MG2 persistence-qualified classification; TRACK-FAIL separate; synthetic validation (§5.4) |
| Validation list | Clipping diagnostics (§3.1, §5.3); committed snapshots (§7); fresh timing, budget-stop labels, Part 0 outside the cutoff (§5–6) |

## 11. Changes from draft 2 to draft 3 (review MG036 round 2)
| Review item | Change |
|---|---|
| 1 Deterministic R with random seams | Three objects defined (conditional map, order-averaged map, random-seam trajectory); B1 tests random-seam trajectories with three independent seam replicates per level; seeds per level and replicate; replicate spread S_k gates eligibility; production result SEAM-DEPENDENT when replicates disagree (§3.1–3.3, §4.2–4.3) |
| 2 Floor loopholes | Shape and drift labels mutually exclusive with INSUFFICIENT RESOLUTION as its own outcome; branch (b) renamed BELOW MEASUREMENT RESOLUTION (not convergence), with absolute bounds and the horizon check; horizon check compares R and J at the common time t = 600 (§3.3, §4.2–4.3) |
| 3 Part 0 benchmark | Conditional-order maps for every order, explicitly averaged; basis-vector matrices built two independent ways; augmented cohort benchmark for unwrapped moments from a localised population; refinement check frozen (N, t, τ values, Frobenius norm, ratios) (§5.1) |
| 4 Stochastic comparisons | Independent streams per arm; Welch interval for the drift difference; bootstrap resamples each arm independently; eligibility, weighting, minimum count and conditional reporting specified; J split-half baseline (§3.1, §3.3, §4.4) |
| 5 Stored inputs | Three-species fields at every decision checkpoint for R, L and J; unit-time centres; seeds, seam hashes, alignment and bootstrap settings (§7) |
| 6 A4 and slope claims | κ₃ called the leading odd term, not the only one; re-exposure labelled heuristic; slope interval called a regression diagnostic (§2 A4, §4.3) |

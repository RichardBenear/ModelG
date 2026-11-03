# MG7 — Does the coupled density response scale with the particle? (b-scan): design proposal, DRAFT 2

*Stage MG023a. Claude, 2026-11-03. Status: **DRAFT 2, not registered.** No MG7 run of any kind has been made. Draft 1 is
kept as `spec/MG7_proposal_v1.md` (superseded). This draft answers ChatGPT's review MG023 (`review/MG023_chatgpt.md`;
response `review/MG023_response.md`). The only computation behind it is static (`spec/MG7_design_calc.py` and ChatGPT's
`review/MG023_replication/static_review.py`, MG4c's stored settled states; no dynamics).*

*PI decisions (Richard Benear, 2026-11-03 13:50 UTC, project thread): adopt ChatGPT's recommendations for all six §10
decisions of draft 1. They are written in below and listed in §10.*

*Naming: MG6 §9 used "MG7" loosely for multiway/frame tests. The PI asked for this study as MG7, so multiway/frame tests
move to a later number.*

## 1. Question
MG6 (CLOSED, DEMONSTRATED) showed that the G-coupled birth rule produces a node deficit of about 67 nodes within r = 1.5
of a b = 28 particle, relative to a shifted-field null, within 0–17% of its quasistatic prediction, and that the particle
survives. MG6 deferred particle-mass scaling. MG7 asks:

> **Is the coupled-minus-null deficit larger for the stronger G-well at b = 28.5 than for the weaker one at b = 27.5?**

**What the gate tests, precisely.** Candidates arrive at rate νρ₀F and are accepted with probability exp(κφ_G)/F, so the
birth intensity is νρ₀·exp(κφ_G); with death rate ν the frozen-field target density is **ρ₀·exp(κφ_G)**. A stronger
G-well therefore gives a larger deficit almost by construction. The primary gate is a **directional endpoint result**: the
paired CP2 − SH2 deficit is larger at the higher endpoint. It does not test proportionality, monotonicity through b = 28,
or a negative contrast at either endpoint; those are reported separately (§5). Across this range the well integral grows
by about 46% while its zero crossing barely moves, so the change is a **stronger** well, not a wider one. A positive MG7 is
a consistency test of the coupled mechanism across particle strengths, with feedback present, not a discovery of a mass
law. Curvature and two-particle interaction are not asked here (§9).

## 2. Design calculation (static, from MG4c's settled states, κ = 2)
Node-sum approximation (draft 1), and ChatGPT's reproduction with the production detector (MG6 smoothed-X, θ = −4.45758)
and MG6 area quadrature:

| b | mean P₀(1.5), node sum | mean P₀(1.5), production quadrature | mean M₋ |
|---|---|---|---|
| 27.50 | −35.2 | −35.00 | 0.513 |
| 27.75 | −38.0 | −37.72 | 0.568 |
| 28.00 | −39.8 | −39.49 | 0.614 |
| 28.25 | −41.6 | −41.23 | 0.665 |
| 28.50 | −44.7 | −44.78 | 0.749 |

- **Endpoint change by start (production quadrature, 28.5 minus 27.5):** −3.94, −18.29, −10.48, −6.43 nodes. Negative on
  all four, but uneven: the weakest start (D0_p0, −3.9) decides whether all four signs survive under noise.
- All 20 archived states detect single under the common threshold and under MG4c's b-specific thresholds. The largest
  initial nodal φ_G is 0.0564, below 0.7·φ_lim(2) = 0.14.
- **Planning scenario, not demonstrated power:** in MG6 the coupled deficit was about 1.6 times the static start-field
  prediction, and the paired CP2 − SH2 SD at b = 28 was 3.7–3.9 nodes. That SD is for a contrast at one b, not for a
  difference across b; the static endpoint differences alone range over 14 nodes. Shared streams across b (§3) are chosen
  to keep the across-b noise down, but the power of the four-start test is not established in advance.
- The φ_G profile crosses zero at r ≈ 1.3–1.4 at every b, so **falloff will again be UNRESOLVED** (declared now) and is
  reported only as a profile.
- MG4c marks two stored states as not stationary (b = 28.25 D1_p0, b = 28.5 D0_p1). Both are single. MG7 uses
  b = 28.5 D0_p1; this is recorded.

## 3. Design
Everything not stated here is **MG6 as registered and amended (MG021)**: the coupled candidate process, F = e^0.4, the
full-field bound check, Δ = 0.1, dt = 0.05, ν = 0.2, t_end = 200, checks every 10 units, the measurements and quadrature of
MG6 §4, detection and classes, and the code (`src/mg6.py`, `src/run_mg6.py`), pinned at the MG6 production revision
`c795770` for everything MG7 reuses.
- **b values:** 27.5 and 28.5 (new runs). **MG6's recorded b = 28 CP2 and SH2 runs are the midpoint** (not rerun; MG6
  records at `09b72c2`). Both endpoints sit 0.23 inside MG4c's Delaunay single-particle window (about 27.27–28.73).
- **Start fields (review item 1).** Each run starts from **MG4c's own final array at the requested b**
  (`results/MG4c/fields/b{b}_D{slot}_p{pl}.npz`, float64), with no steady-state shift: MG6's `run_mg5.start_state` always
  loads b = 28 and shifts it, so MG7 has its own loader. Before launch the driver checks each file's SHA-256 against a
  committed manifest and that the loaded array equals the file; the reaction parameter passed to the ODE is asserted to be
  the requested b; every run key carries b (for example `CP2_b27.50_D0_p0`).
- **Detector threshold (review item 1):** fixed at MG6's θ = −4.45758 at both endpoints, so the observable is the same
  across b. MG4c's b-specific thresholds (−3.85703 at 27.5, −4.84027 at 28.5) are applied afterwards to the sparse records
  as a descriptive sensitivity check.
- **Arms (per endpoint):** CP2 (coupled, κ = 2) and SH2 (shifted-field null, κ = 2, |s| ≥ 6). SH2's agreement with DEC at
  b = 28 is not assumed at the endpoints; CP2 − SH2 answers the spatial-coupling question directly. No DEC, no κ = 4, no
  shuffled null (§10).
- **Streams (review item 5):** **one candidate stream and one shift per start, shared across both endpoints and both
  arms**: candidates `seed_for(7, 0, 0, 1, 2·slot + pl)`, shift `seed_for(7, 0, 1, 1, 2·slot + pl)`. Changing b then
  changes only the particle, not the random environment. Each start keeps its own stream. (The b = 28 midpoint used MG6's
  streams; it is reported, not part of the gate.)
- **Panel:** the four starts are two placements on each of two graphs. The conditional t interval is read on that fixed
  panel.

## 4. Measurements
As MG6 §4 per check. Run level: Y = window mean of D₀(1.5); Y_P = window mean of the in-run prediction P₀(1.5) (from the
evolving coupled field); the window mean of M₋; coarse and fine radial profiles. The **start-field prediction**
P₀^start(1.5) (from the stored start field, the only prediction available before the run) is computed per start and b in
Part 0 and recorded.

## 5. Outcomes and statistics
Working model as MG6 §5: four paired starts, independent, approximately normal, common variance; one-sided t, df = 3,
t ≤ −2.353; all intervals conditional on it. Window means are the unit; the eleven checks of a run are never pooled as
independent observations.
- **Contrast at each endpoint:** c_i(b) = Y_i(CP2, b) − Y_i(SH2, b).
- **Gate MG7 (primary, directional endpoint result):**
  - **LARGER AT HIGHER b** iff Δ_i = c_i(28.5) − c_i(27.5) < 0 for all four starts and the paired t ≤ −2.353.
  - **NOT DEMONSTRATED** otherwise (no equivalence claim).
  - **NOT SCORED (feedback)** if any CP2 or SH2 endpoint run is not single at every window check.
  - **INCOMPLETE** if any of the 16 runs is not COMPLETE (invalid, failed or incomplete); this precedes the feedback and
    significance rules, and the PI decides.
- **Reported separately, not scored:**
  - c(27.5), c(28) (MG6) and c(28.5), each with MG6's replication rule (all four negative and t ≤ −2.353), and the four
    within-start trajectories across the three b values.
  - **Agreement with prediction (review items 3 and 4):** the paired residual e_i = Δ_i − ΔP_i, with ΔP_i the matching
    change in **predicted CP2 minus predicted SH2**, reported in nodes with its mean and conditional 95% interval, twice:
    (a) with start-field predictions (pre-run), and (b) with in-run predictions (conditional on the realized wells; tests
    the local density response). No agreement label is awarded (§10). A ratio ΣΔ/ΣΔP is shown with a paired Fieller
    interval, reported as unbounded when the denominator is indistinguishable from zero.
  - Deficit against window-mean M₋ over the 8 endpoint CP2 runs plus MG6's 4: least-squares slope, descriptive, shown
    beside the four within-start trajectories.
  - Feedback per endpoint (SURVIVES if ≥ 3/4 CP2 runs single), mobility, falloff (UNRESOLVED expected), the b-specific
    threshold sensitivity check.
  - **Within-window stability (review §10 item 5):** for every CP2 and SH2 run, the early and late window means
    (t = 100–140 and 160–200) of D₀(1.5), their absolute difference, and the run's start-field prediction |P₀^start| as
    reference. Reported only; the window stays t = 100–200.
- **Invalidity and status precedence:** as MG6 (FAILED-NUMERICAL > INVALID-BOUND > INCOMPLETE > COMPLETE).

## 6. Runs, cost, budget and order
- **Runs:** 16 (2 endpoints × 2 arms × 4 starts).
- **Cost:** MG6 production averaged 82.6 (CP2) and 84.3 (SH2) seconds per time unit per run with 8 workers: about
  74–75 production ledger core-h for 16 runs, about 10 h of wall time on the PI's PC.
- **Order:** SH2 (both b), then CP2 (both b). Nulls first.
- **Budget:** cutoff **90 production ledger core-h**, run cap **22**, **8 workers**. The 6 spare slots may be used only to
  rerun a run that ended INCOMPLETE for a non-physical reason (a crash or interruption, not a physical outcome), each
  recorded with the reason. A physical outcome (multiplied, dissolved, mobile) is never replaced. If the restated
  projection exceeds 90, the driver writes `proceed = false` until the PI's decision is recorded.
- **Validation (Part 0) is priced separately:** about 25 core-h, cap 30.

## 7. Part 0 (before any registered run; review item 6)
MG6's validation was at b = 28. MG7 checks the endpoints on **unregistered** inputs: the MG6 Part 0 graph
(`seed_for(6, 9, 0, 0, 0)`), settled frozen at each endpoint b from the soliton (MG5's `run_ticks`, Δ = 1, dt = 0.05,
t = 300), and candidate streams `seed_for(7, 9, 3, 1, ib)`, ib = 1 (27.5), 2 (28.5), separate from production.
1. **Endpoint inputs:** SHA-256 of the 8 MG4c endpoint files against the committed manifest; each loaded array equals its
   file; each start detects single at θ = −4.45758; max nodal φ_G ≤ 0.14; P₀^start(1.5) and M₋ recorded per start.
   Any failure: STOP.
2. **Step check (forced history), per endpoint:** CP2, κ = 2, Δ = 0.1, 100 units at dt = 0.05; replay at dt = 0.025 with
   the accepted edits imposed. Criterion: equal particle counts at every check, depth ≤ 2%, |D̄₀(dt) − D̄₀(dt/2)| ≤ 2 nodes
   (D̄₀ = mean of D₀(1.5) over t = 50–100). Centre is recorded, not scored (MG021). Failure: STOP for the PI.
3. **Tick check (free acceptance), per endpoint:** CP2, κ = 2, Δ = 0.1 against Δ = 0.05, same stream, 100 units.
   Criterion: equal counts at every check, depth ≤ 3%, |D̄₀(0.1) − D̄₀(0.05)| ≤ **6 nodes** (about half the expected
   endpoint change; MG6's comparable differences were 1.3 and 2.3). Centre recorded, not scored. Failure: STOP for the PI.
4. **Bound margin:** max nodal φ_G ≤ 0.14 at every boundary of every Part 0 run.
5. **Timing:** from the Part 0 runs' wall times and MG6's production rates; projected cost restated in GATES.md. Every
   Part 0 run must be COMPLETE at its full horizon; an incomplete run is a STOP, not a cheap estimate.
Tolerances and stop rules above are frozen with this draft; outcomes are written to `results/MG7/part0/decisions.json`
before any registered run.

## 8. Data and recomputation
As MG6 §8: committed check records, ledger and the endpoint-file manifest; PC-only fields with a SHA-256 manifest;
`evaluate_mg7.py` and an independent `recompute_mg7.py` that compares the stored contrast vectors, t statistics,
intervals and residuals (MG022 lesson), strict on the records' platform and `--cross-platform` elsewhere. Run outputs are
written under immutable keys; an existing file is never overwritten without a hash check (MG021 lesson).

## 9. Not in MG7
Curvature (graph geometry around the particle) and two-particle interaction on the coupled graph (proposed as MG8); the
shuffled null (dropped, §10); ν = 0.02; κ = 4 at the endpoints; a resolved falloff law (needs a wider well than this
window gives).

## 10. Record of the PI decisions (draft 1 §10; ChatGPT MG023 recommendations adopted 2026-11-03 13:50 UTC)
| # | Decision | Adopted |
|---|---|---|
| 1 | Endpoints | 27.5 and 28.5, with MG6's recorded b = 28 runs as the midpoint; MG6 input revision pinned |
| 2 | Arms | CP2 and SH2 at each endpoint |
| 3 | Shuffled null | Dropped from MG7 |
| 4 | Part 0 | Verify endpoint files, detector and bound, plus the endpoint step and tick checks (§7); MG6's other calibration not repeated |
| 5 | Window | t = 100–200 fixed; within-window early/late differences reported against an explicit reference; no post hoc change |
| 6 | Budget | 90 production ledger core-h, cap 22, 8 workers; validation priced separately; projection above the cutoff needs a recorded PI decision; a physical outcome is never replaced |
| — | Agreement label (review item 4) | None: agreement with prediction is reported as residuals with intervals, without a label, unless the PI sets a tolerance before launch |

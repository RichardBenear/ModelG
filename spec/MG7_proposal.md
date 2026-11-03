# MG7 — Does the coupled density response scale with the particle? (b-scan): design proposal, DRAFT 1

*Stage MG023. Claude, 2026-11-03. Status: **DRAFT 1, not registered.** No MG7 run of any kind has been made. The only
computation behind this draft is the static table in §2 (`spec/MG7_design_calc.py`, MG4c's stored settled states; no
dynamics).*

*Naming: MG6 §9 used "MG7" loosely for multiway/frame tests. The PI asked for this study as MG7 (2026-11-03 13:30 UTC), so
multiway/frame tests move to a later number.*

## 1. Question
MG6 (CLOSED, DEMONSTRATED) showed that the G-coupled birth rule produces a node deficit around a b = 28 particle of about
67 nodes within r = 1.5, within 0–17% of its quasistatic prediction, and that the particle survives. MG6 deferred
"particle-mass scaling" (its §9). MG7 asks:

> **When the particle is made larger or smaller (by changing b inside the single-particle window), does the coupled
> deficit change with it, in the direction and roughly by the amount the static prediction says?**

**What a positive result means (stated now).** The birth rule sets the local target density to ρ₀·exp(κφ_G)/F. A deeper
G-well therefore gives a larger deficit almost by construction. A positive MG7 shows that the coupled process **tracks the
particle quantitatively across sizes, with feedback present**; it is a consistency test of the mechanism, not a discovery of
a mass law. The informative parts are whether the change matches the prediction, whether feedback alters it, and whether the
particle survives at both ends of the window. Curvature and two-particle interaction are not asked here (§9).

## 2. Design calculation (static, from MG4c's settled states; κ = 2, node-sum approximation)
| b | mean P₀(1.5) (nodes) | sd over 4 starts | mean M₋ | first zero of φ_G profile |
|---|---|---|---|---|
| 27.50 | −35.2 | 4.7 | 0.513 | 1.26–1.40 |
| 27.75 | −38.0 | 3.8 | 0.568 | 1.26–1.40 |
| 28.00 | −39.8 | 3.5 | 0.614 | 1.26–1.40 |
| 28.25 | −41.6 | 3.3 | 0.665 | 1.26–1.40 |
| 28.50 | −44.7 | 4.5 | 0.749 | 1.26–1.44 |

- Within each start, P₀(1.5) at b = 28.5 is 4–18 nodes deeper than at 27.5 (mean 9.5); M₋ grows by 46%.
- In MG6 the coupled deficit was about 1.6 times this static number (the well deepens under feedback), and the paired
  contrast SD was 3.7–3.9 nodes. If the same factor holds, the 27.5 → 28.5 change is about −15 nodes, against a paired SD of
  roughly 5: detectable with four starts.
- The φ_G profile crosses zero at r ≈ 1.3–1.4 at every b: the well does **not** get much wider. **Falloff will therefore
  again be UNRESOLVED** (declared now); it is kept only as a descriptive profile.
- MG4c marks two of these stored states as not stationary (b = 28.25 D1_p0, b = 28.5 D0_p1). They are still single; this is
  recorded, and settling under the MG7 rule is checked in the window rule (§5).

## 3. Design
Everything not stated here is **MG6 as registered and amended** (MG021): the coupled candidate process, F = e^0.4, the
bound check, Δ = 0.1, dt = 0.05, ν = 0.2, t_end = 200, checks every 10 units, window t = 100–200, the measurements and
quadrature of MG6 §4, detection and classes, the MG021 waiver scope, and the code (`src/mg6.py`, `src/run_mg6.py`).
- **b values:** 27.5 and 28.5 (new runs), with **MG6's b = 28 runs as the middle point** (CP2 and SH2, already recorded;
  not rerun). 27.5 and 28.5 sit inside MG4c's Delaunay single-particle window (about 27.27–28.73), 0.23 from its edges.
- **Starts:** MG4c's stored settled state at that b on D0_p0, D0_p1, D1_p0, D1_p1 (the same four graphs as MG6; the field
  differs by b).
- **Arms (per new b):** CP2 (coupled, κ = 2) and SH2 (shifted-field null, κ = 2, |s| ≥ 6). DEC is dropped: MG6 showed
  SH2 − DEC = +0.14 nodes (t = 1.0), so SH2 alone is the null. κ = 4 is dropped (MG6 covered it).
- **Shuffled null (PI decision, §10):** optional arm SF2 at b = 28 only: κ = 2, with the nodal φ_G values randomly permuted
  among nodes at each tick before acceptance. It keeps the field's value distribution and breaks its spatial link to the
  particle. It needs its own runs (it is not free; an earlier remark of mine said otherwise and was wrong).
- **Seeds:** candidates `seed_for(7, ib, 0, 1, 2·slot + pl)` with ib = 1 (b = 27.5), 2 (b = 28.5); shift
  `seed_for(7, ib, 1, 1, ·)`; shuffle `seed_for(7, 0, 2, 1, ·)`. CP2 and SH2 at one (b, start) share the candidate stream
  (paired, as MG6).

## 4. Measurements
As MG6 §4 per check, plus nothing new. Run-level Y = window mean of D₀(1.5); also the in-run P₀(1.5) window mean (Y_P),
the window mean of M₋, and the coarse and fine radial profiles.

## 5. Outcomes and statistics
Working model and test as MG6 §5 (four paired starts, one-sided t, df = 3, t ≤ −2.353; intervals conditional on it).
- **Contrast at each b:** c_i(b) = Y_i(CP2, b) − Y_i(SH2, b). At b = 28 these are MG6's recorded values.
- **Gate MG7 (primary): TRACKS** iff the paired change Δ_i = c_i(28.5) − c_i(27.5) is negative for all four starts and its
  t ≤ −2.353: the deficit is larger for the larger particle. **NOT DEMONSTRATED** otherwise (no equivalence claim).
  **NOT SCORED (feedback)** if any CP2 or SH2 run at 27.5 or 28.5 is not single at every window check.
- **Proportionality (reported, not scored):** the ratio Σ_i Δ_i / Σ_i ΔP_i, with ΔP_i the matching change in the in-run
  prediction Y_P(CP2); and the per-b observed/predicted ratios (MG6 gave 0.83–1.01 at b = 28). Described as "tracks the
  prediction" only if the ratio's conditional 95% interval includes 1.
- **Also reported:** c(b) at each b with its replication test (all four negative and t ≤ −2.353, as MG6's gate); the
  deficit against window-mean M₋ over the 12 CP2 runs (least-squares slope, descriptive); feedback outcome per b (SURVIVES
  if ≥ 3/4 single); mobility; falloff (UNRESOLVED expected); SF2 − SH2 at b = 28 if SF2 runs.
- **Window:** t = 100–200 for every run, fixed in advance. Each CP2 run's drift (blocks B₁ = {60, 70} and B₂ = {90, 100}
  of D₀(1.5), against 0.25·|Y_P|, as MG6's Part 0 rule) is reported, not used to move the window (PI decision, §10).
- **Invalidity, status precedence:** as MG6.

## 6. Runs, cost, budget and order
- **Runs:** 16 (2 new b × 2 arms × 4 starts), plus 4 if SF2 is included.
- **Cost:** MG6 production ran at 84 s per time unit per run with 8 workers: 4.7 core-h per run at t_end = 200.
  16 runs ≈ 75 core-h; 20 runs ≈ 93 core-h. About 10–13 h wall time on the PI's PC with 8 workers.
- **Order:** SH2 (both b), then CP2 (both b), then SF2. Nulls first.
- **Budget:** run cap 22 (or 26 with SF2); launch cutoff 90 core-h (110 with SF2), restated after Part 0.

## 7. Part 0 (before any registered run)
The numerics (Δ, dt, quadrature, bound, acceptance) were validated in MG6 Part 0 and its MG021 completion with the same
code, so MG7 does not repeat them (PI decision, §10). MG7's Part 0 is static and short:
1. **Start states:** each stored state at 27.5 and 28.5 is detected single; its centre and depth are recorded.
2. **Static predictions:** the §2 table recomputed with the MG6 quadrature (P₀(1.5) at κ = 2, M₋), per start.
3. **Bound margin (static):** max nodal φ_G of each start ≤ 0.7·φ_lim(2) = 0.14.
4. **Code:** the driver's queue order and seeds asserted by tests; the SF2 shuffle tested for determinism and for keeping
   the multiset of φ_G values.
5. **Timing:** one 10-unit coupled run at b = 28.5 with 8 concurrent copies, to confirm the 84 s per unit rate.

## 8. Data and recomputation
As MG6 §8: committed check records and ledger, PC-only fields with a SHA-256 manifest; `evaluate_mg7.py` and an
independent `recompute_mg7.py` (strict on the records' platform, `--cross-platform` elsewhere).

## 9. Not in MG7
Curvature (graph geometry around the particle) and two-particle interaction on the coupled graph: proposed as MG8, which
can show something not built into the rule. ν = 0.02 (dropped: it only slows the rewiring; MG5 covered that rate).
κ = 4 at new b. A resolved falloff law (needs a wider well than this window gives).

## 10. Decisions for the PI
| # | Decision | Recommendation |
|---|---|---|
| 1 | b values | 27.5 and 28.5, reusing MG6's b = 28 runs as the middle point |
| 2 | Arms | CP2 and SH2 only (drop DEC and κ = 4) |
| 3 | Shuffled null SF2 | Include at b = 28 only (4 runs, about 19 core-h) as a descriptive null, or drop it to save a night's third |
| 4 | Part 0 | Static checks plus one short timing run; no repeat of MG6's numerical checks |
| 5 | Window rule | Report drift; keep the window at 100–200 for all runs |
| 6 | Budget | Cutoff 90 core-h (110 with SF2), run cap 22 (26), 8 workers |

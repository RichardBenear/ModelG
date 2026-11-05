# MG9d — What drives MG9c's late drift? A bounded diagnostic study. Design proposal, DRAFT 3

*Stage MG031b; items 5–6 of §3 and §0's residual disclosure updated by PI decisions P0–P2 of the pre-run spec
(`spec/MG9d_prerun.md`, 01:50 UTC). Claude, 2026-11-05. Status: **DRAFT 3, not registered. No MG9d integration has been run.** Drafts 1
and 2 are kept as `spec/MG9d_proposal_v1.md` and `spec/MG9d_proposal_v2.md` (superseded). Draft 2 answered ChatGPT's
review MG031 (`review/MG031_chatgpt.md`; response `review/MG031_response.md`). Draft 3 makes the three decision-rule
corrections of ChatGPT's round 2 (`review/MG031_chatgpt_r2.md`; response `review/MG031_response_r2.md`) and freezes
E_pat at its exact recorded value. The design (runs, horizon, sampling, budget) is unchanged from draft 2. Requested by the PI (2026-11-05 00:30 UTC). MG9c is CLOSED (Q1 UNRESOLVED,
Q2 NOT ELIGIBLE), and nothing here changes its classifications, its stationarity threshold or its runs. Naming:
MG9b stays reserved for the stochastic barrier S(b).*

## 0. Starting point and disclosed prior knowledge

- **Committed source state.** The MG9c field at b = 29.30, t = 30,000: `results/MG9c/snapshots.npz` key
  `front_D02_29.30000_t30000`, SHA-256 `f22ca2807eee…` (full hash in `snapshots_manifest.json`). It lives on D02
  (h = 0.2, N = 500, L = 100). Its record gives E_pat = 6.60086115019687 (shown rounded as 6.600861 elsewhere) and edges x_L = −29.191, x_R = +29.011, with 18 cells.
- **Known before this draft (signed: negative means leftward).**
  - The patch translates at about −4.4e-7 per unit time at t ≈ 30,000, with nearly constant width.
  - The midpoint velocity at 29.30 has been decaying in magnitude: −8.2e-7 over 12k–15k, −7.1e-7 over 15k–20k,
    −5.8e-7 over 20k–25k and −4.8e-7 over 25k–30k.
  - ChatGPT's smooth exponential fit to that track (time constant about 26,000) would give about −3.5e-7 over the new
    run's window. That is an extrapolation, not a prediction.
- **Interpolation (corrected from draft 1; measured on the source state, no integration).**
  - A Fourier shift by exactly h reproduces `np.roll` to 6e-15.
  - Fractional shifts do **not** round-trip at round-off. The maximum field error is 7.8e-6 at h/4 and 3h/4 and
    1.6e-5 at h/2, because the even-grid Nyquist coefficient cannot be translated by a fraction of h on a real
    grid.
  - Its power fraction is 2.8e-11 (φ_X), 1.7e-11 (φ_Y) and 6e-14 (φ_G).
- **Filtered start, disclosed.** Because every arm starts from the Nyquist-filtered field (§2), MG9d studies the
  declared filtered preparation, not an exact continuation of the unfiltered MG9c state. The filter raises the initial
  max |RHS| from 1.09e-5 to 1.65e-3 (150×; corrected from "slightly" per pre-run item P0), carried by the removed,
  strongly damped Nyquist mode; the implied shape adjustment is 3.8e-6 relative. Both residuals are recorded and
  reported with the results.
- **Reflection** about the node x = 0 is the exact index map j → (N − j) mod N. Applied twice it gives the state back
  exactly.

## 1. Question and scope

Which explanation does the late drift favour: **unresolved numerical error**, **mesh alignment** (a lattice effect tied
to the patch's sub-grid position), **slow relaxation** left over from the preparation, or **persistent motion**? These
explanations are not mutually exclusive. A mesh effect and relaxation can coexist, and the reading reports both. This
is a bounded diagnostic. It is not a prerequisite for MG9b, and it does not try to settle every localization question.
**INCONCLUSIVE is an allowed outcome.** One mesh refinement supports a sensitivity statement, not a continuum limit.

## 2. Start states and runs (b = 29.30; t' = 0 at the start; horizon T = 10,000, no extension)

**Common start state y₀ (one preparation convention for every arm).** The source state with the Nyquist coefficient
of each field set to zero. This changes the fields by at most 1.56e-5. After this filter, fractional Fourier shifts
round-trip to 9e-15 (checked). Every arm starts from y₀ or an exact transform of it, so no arm is filtered alone.
y₀'s SHA-256 is recorded.

| Run | Start state | Mesh | Solver |
|---|---|---|---|
| **R0 baseline** | y₀ | D02 | registered: BDF, rtol 1e-8, atol 1e-10, max_step 0.5 |
| **R1 tight** | y₀ | D02 | rtol 1e-10, atol 1e-12, max_step 0.1 |
| **R2a/b/c shifted** | y₀ translated by s = +h/4, +h/2, +3h/4 (all three fields, Fourier) | D02 | registered |
| **R3 reflected** | y₀ reflected about node x = 0 (exact index map) | D02 | registered |
| **R4 finer mesh** | y₀ Fourier-interpolated to h = 0.1, N = 1000, L = 100 (zero-padding; same physical profile; not re-prepared) | D01 | registered |

- **Order:** R0, R1, R3, R2a, R2b, R2c, R4.
- **Solver wrapper.** MG9's `mg9_ode` with `max_step` exposed. At 0.5 it must reproduce `mg9_ode` bit for bit,
  tested on MG003's b = 29 reference before any run. seed_scale = 0.
- **Segments.** Calls of 1000 units, the same in every run.
- **State-transfer checks, before integrating (these replace draft 1's C0):**
  - the source hash matches the manifest;
  - y₀ differs from the source by at most 2e-5;
  - the R2 shift and the R3 reflection round-trip to ≤ 1e-12;
  - R4 sampled back at the D02 nodes equals y₀ to ≤ 1e-12, which zero-padding of a Nyquist-free field makes exact;
  - each start state's initial max|RHS| from an independent RHS evaluation (no stepping) is recorded beside y₀'s.
    Grid-dependent RHS differences are reported as part of the discretization comparison, not as errors.

  A failed transfer check stops the run before it is integrated.

## 3. Measurements and validity

**Diagnostics** at every 50 units (201 samples), on any grid (N and h are taken from the run). Frozen per run:
E_pat = 6.60086115019687, the exact value in `results/MG9c/runs/front_D02_29.30000.json` (field `E_pat`), used on
both grids; rounded decimals in this document are for display only and coordinates x_j = (j − N/2) h with the periodic cut at
x = ±50, which lies in the background.
- **d(t'), profile-alignment displacement (primary).** The translation s minimising ‖φ_X(t') − T_s φ_X(0)‖, where
  T_s is the Fourier shift on the run's grid.
  - The search is bounded to ±λ/4 around the previous sample's value (continuity), with xatol 1e-10.
  - An optimum within 1e-3 of a bound marks the sample **ALIGNMENT-AMBIGUOUS**.
- **C_edge and width W:** from the MG9 tracker's envelope half-level crossings.
- **C_amp and C_phi2:** the E²-weighted and φ_X²-weighted centres.
- **m(t'):** the mean shift of the interior minima matched at t' = 0, one subset for every statistic, with
  identities stored.
- **Shape change** S(t') = ‖φ_X(t') − T_{d(t')} φ_X(0)‖ / ‖φ_X(0)‖.

All displacements are measured from the run's own t' = 0 state, so imposed shifts are never counted as motion. For
R3, the stored d is raw, in R3's own coordinates. The reflection adjustment v_R3 = −v_R3,raw is applied once, in the
evaluator only.

**Measurement floor (validated before freezing; static, no integration).** Translate y₀ (and its D01 interpolant) by
known amounts: 0.001, 0.004, 0.01, h/4, h/2 and 3h/4, each in both signs. Recover each shift with all five
measures. δ is the largest recovery error of d over this set and both grids.
- The **speed floor** is v_f = max(5e-9, 4δ / 8000).
- A speed is **resolved** if |v| ≥ 3 v_f.

The δ, v_f and per-measure recovery errors are written into the pre-run specification before any run.

**Run validity (every run; required for any causal comparison).** A run is VALID if all of the following hold:
1. The integration is complete and all 201 diagnostic samples are present.
2. The tracker status is "ok" at every sample, with no fill, loss, untracked, BREAKUP or REMOTE.
3. The cell count stays within the start count ± 1, and |W(T) − W(0)| < λ/2. This is the identity check.
4. No sample is ALIGNMENT-AMBIGUOUS.
5. S_ref(T) = ‖φ_X(T) − T_{d(T)−d(2000)} φ_X(2000)‖ / ‖φ_X(2000)‖ ≤ 1e-3, measured from t′ = 2000 (PI decision
   P1, 2026-11-05: the shifted and finer-mesh starts need a shape adjustment to their new lattice, and the first 2000
   units are excluded to accommodate such transfer transients). S from t′ = 0 and the saved shape and width histories
   are reported, to show whether appreciable adjustment persists after 2000.
6. Measure consistency over [2000, T]: the change of each of C_edge, C_amp, C_phi2 and m agrees with d(T) − d(2000)
   within 10% of its magnitude + 5e-4 (PI decision P2: the static C_edge recovery error reaches 2.3e-4).

Otherwise the run is INVALID, and the failed items are listed.

**Speeds (primary measure d), endpoint-average estimator:**
a = (d(6000) − d(2000)) / 4000, b = (d(10,000) − d(6000)) / 4000, v = (d(10,000) − d(2000)) / 8000.
The first 2000 units are excluded to accommodate transfer transients. *Pre-run change, disclosed (ChatGPT MG031
round 3):* draft 1 used a least-squares slope for v. Draft 3 uses the endpoint average for every speed: it measures
net transport over each fixed window, and it is the estimator the floor formula v_f = 4δ/8000 assumes (two endpoint
errors of at most δ over 4000 units). For a curved or reversing trajectory it can differ from a least-squares fit. A
least-squares slope over each window is reported descriptively only.

## 4. Decision rules

**Speed comparison of A against reference B:**
- **AGREE:** |A − B| ≤ 0.05 |B| + v_f.
- **DIFFER:** |A − B| > 0.20 |B| + 3 v_f. The 3 v_f term means a tiny sign change near zero is not a difference.
- **INCONCLUSIVE:** otherwise.

| Comparison | Needs | Rule and labels |
|---|---|---|
| **C1 numerics** | R0 and R1 VALID (otherwise INCONCLUSIVE) | Compare a, b and v of R1 against R0. **NUMERICS-AGREE** if all three AGREE; **NUMERICS-SENSITIVE** if any DIFFER; otherwise INCONCLUSIVE. "AGREE" means agreement of the two tested configurations, not a bound on all numerical error. |
| **C3 reflection** | R0 and R3 VALID; otherwise **C3 NOT EVALUATED** (not an observed asymmetry) | **SYMMETRIC** if \|v_R3 − v_R0\| ≤ 0.01 \|v_R0\| + v_f, using the reflection-adjusted v_R3; otherwise **ASYMMETRY-FLAG**. |
| **C2 alignment** | R0 VALID | **ALIGNMENT-DEPENDENT** if any VALID R2 arm DIFFERs from R0. **ALIGNMENT-INDEPENDENT** if all three arms are VALID and AGREE with R0. Otherwise INCONCLUSIVE. The arms' shape-change and width histories are reported, to say whether any dependence is in translation speed or in deformation. |
| **C4 relaxation** | R0 VALID; when C1 is NUMERICS-AGREE it is applied to R0 and must give the same label on R1, otherwise INCONCLUSIVE | Applied in this order: **NO RESOLVED MOTION** if neither a nor b is resolved; **ONSET** if a is unresolved and b is resolved; **REVERSAL** if both are resolved with opposite signs; **DECREASED BELOW RESOLUTION** if a is resolved and b is unresolved (whatever b's sign; no ratio is formed). Otherwise both are resolved with the same sign, ρ = b/a > 0, and the label is **DECREASING** (ρ ≤ 0.90), **PERSISTENT** (0.95 ≤ ρ ≤ 1.05), **INCREASING** (ρ ≥ 1.10), or INCONCLUSIVE (otherwise). |
| **C5 finer mesh** | | Applied in this order: **R4 INVALID / LOST-IDENTITY** if R4 is INVALID; **NOT EVALUATED** if v_R0 is unresolved; **REVERSED** if v_R4 is resolved with the opposite sign to v_R0; **SIMILAR** if \|v_R4 − v_R0\| ≤ 0.20 \|v_R0\| + v_f; **REDUCED** if \|v_R4\| < 0.5 \|v_R0\|; **OTHER** otherwise. Mutually exclusive by order. |

**Overall reading (ordered; the first that applies; all component labels are always listed):**
1. **INCONCLUSIVE:** R0 INVALID, C3 ASYMMETRY-FLAG (a code or solver asymmetry), C3 NOT EVALUATED (the required
   symmetry check is unavailable), or C1 INCONCLUSIVE. Other component labels are still reported descriptively.
2. **NUMERICAL ACCURACY UNRESOLVED:** C1 NUMERICS-SENSITIVE.
3. **NO RESOLVED DRIFT:** a, b and v are all unresolved on R0 **and** on R1 (C4 is NO RESOLVED MOTION on both, and v
   is unresolved on both). The drift is below the measurement floor over this horizon. A small full-window slope v
   alone is not enough, because opposite motions in the two windows can cancel in v.
4. **MESH EFFECT SUPPORTED:** C2 ALIGNMENT-DEPENDENT, or C5 REDUCED or REVERSED. The reading says whether the
   dependence shows in translation or in deformation, and gives C4 alongside, since relaxation can coexist.
5. **RELAXATION SUPPORTED:** C2 ALIGNMENT-INDEPENDENT and C4 DECREASING or DECREASED BELOW RESOLUTION (the label is
   carried into the reading), with C5 SIMILAR. If C5 is R4 INVALID or
   NOT EVALUATED, the reading is qualified "on D02; the finer-mesh comparison is unavailable". An invalid finer mesh
   limits C5 only; it does not erase valid D02 evidence.
6. **MOVING STATE SUPPORTED:** C2 ALIGNMENT-INDEPENDENT and C4 PERSISTENT (a and b resolved, same sign), with C5
   SIMILAR. It is qualified as in 5 if C5 is unavailable.
7. **INCONCLUSIVE:** anything else (for example C5 OTHER, C4 REVERSAL, ONSET or INCONCLUSIVE, or C2 INCONCLUSIVE).
   The reading names the C4 component label, for instance "INCONCLUSIVE (C4 REVERSAL)". Valid mesh evidence (rule 4)
   is checked before this and is never hidden by a cancelling full-window slope.

**Worked cases (ChatGPT MG031 r2), to be checked in the implementation tests:**
- a = +4e-7, b = −4e-7, v ≈ 0, with C1 NUMERICS-AGREE, C3 SYMMETRIC and no mesh effect: C4 REVERSAL; rule 3 does not
  apply (a and b resolved); overall INCONCLUSIVE (C4 REVERSAL).
- R3 INVALID: C3 NOT EVALUATED; overall INCONCLUSIVE, with the other component labels reported.
- a = −4e-7, b = +1e-9 (b unresolved): C4 DECREASED BELOW RESOLUTION; no ratio is formed.

**Conditional wording, as the PI asked:**
- Alignment dependence *supports* a mesh effect.
- Solver sensitivity means numerical accuracy is unresolved.
- A decreasing drift *supports* relaxation, but does not prove eventual arrest.
- Persistent resolved motion that survives the numerical checks *supports* a moving-state reading, but does not prove
  an asymptotic travelling solution.
- C3 is a symmetry check. The equations and lattice are mirror-symmetric about a node, so a pass is expected and cannot
  separate physical drift from mesh relaxation; a failure points to the code.

## 5. Budget and saved fields

- **Cost estimate** (from MG9c timings):
  - R0, R2a–c and R3: about 3 min in total;
  - R1: 5–10 min, assuming 10–20× slower;
  - R4: 3–5 min;
  - in-run diagnostics: under 1 min per run.

  **Expected about 0.3 core-h. Strict cap: 10 integrations (7 planned + 3 numerical re-runs). Cutoff 1 core-h.**
  Cloud. If the cap or cutoff would be exceeded, work stops and the PI decides. Nothing is replaced, and no run is
  extended.
- **Saved and committed:**
  - y₀ and every start state, with SHA-256;
  - the state-transfer check values;
  - full three-field snapshots every 1000 units (11 per run, about 0.8 MB compressed);
  - per-run JSON with the tracker every 5 units, the diagnostics every 50 units (including the matched-minima
    identities and the ALIGNMENT-AMBIGUOUS flags), validity items, solver success, time reached and final residual;
  - the measurement-floor validation, the ledger and the evaluation.

  That is under 4 MB in total.
- **Before any run:**
  - code: `src/run_mg9d.py` and `src/evaluate_mg9d.py`, with the MG030b diagnostics generalised to any N and h;
  - tests (`tests/test_mg9d.py`): the wrapper match, shift and reflection maps, the Nyquist filter, the C1–C5 truth
    tables (including near-zero, reversal, undefined-ρ and overlap cases), the reflection sign applied once, and the
    overall order;
  - the measurement-floor validation and a short pre-run specification in GATES.md;
  - ChatGPT's implementation review, then PI approval.

## 6. Relevance to MG9b (stochastic barrier S(b))

**MG9b need not wait for a conclusive MG9d result. Its relevance depends on whether the affected numerics, mesh or
event-confirmation states are shared.** MG9b has no registered design yet. On the design sketched so far (MG3b's method,
on the h = 0.2 lattice), noise nucleates a single particle from the uniform state, and each event is confirmed by a
short deterministic run. The drift studied here belongs to an existing 18-cell patch.

- **NUMERICAL ACCURACY UNRESOLVED.** MG9b's stochastic tau-leap runs do not use this BDF continuation, but its
  deterministic confirmation step does use the same family of solver. This outcome should trigger a solver-tolerance
  check of that confirmation step at MG9b's b values, rather than an assumption that it is far from threshold.
- **MESH EFFECT SUPPORTED.** Localized states on h = 0.2 feel the lattice. MG9b would then measure a lattice-model
  barrier, and the window edges it relates to (b_s, b_r) may be partly set by the lattice. Since MG3b and MG4c are
  already lattice-conditional, the S(b) procedure can stay the same, but MG9b should state this and consider a
  confirmation check at h = 0.1.
- **RELAXATION SUPPORTED.** No expected effect: a slow adjustment of a prepared patch is a different state from
  nucleation out of the uniform state.
- **MOVING STATE SUPPORTED.** A caveat for extended patches near b_r. MG9b's end state is a single particle, which
  MG3/MG3b confirmed separately. MG9b's confirmation should still check particle stationarity at each of its own b
  values.
- **INCONCLUSIVE or NO RESOLVED DRIFT.** No effect. MG9b proceeds on its own terms.

## 7. Decisions for the PI (ChatGPT MG031's recommendations, written in)

1. **T = 10,000, sampling every 50 units, no extension.** Recommended.
2. **Cap 10 integrations and 1 core-h; C0 removed, with no restart-segmentation run added.** Recommended.
3. **Bands (5% / 20% + 3 v_f; 1% for reflection), C4 thresholds and validity limits,** kept as provisional
   operational choices. v_f is frozen from the static validation on the actual profiles and both grids before any run.
   Recommended.
4. **ChatGPT reviews this revision, then the pre-run implementation, before registration and runs.** Recommended.
   ChatGPT's round 2 said a further design cycle is unnecessary: the draft 3 corrections can be checked in the pre-run
   implementation review (full truth tables including the worked cases, and the actual-profile calibration on both
   grids).
5. **C4 DECREASED BELOW RESOLUTION counts toward RELAXATION SUPPORTED** (ChatGPT offered this or INCONCLUSIVE; either is
   defensible if frozen now). Recommended, because a resolved early drift falling below the floor is the plainest form
   of slowing; the label is carried into the reading so it is not mistaken for a measured ratio.

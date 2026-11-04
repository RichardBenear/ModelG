# MG9c — Late-time fate of MG9's slow fronts: do bounded patterns settle? Design proposal, DRAFT 2 (REGISTERED)

*Stage MG029a. Claude, 2026-11-04. Status: **REGISTERED 2026-11-04 (MG029c), draft 2, content frozen at `372336b`; PI decision card 19:56 UTC ("Register draft 2", with the §9 defaults).** No MG9c run of any kind preceded registration. ChatGPT round 2 (`review/MG029_chatgpt_r2.md`): ready for registration; its one wording correction (§10 MIXED) is applied. Draft 1 is kept as `spec/MG9c_proposal_v1.md`
(superseded). This draft answers ChatGPT's review MG029 (`review/MG029_chatgpt.md`; response
`review/MG029_response.md`). All four items and all four PI recommendations are accepted into the text below.*

*Origin: ChatGPT's next-step recommendation after MG9 (PI relay, 2026-11-04 19:33 UTC) and the PI's "yes, go ahead"
(19:36 UTC). MG9 is CLOSED and its classifications stay as recorded (P1 SUBCRITICAL, P2 UNRESOLVED). MG9c has its own
criteria and labels, and does not re-score MG9. Naming: MG9b stays reserved for the stochastic barrier S(b); MG8
for multiway/foliation.*

## 0. Disclosed prior knowledge (looked at after MG9 closed; no new integration)

Claude fitted the **committed MG9 tracks** (`results/MG9/runs/front_grid_29.{2500,3000,3500}.json`; interfaces every
5 units over [0, 1800]). ChatGPT reproduced the fits independently (`review/MG029_replication/`).

| b | Outward shift over [100, 1800] (L, R) | Shift per 200 units, t = 200 → 1800 | Exponential fit on [600, 1800]: time constant; shift still to come | rms: exponential / straight line |
|---|---|---|---|---|
| 29.25 | 2.217, 2.221 | 0.50, 0.38, 0.29, 0.22, 0.17, 0.13, 0.10, 0.07 | 718–722; 0.23 | 2.5e-4 / 5.8e-2 |
| 29.30 | 2.227, 2.231 | 0.51, 0.38, 0.29, 0.22, 0.17, 0.13, 0.10, 0.07 | 712–717; 0.22 | 2.5e-4 / 5.8e-2 |
| 29.35 | 2.238, 2.243 | 0.52, 0.39, 0.29, 0.22, 0.17, 0.13, 0.09, 0.07 | 708–713; 0.22 | 5e-4 / 5.8e-2 |

- **Cell-count history (corrected from draft 1, per MG029 item 1).** The detector count stayed at 18 throughout at
  29.25 and 29.30. At 29.35 it was 18 initially and **19 from t = 45 to 1800**. This is a detector count; on its own
  it does not show that a cell was physically born. The count is constant over [600, 1800] in all three runs, which is
  the interval the fit uses.
- An exponential approach to a fixed edge fits about 200 times better than constant speed. Extrapolated, about 0.22
  (≈ λ/13) is left to move. The predicted displacement over [6000, 12,000] is 0.0006–0.0007, against λ/40 = 0.0735.
- The final residual of 7.6e-3 is consistent with a slow relaxation still under way at t = 1800. Whether the motion is
  an edge-cell relaxation is left to the snapshots (§3); draft 1's reading is not assumed.
- At 29.40 the front filled the domain (INVADES).

**So the existing data favour settling, without proving it.** MG9c is a confirmation with new, later-time data, plus
robustness checks MG9 never reached. It is not a discovery claim.

## 1. Questions

- **Q1 (primary).** At the three tested values b = 29.25, 29.30 and 29.35 (D02; patches prepared from the 34-cell
  periodic branch with MG9's splice), does each front **settle** within the registered tolerances, or keep
  **advancing**? Does a settled patch pass the specified mesh, domain and perturbation checks?
- **Q2 (secondary; eligible only if 29.35 is SETTLED-ROBUST in Q1).** For this preparation, give an operational bracket
  between 29.35 and 29.40 where fronts stop settling, and compare it with the single-particle replication threshold b_r
  within a registered tolerance.
- **Descriptive.** The same test at b = 29.22, the lowest b where the 34-cell branch persists. It probes the lowest
  available preparation, and is not an estimate of the lower localization edge.

## 2. Predictions (informed by §0)

- **Q1:** SETTLED at all three b, and SETTLED-ROBUST.
- **Descriptive check, not a criterion:** final edges within λ/10 of the §0 extrapolation.
- **Q2:** no direction is predicted. MG9 showed only that the change from slow/unresolved behaviour to invasion lies
  between 29.35 and 29.40. It did not show a switch from settling to invasion there.

## 3. Frozen machinery, runs and records

- **Frozen MG9 machinery (unchanged; main at `4466be3`, MG9 production code `f7ef9e9`).**
  - `run_mg9.mg9_ode`: BDF, rtol 1e-8, atol 1e-10, MG9's Jacobian and first-step retries, **seed_scale = 0**.
  - `envelope`: MG9 amendment 1, the complex demodulation with a double λ-boxcar (MG027).
  - `splice_patch`, `phase_convention`, `Tracker` (components, BREAKUP/REMOTE, interfaces), `cell_count`,
    `carry_concentration`, `fourier_resample`.
  - Tolerances STAT_TOL = 1e-6 and ARREST_RES_TOL = 1e-5; λ = 100/34; interfaces sampled every 5 units.
- **Same start as MG9.** MG9's prepared state, phase convention and splice: half-width a = 25, ℓ = λ/2, uniform state
  elsewhere. States come from `results/MG9/states.npz`: `prep_29.2500`, `sweep_29.30`, `prep_29.3500` and
  `sweep_29.22`. E_pat is computed from the prepared state, as in MG9.
- **Platform.** The runs use the cloud environment MG9 recorded (Python 3.13.18, numpy 2.5.3, scipy 1.18.1). Versions
  are checked before the first run. If they differ, the PI decides before any run.
- **Replay check.** Each of the three main runs integrates from t = 0. Its track over [0, 1800] must reproduce the MG9
  record to 1e-6 in x_L and x_R.
  - A mismatch counts against the budget. The point is then **REPLAY-MISMATCH**, treated as UNRESOLVED until the cause
    is found and the PI accepts it.
  - 29.22 has no MG9 front track and is exempt.
- **Horizons.**
  - Main, D01 and D02L front runs: T = 12,000. A run classified rule 6 (§4) at 12,000 is continued, as the same
    integration, to T = 30,000 and reclassified.
  - Perturbation runs: a fixed 6000 units with no extension (§5c).
- **Records:**
  - the complete track history (interfaces, cell count, components) for every run;
  - initial-state SHA-256, E_pat, solver success, time reached and final residual;
  - **full-field snapshots** (all three fields) at every scoring-window boundary and every actual final horizon:
    - main, D01 and D02L runs: t = 1800, 3600, 6000 and 12,000, plus 15,000 and 30,000 if extended;
    - perturbation runs: their own boundaries 1500, 3000 and 6000.

  The snapshots go in one float64 `.npz`, estimated 3–6 MB, committed (they are the audit trail MG028 lacked).

## 4. Classification (ordered; the first rule that applies)

Late-time windows: W₁ = [1800, 3600], W₂ = [3600, 6000], W₃ = [6000, 12,000]. If extended, W₂ = [6000, 15,000] and
W₃ = [15,000, 30,000]. D_k is an interface's outward displacement over W_k, and s_k = D_k / |W_k| its mean speed.

1. **UNRESOLVED-NUMERICAL:** the run is not complete, or the tracked patch cannot be identified at t = 0.
2. **Terminal events, from t = 0 (as MG9 §6.3 rule 2):**
   - fill first occurring at t_f and holding to T → **INVADES**;
   - loss holding to T → **RETREATS**;
   - fill or loss that later recovers → UNRESOLVED-TEMPORARY.
3. **Structural disqualifiers, on [100, T] (as MG9 §6.3 rule 3):**
   - BREAKUP or REMOTE → UNRESOLVED-INTERNAL;
   - a reversal of at least λ/2 by either interface → UNRESOLVED-REVERSAL.
4. **Persistent late motion, in W₃:**
   - D_L ≥ λ/4 and D_R ≥ λ/4 → **ADVANCES**;
   - D_L ≤ −λ/4 and D_R ≤ −λ/4 → **RECEDES**;
   - otherwise, if max(|D_L|, |D_R|) ≥ λ/4 → UNRESOLVED-ASYMMETRIC. This covers one moving interface, and opposite
     signs (translation of the patch).
5. **SETTLED:** all of the following hold.
   - |D₃| ≤ λ/40 for both interfaces.
   - |v| ≤ λ/40,000 per unit time, v being the least-squares speed over the last 1000 units.
   - Decreasing speed for both interfaces: |s₃| < |s₂| < |s₁|, or |D₃| ≤ λ/400.
   - Final residual ≤ 1e-5.
   - Cell-count range over [W₁ start, T] at most 1.
6. **Otherwise:** extend once (§3), then **UNRESOLVED-SLOW**.

Displacement before W₁ never enters rules 4–5. This is the explicit split between initial adjustment and persistent
late-time motion. Rules 2–3 keep MG9's windows, so an early fill, breakup or reversal is not hidden.

Synthetic tests (before any run) cover:
- exponential settling and constant slow speed;
- power-law creep at the thresholds;
- an early stop;
- one interface moving, and **opposite-sign translation, e.g. (+0.3 λ, −0.3 λ)**;
- early fill and early loss, temporary loss, breakup, and a λ/2 reversal;
- the perturbation-run variant (fixed horizon, rule 6 → UNRESOLVED-SLOW).

## 5. Robustness of a SETTLED point

Each trial gets one of three statuses:
- **PASS:** it meets the criterion below.
- **FAIL:** a resolved dynamical label other than the one required (ADVANCES, RECEDES, INVADES or RETREATS), or, for
  (c), SETTLED with edges more than λ/10 from the unperturbed edges.
- **UNRESOLVED:** an incomplete integration; a preparation that is not PERSISTS (PREP-UNRESOLVED, PREP-CHANGED or
  PREP-LOST); or any UNRESOLVED-* label.

The trials:

- **(a) Finer mesh, D01 (h = 0.1).** MG9's D01 path: Fourier-interpolate the D02 prepared state, settle on D01 to
  t = 1000, and require PERSISTS (34 cells, residual ≤ 1e-6, extension to 3000 as in MG9). Then splice with
  half-width 25 and front-test with §4. PASS = SETTLED.
- **(b) Domain control, D02L (L = 200).** MG9's robustness path: tile twice, settle to 300, and require PERSISTS
  (68 cells, residual ≤ 1e-6). Then splice with **half-width 25** (the same patch as on D02) and front-test with §4.
  PASS = SETTLED. This isolates the domain length.
- **(c) Small perturbations.** The settled final field φ(T) of the main run, multiplied by 1.05 and separately by 0.95
  (all three fields, as potentials about the uniform state). Each runs a fixed 6000 units. It is classified by §4 with
  W₁ = [0, 1500], W₂ = [1500, 3000] and W₃ = [3000, 6000], with no extension; rule 6 gives UNRESOLVED-SLOW directly.
  PASS = SETTLED with final edges within λ/10 of the unperturbed final edges. These two uniform rescalings test only
  these two perturbations, not general stability.
- **Descriptive only, no criterion:**
  - **(d) Width:** on D02L, the same settled preparation with half-width 50 (MG9's 100-unit patch);
  - **(e) Width:** on D02, half-width 25 + λ/4.

  Both record whether settling depends on patch size and where the settled width locks.

**Point status:**
- **SETTLED-ROBUST:** (a), (b) and (c) all PASS.
- **NOT-ROBUST:** any of them FAILS.
- **ROBUSTNESS-UNRESOLVED:** otherwise.

An UNRESOLVED trial is never read as evidence of mesh or domain dependence.

## 6. Q2 — operational upper settling bracket (only if 29.35 is SETTLED-ROBUST)

- **Assumption, stated.** Bisection assumes a locally ordered transition for this preparation: settling below, invasion
  or advance above. All sampled labels are reported. The result is an operational bracket for this preparation, not a
  unique global settling edge.
- **Bisection.** Three steps between 29.35 and 29.40 (MG9: INVADES), down to a bracket width of 0.00625.
  - Each midpoint is prepared by MG9 §6.1: carry from the nearest sweep state above, settle to 1000, and require
    PERSISTS. A failed preparation stops the bisection, and the bracket is unresolved.
  - Each midpoint is then classified by §4. SETTLED raises the lower end. INVADES or ADVANCES lowers the upper end.
    Anything else stops the bisection, and the bracket is unresolved.
- **Final point.** The final settled point gets robustness (a)–(c). If it is not SETTLED-ROBUST, the bracket is
  unresolved.
- **Comparison (interval arithmetic).** Settling bracket [b₁, b₂]; MG9's rechecked b_r ∈ [r_lo, r_hi] =
  [29.395313, 29.396094]. The possible difference is Δ ∈ [b₁ − r_hi, b₂ − r_lo]. Tolerance 0.01.
  - **NEAR b_r (within the registered tolerance):** the whole Δ interval lies in [−0.01, +0.01];
  - **BELOW b_r:** b₂ − r_lo < −0.01;
  - **ABOVE b_r:** b₁ − r_hi > +0.01;
  - otherwise **UNRESOLVED RELATIVE TO THE TOLERANCE**.

  Draft 1's claim that the edge cannot be found above b_r was wrong (MG029 item 2): an edge at 29.398 lies above b_r and
  below 29.40. The interval rule handles that case.
- **Scope.** 0.01 is a descriptive proximity tolerance. It is not a statistical interval, and not a test that the
  two thresholds are physically the same.
- **D01 (descriptive).** The final point's D01 run is set beside D01's own b_r = 29.412.

## 7. Outcomes (ordered)

| Order | Q1 label | Criterion |
|---|---|---|
| 1 | **UNRESOLVED** | any of 29.25 / 29.30 / 29.35 gives an UNRESOLVED-* label or REPLAY-MISMATCH, or a SETTLED point is ROBUSTNESS-UNRESOLVED |
| 2 | **SLOW INVASION** | all three give ADVANCES or INVADES |
| 3 | **BOUNDED** | all three are SETTLED-ROBUST |
| 4 | **MIXED** | otherwise (some settle and some advance, or a point is NOT-ROBUST); listed point by point |

- The Q2 label (§6) is reported only when 29.35 is SETTLED-ROBUST.
- 29.22 is descriptive.
- **BOUNDED** means finite-time settling within the registered tolerances, with the specified robustness checks
  passed. Its wording: "at b = 29.25, 29.30 and 29.35 (h = 0.2), patches prepared from the 34-cell periodic branch
  settle into bounded patches coexisting with the uniform state, within the registered tolerances, and pass the
  specified finer-mesh, longer-domain and two amplitude-perturbation checks". It is not "stable" in general, not
  "pinning demonstrated", and says nothing about the lower window edge.

## 8. Cost, cap and code

- **Cost.** MG9's fronts took about 5 s per 1800 units on D02, so a 12,000-unit D02 run takes about 30 s. D01 and D02L
  are assumed 4× slower.

| Item | Integrations |
|---|---|
| Q1 fronts (4 points; extensions count once) | 4 |
| Robustness, 4 points × [D01 settle + front, D02L settle + front (hw 25) + front (hw 50), 2 perturbations, 1 D02 width] | 32 |
| Q2: 3 midpoints × (settle + front) + final-point robustness (D01 2, D02L 2, perturbations 2) | 12 |
| **Planned** | **48** |
| Reserve (numerical re-runs) | 12 |
| **Run cap** | **60** |

  **Expected under 0.5 core-h. Cutoff 2 core-h.** It runs in the cloud, resumably. No PC time.
- **Stop rule.** If the cap or cutoff would be exceeded, the work stops and the PI decides. The order of work is Q1
  fronts, then robustness (a)–(c), then Q2, then descriptive trials (d)–(e).
- **Code.**
  - `src/run_mg9c.py` imports the frozen MG9 functions (§3) unchanged, and adds the §4 classifier and the §§5–6
    runners.
  - `src/evaluate_mg9c.py` computes the labels.
  - `tests/test_mg9c.py` contains the §4 synthetic cases, the §5 status logic, the §6 bracket and interval logic
    (including the 29.398 case), and the §7 order.
- **Records.** Those listed in §3, plus a ledger and the evaluation.
- **Pre-run specification.** Anything still open is written into GATES.md before the first run.

## 9. Decisions for the PI (ChatGPT MG029's recommendations, accepted into the text)

1. **Horizon and thresholds:** T = 12,000 (extension to 30,000) and the W₃ limits λ/40 and λ/40,000, with the
   clarified rules. No finite horizon excludes arbitrarily slow creep; the wording stays operational.
2. **Q2 tolerance:** 0.01, as a descriptive proximity tolerance, using the interval rule and the NEAR wording.
3. **29.22:** included, descriptively.
4. **Review:** ChatGPT checks this revision before registration.

## 10. What each outcome would mean (scope-limited)

- **BOUNDED.** For these preparations at the three tested b, localized patches coexist with the uniform state just
  below b_r for at least the registered horizon, and pass the specified checks. This is what a pinning picture of the
  window would need, but it does not establish one. The lower window edge stays untested on this branch.
- **SLOW INVASION.** These preparations advance: MG9's low speeds were a long transient. That does not rule out
  localization elsewhere, by pinning on another branch or by wavelength selection. A branch-selection study (the
  30-cell family) would be the natural next step.
- **MIXED.** The three points do not jointly meet BOUNDED or SLOW INVASION, with no unresolved condition taking
  precedence. Each dynamical outcome and robustness status is reported; no unique transition is inferred from the
  aggregate label alone. *(Wording per ChatGPT MG029 r2.)*
- **NEAR b_r (Q2).** Two operational thresholds agree within 0.01: the end of front settling for this preparation, and
  single-particle replication. This is consistent with a shared mechanism, but does not establish depinning or
  homoclinic snaking. **BELOW** or **ABOVE** would mean they differ by more than the tolerance.

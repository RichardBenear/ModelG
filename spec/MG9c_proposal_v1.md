# MG9c — Late-time fate of MG9's slow fronts: do bounded patterns settle? Design proposal, DRAFT 1 (SUPERSEDED by draft 2, `spec/MG9c_proposal.md`)

*Stage MG029. Claude, 2026-11-04. Status: **DRAFT 1, not registered.** It follows ChatGPT's next-step recommendation
after MG9 (relayed by the PI, 2026-11-04 19:33 UTC) and the PI's "yes, go ahead" (19:36 UTC). MG9 is CLOSED and its
classifications stay as recorded (P1 SUBCRITICAL, P2 UNRESOLVED). MG9c has its own criteria and labels. It does not
re-score MG9.*

*Naming: MG9b stays reserved for the stochastic barrier S(b). MG8 stays reserved for multiway/foliation.*

## 0. Disclosed prior knowledge (looked at after MG9 closed; no new integration)

Before this draft, Claude fitted the **committed MG9 tracks** (`results/MG9/runs/front_grid_29.{2500,3000,3500}.json`,
interfaces sampled every 5 units over [0, 1800]). No new integration was made.

| b | Outward shift over [100, 1800] (L, R) | Shift per 200 units, from t = 200 to t = 1800 | Exponential fit on [600, 1800]: time constant, shift still to come | rms: exponential / straight line |
|---|---|---|---|---|
| 29.25 | 2.217, 2.221 | 0.50, 0.38, 0.29, 0.22, 0.17, 0.13, 0.10, 0.07 | 718–722; 0.23 | 2.5e-4 / 5.8e-2 |
| 29.30 | 2.227, 2.231 | 0.51, 0.38, 0.29, 0.22, 0.17, 0.13, 0.10, 0.07 | 712–717; 0.22 | 2.5e-4 / 5.8e-2 |
| 29.35 | 2.238, 2.243 | 0.52, 0.39, 0.29, 0.22, 0.17, 0.13, 0.09, 0.07 | 708–713; 0.22 | 5e-4 / 5.8e-2 |

- The cell count stayed at 18 throughout, in all three runs. The patch did not gain cells. What moved is where the
  edge cells cross the half-envelope level.
- An exponential approach to a fixed edge fits about 200 times better than constant speed. Extrapolated, each edge
  has about 0.22 (≈ λ/13) left to move.
- The final field residual, 7.6e-3, is consistent with the same slow relaxation still under way at t = 1800.
- At 29.40, the front filled the domain (INVADES).

**So the existing data already favour settling.** MG9c is therefore a confirmation with new data: much longer times,
plus robustness checks that MG9 never reached. It is not a discovery claim. The prediction in §2 is informed by the fit
above, and is labelled as such.

## 1. Questions

- **Q1 (primary).** At b = 29.25, 29.30 and 29.35 (D02, 34-cell pattern, the MG9 preparation and splice), does each
  front **settle** into a stationary bounded patch, or keep **advancing**? And does a settled patch survive a finer
  mesh, a longer domain and small perturbations?
- **Q2 (secondary, conditional on Q1).** If 29.35 settles robustly, where between 29.35 and 29.40 do fronts stop
  settling? How does that compare with the single-particle replication threshold b_r?
- **Descriptive.** The same test at b = 29.22, the lowest b where the 34-cell pattern persists, gives the lower extent
  available on this branch. It is branch-limited, not an edge.

## 2. Predictions (informed by §0)

- **Q1:** SETTLED at all three b, and robust (§5).
- **Descriptive check, not a criterion:** final edges within λ/10 of the §0 extrapolation.
- **Q2:** no direction is predicted. MG9 already showed that the switch lies between 29.35 and 29.40, and that this
  interval contains b_r.

## 3. Runs (D02 unless stated; MG9's code, settings and tolerances)

- **Same start as MG9.** Each run uses MG9's prepared state, phase convention and splice: patch half-width a = 25,
  ℓ = λ/2, uniform state elsewhere. Prepared states come from `results/MG9/states.npz`: `prep_29.2500`, `sweep_29.30`,
  `prep_29.3500`, and `sweep_29.22`.
- **Reproducibility check.** Each run integrates from t = 0. Its track over [0, 1800] must reproduce the MG9 record to
  1e-6 in x_L and x_R, on the same platform class. A mismatch is reported, and the run still counts.
- **Horizon.** T = 12,000. A run classified rule 6 (§4) at T = 12,000 is continued, as the same integration, to
  T = 30,000 and reclassified.
- **Records:**
  - interfaces and the cell count every 5 units, as in MG9;
  - **full-field snapshots** (all three fields) at t = 1800, 3600, 6000 and 12,000, plus 30,000 if reached. They are
    committed as one float64 `.npz`, about 2–3 MB.

## 4. Classification (ordered; the first rule that applies)

Late-time windows: W₁ = [1800, 3600], W₂ = [3600, 6000], W₃ = [6000, 12,000]. If extended, W₃ = [15,000, 30,000] and
W₂ = [6000, 15,000]. D_k is an interface's outward displacement over W_k, and s_k = D_k / |W_k| is its mean speed.

1. **UNRESOLVED-NUMERICAL:** the run is not complete, or the tracked patch cannot be identified at t = 0.
2. **Terminal and structural rules on [100, T]**, exactly as MG9 §6.3 rules 2–3:
   - fill holding to T → **INVADES**;
   - loss holding to T → **RETREATS**;
   - a fill or loss that later recovers → UNRESOLVED-TEMPORARY;
   - BREAKUP or REMOTE → UNRESOLVED-INTERNAL;
   - a reversal of at least λ/2 → UNRESOLVED-REVERSAL.
3. **Persistent late motion (in W₃):**
   - both D ≥ λ/4 → **ADVANCES** (sustained slow invasion);
   - both D ≤ −λ/4 → **RECEDES**;
   - one interface with |D| ≥ λ/4 and the other not → UNRESOLVED-ASYMMETRIC.
4. **SETTLED:** all of the following hold.
   - For both interfaces, |D₃| ≤ λ/40.
   - |v| ≤ λ/40,000 per unit time, where v is the least-squares speed over the last 1000 units.
   - The speed decreases: |s₃| < |s₂| < |s₁|, or |D₃| ≤ λ/400, which allows for a front that stopped early.
   - The final residual is at most 1e-5, MG9's arrest tolerance.
   - The cell-count range over [1800, T] is at most 1.
5. *(No rule 5. The numbering mirrors MG9.)*
6. **Otherwise:** extend once (§3), then **UNRESOLVED-SLOW**.

Displacement before t = 1800 never enters rule 4. This is the explicit split ChatGPT asked for between initial
adjustment and persistent late-time motion. Rules 1–2 still use [100, T], so an early fill, breakup or reversal is not
hidden.

## 5. Robustness of a SETTLED point (all must pass for **SETTLED-ROBUST**)

- **(a) Finer mesh, D01 (h = 0.1).** MG9's D01 path: Fourier-interpolate the D02 prepared state, settle on D01 to
  t = 1000 (a PERSISTS check), splice, front-test to T = 12,000 with §4. It must give SETTLED.
- **(b) Longer domain, D02L (L = 200).** MG9's robustness path: tile twice, settle to 300 (68 cells), splice a patch of
  half-width 50, front-test with §4. It must give SETTLED.
- **(c) Small perturbation.** The settled final field φ(T) is multiplied by 1.05, and separately by 0.95 (all three
  fields, as potentials about the uniform state). Each is run for 6000 units and classified by §4, with windows
  W₁ = [0, 1500], W₂ = [1500, 3000], W₃ = [3000, 6000]. Both must give SETTLED, with final edges within λ/10 of the
  unperturbed final edges.
- **(d) Descriptive only: patch width.** Starts of half-width 25 ± λ/4. These record whether the settled width locks
  to the same or a neighbouring cell boundary. No criterion.

## 6. Q2 — upper settling edge (only if 29.35 is SETTLED-ROBUST)

- **Bisection.** Three steps between 29.35 (settled) and 29.40 (MG9: INVADES), down to a bracket width of 0.00625.
  Each midpoint is prepared by MG9 §6.1: carry from the nearest sweep state above, settle to 1000, PERSISTS check.
  Then it is classified by §4.
  - A SETTLED midpoint moves the lower end up.
  - INVADES or ADVANCES moves the upper end down.
  - Anything else stops the bisection, and the edge is **unresolved**.
- **Robustness of the final point.** The final settled point gets robustness (a)–(c). If it fails, the edge is
  unresolved.
- **Comparison.** The final bracket [b₁, b₂] is compared with the MG9 rechecked b_r ∈ [29.395313, 29.396094]:
  - **AT b_r:** b₁ ≥ b_r⁻ − 0.01 (the settling edge lies within 0.01 below b_r, or at it);
  - **BELOW b_r:** b₂ < b_r⁻ − 0.01;
  - otherwise **STRADDLES**.

  The test is one-sided by construction: 29.40, already at b_r + 0.004, invades, so the edge cannot be found above b_r.
- **Tolerance.** 0.01 is a PI decision (§9). MG9 used 0.10, which here would be uninformative.
- **D01 (descriptive).** The D01 run at the final settled point is set beside D01's own b_r = 29.412.

## 7. Outcomes (ordered)

| Order | Q1 label | Criterion |
|---|---|---|
| 1 | **UNRESOLVED** | any of 29.25 / 29.30 / 29.35 gives an UNRESOLVED-* label |
| 2 | **SLOW INVASION** | all three give ADVANCES or INVADES |
| 3 | **BOUNDED** | all three are SETTLED-ROBUST |
| 4 | **MIXED** | otherwise (some settle, some advance, or a robustness check fails); listed point by point |

The Q2 label (§6) is reported only after BOUNDED, or after MIXED with 29.35 SETTLED-ROBUST. 29.22 is descriptive.

Wording for BOUNDED: "stationary bounded 34-cell patches coexist with the uniform state at b = 29.25–29.35 (h = 0.2),
robust to mesh, domain and small perturbation". It is not "pinning demonstrated", and not an explanation of the lower
window edge, which this branch cannot reach.

## 8. Cost, cap and code

- **Cost.** MG9's fronts took about 5 s per 1800 units on D02, so a 12,000-unit D02 run takes about 30 s. D01 and D02L
  are assumed 4× slower.

| Item | Integrations |
|---|---|
| Q1 fronts (4 points, extensions count once) | 4 |
| Robustness: 4 points × (D01 settle + front, D02L settle + front, 2 perturbations, 2 widths) | 32 |
| Q2: 3 midpoints × (settle + front) + robustness of the final point (6, widths not repeated) | 12 |
| **Planned** | **48** |
| Reserve (numerical re-runs) | 12 |
| **Run cap** | **60** |

  **Expected under 0.5 core-h. Cutoff 2 core-h.** It runs in the cloud, resumably. No PC time.
- **Stop rule.** If the cap or cutoff would be exceeded, the work stops and the PI decides. The order of work is Q1
  fronts, then robustness, then Q2.
- **Code.**
  - `src/run_mg9c.py` imports MG9's `mg9_ode`, splice, envelope, tracker and preparation unchanged, and adds the §4
    classifier and the runners for §§5–6.
  - `src/evaluate_mg9c.py` computes the labels.
  - `tests/test_mg9c.py` covers the §4 rules on synthetic tracks (exponential settling, constant slow speed, power-law
    creep at the threshold, early stop, asymmetry), the §6 bracket logic and the §7 outcome order.
- **Records.** Per-run JSON as in MG9, a ledger, the evaluation, and the snapshots `.npz`.
- **Pre-run specification.** Anything left open is written into GATES.md before the first run.

## 9. Decisions for the PI

1. **Horizon and thresholds:** T = 12,000 (extension to 30,000); W₃ at λ/40 and λ/40,000. These are generous to a
   pure exponential (with MG9's fitted time constant, the predicted D₃ is about 6e-4 ≈ λ/5000). They reject any
   constant creep faster than about λ/240,000 per unit time (λ/40 over 6000 units).
   Recommended: as written.
2. **Q2 tolerance:** 0.01 (recommended), or none, in which case Q2 is purely descriptive.
3. **29.22 included as descriptive:** recommended yes; it costs nearly nothing.
4. **Review:** a quick ChatGPT pass before registration (recommended), or register directly.

## 10. What each outcome would mean

- **BOUNDED.** Localized patterned states coexist stably with the uniform state just below b_r. This is what a pinning
  (homoclinic snaking) picture of the window would require. The lower window edge stays untested on this branch. The
  next study would be branch selection (the 30-cell family) to reach it; S(b) (MG9b) becomes better motivated.
- **SLOW INVASION.** MG9's low speeds were a long transient toward invasion. Localization would then have to come from
  wavelength and branch selection, not front pinning; that study comes next.
- **MIXED.** The settling region ends inside 29.25–29.35, or depends on mesh or domain. The point-by-point list says
  which.
- **AT b_r (Q2).** The settling of extended patterns stops where single particles start to replicate. That links the
  upper window edge to front depinning. **BELOW b_r** would mean the two thresholds are distinct.

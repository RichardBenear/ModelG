# MG9d — What drives MG9c's late drift? A bounded diagnostic study. Design proposal, DRAFT 1

*Stage MG031. Claude, 2026-11-05. Status: **DRAFT 1, not registered. No MG9d integration has been run.** Requested
by the PI (2026-11-05 00:30 UTC) after discussion with ChatGPT. MG9c is CLOSED (Q1 UNRESOLVED, Q2 NOT ELIGIBLE), and
nothing here changes its classifications, its stationarity threshold, or its runs. Naming: MG9b stays reserved for
the stochastic barrier S(b).*

## 0. Starting point and disclosed prior knowledge

- **Primary start state.** The committed MG9c field at b = 29.30, t = 30,000:
  `results/MG9c/snapshots.npz` key `front_D02_29.30000_t30000`, SHA-256 `f22ca2807eee…` (full hash in
  `snapshots_manifest.json`). This is D02 (h = 0.2, N = 500, L = 100), E_pat = 6.6009, edges at x_L = −29.191 and
  x_R = +29.011, 18 cells.
- **Known before this draft (MG9c and its diagnostics):**
  - The patch translates leftward at about 4.4e-7 per unit time at t ≈ 30,000.
  - The width is nearly constant.
  - The speed at 29.30 has been falling: about 8.2e-7 over 12k–15k, 7.1e-7 over 15k–20k, 5.8e-7 over 20k–25k, and
    4.8e-7 over 25k–30k. That is roughly 18% per 5000 units.
  - The final residual, 1.09e-5, is consistent with that drift.
- **Static checks on the start state (no integration).**
  - Spectral content: the power fraction in the Nyquist mode is 2.8e-11 (φ_X), 1.7e-11 (φ_Y) and 6e-14 (φ_G). In the
    top 10% of wavenumbers it is 4.2e-7 (φ_X) or less. So Fourier translation by a fraction of h is well defined
    here.
  - A Fourier shift by exactly h reproduces a one-node `np.roll` to 6e-15.
  - Reflection about the node x = 0 is the exact index map j → (N − j) mod N. No interpolation is involved.

## 1. Question and scope

Which of four explanations does the late drift favour: **unresolved numerical error**, **mesh alignment** (a lattice
effect tied to the patch's sub-grid position), **slow relaxation** left over from the preparation, or **persistent
motion**? This is a bounded diagnostic. It is not a prerequisite for MG9b, and it does not try to settle every
localization question. **INCONCLUSIVE is an allowed outcome.**

## 2. Runs (b = 29.30; all start at t' = 0 from the state named; horizon T = 10,000)

| Run | Start state | Mesh | Solver | Purpose |
|---|---|---|---|---|
| **R0 baseline** | the start state | D02 | registered: BDF, rtol 1e-8, atol 1e-10, max_step 0.5 | reference continuation |
| **R1 tight** | the start state | D02 | rtol 1e-10, atol 1e-12, max_step 0.1 | numerical accuracy |
| **R2a/b/c shifted** | start state translated by s = h/4, h/2, 3h/4 (all three fields, periodic Fourier interpolation) | D02 | registered | mesh alignment |
| **R3 reflected** | start state reflected about the node x = 0 (exact index map) | D02 | registered | symmetry check |
| **R4 finer mesh** | start state Fourier-interpolated to h = 0.1, N = 1000, L = 100 (same physical profile; no re-preparation) | D01 | registered | mesh refinement of the same patch |

- **Order of runs:** R0, R1, R3, R2a, R2b, R2c, R4. The reflection check comes before the shifts, because a failed
  symmetry check makes everything else INCONCLUSIVE.
- **Fixed horizon.** T = 10,000 for every run, with no extension. A run is never continued until it passes.
- **The solver wrapper** is MG9's `mg9_ode` with `max_step` exposed as a parameter. This is the only change. At
  max_step 0.5 it must reproduce `mg9_ode` bit for bit (tested before any run). seed_scale = 0.
- **Segments** are 1000 units long, each a separate call starting from the previous final state. This is the same in
  every run, so restarts do not differ between runs.

### Interpolation error (documented per run, before integrating)
- **R2:**
  - the spectral power fractions of the start state above;
  - the round-trip check: shift by s and back by −s, and the maximum difference (expected at round-off);
  - the one-node check (s = h against `np.roll`);
  - the **initial residual** max|RHS| of each shifted state against the unshifted state's 1.09e-5, computed from an
    independent RHS evaluation with no time stepping.
- **R4:** the same initial-residual comparison on D01, and the maximum difference between the D01 state sampled back
  at the D02 nodes and the original.
- These numbers describe how far each start departs from "the same physical state". They are reported, not used as
  criteria.

## 3. Measurements (from each run's own t' = 0 state)

At every 50 units (201 samples), using the MG9 tracker and the MG030b diagnostic code (with its corrected
matched-minima subset):
- boundary midpoint C_edge and width W from the tracked envelope half-level crossings;
- the E²-weighted and φ_X²-weighted centres;
- the mean shift of the interior minima matched at t' = 0 (one subset for every statistic);
- the **profile-alignment displacement** d(t'): the translation that best maps the run's own t' = 0 φ_X onto φ_X(t').
  It is the primary measure;
- the translation-adjusted shape change ‖φ_X(t') − shifted φ_X(0)‖ / ‖φ_X(0)‖.

All displacements are measured from the run's own start, so an imposed translation is never counted as motion. For
R3, displacements are multiplied by −1 to compare with R0. R4 uses the same physical coordinates.

**Derived quantities (primary measure d):**
- v_X: the least-squares slope of d over [2000, 10,000]. The first 2000 units are excluded so that any
  interpolation transient settles.
- ρ_X: the slope over [6000, 10,000] divided by the slope over [2000, 6000] (deceleration ratio).
- **Noise floor** v_f = 5e-9, about 1% of the expected speed and about 40 times the speed resolution implied by the
  ~1e-6 alignment error over 8000 units. The other measures are reported beside d; a disagreement among them
  larger than 10% of |d(10,000)| is flagged.

## 4. Decision rules (ordered; each comparison gets a label)

Agreement between two speeds: **AGREE** if |v_A − v_B| ≤ 0.05 |v_B| + v_f; **DIFFER** if
|v_A − v_B| > 0.20 |v_B| + v_f; otherwise **INCONCLUSIVE**.

| Comparison | Rule | Labels |
|---|---|---|
| **C0 restart** | v_R0 against MG9c's 29.30 speed over [25,000, 30,000] (4.8e-7, from the committed track) with the AGREE/DIFFER rule, at a 25% band | CONTINUOUS / RESTART-SENSITIVE |
| **C1 numerics** | v_R1 against v_R0; also \|ρ_R1 − ρ_R0\| ≤ 0.05 | NUMERICS-RESOLVED (both agree) / NUMERICS-SENSITIVE (DIFFER) / INCONCLUSIVE |
| **C3 reflection** | −v_R3 against v_R0, with a 1% band: \|v_R3 + v_R0\| ≤ 0.01 \|v_R0\| + v_f | SYMMETRIC / ASYMMETRY-FLAG |
| **C2 alignment** | v for s = 0 (R0), h/4, h/2, 3h/4 | ALIGNMENT-INDEPENDENT (all AGREE with R0) / ALIGNMENT-DEPENDENT (any DIFFER, including a sign change) / INCONCLUSIVE |
| **C4 relaxation** | ρ_R0 | DECREASING (ρ ≤ 0.90) / PERSISTENT (0.95 ≤ ρ ≤ 1.05) / INCREASING (ρ ≥ 1.10) / INCONCLUSIVE (otherwise) |
| **C5 finer mesh** | Identity first: the tracked patch is present at every sample, with no fill, loss, breakup or remote component; the cell count stays within 18 ± 1; \|W(10,000) − W(0)\| < λ/2. Then v_R4 against v_R0 | LOST-IDENTITY / SIMILAR (\|v_R4 − v_R0\| ≤ 0.20 \|v_R0\| + v_f) / REDUCED (\|v_R4\| < 0.5 \|v_R0\|, or the opposite sign) / OTHER |

**Overall reading (ordered; the first that applies):**
1. **INCONCLUSIVE:** C3 ASYMMETRY-FLAG (a code or solver asymmetry), C0 RESTART-SENSITIVE, or C1 INCONCLUSIVE.
2. **NUMERICAL ACCURACY UNRESOLVED:** C1 NUMERICS-SENSITIVE.
3. **MESH EFFECT SUPPORTED:** C2 ALIGNMENT-DEPENDENT, or C5 REDUCED. The C4 label is reported alongside, since a mesh
   effect and relaxation can coexist.
4. **RELAXATION SUPPORTED:** C2 ALIGNMENT-INDEPENDENT, C5 SIMILAR, C4 DECREASING.
5. **MOVING STATE SUPPORTED:** C2 ALIGNMENT-INDEPENDENT, C5 SIMILAR, C4 PERSISTENT.
6. **INCONCLUSIVE:** anything else (for example, C5 LOST-IDENTITY or OTHER, or a C2 or C4 INCONCLUSIVE). The
   component labels are listed.

**Wording is conditional, as the PI asked:**
- Alignment dependence *supports* a mesh effect.
- Solver sensitivity means numerical accuracy is unresolved.
- A decreasing drift *supports* relaxation, but does not prove eventual arrest.
- Persistent motion that survives the numerical checks *supports* a moving-state reading, but does not prove an
  asymptotic travelling solution.
- C3 is a symmetry check. The equations and lattice are mirror-symmetric about a node, so R3 should reproduce R0
  mirrored to round-off. A pass cannot distinguish physical drift from mesh relaxation; a failure points to the code.

## 5. Budget and saved fields

- **Cost estimate** (MG9c timings: about 30 s per 10,000 units on D02 with the registered solver, including
  tracking; about 60 s on N = 1000 domains):
  - R0, R2a–c, R3: about 3 min in total;
  - R1: assumed 10–20× slower, 5–10 min;
  - R4 (D01, stiffer): assumed 3–5 min;
  - in-run diagnostics every 50 units add under 1 min per run.

  **Expected about 0.3 core-h. Strict cap: 10 integrations (7 planned + 3 numerical re-runs). Cutoff 1 core-h.**
  Cloud. If the cap or cutoff would be exceeded, work stops and the PI decides. Nothing not reached is replaced.
- **Saved and committed:**
  - each run's t' = 0 state (with SHA-256);
  - full three-field snapshots every 1000 units (11 per run: about 0.8 MB float64, compressed);
  - per-run JSON with the tracker every 5 units, the diagnostics every 50 units, solver success, time reached, and the
    final residual;
  - the interpolation-error numbers of §2, the ledger and the evaluation.

  That is under 4 MB in total.
- **Before any run:** code (`src/run_mg9d.py`, `src/evaluate_mg9d.py`) and tests (`tests/test_mg9d.py`). The tests
  cover:
  - the wrapper's bit-for-bit match at max_step 0.5 (on MG003's b = 29 reference, as MG9's tests do);
  - shift and reflection maps on synthetic fields;
  - the decision rules on synthetic speed sets;
  - the overall reading order.

  A short pre-run specification goes in GATES.md, then a ChatGPT review, then PI approval.

## 6. What each outcome would mean for MG9b (stochastic barrier S(b))

MG9b would measure the rate at which noise nucleates a single particle from the uniform state at a few b values
(MG3b's method on the h = 0.2 lattice). The drift studied here concerns extended 18-cell patches that already exist.
So **no outcome changes MG9b's design materially.** The effects are on interpretation and wording:

- **NUMERICAL ACCURACY UNRESOLVED:** little effect. MG9b's stochastic tau-leap runs do not use this BDF continuation,
  and its event confirmation uses short deterministic runs whose outcome (particle or not) is far from tolerance
  level. It would argue for one solver-tolerance spot check in MG9b's confirmation step.
- **MESH EFFECT SUPPORTED:** moderate effect on interpretation only. It would mean localized states on h = 0.2 feel the
  lattice. MG9b, also on h = 0.2, would then be measuring a lattice-model barrier, and the particle-window edges
  (b_s, b_r) it relates to may themselves be partly lattice-set. Since MG3b and MG4c already treat everything as
  lattice-conditional, the S(b) procedure is unchanged, but MG9b should say this.
- **RELAXATION SUPPORTED:** no effect. A slow leftover adjustment of a prepared patch says nothing about nucleation from
  the uniform state.
- **MOVING STATE SUPPORTED:** small effect. It would add a caveat that "stationary localized state" is not guaranteed
  for extended patches near b_r. MG9b's end state is a single particle, whose stationarity MG3/MG3b confirmed
  separately, so S(b) is unaffected.
- **INCONCLUSIVE:** no effect. MG9b can proceed on its own terms.

## 7. Decisions for the PI

1. **Horizon T = 10,000 and the 50-unit diagnostic sampling.** Recommended: as written. It is long enough to resolve
   the expected ~0.004 displacement and a deceleration ratio near 0.85, if MG9c's trend continues.
2. **Agreement bands** (5% AGREE / 20% DIFFER; 1% for reflection), **noise floor** v_f = 5e-9, and the **C4
   thresholds** 0.90 / 0.95–1.05 / 1.10. Recommended: as written.
3. **The cap of 10 integrations and the 1 core-h cutoff.**
4. **A ChatGPT review before registration.** Recommended.

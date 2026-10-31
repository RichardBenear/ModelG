# MG4c — Is there a single-particle window on unweighted Delaunay graphs? (proposal, draft 2)

*Stage MG011a. Claude, 2026-10-31. Status: **DRAFT 2 for ChatGPT review; not registered; no MG4c run has been made.***
*Draft 1 is preserved as `spec/MG4c_proposal_v1.md` (superseded). Review: `review/MG011_chatgpt.md`. Response:
`review/MG011_response.md`. MG4c follows MG4 (CLOSED, FAIL) and does not re-score it. The name MG4b stays reserved for
the Wolfram hypergraph snapshots.*

## 1. Question
At b = 29, every unweighted Delaunay (D) and RGG (R) run in MG4 *multiplied*. The triangular comparator T stayed single,
and on T b = 29 lies 0.2 below the replication edge (T's window is 27.33 < b < 29.20, from MG4 Part 0). MG4c asks two
narrow questions, on the D panel graphs D0 and D1 only:

1. **Window:** at lower b, does a D graph keep a single particle, as MG4's detector sees it, from MG4's initial
   condition?
2. **Operator:** at b = 29, does replacing the unweighted Laplacian with a geometry-weighted one change D0 p0's outcome?

MG4c is an **estimation study with pre-registered outcome categories** (like MG3b). It has no PASS/FAIL verdict.

## 2. Fixed from MG4 (unchanged)
- Rules v1.1 and the eq. 17 parameters except b. The 20 × 20 commensurate torus (122 × 142 at h_t = 0.2).
- Panel graphs D0 and D1, rebuilt from MG4's Part 0 records, with their calibrated h. Placements p0 and p1 of each, as
  in MG4 Part A.
- The initial condition at every b is the continuum soliton (MG1, b = 29).
- The solver: BDF with rtol 1e-8, atol 1e-10, max_step 0.5 and nested-dissection ordering, from t = 0 to 300.
- Checks at t = 200, 210, …, 300.
- MG4's smoothing, particle detection and outcome classes (multiplied > decayed > single > unresolved).
- MG4's node-sum balance check, with the change in §7.

**"The four D runs" at a value of b** means D0 p0, D0 p1, D1 p0 and D1 p1. **"All-single at b"** means all four are
classed *single*. An *unresolved* run counts as not single.

## 3. Detection threshold (fixed bootstrap, no circularity)
- At every b that is tested, one T run is made, from the centre placement (MG4's T p0).
- **θ_b = ½ × the global minimum of T's smoothed φ_X at t = 300.** That minimum is a number read off the field, with no
  classification. T is then classified with θ_b.
- **If T is single at b,** the four D runs at b are classified with θ_b.
- **If T is not single at b** (expected below about 27.33), the D runs at b are classified with the **fixed reference
  θ_29 = −5.117354**, MG4's registered threshold. No "nearest" rule is used. Classes, once assigned, are never revised
  by later runs.
- For every D run, the global minimum of its smoothed φ_X at t = 300 and its class under θ_29 are also reported. A
  shallow structure that the threshold misses therefore stays visible.
- Any negative finding is conditional on this detector.

## 4. Design
**Part 1 — coarse scan.** b ∈ {28.75, 28.5, 28.25, 28.0, 27.75, 27.5, 27.25, 27.0}. At each b, T × 1 and D × 4: 40 runs.
MG4's Part A runs give D0/D1 × p0/p1 at **b = 29**, where all four multiplied. They serve as the b = 29 grid point and
are not rerun.

**Part 2 — refinement.** Each step tests one new b with T × 1 and the four D runs (5 runs). An edge gets at most 6
steps. Each step halves a bracket whose two ends differ in the predicate being tracked.
- **Case A: some grid b (Part 1, or b = 29) is all-single.**
  - Take the contiguous block of all-single grid points that contains the highest all-single b. Other all-single
    blocks are reported but not refined.
  - **Upper edge:** bisect between the block's top b and the next grid b above it, on the predicate "all-single".
  - **Lower edge:** bisect between the block's bottom b and the next grid b below it, on the same predicate.
  - An edge with no grid point beyond it is reported as unbracketed: "upper edge > the top b" or "lower edge < 27.0".
    No steps are spent on it, and the scan is not extended.
- **Case B: no grid b is all-single.** The narrow-window search follows the four-run count of single runs.
  - Take the highest grid b at which **no D run is multiplied**, and the grid b just above it. Bisect between them on
    the predicate "at least one of the four D runs is multiplied".
  - The search ends early if a step is all-single. That b then counts as an all-single point (§5).
  - If every grid b down to 27.0 has a multiplied run, there is no bracket. Case B then runs no refinement and reports
    "multiplied at every tested b". The scan is not extended below 27.0.
- **Resolution:** 6 steps reduce a 0.25 bracket to 0.004. This locates an edge inside a valid bracket. **It does not
  guarantee finding a window narrower than the grid**, because Case B finds a single-particle region only if a
  bisection point lands in it.
- A disconnected single-particle region, meaning non-adjacent all-single points, is reported with its points. Its
  edges are not refined.

**Part 3 — operator diagnostic (a diagnostic amendment to rules v1.1; not a change of MG4's rules).**
- This part uses the **cotangent (P1 finite-element, lumped-mass) Laplacian** on the same periodic Delaunay triangulation:
  (L_w f)_i = (1/A_i) Σ_j w_ij (f_j − f_i), with w_ij = ½(cot α_ij + cot β_ij), where α and β are the angles opposite
  edge ij. A_i is node i's barycentric area, one third of the areas of its triangles.
- It replaces only the diffusion operator. The reactions and all the settings are as in §2. No h calibration is used.
- Two runs at b = 29: D0 p0 with L_w, and T p0 with L_w as the control.
- **Construction checks (before any Part 3 run; any failure stops Part 3):**
  - **Periodic faces:** triangles come from the 3 × 3 tiling. A face is kept once, if and only if its centroid lies in
    the central tile. The sum of the kept faces' areas equals Lx·Ly to relative 1e-12, and every edge of the D graph lies
    in exactly two kept faces.
  - **Weights and areas:** every w_ij ≥ −1e-12 (Delaunay), every A_i > 0, and Σ A_i = Lx·Ly.
  - **Weighted symmetry:** K = A·L_w is symmetric to 1e-12, L_w annihilates constants, and Σ_i A_i (L_w f)_i = 0 for random
    f to 1e-12·Σ|terms|.
  - **On T:** L_w equals MG4's L_G (h² = 1.5a²) entrywise to 1e-12.
  - The Jacobian with L_w matches finite differences.
- **Control:** T with L_w must reproduce MG4's T p0 particle: the same class (single), and a final φ_X field equal to
  MG4's to 1e-6 (the operators are equal on T). If not, Part 3 has status **O-invalid**.
- **Conservation** for L_w uses A-weighted node sums, I_S = Σ_i A_i φ_S,i, with the same identities and tolerances as MG4.

## 5. Outcomes (pre-registered statuses; no PASS/FAIL)
**Window (Parts 1–2):**
- **W1 — single-particle points found:** at least one b, on the grid or in refinement, is all-single. Report:
  - every all-single b;
  - the refined edges with their brackets, or "unbracketed";
  - depth and λ against T at the same b.
- **W2 — none detected:** no tested b is all-single, and no D run at any tested b is *single*.
- **W3 — partial:** some D runs are single, but no tested b is all-single. Report the counts by b and run.
- **INCOMPLETE:** Part 1 did not finish (run cap or ceiling). Report whatever data were completed.
- **NOT SCORED:** an implementation check failed (§7) and is not resolved by the recorded precision recheck.

**Operator (Part 3):**
- **O1:** D0 p0 with L_w is *single*.
- **O2:** it is *multiplied* or *decayed*.
- **O3:** it is *unresolved*. This does not support either reading.
- **O-invalid:** the T control or a construction check failed.

**How they read (stated before any run; limited to what was tested):**
- **W1:** at the b values found, D graphs keep a single particle under MG4's detector and initial condition, so b = 29 is
  outside D's single-particle range. MG5 design may consider those b values.
- **W2:** no single-particle outcome was detected at the tested b values from this initial condition. This does *not*
  show that no window exists.
- **O1:** for D0 p0 at b = 29, changing the operator changes the outcome. That shows the unweighted operator matters
  there. It does *not* show that every MG5 mapping needs weights.
- **O2:** the weighted operator did not change the outcome there.
- **Programme policy, not a scientific conclusion:** if the result is W2 with O2, the PI and Claude stop and reconsider
  before designing MG5.

## 6. Budget
| Item | Runs (cap) |
|---|---|
| Part 1: 8 × (T 1 + D 4) | 40 |
| Part 2: 2 edges × 6 steps × (T 1 + D 4) | 60 |
| Part 3: D weighted + T weighted | 2 |
| **Total** | **102** |

- **Cost:** D runs took about 13 min each on the PI's PC in MG4, and T runs about 4 min. Projected about 10 core-h,
  with an upper bound near 18 if every refinement step is used.
- **Ceiling:** 20 core-h, a **launch cutoff**. No job starts once the ledger reaches 20. Jobs already running finish and
  are reported. Since at most 9 jobs (one batch) run at once, the hard maximum is about 20 + 9 × 0.25 ≈ 22.5 core-h.
  This hard maximum is registered as part of the budget, not treated as an overrun.
- **Run cap:** 102. A run-cap or ceiling stop before Part 1 finishes gives **INCOMPLETE**. A stop during Part 2 or Part 3
  is reported with the brackets reached so far.

## 7. Implementation (fixed in GATES.md before any run)
- **Driver:** a new `src/run_mg4c.py`. It is resumable, enforces the run cap and the ceiling, logs in UTF-8, and runs on
  the PI's PC with `MODELG_WORKERS=9`. Each Part 2 step waits for its 5 runs before the next b is chosen. The rules in §4
  are coded as pure functions of the recorded classes and unit-tested on constructed class tables, including both
  cases, the unbracketed edges, disconnected regions and unresolved runs.
- **Operator:** `graph_ode.make_system` gets an optional operator argument (L_w and the area vector A), with the §4
  construction checks as unit tests.
- **Balance check and its recheck:**
  - The registered float64 check is unchanged.
  - A flag is resolved only by a recheck of that run's own b and operator, in genuine higher precision. That means
    `mpmath` at 40 digits, which is installed on the PI's machine; otherwise exact `fractions.Fraction` arithmetic.
    Not `np.longdouble`, which on the PI's Windows build has float64 precision.
  - The recheck takes b, the operator and the run directory from the run record, not from fixed MG4 paths.
  - Resolutions are recorded in a file before any outcome is read, as in MG4.
- **Random streams:** none are new. The graphs and placements are MG4's.

## 8. Not in MG4c
R graphs, weighted panels, noisy runs, Wolfram hypergraph snapshots (MG4b), and any change of the initial condition
with b.

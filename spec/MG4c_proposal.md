# MG4c — Is there a single-particle window on unweighted Delaunay graphs? (proposal, draft 1)

*Stage MG011. Claude, 2026-10-31. Status: **DRAFT for ChatGPT review; not registered; no MG4c run has been made.***
*Follows MG4 (CLOSED, FAIL). It does not re-score MG4. The name MG4b stays reserved for the Wolfram hypergraph snapshots
(MG4 proposal §11).*

## 1. Question and why it decides what comes next
At b = 29, every unweighted Delaunay (D) and RGG (R) run in MG4 *multiplied*: the particle copied itself until 33–37
particles filled the torus. The triangular comparator T stayed single, and on T b = 29 lies only 0.2 below the
replication edge (T's window is 27.33 < b < 29.20, MG4 Part 0). Two questions follow, and MG4 answered neither:

1. **Window:** at lower b, is there a range in which an unweighted D graph keeps a single particle?
2. **Operator:** is MG4's outcome caused by the unweighted graph Laplacian? Its local action on a smooth profile follows
   the raw edge moments, not the calibrated h (MG4's C3(b) finding). A geometry-weighted Laplacian on the same graph
   would show this.

Why it matters: a Wolfram-style network (MG5) is at least as irregular as D. If D has a window, MG5 can be designed at a b
inside it. If D has no window but the weighted operator keeps a single particle, the mapping needs rules that carry edge
weights or local geometry. If neither works, the programme stops and rethinks before MG5.

This is an **estimation study with pre-registered outcome categories** (like MG3b). It has no PASS/FAIL verdict.

## 2. Fixed from MG4 (unchanged)
- Rules v1.1, eq. 17 parameters except b. The 20 × 20 commensurate torus (122 × 142 at h_t = 0.2).
- Panel graphs D0 and D1, rebuilt from MG4's Part 0 records, with their calibrated h.
- Placements p0 and p1 of each graph, as in MG4 Part A.
- The continuum soliton initial condition (MG1, b = 29), the solver (BDF, rtol 1e-8, atol 1e-10, max_step 0.5,
  nested-dissection ordering), t = 0–300, and checks at t = 200, 210, …, 300.
- The smoothing, particle detection and outcome classes (multiplied > decayed > single > unresolved) of MG4 §5.1–5.2.
  The node-sum balance check also stays, with one addition: a flag is resolved by the extended-precision recheck
  (`src/check_mg4_balances.py`), recorded before the outcome is read (MG4 lesson).

The initial condition stays the b = 29 soliton at every b. This tests whether a particle-like start settles to a single
particle at that b. It is the same choice as MG4 Part 0's T bisection.

## 3. Threshold at each b
MG4 set θ = ½ × T's mean smoothed depth at b = 29. The particle's depth changes with b, so MG4c runs **T at every b
tested** (one run from the centre placement) and uses θ_b = ½ × T's smoothed depth at that b.
- If T is not *single* at some b, that b has no threshold. Its D runs are still made and classified with the nearest
  valid θ, and they are reported as such.
- T's own class at each b also re-checks T's window (expected single for 27.33 < b < 29.20).

## 4. Design
**Part 1 — coarse scan.** b ∈ {28.75, 28.5, 28.25, 28.0, 27.75, 27.5, 27.25, 27.0}, 8 values with a step of 0.25. At each b:
T × 1 and D × 4 (D0 and D1 × p0 and p1). That is 40 runs.

**Part 2 — refinement** (on D0 p0 only, by bisection, 6 steps per edge, resolution 0.004):
- If some grid b has D 4/4 single: bisect the upper edge (between the highest all-single b and the next b above it) and
  the lower edge (between the lowest all-single b and the next b below it).
- If no grid b has any single D run, and the scan goes from *multiplied* straight to *decayed*: bisect between the
  lowest multiplied b and the highest decayed b, to look for a window narrower than 0.25.
- Each bisection step that lands on a single run is confirmed on the other three D runs at that b (4 runs per step at
  most).

**Part 3 — operator diagnostic** (one run plus its comparator): D0 p0 at b = 29 with the **cotangent (P1 finite-element,
lumped-mass) Laplacian** on the same Delaunay triangulation:
- (L_w f)_i = (1/A_i) Σ_j ½(cot α_ij + cot β_ij)(f_j − f_i), where α and β are the angles opposite edge ij and A_i is the
  node's barycentric area (one third of its triangles' areas). It uses the embedding and needs no h calibration.
- On T it reduces exactly to MG4's operator: weight 1/√3 and A = (√3/2)a², so (2/3a²) Σ (f_j − f_i), i.e. h² = 1.5a².
- One T run with L_w at b = 29 checks that the implementation reproduces MG4's T particle (same class, depth within 1%).
- Conservation: the weighted sums Σ_i A_i φ_S,i obey MG4's balance identities. The check uses A-weighted node sums.

## 5. Outcomes (pre-registered categories; no PASS/FAIL)
- **W1 — window exists:** some grid b has D at 4/4 single. Report the window edges from Part 2 next to T's (27.33, 29.20),
  plus depth and λ against T at the same b.
- **W2 — no window at 0.25 resolution, none found in refinement:** no grid b and no refined b has any single D run.
- **W3 — partial:** single runs occur, but at no b on all four D runs. Report the counts.
- **Operator:** the weighted D run is *single* (O1) or *multiplied/decayed/unresolved* (O2). It is reported with T's
  weighted control.

How they read (stated before any run):
- W1: the particle survives on unweighted irregular graphs with a shifted b.
- W2 with O1: the unweighted operator is the problem; MG5 would need weights or geometry.
- W2 with O2: a deeper problem; stop and rethink before MG5.
- W3: graph-to-graph or placement sensitivity; reported as such. No further runs without a new registration.

## 6. Budget (run caps bind; core-hour ceiling is a stop trigger, checked before each job starts)
| Item | Runs (cap) | Projected core-h (PI's PC) |
|---|---|---|
| Part 1: 8 × (T 1 + D 4) | 40 | about 7 |
| Part 2: up to 2 edges × 6 steps × up to 4 runs | 48 | about 2–8 |
| Part 3: D weighted + T weighted | 2 | about 0.3 |
| **Total** | **90** | **about 10–15; ceiling 20** |

Lesson from MG4: runs already in flight when the ledger passes the ceiling are allowed to finish and are reported. Their
cost is part of the registered ceiling, which therefore includes a margin of one batch.

## 7. Implementation (to be fixed in GATES.md before any run)
- New driver `src/run_mg4c.py`: resumable, enforces the run cap and the ceiling, logs as UTF-8, and runs on the PI's PC
  with `MODELG_WORKERS=9`.
- `graph_ode.make_system` gets an optional operator argument (the weighted L_w), with unit tests: L_w on T equals MG4's
  L; L_w annihilates constants; its A-weighted sum of L_w f is 0; the Jacobian matches finite differences.
- Random streams: none are new. The graphs and placements are MG4's.

## 8. Not in MG4c
R graphs (about 1.5 h per run), weighted panels, noisy runs, Wolfram hypergraph snapshots (MG4b), and any change of the
initial condition with b.

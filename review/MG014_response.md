# Response to review MG014 (MG5 draft 1) — Claude, 2026-10-31

All five items are accepted. Draft 2 is `spec/MG5_proposal.md`; draft 1 is preserved as `spec/MG5_proposal_v1.md`.
The PI's three decisions (interpolating transfer, Part B, unweighted operator) are unchanged.

**1. Event rule vs cost — Accepted with modification.** The reviewer is right: draft 1's budget allowed one rebuild per
epoch, while its rule had each birth interpolate on the mesh after earlier events. Rather than build local Delaunay
updates (unmeasured, and a large new component), draft 2 makes the tick-synchronous rule the registered rule (§3.2): all
births of a tick interpolate on the tick-start mesh and field, deletions are applied, then one rebuild. Δ becomes part of
the rule, the Δ → 0 limit is the continuous-time rule, and a tick check (Δ = 0.2 against 0.1 on the same event list, §6)
measures the dependence, with a pre-set fallback (halve Δ, repeat once, then stop for the PI). The fastest rate's tick is
now 0.2 (νΔ = 0.04), and the budget counts one rebuild and one LU per tick: about 20 core-h in all (§10).

**2. Frozen baseline and Part B initial condition — Accepted.** I confirmed both figures by evaluating the stored states:
max|dφ/dt| is 6.3e-13 to 8.5e-13 under the panel h and 0.049, 0.049, 0.018, 0.018 under h = 0.2. Draft 2 keeps each panel's
calibrated h for the whole run (§3.5) **and** adds four matched frozen controls in Part A with MG5's solver, h and horizon
(§5), with an INVALID-CONTROL rule (§8). Part B now keeps absolute concentrations: φ_Y is shifted by −0.0642857 (= −0.9/14;
G₀ and X₀ do not depend on b), as the reviewer computed (§5).

**3. ROS2 validation — Accepted.** Draft 2 separates a step check (dt 0.05 against 0.025, same events and ticks) from the
tick check (Δ 0.2 against 0.1, same event list), both exact replays of one saved event list (§6). It defines the
comparison when particles multiply or dissolve (same class, first qualifying check at most one check apart). Alignment:
every Δ (1, 0.2, 0.1) is an integer multiple of every dt (0.05, 0.025) and divides the check interval, so no partial
steps occur; there is one LU per tick (§4).

**4. A1 wording — Accepted.** A1 is now "recovery under rewriting plus transfer". "An averaged, more regular medium" and
"smoothing by repeated interpolation" are both listed as possible explanations that MG5 does not separate (§3.4, §8). No
weighted arm is added.

**5. Auditable data and statuses — Accepted.** At every check the driver saves the positions, ids and smoothed φ_X of all
nodes below −1.0. Detection at any θ ≤ −1.0 depends only on those nodes, so every class can be recomputed exactly from the
saved records. Full fields, positions and ids are saved every 50 units. Per-tick jump accounting (I_after = I_before − R +
B per species) and an independent-interpolator audit are added. COMPLETE / FAILED-NUMERICAL / INCOMPLETE statuses are
applied before any class, with NOT SCORED and INCOMPLETE rules at the rate and gate level (§7, §8).

**Reviewer answers to draft 1's questions** are recorded in §11 (1.5 radius acceptable as declared; four runs per rate
adequate for the fixed panel, without a general probability claim; solver acceptable conditionally).

One new PI question: whether to commit all full-field snapshots (about 0.5 GB) or only the check records and accounting,
keeping full fields on the PC with checksums (§7).

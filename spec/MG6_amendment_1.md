# MG6 amendment 1 (MG021): tick Δ = 0.1 after the Part 0 STOP — draft 2

Status: **draft 2, 2026-11-02, after ChatGPT review MG021 (`review/MG021_chatgpt.md`; response
`review/MG021_response.md`). PI decision taken; nothing has run under this amendment.** Draft 1 is kept as
`spec/MG6_amendment_1_v1.md` (superseded). This amends the registered MG6 design (`spec/MG6_proposal.md`, draft 3,
frozen at `106b3e2`) and its pre-run implementation specification (GATES.md, MG6, items 1–13). Everything not named
here is unchanged.

## 1. What happened (the registered result stands)
MG6 Part 0 ran on the PI's PC on 2026-11-02, 08:43–11:58 local (Windows 11, Python 3.13.2, numpy 2.2.4, code `410696d`,
8 workers). Records: `results/MG6/part0/decisions.json`, `check_*.npz`, `accepted_*.npz`, `results/MG6/part0_stdout.txt`
(commit `5d351dc`). **Registered outcome: STOP (step or tick check failed twice). It stays a STOP and is not
re-scored.** ChatGPT independently reproduced these records (MG021 §1).

| Check | Result |
|---|---|
| 1. Acceptance and determinism | ok (z = −0.45; identical bitmaps) |
| 2. Density equilibration | ok (all \|z\| ≤ 1.87) |
| 3. Step check, dt 0.05 vs 0.025 (fixed edits) | ok at Δ 0.2 (centre 6.0e-4, depth 0.09%) and Δ 0.1 (4.4e-4, 0.09%) |
| 4. Tick check, Δ 0.2 vs 0.1 | **FAIL, part (a)**: centre 0.27, 0.40, 0.54, 0.17 at t = 40, 50, 60, 90 (limit 0.15); depth up to 3.76% (limit 3%). Part (b) ok: \|ΔD̄₀\| = 1.33 nodes, limit 37.2 |
| 4. Tick check, Δ 0.1 vs 0.05 (registered retry) | **FAIL, part (a), centre clause only**: centre 0.36 at t = 60 (limit 0.15). Counts equal at every check; depth max 2.79% (limit 3%). Part (b) ok: \|ΔD̄₀\| = 2.33 nodes, limit 36.7 |
| 5–8 | not reached |

Record note: the first attempt's Δ 0.1 comparison run and the retry's Δ 0.1 main run have the same configuration and seed,
and share a file name (`check_d0.1_dt0.05.npz`); the second overwrote the first. Byte identity cannot be shown, because the
first file's hash was not recorded. Both give D̄₀ = −148.10758250291826 to every digit, and ChatGPT reproduced the stored
comparison statistics of both attempts from the available files. Future Part 0 code must use immutable attempt IDs or
check hashes before replacing an output.

## 2. Diagnosis (unregistered, after the data were seen)
- **The centre estimate is the MG2–MG5 centroid:** the deepest smoothed node, then a (θ − s)-weighted centroid of the
  nodes within 0.6 of it. Node spacing is about 0.2, and at ν = 0.2 nodes live about 5 time units. Nodes crossing the hard
  0.6 disc and the threshold make the estimate jump.
- **The estimate is stable under the tested time-step refinement on a fixed edit history:** the step check imposes the
  same edits and the centres agree to 6e-4. (This is agreement between two numerical solutions, not a statement about the
  estimator's absolute accuracy.)
- **The paths are common:** in all three Δ runs the particle moves about 1.3–1.4 along the same path. At most checks the
  runs agree to 0.01–0.03, and they rejoin at t = 70 after the t = 40–60 separations.
- **Second estimator** (`review/MG020_centre_diagnostic.py`, repeated independently by ChatGPT): a quadratic fit to the
  smoothed field. On the same field it differs from the centroid by 0.02–0.21. With the fit-based centre:
  - Δ 0.1 vs 0.05: every check ≤ 0.116 (radius-1 fits 0.110, radius-1.5 fits 0.116).
  - Δ 0.2 vs 0.1: t = 40–60 stay at 0.26–0.38.
- **Density robustness (ChatGPT):** with either fitted centre, the two Δ 0.1/0.05 window means of D₀(1.5) differ by 1 node.
  Alternative centres move single snapshots by up to 7–8 nodes but the window means by at most 1.5 nodes.
- **Reading:** estimator sensitivity is a **plausible contributor** to the Δ 0.1 centre failure; it is not a measured
  noise variance, and it is not shown to account for the whole discrepancy (different ticks also change accepted edits and
  the field history). At Δ 0.2 part of the separation is a real difference in where the particle went. The single-snapshot
  centre clause is recorded as a weakness of the tick-check design.

## 3. Amendment (PI decision, 2026-11-02 21:18 UTC; narrowed in draft 2 as ChatGPT recommended)
1. **Tick Δ = 0.1, dt = 0.05** for all MG6 runs (proposal §7 check 4's own fallback value; the step check passed at
   Δ 0.1). Δ 0.2 is not used.
2. **Numerical validation accepted by exception, with a narrow waiver.** Only the **instantaneous-centre clause** of tick
   part (a) for Δ 0.1 vs 0.05 is waived. Still required, and met by the records: equal particle counts at every check,
   depth ≤ 3%, and part (b). The code checks this condition (`waiver()`) and stops if it does not hold.
   This is **not** a passed original tick check, and it does not demonstrate pathwise convergence. MG6's motion
   measures (tracks, D_p) remain descriptive, and this amendment does not establish their tick-size accuracy.
3. **Replacement for proposal §7.6 (transient):** the production transient uses the **Δ = 0.1 run's** D₀(1.5) series
   (`check_d0.1_dt0.05.npz`), with P̄₀ from the **Δ = 0.05 run** (the finer reference of the amended part (b) comparison,
   −146.89) for the 0.25·|P̄₀| tolerance. Blocks and the window rule are otherwise as registered.
4. **Checks 5–8 are completed** by `python src/run_mg6.py part0-mg021` on the PI's PC, with the registered criteria and
   these guards:
   - **Provenance:** the run hashes its inputs (`decisions.json` and the five check records) into the MG021 record. The
     final state `final_d0.1_dt0.05.npz` must have a row in `local_manifest.csv` and a matching SHA-256; a missing row or
     file is a STOP, not a pass.
   - **5. Bound margin:** max nodal φ_G ≤ 0.7·φ_lim(4) = 0.07 over the five committed step/tick runs (per-tick maxima),
     the terminal boundary of the two runs whose final states are kept (hash-verified), and every timing run (all
     boundaries, terminal included). The committed per-tick tables of the five runs stop before the terminal boundary;
     the driver now records the maximum over every boundary (`max_phiG`) for all later runs.
   - **6. Transient:** as in item 3.
   - **7. Quadrature** (i)–(iv) on the verified final state.
   - **8. Timing:** the settle is repeated (deterministic) and must be COMPLETE at t = 300; then 8 concurrent coupled
     runs of 20 units at Δ 0.1, seeds `seed_for(6, 9, 3, 1, 10 + i)`, i = 0–7, **on the grid selected by check 7**.
     Every timing run must be COMPLETE at t = 20.
   - **Authorisation:** `proceed = true` only if all of the above hold and the projected cost is ≤ 80 core-h. A projection
     above 80 writes `proceed = false` with "PI decision required"; the driver then refuses to launch until the PI's
     further decision is recorded.
   - Output `results/MG6/part0/decisions_mg021.json`, which links the original STOP (kept unchanged in
     `decisions.json`), the waiver, the input hashes, checks 5–8, the production configuration and the budget. The driver
     reads it when present.
5. **Budget:** provisional launch cutoff 80 core-h (was 40); run cap 24. Restated in GATES.md from check 8 before any
   registered run.
6. **Not changed:** the gate, its contrasts and thresholds, the arms, seeds, window rule, measurements, estimators and
   the evaluator. A tick-check redesign (for example, estimator calibration on a prescribed continuous field with a known
   centre sampled on varying node sets) is deferred to the next study that needs one.

## 4. Disclosure
Checks 5 and 6 depend only on committed files, and Claude computed them in the cloud while writing the code (bound ok on
the records, no drift: window t = 100–200). ChatGPT reproduced them (B₁ = −147.44, B₂ = −136.44, drift 11 nodes against a
36.7-node limit; largest recorded per-tick φ_G 0.0534). The PC run is the record.

## 5. Review
ChatGPT review MG021 supported the direction as a disclosed exception and asked for four implementation corrections
(provenance, complete settle and timing runs with the margin including timing and terminal boundaries, a budget stop
above 80 core-h, and timing on the selected grid) and narrower wording. All are adopted here; see
`review/MG021_response.md`.

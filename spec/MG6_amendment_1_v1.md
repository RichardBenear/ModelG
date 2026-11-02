# MG6 amendment 1 (MG021): tick Δ = 0.1 after the Part 0 STOP

Status: **drafted 2026-11-02 for ChatGPT review. PI decision taken; nothing has run under this amendment.** It
amends the registered MG6 design (`spec/MG6_proposal.md`, draft 3, frozen at `106b3e2`) and its pre-run implementation
specification (GATES.md, MG6, items 1–13). Everything not named here is unchanged.

## 1. What happened (the registered result stands)
MG6 Part 0 ran on the PI's PC on 2026-11-02, 08:43–11:58 local (Windows 11, Python 3.13.2, numpy 2.2.4, code `410696d`,
8 workers). Records: `results/MG6/part0/decisions.json`, `check_*.npz`, `accepted_*.npz`, `results/MG6/part0_stdout.txt`
(commit `5d351dc`). **Registered outcome: STOP (step or tick check failed twice). It is not re-scored.**

| Check | Result |
|---|---|
| 1. Acceptance and determinism | ok (z = −0.45; identical bitmaps) |
| 2. Density equilibration | ok (all \|z\| ≤ 1.87) |
| 3. Step check, dt 0.05 vs 0.025 (fixed edits) | ok at Δ 0.2 (centre 6.0e-4, depth 0.09%) and Δ 0.1 (4.4e-4, 0.09%) |
| 4. Tick check, Δ 0.2 vs 0.1 | **FAIL, part (a)**: centre 0.27, 0.40, 0.54, 0.17 at t = 40, 50, 60, 90 (limit 0.15); depth up to 3.8% (limit 3%). Part (b) ok: D̄₀ −149.44 vs −148.11, \|P̄₀\| = 148.73 (limit 37.2) |
| 4. Tick check, Δ 0.1 vs 0.05 (registered retry) | **FAIL, part (a)**: centre 0.36 at t = 60 only. Part (b) ok: D̄₀ −148.11 vs −145.77, \|P̄₀\| = 146.89 (limit 36.7) |
| 5–8 | not reached |

Record note: the first attempt's Δ 0.1 comparison run and the retry's Δ 0.1 main run have the same configuration and seed,
so they share a file name (`check_d0.1_dt0.05.npz`) and the second overwrote the first. Both give D̄₀ =
−148.10758250291826 to every digit, so nothing was lost.

## 2. Diagnosis (unregistered, after the data were seen)
- **The centre estimate is the MG2–MG5 centroid:** the deepest smoothed node, then a (θ − s)-weighted centroid of the nodes
  within 0.6 of it. Node spacing is about 0.2, and at ν = 0.2 the nodes turn over every ~5 time units. Nodes crossing the
  hard 0.6 disc and the threshold make the estimate jump.
- **With the node set fixed, the estimate is exact:** the step check imposes the same edits and the centres agree to 6e-4.
- **The paths are common:** in all three Δ runs the particle moves about 1.3–1.4 along the same path. At most checks the
  runs agree to 0.01–0.03, and they rejoin at t = 70 after the t = 40–60 separations.
- **Second estimator** (`review/MG020_centre_diagnostic.py`): a quadratic fit to the smoothed field within radius 1.0. On
  the same field it differs from the centroid by 0.02–0.21, typically about 0.1: half to one node spacing. With the
  fit-based centre:
  - Δ 0.1 vs 0.05: every check ≤ 0.116; the t = 60 failure falls from 0.362 to 0.110.
  - Δ 0.2 vs 0.1: t = 40–60 stay at 0.26–0.38.
- **Reading:** part (a) as built compares single-snapshot centres whose estimator noise (~0.1) is close to its tolerance
  (0.15). At Δ 0.1 the failure is consistent with estimator noise. At Δ 0.2 part of the separation is a real difference in
  where the particle went. This is a design flaw in the tick check, recorded as such. Part (b), the window mean of the
  gate quantity D₀(1.5), agreed within 1–2 nodes (about 1–1.6% of \|P̄₀\|) at both tick sizes, against a 25% tolerance.

## 3. Amendment (PI decision, 2026-11-02 21:18 UTC, decision card in the project thread)
1. **Tick Δ = 0.1** for all MG6 runs (proposal §7 check 4's own fallback value); dt = 0.05 (the step check passed at
   Δ 0.1). Δ 0.2 is not used: its tick comparison shows real divergence in part (a).
2. **Part (b) governs the tick check for this campaign.** Part (a) at Δ 0.1 vs 0.05 is recorded as failed and attributed
   to estimator noise (§2); it does not stop the campaign. This is an amendment after seeing data, disclosed as such.
3. **Checks 5–8 are completed for Δ 0.1** by `python src/run_mg6.py part0-mg021` on the PI's PC, with the registered
   criteria unchanged:
   - **5. Bound margin:** max nodal φ_G over every committed Part 0 coupled run (the five tick/step runs, from their
     per-tick `max_phiG`) ≤ 0.7·φ_lim(4) = 0.07; else STOP.
   - **6. Transient:** from the accepted configuration's run, `check_d0.1_dt0.05.npz` (Δ 0.1, dt 0.05), blocks
     B₁ = {60, 70}, B₂ = {90, 100}, with P̄₀ from the passing comparison (the Δ 0.05 run, −146.89). The proposal's
     "finer run of the passing pair" is read as the run with the accepted configuration, as §7's heading says.
   - **7. Quadrature** (i)–(iv) on `final_d0.1_dt0.05.npz`, whose SHA-256 must match `local_manifest.csv`.
   - **8. Timing:** the settle is repeated (deterministic), then 8 concurrent coupled runs of 20 units at Δ 0.1, seeds
     `seed_for(6, 9, 3, 1, 10 + i)`, i = 0–7, as registered.
   - Output `results/MG6/part0/decisions_mg021.json`; `decisions.json` (the STOP) is kept unchanged. The driver reads the
     MG021 file when present.
4. **Budget:** Δ 0.1 about doubles the cost. The provisional launch cutoff rises from 40 to **80 core-h**; the run cap
   stays 24. The budget is restated from check 8 in GATES.md before any registered run, and the PI decides if the
   projection exceeds 80.
5. **Not changed:** the gate, its contrasts and thresholds, the arms, seeds, window rule, measurements and the evaluator.
   The tick-check redesign (a fit-based or smoothly weighted centre, a tolerance set from a measured noise floor, or
   time-averaged positions) is deferred to the next study that needs a tick check; it is not applied retroactively.

## 4. Disclosure
Checks 5 and 6 depend only on committed files. While writing and testing the MG021 code in the cloud, Claude ran them on
those files: bound margin ok, no drift (window stays t = 100–200). These are the registered rules applied mechanically;
the PC run recomputes them and is the record.

## 5. Review
ChatGPT reviews this amendment before any registered MG6 run (README rule 4). Questions for the reviewer:
1. Is Δ = 0.1 with part (b) governing defensible, given that part (a) failed twice?
2. Is the centroid-noise diagnosis adequate, or is a further check needed before the runs (for example, a noise floor
   from two independent Δ 0.1 runs that differ only in candidate seed)?
3. Is the reading of "finer run of the passing pair" for check 6 right?

# MG9b — pre-run implementation specification and calibration result

*Claude, 2026-11-05. Written after registration (draft 3 frozen at `185aceb`, PI card 03:01 UTC) and **before any MG9b
production run**. For ChatGPT's implementation review. No rule of draft 3 is changed here; §C lists implementation
choices for review.*

Code: `src/mg9b.py` (rules and fits), `src/run_mg9b.py` (driver), `src/evaluate_mg9b.py` (evaluation from records),
`src/calib_mg9b.py` (synthetic calibration), `tests/test_mg9b.py` (24 tests after amendment 1, all pass). Changed shared code:
`src/confirm.py` (optional `checks`, `rtol`, `atol`; completion flag) and `src/mg3.py` (records the completion flag in
the confirmation cache). `tests/test_mg3.py` (14) and `tests/test_mg3b.py` (7) still pass.

`python src/run_mg9b.py runs` refuses to start unless `results/MG9b/GO` exists. That file is written only at the PI's
go after the implementation review.

## A. Calibration result (§4.3; cloud, 0.36 core-h of the 4 core-h cap)

`results/MG9b/calibration/calibration.json`. 4000 replicates generated per case. About 83% are retained after the
simulated MG3b Part 1 selection (EQUIVALENT), which matches MG3b's own P(EQUIVALENT | ratio 1) ≈ 0.85.

**Labels per branch: B3 calibrated, B2-20 calibrated, B2-35 not calibrated (not supported: 389 and 178 retained
replicates in the two Q cases, below 400).** So a study that ends with only 29.00 and 29.35 would report any
directional result as EXPLORATORY.

| Case | Retained | Branch B3: n / coverage / neg / pos | B2-20: n / coverage / neg / pos | No trend |
|---|---:|---|---|---:|
| O, S₁ = −0.10 | 3331 | 3331 / 0.952 / 0.732 / 0 | — | 0 |
| O, S₁ = −0.05 | 3350 | 3349 / 0.951 / 0.250 / 0.001 | — | 0 |
| O, S₁ = 0 | 3348 | 3347 / 0.949 / 0.023 / 0.028 | — | 0 |
| O, S₁ = +0.10 | 3331 | 3331 / 0.946 / 0 / 0.766 | — | 0 |
| P, S₁ = −0.10 | 3335 | 3333 / 0.946 / 0.734 / 0 | — | 0 |
| P, S₁ = 0 | 3305 | 3305 / 0.955 / 0.021 / 0.024 | — | 0 |
| Q, S₁ = −0.10 | 3335 | 265 / 0.970 / 0.675 / 0 | 1100 / 0.955 / 0.308 / 0.001 | 0.474 |
| Q, S₁ = 0 | 3346 | 101 / 0.970 / 0.020 / 0.010 | 1203 / 0.952 / 0.025 / 0.023 | 0.557 |

Monte Carlo SE of coverage is about 0.004 at n ≈ 3300 and 0.006 at n ≈ 1100. No numerical failures occurred in any case.

- **Power (scenario O):** 0.73 for S₁ = −0.10 and 0.25 for S₁ = −0.05.
- **Linear-model check:** rejects in 5% of linear B3 replicates (0.049–0.055), as it should, and in 91% of the
  non-monotone case N. In N the trend label was NEGATIVE in 76% of replicates, which is why §5's LINEAR MODEL MISFIT
  rule matters.
- **Design load:** in O the two new b use about 300 runs on average (95th percentile 368), against the cap of 534.
  Capped cells average 0.3–1.0 per replicate. In Q, π_f = 0.2 per screening decision often stops Part A early (29.20
  not placed in about 60% of replicates), which is why its branches are thin.
- Calibration is conditional on screens passing (or failing at random) and on the first-passage model (§4.3).

## B. What the code does (map to draft 3)

- **§2 run protocol:** `mg3.new_state(b, Ω, τ, key, b_det=b)`; the driver asserts b_det = b in every record.
  Streams: `seed_for(9, b_idx, Ω, stage, rep)` with b_idx 1 (29.20) and 2 (29.35), and stages Part A 0, clip 1,
  screen 2, production 3. Gate 9 is unused by earlier stages.
- **§2.1 Part A:** `mg9b.parta_next` (tested on stop, top-cell, screen-failure and end-of-grid cases). The clip pilot is
  rep 0 of the "CLIP" cell (T = 20, no τ halving). Crowding uses the 4 runs.
- **§2.2 Part B:** `mg9b.partb_cells` gives the unique candidate list (the review's 150 → 135, 120 failing example ends
  with 150, 165, 180). Screen = clip pilot (rep 0) plus 4 crowding runs to t = 50 (reps 1–4) in the "SCR" cell.
  Production uses `mg9b.stop20_need` with whole batches.
- **§2.2 identifiability:** `mg9b.identifiable` (two distinct Ω with ≥ 1 event).
- **§3 V-b1:** `run_mg9b.vb1`: the reference continued to 2000, one deep minimum at each of the checks 1900–2000,
  linked; parabolic-vertex drift ≤ λ*/4; depth within 5% of A*. It is run as the first step of the PC run.
- **§3 V-b2:** `run_mg9b.vb2`: the first 10 unique (run, snapshot) fields in the §3 order, a late window 1900–2000 via
  `confirm_field(checks=...)`, completion of all 20 integrations required first. A test shows that nine agreements
  plus one incomplete pair gives CONTROL-NUMERICAL-FAILURE.
- **§3 V-b3:** `run_mg9b.vb3`: all confirmations in the Ω_T cell repeated at rtol 1e-10, atol 1e-12.
- **§4:** `mg9b.fit_b` reproduces MG3b's recorded slope, profile interval and deviance to 1e-8 (test).
  `mg9b.fit_trend` gives the profile CI for S₁ in the bracket [−2, 2] (unbounded sides reported as ±∞) and the 1-df LR
  check with three b.
- **§5:** `mg9b.primary` (insufficient, NUMERICAL, directional, EXPLORATORY for an uncalibrated branch or a
  BUDGET-SHORT cell). `evaluate_mg9b.py` adds the qualifiers.
- **§6 budget:** one ledger (`results/MG9b/ledger.json`), cutoff 50 core-h, checked before every pilot, screen, batch
  and control set. A budget stop writes `budget_stop`, marks the open production cell BUDGET-SHORT, and blocks
  restarts until the PI clears it (tested).

## C. Implementation choices for review (no rule changes)

1. **Confirmation completion.** `confirm_field` now returns `cls = "incomplete"` when the solver returns fewer check
   times than asked. Before, an incomplete integration would have been classified on the shorter list. Default calls are
   otherwise unchanged (test: default and explicit default tolerances give identical results). This affects MG3/MG3b
   replays only if a solver failure occurred, and none was recorded in those studies.
2. **Event-shape "single minimum"** = one deep minimum at t = 300 in the first event's confirmation (`n_deep_end == 1`),
   reported with the n_persistent histogram.
3. **V-b2 timing.** V-b2 runs after a b's Part B (the latest time §6 allows), so its sample can include Part B fields
   if Part A has fewer than ten.
4. **V-b1 platform.** V-b1 runs on the PC as the first step of the production run (§6 order). A cloud preview, not part
   of the record, is reported in §D.
5. **Reference cache.** `runner.reference` for b = 29.20 and 29.35 was computed in the cloud and added to
   `results/MG003/ref_profiles.npz` (committed), so the PC run uses exactly these references. λ* = 3.042 (29.20) and
   3.027 (29.35); A* = −9.78 and −10.15.
6. **Records.** Per-run JSON (committed); confirmation fields per run in `results/MG9b/fields/` (gitignored, on the PC);
   the V-b2 sample's fields are committed as `vb2_<b>_fields.npz`; the V-b1 final field as `vb1_<b>_final.npz`.
7. **Cost projection for the PC.** Calibration suggests about 300 stochastic runs in total. With MG3b's measured
   cost (about 280 core-s per run at τ, less when the run stops 50 after an event), that is about 25 core-h, within
   the projection of about 35.

## D. Cloud preview of V-b1 (not the registered record)

Run in the cloud into a scratch folder (not committed as the record; the PC run repeats V-b1 as its first step).
Both references **pass** with wide margins:

| b | Minima per check (1900–2000) | Linked | Drift (limit λ*/4) | Depth change | max \|RHS\| at 0 / 2000 | Wall |
|---|---|---|---|---|---|---|
| 29.20 | 1 at all 11 | yes | 1.2e-15 (0.76) | 7e-16 relative | 2.8e-13 / 2.3e-13 | 36 s |
| 29.35 | 1 at all 11 | yes | 2.3e-14 (0.76) | 1e-15 relative | 2.3e-13 / 3.1e-13 | 35 s |

The preview found one driver bug before any record existed: V-b1 failed to save its final field when the output
folder did not yet exist. It is fixed, and a test now runs V-b1 on a fresh folder.

## E. Amendment 1 (after ChatGPT MG034 round 3; ADOPTED by the PI, card 03:27 UTC)

Budget truncation is now a study-level condition (`spec/MG9b_amendment_1.md`). Code: `run_mg9b.record_truncation`
(writes `rate_complete` / `budget_truncated` into the ledger's `budget_stop`), `evaluate_mg9b.evaluate` (study-level
flag), `mg9b.primary` (reason "BUDGET-TRUNCATED"). Wording: SLOPE-ELIGIBLE, the bracket note, branch B2-new. Four tests
added (24 in total). The calibration is unaffected: it does not model the budget stop, and B2-new had 0 replicates.

## F. Fixes after the amendment-1 implementation review (MG034)

`review/MG034_amendment1_response.md`: a failed reference now counts as finished rate collection; unstarted accepted
cells are distinguished from BUDGET-SHORT cells; tests use an in-process pool and check that the real output folder is
untouched; event shapes are matched within each run; V-b3 with no pairs is CONTROL-INSUFFICIENT. 27 tests, all pass.

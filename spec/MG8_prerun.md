# MG8 pre-run implementation specification

*Stage MG037. Claude, 2026-11-05. Written after registration (draft 4 frozen at `3ab4ec8`) and before any production
run. It records how the code implements points the proposal leaves open. It changes no registered criterion. Code:
`src/mg8.py`, `src/run_mg8.py`, `src/evaluate_mg8.py`, `src/part0_mg8.py`; tests `tests/test_mg8.py` (18; 24 after §7).*

## 1. Schedulers and random numbers
- `mg8.run_stoch` and `mg8.run_det` implement §3.1 with one numba generator per run. Per step the order draw comes first:
  `np.random.randint(0, N)` for the seam (R and L), `np.random.permutation(N)` for P, nothing for J. The node updates
  follow in processing order, with the per-node draw sequence of `lattice.tau_leap_euler_1d` (G, X and Y exits via
  `_exits`, then the R1 Poisson source).
- **L:** seam s′ = N − 1 − s and order s′, s′ − 1, …; at each node the first hop draw goes to i + 1 and the second to
  i − 1. With the start mirrored about the node index (j → N − 1 − j) and R's seed, L reproduces R mirrored exactly
  (test). Mirror checks: deterministic, every 50-unit checkpoint to 1e-12 relative; stochastic, the SHA-256 of L's
  recorded counts (mapped back to R's orientation) equals R's.
- **Segmented runs:** a seed below 0 continues the generator state. A segmented run is bit-identical to one call
  (test). Deterministic runs advance one time unit per call, and the seed is given only on the first call.
- **Time unit:** one "unit" is round(1/τ) steps (1500 · 2^k), that is 1.00005 time units at every level. Recorded times
  are steps × τ.
- **Seeds:** B1 `seed_for(8, 1, k, 0, j)` (J uses j = 0 but draws nothing; L uses R's j = 0 seed); B2
  `seed_for(8, 0, a, 0, rep)` with a = 0, 1, 2 for J, R, P; the B2 L mirror run uses R's rep-0 seed and the mirrored start;
  Part 0 `seed_for(8, 9, …)`.

## 2. Runs
- **The two continuations (§3.2)** are implemented by running the k = 6 J and R (j = 0) records straight to T = 600, with
  field checkpoints every 50 units. Because the runs are deterministic, their state at t = 300 is identical to that of a
  run stopped at 300 and continued. So there are 55 run records for the 57 registered runs. The T = 300 decisions at
  k = 6 use the t = 250 and 300 checkpoints of those records; the horizon check uses t = 550 and 600.
- **Driver order:** B1 with k = 0 first, then k = 1–6, then B2. Batches of `MODELG_WORKERS` runs; the ledger is checked
  before each batch (a batch already started finishes; §6 "runs in progress finish").

## 3. Measurements
- **Shape alignment:** before alignment, the Nyquist Fourier mode (N = 500 is even) is removed from both fields. Without
  this step a sub-grid Fourier shift is not exactly invertible: a shift and its inverse leave a relative error of about
  2e-5 to 6e-5 on the reference profile, which would set an artificial floor near the BELOW MEASUREMENT RESOLUTION
  bound. The reference's Nyquist amplitudes are 3e-5, 3e-2 and 2.5e-3 (G, X, Y) in FFT units. With the filter, a known
  sub-grid shift is recovered to 1e-10 with E ≈ 5e-11. Part 0's brute-force reference uses the same filtered fields.
- **Centre:** parabolic minimum of smoothed φ_X (σ = 0.4) at its global minimum, wrapped to [−50, 50), then unwrapped in
  time.
- **Stochastic records:** counts are recorded every unit. `measure.analyze_run` runs with t_discard = 150, T = 500, K = 7
  blocks of 50 units, and min_snap_per_block = 1 (its own eligibility flag is not used; §3.3 eligibility is computed
  from its events plus TRACK-FAIL). The per-run profile is the n_per_block-weighted mean of the block profiles. TRACK-FAIL
  (`run_mg8.track_fail`) is the first snapshot at which the x0 track is unmatched while a detection lies within 3 units
  of its last matched position.
- **D_p:** MG2's method (MSD at a lag of 10 snapshots, i.e. 10 units, over eligible runs), descriptive.
- **Fate CIs:** exact Clopper–Pearson 95% intervals for the number of runs with each event type.

## 4. Storage
Committed: `results/MG8/runs/*.json` (records), `results/MG8/fields/*.npz` (B1 three-species fields at every 50-unit
checkpoint), `results/MG8/snap/*.npz` (B2 smoothed φ_X every 5 units, float32), the ledger, `part0/`. PC only:
`results/MG8/full/` (B2 full counts; SHA-256 in each record), gitignored.

## 5. Disclosure
While writing and timing the code (before Part 0), I ran short trajectories on unregistered seeds: deterministic and
stochastic R and P for up to 10 units at τ₀. One deterministic R check moved max |X − X_ref| by about 6.6 over 10 units.
I recorded nothing else from these runs. They are timing and regression checks of the kind Part 0 itself makes, and they
change nothing in the registered design.

## 6. Part 0 record (cloud, 2026-11-05; Linux, Python 3.13, numpy 2.5.3)
- **Attempt 1** (`results/MG8/part0/attempt1/`): **STOP.** Check 1(b) failed on one replicate (R, N = 12, r = 0.05:
  chi-square p = 1.8e-8 < 1e-6). Cause: an implementation fault in the test, not in the simulator. All low-expectation
  cells were pooled into one cell that could itself have an expected count below 5 (0.30 here, with 4 tokens observed;
  the exact Poisson tail probability of that is about 3e-4). Every other part of 1(b) passed: KS p-values 0.12–0.91, all
  moment z-scores with mean within ±0.09 and SD 0.97–1.04, and every across-replicate moment within 4 SE. 1(a), 1(c),
  regression, mirror and synthetic checks passed.
- **Fix:** `part0_mg8.pooled_chisq` merges cells from each tail inward until every pooled cell has expected count ≥ 5,
  the standard reading of §5.1(b)'s "cells with expected count < 5 pooled". It has a test reproducing the attempt-1 case.
  Nothing else changed.
- **Attempt 2** (`results/MG8/part0/`): all checks pass. 1(b) min p = 3.0e-5, KS p 0.15–0.91. 1(c) last ratios R 0.4967,
  0.4984 and P 0.5039, 0.5018, matching ChatGPT's static reconstruction. No clipping in any scheduler. Brute-force
  alignment matches within 0.02%. Projected production cost **3.9 core-h** (cloud timing: 0.78–0.86 core-s per unit
  stochastic, 0.020 core-s per τ₀-unit deterministic), against the 20 core-h cutoff.
- Registered rule: a §5.1 failure is a STOP for the PI. `decisions.json` therefore keeps `proceed = false`, with
  `proceed_if_pi_accepts_attempt2 = true`, until the PI rules.
- **PI ruling (2026-11-05 22:00 UTC):** attempt 2 accepted; `proceed = true`.


## 7. Changes after ChatGPT's implementation review (MG037a, 2026-11-05; before any production run)
Review: `review/MG037_chatgpt.md`; response: `review/MG037_response.md`. No registered criterion changes.
- **Checkpoints by unit, not by time (item 1).** Each recorded unit spans 1.00005 time units, so nominal unit 300 sits at
  t = 300.015 and the old time-bound windows dropped it. Deterministic records now carry `unit` (1, 2, …) and
  `finite_units` (all three fields finite after that unit). The evaluator selects every window and checkpoint by unit,
  with both ends included: SINGLE through U means one detection at every unit 50…U, every unit present and finite up to
  U; v over [U − 150, U] uses units U − 150…U (slope against the physical times). A failure after 300 in a run continued
  to 600 does not affect its assessment at 300. In B2 the `analyze_run` windows already include their end snapshots
  (blocks hold units 150–199, …, 450–500).
- **Implementation checks gate scoring (item 2).** `evaluate_mg8.impl_checks` gives each mandatory check (Part 0
  decision, deterministic mirror at k = 0 and k = 6, stochastic mirror) the status PASS, FAIL or INCOMPLETE. Any FAIL:
  every label is WITHHELD (§4.1 STOP, PI decides). The driver evaluates the checks at startup and after every batch and
  stops dispatching (exit 3, `status.json` state STOP-IMPLEMENTATION) on a FAIL. The production label needs all k = 0
  records and a passed k = 0 mirror; the refinement label needs all B1 records and both deterministic mirrors; B2 is
  marked NOT FINAL until the stochastic mirror passes.
- **Recorded deviation D1: the scored E is computed on Nyquist-free fields (item 3).** The registered text defines e_s on
  the fields themselves; the code removes the Nyquist mode (the alternating component (−1)^j) from both fields first,
  for the reason in §3. A difference confined to that one mode is therefore invisible to the scored E. To show what is
  removed, every E (B1 replicates, horizon, B2 E_stoch) is reported with two companions, which do not enter any label:
  **E_raw**, the registered definition applied literally to the unfiltered fields (own σ*), whose numerical floor is about
  2e-5 to 6e-5; and **E_nyq**, the removed-mode discrepancy max_s ‖Nyq(A) moved by the integer node shift nearest σ* −
  Nyq(B)‖ / ‖φ_B‖ (an integer shift moves that mode exactly; a sub-grid shift of it is undefined on the grid). **PI
  sign-off 2026-11-05 23:37 UTC:** D1 accepted (filtered E scored, E_raw and E_nyq reported); production go given.
- **TRACK-FAIL scans every end of the track (item 4).** After an unmatched snapshot that does not qualify, the scan
  continues; each later end of the track (first unmatched snapshot after a matched one) is tested again.
- **Accounting (item 5).** Run records, the ledger, `status.json` and npz files are written to a temporary file and
  renamed into place; the record is written last and marks completion. On startup the ledger is reconciled with the
  completed records (a missing cost is taken from the record's `wall`, listed under `reconciled`). `status.json` holds
  running / complete / BUDGET-SHORT / STOP-IMPLEMENTATION, and the evaluator applies §6: under BUDGET-SHORT, 4.2 is
  scored only if every k = 0 record exists (else NOT SCORED (BUDGET-SHORT)), 4.3 is NOT RESOLVED BY THIS LADDER
  (BUDGET-SHORT) unless all of B1 completed, and B2 is descriptive on the completed runs. Any other incomplete state gives
  NOT FINAL.
- **Stored inputs and production-kernel checks (item 6).** Each deterministic record stores `seam_sha256`, the SHA-256 of
  its realised seam sequence (int64 little-endian, as drawn; regenerated by `mg8.seam_sha256` after the run; empty for J).
  Each B2 snapshot file stores the run's K = 7 block profiles and `n_per_block`. A supplementary Part 0 check
  (`python src/part0_mg8.py --supplement` → `results/MG8/part0/supplement_MG037.json`, not a rerun of Part 0) runs the
  production kernels diffusion-only: `run_det` for J, R and L over 20 steps matches the product of the conditional maps
  for the regenerated seams to 1.0e-15 (this also confirms the seam regeneration), and one `run_stoch` step from 10⁴
  tokens follows Multinomial(n₀, T_o e₀) given the drawn order (500 replicates each for J and R at N = 12, P at N = 6,
  r = 0.05 and 0.2: KS p 0.17–0.99, minimum p 2.6e-4, aggregate |z| ≤ 2.93). Pass.
- Tests: 24 (six new: endpoints, check gating and budget states, TRACK-FAIL rescan, Nyquist report, production-kernel
  cross-check, seam hash and ledger reconciliation).

# Response to ChatGPT's MG8 implementation review (MG037)

*Claude, 2026-11-05. Review: `review/MG037_chatgpt.md` (commit `e050fbe`). All six items accepted and fixed before any
production run. Details are in `spec/MG8_prerun.md` §7. Tests: 24, all pass (six new).*

**1. [High] Decision checkpoint excluded. Accepted; fixed.** You were right: nominal unit 300 is recorded at
t = 300.015 and the time-bound windows dropped it, so a failure exactly at 300 or 600 went unseen. Deterministic records
now store the unit index and per-unit finiteness, and the evaluator selects every window and checkpoint by unit with both
ends included. SINGLE through U also requires every unit present and finite up to U, so a failure after 300 in a continued
run no longer affects the assessment at 300. Regression test `test_evaluator_endpoints_included` covers zero detections
exactly at 300 and at 600, a non-finite field after 450 in the continued run, and a velocity window whose only change is
at its end unit. B2 was already safe (`analyze_run` blocks include their end snapshots).

**2. [High] Failed mirror checks do not stop scoring. Accepted; fixed.** `evaluate_mg8.impl_checks` gives the Part 0
decision, the k = 0 and k = 6 deterministic mirrors and the stochastic mirror the status PASS, FAIL or INCOMPLETE. Any
FAIL withholds every label (§4.1 STOP). INCOMPLETE gives NOT FINAL (or NOT SCORED where §6 says so), never a label. The
driver runs the same checks at startup and after every batch, and stops dispatching on a FAIL (exit 3, state
STOP-IMPLEMENTATION). Test: `test_evaluator_gates_on_implementation_checks`.

**3. [Medium] Nyquist filtering changes the observable. Accepted; recorded as deviation D1, pending the PI's sign-off.**
The scored E stays on Nyquist-free fields, because the literal definition has a numerical floor of about 2e-5 to 6e-5
from the non-invertible sub-grid shift of that mode, which would bend the refinement slope near the bottom of the ladder.
To make the change auditable, every E (B1 replicates, horizon, B2 E_stoch) is now reported with E_raw (the registered
definition applied literally to the unfiltered fields) and E_nyq (the removed-mode discrepancy at the integer shift
nearest σ*). Neither enters a label. Test `test_nyquist_report`: a difference confined to the Nyquist mode gives scored
E < 1e-9, while E_nyq equals the injected amount and E_raw sees it.

**4. [Medium] TRACK-FAIL stops after the first gap. Accepted; fixed.** The scan now tests every end of the x0 track (the
first unmatched snapshot after a matched one) and continues past ends that do not qualify. Test
`test_track_fail_rescans_after_gap` uses your sequence: an initial empty gap, reconnection, a gap with no nearby
detection, reconnection, then an over-link jump, which is flagged.

**5. [Medium] Resume and budget accounting. Accepted; fixed.** Records, ledger, status and npz files are written atomically
(temporary file, then rename), with the record written last. On startup the ledger is reconciled against the completed
records from their stored wall times. `results/MG8/status.json` persists the state (running, complete, BUDGET-SHORT,
STOP-IMPLEMENTATION), and the evaluator applies the §6 partial-result rules under BUDGET-SHORT. Tests:
`test_seam_hash_and_ledger_reconcile` and the budget cases in `test_evaluator_gates_on_implementation_checks`.

**6. [Medium] Validation gap and missing stored inputs. Accepted; fixed.**
- Deterministic records now store `seam_sha256`, the SHA-256 of the realised seam sequence.
- B2 snapshot files now store the K = 7 block profiles with `n_per_block`.
- A supplementary Part 0 check (`results/MG8/part0/supplement_MG037.json`; Part 0 itself is not rerun) runs the production
  kernels with only diffusion switched on:
  - `run_det` (J, R, L; 20 steps) matches the product of the conditional maps for the regenerated seams to 1.0e-15. That
    also confirms the seam regeneration behind the hashes.
  - One `run_stoch` step follows the exact Multinomial(n₀, T_o e₀) given the drawn order: J and R at N = 12, P at N = 6,
    r = 0.05 and 0.2, 500 replicates each; KS p 0.17–0.99, minimum p 2.6e-4. Pass.

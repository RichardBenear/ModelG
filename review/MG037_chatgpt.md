*ChatGPT implementation review of MG8 (code, tests, pre-run spec and Part 0 at `e050fbe`), pasted by the PI 2026-11-05 22:06 UTC.
Saved verbatim.*

Recommendation: fix several implementation issues before production. The update kernels look consistent with the agreed scheduler semantics, but the evaluator can miss endpoint failures, and failed implementation checks do not currently prevent scientific scoring.
I reviewed commit `e050fbe`, including the registered Draft 4, pre-run specification, kernels, driver, evaluator, tests, and both Part 0 attempts. I ran 17 of the 18 supplied tests in an isolated scratch copy; all passed. I also reproduced the endpoint-classification bug below. I did not run production or modify repository files.

1. [High] Time-window selection excludes the decision checkpoint.
The driver records nominal checkpoint 300 at actual time `300.015`, because each recorded unit spans `1.00005`. But `single_through(..., 300)` and `window_v(..., 150, 300)` use an upper bound of `300.000001`.
Consequently, the evaluator excludes checkpoint 300. I reproduced both failures:
   * A zero-detection count at checkpoint 300 still returns `single=True`.
   * A zero-detection count at checkpoint 600 also returns `single=True`.
Velocity windows likewise omit their nominal final checkpoint. In addition, the single `finite` flag covers the entire continued run, allowing a failure after 300 to contaminate the assessment at 300.
Fix: distinguish nominal checkpoint indices from physical times consistently, include the intended endpoints, and store/check numerical validity by checkpoint. Add regression tests with failures exactly at 300 and 600. [Evaluator, lines 28–37](https://github.com/RichardBenear/ModelG/blob/e050fbe9634964dc60272b28f5f6fa3b8d770a4c/src/evaluate_mg8.py#L28)
2. [High] Failed mirror checks are reported but do not stop scoring.
The evaluator calculates deterministic mirror results and the stochastic mirror hash comparison, but neither controls whether production/refinement labels are emitted. The driver also proceeds through later batches without checking completed production mirror runs.
Thus a failed registered implementation check can coexist with a scientific verdict. Missing mirror records also do not prevent scoring.
Fix: require all applicable implementation checks to pass before issuing final scientific labels. Distinguish an incomplete check from a failed check, and stop further dispatch when a completed mandatory check fails. [Mirror evaluation](https://github.com/RichardBenear/ModelG/blob/e050fbe9634964dc60272b28f5f6fa3b8d770a4c/src/evaluate_mg8.py#L88), [driver dispatch](https://github.com/RichardBenear/ModelG/blob/e050fbe9634964dc60272b28f5f6fa3b8d770a4c/src/run_mg8.py#L169)
3. [Medium] Nyquist filtering changes the registered shape observable.
`shape_error` removes the Nyquist mode from both fields before computing the discrepancy. The pre-run document explains the numerical motivation, but the frozen proposal defines the discrepancy on the fields themselves.
This is a defensible measurement choice, but it changes what E can detect. A difference confined to the removed mode becomes invisible. The filtered brute-force benchmark validates agreement with that modified observable; it does not establish equivalence to the registered one.
Fix: document this explicitly as a measurement amendment or recorded deviation before production, including its application to stochastic profiles. Preferably report the removed-mode discrepancy separately. [Filtering implementation](https://github.com/RichardBenear/ModelG/blob/e050fbe9634964dc60272b28f5f6fa3b8d770a4c/src/mg8.py#L342), [pre-run rationale](https://github.com/RichardBenear/ModelG/blob/e050fbe9634964dc60272b28f5f6fa3b8d770a4c/spec/MG8_prerun.md)
4. [Medium] TRACK-FAIL stops searching after the first unmatched snapshot.
`track_fail` returns `False` immediately when the first unmatched snapshot lacks a nearby detection. It never checks later snapshots.
An original track can reconnect after a temporary gap and subsequently suffer a genuine linking failure. That later failure would be missed, potentially making the run eligible for drift and shape estimates.
Fix: continue scanning after an unmatched snapshot that does not meet TRACK-FAIL criteria. Add a synthetic sequence containing an initial empty gap, reconnection, and a later over-link jump. [Tracking check](https://github.com/RichardBenear/ModelG/blob/e050fbe9634964dc60272b28f5f6fa3b8d770a4c/src/run_mg8.py#L92)
5. [Medium] Resume and budget-stop behavior need to preserve the registered accounting.
A worker writes its completed run record before the parent updates the ledger. If execution stops between those operations, restarting skips that record but does not recover its cost from the stored `wall` value. The ledger can therefore undercount completed work.
Separately, budget stopping only prints `BUDGET-SHORT`; the evaluator does not implement the specified partial-result labels. An incomplete ladder can simply lack a refinement result.
Fix: reconcile the ledger with completed run records on startup, write records atomically, and persist an explicit completion/budget status that the evaluator uses for partial scoring. [Driver, lines 157–188](https://github.com/RichardBenear/ModelG/blob/e050fbe9634964dc60272b28f5f6fa3b8d770a4c/src/run_mg8.py#L157)
6. [Medium] Some promised validation and reproducibility inputs are missing.
The deterministic matrix benchmark exercises the separate `diffusion_map` helper, rather than the production `run_det` path. The stochastic cohort benchmark likewise uses a separate benchmark kernel. These validate the intended algorithm, but leave a gap between the exact benchmark and production code.
Also, production records omit the promised realised seam-sequence hashes and stochastic block profiles. They retain the weighted profile and block counts, which cannot independently reconstruct the weighting operation.
Fix: add a small diffusion-only cross-check against the actual production kernels, and save the specified hashes and block profiles. [Matrix helper](https://github.com/RichardBenear/ModelG/blob/e050fbe9634964dc60272b28f5f6fa3b8d770a4c/src/mg8.py#L210), [stochastic record construction](https://github.com/RichardBenear/ModelG/blob/e050fbe9634964dc60272b28f5f6fa3b8d770a4c/src/run_mg8.py#L144)

The Part 0 pooling correction appears reasonable, and the failed first attempt and PI acceptance are transparently recorded. Its operator-refinement ratios agree with my earlier independent calculation. The passing kernel regression, mirror, conservation, and synthetic tests are useful evidence.
Disposition: retain the registered scientific design and accepted Part 0 record. Correct the evaluator and driver issues, resolve the filtered-observable discrepancy, and run focused regression checks before the production go.

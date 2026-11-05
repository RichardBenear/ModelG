# MG032 — MG9d Draft 3 implementation review

Reviewed revision: `27e0d4cf81d0235daf9c60b4bfcd5962cfc0fe4d`. The proposal identifies its registered scientific content as frozen at `69ad83e` (PI card 2026-11-05 02:00 UTC).

**Verdict: the scientific implementation and static calibration pass this review, but fix two failure/budget-handling issues before launching production. No new design cycle is needed.** These are implementation defects, not objections to the PI's accepted residual risks or requests to change the registered thresholds.

## 1. Checks that pass

- All 16 repository tests pass locally, including the short MG003-state solver-wrapper equivalence check, the synthetic end-to-end run, reflection sign handling, and the Draft 3 worked classification cases. I ran no MG9d production trajectory.
- Independently reconstructed the common state by explicitly subtracting each field's alternating Nyquist mode; generated fractional shifts using a full complex FFT and the finer state using SciPy resampling. All seven starting states agree with production within `8.9e-15`. The original source SHA-256 matches the committed manifest. Filtering changes the source by `1.5622865851e-5`, below the registered `2e-5` limit.
- Independently expanded the reaction polynomial and used array-based finite differences for diffusion. All recorded initial maximum residuals reproduce within `2.42e-13`. The source, R0, R2b and R4 residuals are respectively approximately `1.09217e-5`, `0.001654685`, `0.36680786`, and `4.85093386`. The substantial interpolation-induced residuals are real and already disclosed; transfer accuracy does not imply dynamical equilibrium.
- Checked the analytic Jacobian against directional finite differences of that independent RHS on R0 and R4; relative errors are below `1.5e-10`.
- Reran all 24 production static calibration cases and independently recovered known full-FFT shifts by solving the derivative of the alignment objective. Independent recovery error is below `3e-16`; production recovery error is `1.1104158704e-9`. All cases match all 16 interior minima and have no alignment ambiguity. The registered floor remains `v_f = 5e-9`, with resolved speed threshold `1.5e-8`.
- The reproduced largest edge-midpoint error is `2.28574e-4`, and the minima-mean error is `7.56416e-5`. These support the adopted `5e-4` absolute consistency allowance for the tested profiles; they are not universal error bounds for deformed profiles.
- The implementation uses the adopted reference time 2000 for validity items 5–6, while retaining whole-run identity and width checks. Matched-minima displacement uses one common identity subset. Endpoint slopes, descriptive least-squares slopes, C1–C5, and overall classification order agree with the registered text. The previously flagged endpoint-versus-regression ambiguity is resolved.

## 2. Required fix: preserve early numerical failures as incomplete records

Location: `src/run_mg9d.py:103–117`, with persistence at `412–454`.

When all solver attempts fail before the first requested output time, there are no saved solver samples. `mg9d_ode` assumes `sol.t` has a `.size` attribute and then unconditionally indexes `sol.y[:, -1]`. Mocked empty-list outputs produce `AttributeError`; normalized empty arrays produce `IndexError`. Both probes perform the three configured attempts and then crash instead of returning `complete=False`.

The exception escapes `one_run` before `_ledger_finish` and `save_rec`, leaving neither an incomplete run record nor a completed runtime entry. A later invocation can start the same arm again because its record is absent. This loses the intended failure provenance and prevents the evaluator from treating that attempted arm explicitly as INVALID.

**Requested correction:** handle zero returned samples explicitly, preserve the last available state (or input state when none was returned), and record failure status, message and last saved time without inventing a successful endpoint. Finalize elapsed-time accounting and preserve an incomplete record for handled solver exceptions as well. Retain attempt history if numerical reruns are permitted. Add tests for both zero-sample failure and failure after partial output, checking the saved record, INVALID classification and ledger accounting.

## 3. Required fix: make the runtime cutoff effective across interruption and long runs

Location: `src/run_mg9d.py:367–385`, `412–445`; the evaluator uses the same sum of recorded runtimes.

The run-count cap is correctly counted at start. Runtime accounting is weaker: `wall_s=None` is summed as zero, and elapsed time is finalized only after the run loop returns normally. My mocked ledger containing an unresolved prior attempt permits another start. A ledger at 3599 seconds also permits another complete run; there is no runtime check within that run. Thus the registered one-core-hour cutoff is currently only a check of completed earlier runs, not an enforced total budget.

**Requested correction:** do not silently treat an unfinished attempt as free. Reconcile its consumed runtime, or stop for PI reconciliation before resuming. Persist elapsed runtime during execution and define/enforce a stopping mechanism for the remaining budget. If only segment-boundary checking is intended, disclose the possible overshoot and obtain a PI-approved operational policy rather than claiming a strict cutoff. Add interruption/resume and near-cutoff tests. A budget stop should leave auditable partial output and a clear stop reason.

## 4. Small cleanup, not a scientific blocker

`stage_runs()` calls `stage_checks(with_relax=False)`, which overwrites the committed `transfer.json` and removes its existing descriptive `relax_estimate` fields. Preserve the frozen pre-run calibration or write runtime verification separately. Also update the stale “not yet registered” pre-run introduction and `PROPOSED` comments to distinguish historical proposals from adopted settings.

## 5. Recommendation to the PI

Keep the registered science unchanged. Ask Claude to fix the two operational issues and demonstrate the added failure-path tests; then proceed with the seven planned arms after the implementation check and your launch approval. The starting-state residuals, finite horizon, and possible R4 deformation remain accepted interpretive limitations, not newly discovered errors. I have not independently recomputed the descriptive pseudoinverse relaxation estimates or asserted that transients finish by time 2000.

Audit script and machine-readable evidence: `review/MG032_replication/audit.py` and `audit.json`. No MG9d production integrations or repository source edits were made in this review.

# MG012: independent review of MG4c results

Reviewed commit: `37c9c34b6574c77b42d1a2d1e65957754b518d7a`, 2026-10-31.

**Decision:** W1 and the observed weighted single-particle outcome are reproduced. The T operator's tolerance convention/deviation needs explicit documentation before an unqualified claim that every registered construction check passed. No new trajectory campaign is indicated. See [review.md](review.md).

Reproduce using Python with NumPy, SciPy, Numba, NetworkX and mpmath:

```
python review/MG012_replication/download_inputs.py
python review/MG012_replication/replicate.py
python review/MG012_replication/audit_tolerance.py
python review/MG012_replication/inputs/tests/test_mg4.py
python review/MG012_replication/verify_tests.py
```

The download is about 117 MB, pinned to a commit and verified by Git blob hashes. `inputs/`, temporary test outputs and caches are excluded from this package. Original input bytes remain in the public repository.

- `replication.json`: independently reconstructed classes, thresholds, brackets, shapes, RHS and balances for all 102 runs and 1,122 saved checks. Graph positions/edges use registered builders. Measurements and ODE/finite-element calculations are separately written. Weighted geometry uses a different periodic-face selection and gradient-based FEM stiffness assembly, rather than the source cotangent formula.
- `tolerance_audit.json`: registered entrywise 1e-12 versus implemented scaled tolerance, with measured matrix and T-field differences.
- `test_results.json`: 10 MG4c tests passed, 2 upstream fork-dependent tests skipped on Windows. The separate MG4 test command passed all 11 tests. The wrapper places temporary files in this workspace because the initial direct invocation could not create test subdirectories in the sandbox's default temporary directory.

Registered full-size trajectories were not rerun. Upstream validation includes small unregistered-graph integrations. No upstream proposal, source, result or gate file was changed. This review package has not been committed or pushed.

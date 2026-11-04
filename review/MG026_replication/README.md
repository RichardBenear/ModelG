# MG026 static design checks

Review target is MG9 Draft 1 at 99229bc0beb849f877fae17054aba50cddcb770e.

Run with Python, NumPy and SciPy.

```text
python review/MG026_replication/design_checks.py
```

The script constructs the Model G Jacobian directly, finds the stationary onset over continuous Laplacian eigenvalues and over finite periodic lattice modes, and computes finite-mode growth-rate derivatives. It also totals the integrations listed in the draft.

design_checks.json contains the output. No nonlinear integration, Landau cubic calculation, branch continuation or front-pinning test was performed.

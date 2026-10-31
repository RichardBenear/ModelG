# MG010: independent review of MG4 results

Reviewed 2026-10-31 at ModelG commit `41ee3de73d30e80fc0aeb95350dc00070c1c6f58`.

**Decision:** accept the scientific FAIL for the registered panels. Reconcile the final record with the evaluator's unresolved status. See [review.md](review.md).

Reproduce with Python, NumPy, SciPy, Numba, and NetworkX:

```
python review/MG010_replication/download_inputs.py
python review/MG010_replication/replicate.py
python review/MG010_replication/inputs/tests/test_mg4.py
```

The manifest pins inputs to the reviewed commit and Git blob hashes. Downloads are about 103 MB, mostly already-published fields. Inputs are ignored here to avoid duplicating repository data. `verified_inputs.json` records the initial download; the test source was subsequently fetched at the same pinned commit and added to the manifest.

`replicate.py` reconstructs graphs using registered builders, then uses independently written detection, centroid/link checks, Gaussian smoothing, graph Laplacian, and 40-digit Decimal balance arithmetic. It checks all 80 runs at all 11 saved times, and C3(b) at all 68 main-panel placements. It reproduces T's spectral count analytically and checks saved irregular-graph spectral arithmetic. It does not repeat full trajectories or solve full D/R/N3 eigensystems. Results: [replication.json](replication.json).

All 10 upstream MG4 tests passed locally. No upstream sources, criteria, records, or gate entries were changed. Downloaded inputs and Python caches are excluded from this replication package.

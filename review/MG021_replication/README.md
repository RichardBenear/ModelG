# MG021 amendment review: MG6 Part 0 STOP and proposed exception

Reviewed at `ba765876dfde029388d5769a0cd690d25b38aeee`, amendment blob `917a5c82bc32ce54a1902d35a63a9494c2473e6f`.
Full review: `../MG021_chatgpt.md`.

```text
python review/MG021_replication/download_inputs.py
python review/MG021_replication/record_checks.py
python review/MG021_replication/guard_probes.py
```

The manifest pins 81 public files (1,816,504 bytes); all Git blob hashes were verified. Requires Python, NumPy, SciPy, Numba and mpmath. Inputs and caches are excluded from publication.

`record_checks.py` does not import the simulator. It separately generates the prescribed candidate/lifetime stream using NumPy and its specified seed/draw order, rebuilds full node sets at all 55 recorded checks, and compares sparse node positions/IDs, cumulative counts, coarse annular counts, sparse detections and the registered step/tick statistics. The generated candidate hash is reported for reproducibility; the Part 0 JSON does not expose a reference candidate hash against which to compare it. The reconstructed positions and counts match.

It also independently fits local quadratics with radii 1 and 1.5, verifies positive-definite Hessians and stationary points inside the fitting discs, and recomputes counts about those centres and common midpoint centres. These are exploratory sensitivity calculations only. They do not rescore Part 0, establish a noise variance or prove an estimator is unbiased. No ODE trajectory or new physical run is generated.

`guard_probes.py` calls the reviewed `main_amended` with physics, timing and file fixtures replaced by mocks. It demonstrates that this pinned code writes `proceed: true` for missing manifest verification, failed timing/settling fixtures, timing margins exceeding 0.07, and a projection above 80 hours. It also exposes that a refined selected grid is not passed into timing jobs. These are control-flow tests of invalid hypothetical inputs, not failures observed in MG6 results. Its assertions intentionally verify the holes present in the pinned revision. Temporary fixture directories are created and removed only within this review's ignored inputs directory.

No production campaign was run, no STOP was changed, and no upstream simulator, proposal, result or gate was edited. PC-only final fields were unavailable; quadrature check 7 and hardware timing check 8 were not independently reproduced.

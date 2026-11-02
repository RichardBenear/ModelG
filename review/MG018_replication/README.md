# MG018: MG5 results review and replication

Reviewed snapshot: `aa237e0f892cc4709ecc87cf72cd483f8ec8d536`.
The authoritative review is [`../MG018_chatgpt.md`](../MG018_chatgpt.md).

```text
python review/MG018_replication/download_inputs.py
python review/MG018_replication/review_checks.py
python review/MG018_replication/replay_checks.py
python review/MG018_replication/inputs/src/recompute_mg5.py
```

Requires Python, NumPy, SciPy, Numba and mpmath; the tested environment is recorded in `review_checks.json`. The download manifest pins 119 files (17,555,939 bytes), including all MG5 checks, run records, Part 0 records, relevant source and specification. Every input was checked against its Git blob SHA-1; hashes identify repository content, not experimental correctness. Downloaded inputs and caches are excluded from publication.

`review_checks.py` does not import the MG5 simulator or either evaluator. It implements detection with periodic KD-tree neighborhoods and component traversal, then independently calculates every count sequence, multiplied subtype, single-particle track, mobility status, MSD fit, aggregate outcome, paired Part 0 comparison, jump identity, population accounting, depth comparison and ledger total. Detection starts from saved smoothed X values; smoothing itself and the ODE solution are not reproduced.

`replay_checks.py` uses the pinned registered event generator to reproduce all 16 event lists, checks their SHA-256 values with separately written serialization, then independently checks tick counts, ephemeral nodes, sparse node lifetimes and newborn positions. Event lists remain in memory and are not saved. This is canonical event replay, not an independent stochastic generator or a re-integration of the dynamics. `code_provenance.json` records that the simulator source blob is identical at both production revisions.

The repository's separately written `recompute_mg5.py` also agreed on all 24 run classes and all single-run mobility labels/displacements. Its script currently prints differences rather than returning a failing exit status, and does not independently aggregate the gate or recompute D_p; the review script checks these explicitly.

No registered trajectory was rerun. The local full-field files listed in the manifests were not available, so their checksums, ODE residuals and original smoothing/interpolation calculations could not be independently verified. No proposal, simulator, result, threshold, status or gate file was modified.

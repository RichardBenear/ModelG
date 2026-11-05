# MG033 — independent MG9d result replication

Reviewed revision: `d3bcfd35ad5b85a88ac59a6cdfd1e4cd78f06d98`.

```text
python review/MG033_replication/audit.py /path/to/ModelG
```

Requires NumPy and SciPy. Reads committed MG9c/MG9d data, imports no ModelG production modules, runs no integrations, and writes `audit.json` beside the script. Input archives are not duplicated here.

Independently reconstructs starts, snapshot hashes, Fourier alignment through analytic objective-derivative root finding, envelope crossings, weighted centres, parabolic minima, and expanded reaction/diffusion residuals. Rechecks complete record schedules, validity, solver settings, segment attempts, the actual registered decision branches, and ledger totals. It is a result-specific audit, not a replacement generic evaluator.

All 77 saved field states are checked. Intermediate diagnostic/tracker records are checked as recorded data because full fields were saved only every 1000 units. No independent production trajectory replication is claimed.

Separately ran `python tests/test_mg9d.py`: all 22 repository tests pass. Those tests include short non-MG9d solver checks and synthetic failure/run records, not production MG9d continuations.

Full review: `review/MG033_chatgpt.md`.

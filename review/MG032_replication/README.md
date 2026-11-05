# MG032 implementation audit

Reviewed ModelG revision: `27e0d4cf81d0235daf9c60b4bfcd5962cfc0fe4d`.

Run from a checkout with the committed MG9c snapshot archive and MG9d source:

```text
python review/MG032_replication/audit.py /path/to/ModelG
```

Requires NumPy and SciPy. Writes `audit.json` beside the script. It does not modify the supplied checkout or perform any time integration. The source snapshot is read from the checkout; no input archive is duplicated here.

Checks include independently constructed states, expanded RHS, Jacobian directional differences, production calibration rerun, independent alignment recovery, and mocked empty-solver-output and ledger failure paths. The latter record defects present at the reviewed revision; they are not successful-behavior assertions. After a fix, changed failure-probe output is expected.

The repository's separate `python tests/test_mg9d.py` suite also passed (16 tests). That suite includes a short integration of its MG003 reference state and synthetic MG9d records, not production continuation of the MG9d starting state.

Full findings: `review/MG032_chatgpt.md`.

# MG035 result replication

Reviewed ModelG revision `5e476adaeb3b57b231c88ebb3872c762df10becc`. Run with Python, NumPy and SciPy from a checkout with the committed MG003 reference archive and MG3b/MG9b results.

```text
python review/MG035_replication/audit.py /path/to/ModelG
python review/MG035_replication/replay_controls.py /path/to/ModelG
```

Each script writes JSON beside itself. The first reconstructs first-passage records, placement, stopping, rate fits, control pairing and ledger totals. The second integrates the 20 archived control fields with an independent deterministic solver and checks both classification windows. Neither script imports ModelG production functions. The full review distinguishes field replay from verification of recorded control pairs.

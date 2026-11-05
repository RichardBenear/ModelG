# MG031 Draft 1 static checks

Reviewed revision `02e3a37ee4328ef9900673f1cec9ccd91010db50`.

Run `python review/MG031_replication/design_checks.py CHECKOUT` with NumPy and SciPy. CHECKOUT must contain `results/MG9c/snapshots.npz` and the b = 29.30 front record. Output is `design_checks.json` alongside the script.

Checks: actual starting-state Fourier shift/round-trip errors; a smooth relaxation fit to existing midpoint tracks demonstrating the C0 time-window confound; and explicit near-zero/ratio logic examples. No production modules are imported and no ODEs are integrated. Extrapolation is a counterexample to the decision rule, not a forecast gate or newly generated physics result.

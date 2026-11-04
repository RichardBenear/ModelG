# MG028 replication

Reviewed ModelG revision 6e3234b66d92e48e116cbc5e8c81211660a06ab7.

Run with Python, NumPy and SciPy.

```text
python review/MG028_replication/audit.py CHECKOUT
```

The checkout must contain results/MG9. audit.json stores the output.

The script imports no production Model G modules. It reconstructs the cubic coefficient using phase quadrature and linear solves, refits the recorded amplitude series, rescores the actual front trajectories relevant to the observed outcomes, verifies saved field hashes and preparation residuals, and sums the ledger.

It does not rerun nonlinear evolution or reconstruct front tracks from time-resolved fields, which are absent from these committed records.

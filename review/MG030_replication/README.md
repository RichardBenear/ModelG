# MG030 independent replication

MG030 is the MG9c late-time-front result; MG9b remains reserved for the stochastic barrier study.

Reviewed revision `810ad932de904d633c1411c5a5c6593e8cc4e47a`.

Run with Python, NumPy and SciPy:

```text
python review/MG030_replication/audit.py CHECKOUT
```

The checkout must contain `results/MG9` and `results/MG9c`, including the committed snapshots. Output is `audit.json` beside the script. Inputs are not duplicated here.

The script imports no production modules and runs no ODE integration. It independently checks all nine observed trajectory classifications, the three replay records, 50 snapshot hashes, the snapshot schedule, final residuals from the expanded reaction polynomial, FFT-convolution envelope crossings and cell counts at saved times, preparation statuses, perturbation/width differences, and the ledger. Its checks for earlier terminal/structural rules are assertions tailored to these observed event-free tracks, not a general replacement classifier.

The translation projection is a post-result descriptive calculation: project the deterministic RHS onto minus the Fourier spatial derivative in the unweighted Euclidean norm of the three potential fields. The reported relative remainder measures departure from pure instantaneous translation under that diagnostic; it is not a registered gate, an energy norm, a traveling-wave proof or a continuum-limit test.

Separately, all 20 repository synthetic tests in `tests/test_mg9c.py` passed using that file's standalone runner. The local environment is Python 3.13 on Windows with NumPy 2.2.4 and SciPy 1.16.1; the production trajectories were generated in the registered Linux environment. No nonlinear replay was attempted locally.

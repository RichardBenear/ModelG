# MG025 replication

Reviewed ModelG revision 0b86bf7 on 2026-11-04, using Python 3.13 and NumPy 2.2.4 on Windows.

The checkout needs src, results/MG7, results/MG4/part0, and the eight endpoint fields named in results/MG7/endpoint_manifest.csv. The endpoint fields used here came from the earlier local MG4c archive and were hash-verified by the strict replay against that manifest.

Run from the review workspace, replacing CHECKOUT with the path to that revision.

```text
python CHECKOUT/src/recompute_mg7.py
python review/MG025_replication/audit.py CHECKOUT
```

strict_replay.txt records the first command. audit.json records the second. The supplemental audit regenerates candidate hashes and every edited-boundary population independently, computes count-based window means and timing totals, and checks Part 0 count and depth differences. It reuses the repository's independently written MG6 replay detector for the alternative-threshold check, as disclosed in the review.

The ODE evolution and full-field quadratures were not rerun. Inputs remain outside this package.

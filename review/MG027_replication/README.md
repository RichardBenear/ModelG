# MG027 envelope replication

Reviewed revision f7ef9e9, MG9 amendment 1, on 2026-11-04.

Run with Python, NumPy and SciPy, supplying a checkout containing src and tests.

```text
python review/MG027_replication/envelope_checks.py CHECKOUT
```

envelope_checks.json records independent convolution and edge-location checks on synthetic cosine patches, and five selected repository synthetic tests. The tests use no Model G evolution or primary MG9 outcome.

The independent envelope uses fractional physical cell overlaps and FFT circular convolution. Its crossing measurements do not use the production tracker. The selected repository tests separately exercise that tracker and its classifications.

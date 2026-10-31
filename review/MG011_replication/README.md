# MG011: MG4c proposal reviews and design checks

Reviewed 2026-10-31 at repository commit `7d2b174800adb5f12809a3e1c5867ba35709b4f8`.
Proposal blob: `00ee82eb2add64074f7401f0152a15d4bbbdfadf`.

Draft 1 review: [review.md](review.md). Draft 2 review: [review_r2.md](review_r2.md), at merged PR #19, commit `877e9a6858aa69ff6678c96af1cc542bf8a3208a`, proposal blob `5ca1f4b30dffe88a2d769675b57256e3f3c9aa44`.

Draft 3 review: [review_r3.md](review_r3.md), commit `e9c5e81cbc825512fe9cb1887d24929a862c5c38`, proposal blob `cee43bed8f7cbae7ad728e9d5dcc473e4227c309`.

Current recommendation: support registration of Draft 3. All round-2 design blockers are resolved; no additional proposal revision or expanded experiment is requested.

`design_checks.py` checks the trial arithmetic, an invalid panel-edge bracket, a gap in the W categories, a narrow interval missed by six binary subdivisions, and T's weighted normalization. Run:

```
python review/MG011_replication/design_checks.py
python review/MG011_replication/design_checks_r2.py
python review/MG011_replication/design_checks_r3.py
```

Outputs are [design_checks.json](design_checks.json) and [design_checks_r2.json](design_checks_r2.json). The latter demonstrates the self-normalized threshold's amplitude invariance on constructed decaying fields, verifies the corrected 102-run cap and refined W1 coverage, and supplies O2 and launch-cutoff counterexamples. These scripts need only the Python standard library.

The round-3 check [design_checks_r3.json](design_checks_r3.json) confirms that the new fixed eligibility floors reject the round-2 vanishing-amplitude counterexamples and admit the finite-amplitude example. It does not certify a stationary particle or predict the actual T scan.

These are arithmetic and decision-logic examples, not Model G simulations or predictions. No MG4c dynamics were run. Reference inputs and Python caches are excluded from this package. No upstream proposal, source, result or gate file was changed by these reviews.

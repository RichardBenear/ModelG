# MG014: MG5 proposal reviews and static checks

Latest: **[Draft 2 review](review_r2.md)**, at commit `a1c4b14c81cd5f8c72956b4ac7da9cc05153d491`, proposal blob `470d5621e12fb689e1444d83d5011266afaf03d0`.

Draft 2 resolves the main scientific-design objections. Two corrections remain before registration: the event-list storage commitment and the disposition of interpolation-audit failures. Implementation and Part 0 can proceed once those are fixed; production remains conditional on validation and controls. No MG5 trajectory or event stream was generated.

## Draft 2 checks

After downloading the default inputs below:

```text
python review/MG014_replication/draft2_checks.py
python review/MG014_replication/download_inputs.py --manifest download_manifest_r2.json --output inputs_draft2 --verification verified_inputs_r2.json
```

The second command retrieves the pinned Draft 2 documents. Locally these were fetched through the GitHub connector and blob-hash verified after a direct shell download failed due to unavailable DNS/network access. `draft2_checks.*` adds storage arithmetic, tick/step alignment including fallbacks, batch bookkeeping, and full/sparse detection comparisons on four stored states and nine synthetic cases. All passed. It uses the inspected existing detector, not an independent implementation. Requirements: Python, NumPy, SciPy, Numba and mpmath.

## Historical Draft 1 review and checks

Reviewed 2026-10-31 at commit `3179342a9b559af06beabdfea98890c627a97e12`.
Proposal blob: `9d6dd2b6800ffcbcef61923a232180c0f341c944`.

Recommendation: support the study, but revise the event/transfer implementation, frozen controls, numerical validation and audit records before registration. See [review.md](review.md).

```
python review/MG014_replication/download_inputs.py
python review/MG014_replication/design_checks.py
```

The input manifest pins 59 public files (about 5.3 MB), including the four MG4c b=28 states. Inputs are excluded from the review package. The local run reused matching cached files and verified all 59 Git blob hashes. Downloading from the pinned public commit reproduces them.

`design_checks.json` contains event-work arithmetic, an interpolation-bias example, the Part B background conversion, and static RHS evaluations of the four stored states under their original h and MG5's proposed h=0.2. No MG5 trajectory, event stream or particle outcome was simulated. The static RHS calculation reuses the inspected registered graph builder and equation implementation; it is not an independent ODE implementation.

No upstream proposal, source, result or gate file was edited. Publication includes only this review directory's text files, excluding reference inputs and caches. `review.md` and `design_checks.*` remain Draft 1 records; their h=0.2 and Delta=0.25 objections are resolved by Draft 2.

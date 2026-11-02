# MG019: MG6 Draft 1 proposal checks

Proposal reviewed at `53a2adc89e330d2ca2686af9bde1665b0f57003e`; proposal blob `ea73629494a7cd9d4950c1bb29dae17e1cc4a319`. `proposal_snapshot.json` preserves the reviewed proposal, design calculation, roadmap and glossary as returned by the GitHub connector, with their blob IDs.

These are static design checks, not MG6 simulations. The four MG4c b=28 float64 states and source are pinned separately in `download_manifest.json` at `3179342a9b559af06beabdfea98890c627a97e12`. All referenced files were Git-blob-hash verified. Locally, the script reused the previously verified MG014 input cache. To reproduce independently from the repository root:

```text
python review/MG019_replication/download_inputs.py
python review/MG019_replication/static_checks.py
```

Requires Python, NumPy, SciPy, Numba and mpmath. Inputs and caches are ignored. The script uses the inspected MG4 graph builder and detector to place the four stored particles; it independently evaluates the proposed coupling functions, clipping, count-response arithmetic, global-baseline correction, and the Student-t tail probability for a statistic of 3 from four independent Gaussian run summaries. This last calculation illustrates a statistical assumption; it does not establish that real run summaries are Gaussian or identically distributed.

The design-level equal-area quadrature is used only on the initial near-uniform point sets. It is explicitly not offered as valid integration over the later spatially biased MG6 node set. No Voronoi-area reconstruction, MG6 trajectory, event stream or particle outcome was generated. No upstream file was modified or published.

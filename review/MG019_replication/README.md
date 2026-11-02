# MG019: MG6 Draft 1 proposal checks

Round 2: `../MG019_chatgpt_r2.md` reviews Draft 2 at `92d8f2b8505526f1e6e559d450352682f574d6a3`. `proposal_snapshot_r2.json` preserves and hash-identifies Draft 2, the archived Draft 1 and Claude's response. `draft2_checks.py` / `draft2_checks.json` add exact-torus grid arithmetic and periodic barycentric quadrature of the four stored starting G profiles. Run after obtaining the original manifest inputs:

```text
python review/MG019_replication/draft2_checks.py
```

The check uses a commensurate midpoint grid near spacing 0.05 and 0.1-wide radial bins to estimate the first G zero. It finds no complete eligible falloff bins under Draft 2's specified fit range for any of the four stored states. It also calculates the fraction of ephemeral candidates skipped before acceptance. This is a static design check, not an evolving-profile result or a certified grid-convergence study. Round-2 snapshots were verified against all three Git blob hashes; the original Draft 1 files remain historical records.

Proposal reviewed at `53a2adc89e330d2ca2686af9bde1665b0f57003e`; proposal blob `ea73629494a7cd9d4950c1bb29dae17e1cc4a319`. `proposal_snapshot.json` preserves the reviewed proposal, design calculation, roadmap and glossary as returned by the GitHub connector, with their blob IDs.

These are static design checks, not MG6 simulations. The four MG4c b=28 float64 states and source are pinned separately in `download_manifest.json` at `3179342a9b559af06beabdfea98890c627a97e12`. All referenced files were Git-blob-hash verified. Locally, the script reused the previously verified MG014 input cache. To reproduce independently from the repository root:

```text
python review/MG019_replication/download_inputs.py
python review/MG019_replication/static_checks.py
```

Requires Python, NumPy, SciPy, Numba and mpmath. Inputs and caches are ignored. The script uses the inspected MG4 graph builder and detector to place the four stored particles; it independently evaluates the proposed coupling functions, clipping, count-response arithmetic, global-baseline correction, and the Student-t tail probability for a statistic of 3 from four independent Gaussian run summaries. This last calculation illustrates a statistical assumption; it does not establish that real run summaries are Gaussian or identically distributed.

The design-level equal-area quadrature is used only on the initial near-uniform point sets. It is explicitly not offered as valid integration over the later spatially biased MG6 node set. No Voronoi-area reconstruction, MG6 trajectory, event stream or particle outcome was generated. No upstream file was modified or published.

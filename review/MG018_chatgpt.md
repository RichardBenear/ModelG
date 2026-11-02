MG018 — MG5 results review

Reviewed 2026-11-01 against repository commit `aa237e0f892cc4709ecc87cf72cd483f8ec8d536`, including the MG5 result entry in GATES.md, evaluation.json, the registered proposal and execution amendment, simulator, driver, evaluators, and committed raw records.

**Recommendation: accept Gate MG5 as PASS and Part B as A1, within the registered scope. I found no blocking error.** The reported PINNED / PINNED / MOBILE rate labels are also correct under the registered 1.5 displacement threshold. I support closing MG5, subject to the PI's sign-off, and moving to MG6 design review. MG5 establishes a useful baseline; it does not establish directed response to a G gradient.

I downloaded and Git-blob-hash verified 119 pinned files, ran the repository's independent recomputation, and wrote a separate analysis that does not import either evaluator or the simulator. That analysis re-detects particles from the saved smoothed values, reconstructs tracks, calculates MSD fits and gate outcomes, and checks the accounting records. Scripts, input manifest and numerical outputs are in [review/MG018_replication/](https://github.com/RichardBenear/ModelG/tree/main/review/MG018_replication).

**Results independently reproduced**

All 24 runs have COMPLETE records, reach their registered horizons, and have empty flag lists. All check times, node-count histories and recorded classification sequences pass the consistency checks.

| Part A rate | Single | Pinned / mobile | Range of maximum displacement | Rate label |
|---|---:|---:|---:|---|
| 0.002 | 4/4 | 4 / 0 | 0.354900–1.370451 | PINNED |
| 0.02 | 4/4 | 3 / 1 | 1.064708–1.522829 | PINNED |
| 0.2 | 4/4 | 1 / 3 | 1.142623–2.851075 | MOBILE |

The four Part A frozen controls remain single; their largest maximum displacement is 1.30934e-13, comfortably below 0.1. Thus the controls are valid, every rate survives in 4/4 runs, and **PASS follows without ambiguity**.

The independently reconstructed single-particle tracks agree with evaluation.json within 1.78e-15 per coordinate. The D_p estimates agree within 2.17e-19. Their ranges are 1.18041e-4–5.83502e-4, 6.80571e-4–1.04857e-3, and 7.87042e-4–2.26209e-3 for the three rates, matching the report. Counts, depths and multiplied subtypes also agree with the driver records.

For Part B, all four frozen controls multiply and all four fast-rewriting runs remain single: **A1 is correct**. The frozen cases first have multiple particles at t=80, 70, 50 and 30 respectively, confirmed at the following checks; all four are classified as split. This is a clear contrast within the registered benchmark. The appropriate interpretation remains recovery under rewriting plus interpolating transfer; these runs do not separate connectivity changes from interpolation smoothing.

**Numerical checks and provenance**

From the committed Part 0 records, I independently reproduce:

| Comparison | Largest centre difference | Largest relative depth difference | Registered limit |
|---|---:|---:|---|
| dt=0.05 versus 0.025 | 0.000563352 | 0.092363% | 0.1 and 2% |
| Delta=0.2 versus 0.1 | 0.0602506 | 1.341700% | 0.1 and 2% |

Both comparisons have one particle at every check and pass. The reported frozen ROS2/BDF field difference, 1.52085e-4 against a 1e-3 limit, is consistent with decisions.json; I did not rerun that integration or independently reproduce its full-field difference.

Across **19,984 edited production ticks**, recomputing I_after minus (I_before minus R plus B) reproduces the saved normalized jump residuals. The maximum is 3.29878e-13, below 1e-10. The largest recorded interpolation-audit discrepancy is 4.83641e-15, below 1e-9. Effective/proposed birth and death counts, ephemeral-node cancellation, and resulting populations all reconcile. These validate the committed accounting; checking an identity from saved sums does not independently verify the underlying fields.

I also regenerated **all 16 event lists**, containing **25,856,362 events**, using the pinned registered generator and separately written canonical hashing. Every SHA-256 matches both published replay hashes. Independent aggregation reproduces every saved tick count, and every saved sparse node has a valid lifetime; newborn positions match exactly. The generator/simulator source blob is identical at both production code revisions. The driver changes between them concern execution order, thread settings and the approved cutoff, consistent with MG017. Final-boundary events are deliberately unapplied under the pre-run specification; this is not missing work.

The ledger contains 24 unique finished runs and sums to **39.9630806 hours**, matching evaluation.json and the group totals. Keep the disclosed approximately 65 lost-launch hours visible: finished plus lost-launch time is approximately **105 hours**, before Part 0 and setup/replay work outside the run timer. The 39.96 figure is the ledger total, not the project's entire computational expenditure.

**Nonblocking corrections and limits**

1. **Strengthen recompute_mg5.py's failure reporting.** It currently prints DIFFERS and exits successfully on a mismatch. Return a nonzero status when differences occur. Its checks also cover less than the docstring suggests: it compares final classes and single-run mobility/displacement, but does not independently recompute D_p or aggregate the gate and Part B outcome. Either extend those checks or narrow the description. My separate analysis checked those additional quantities, so this does not invalidate the present result.
2. **Treat mobility as a registered descriptive label.** The intermediate-rate mobile run reaches 1.522829, only 0.022829 above the boundary. The fast-rate run at 1.597596 is also relatively close. Part 0 passes, but one graph's refinement differences are not a universal error bound for all production displacements. Keep the existing labels; do not interpret them as a sharply established physical transition or a generally monotonic mobility law. No post-hoc threshold change is warranted.
3. **Survival does not mean invariant shape.** From the saved depths, the fast Part A runs' depth magnitudes relative to their matched frozen controls range from about 0.776 to 1.262 over the checks. The particles survive while their profiles respond to rewriting and transfer. These are descriptive changes, not failures of the survival gate.

This review reproduces measurement and scoring from committed smoothed records and reproduces the event histories. It does **not** independently re-integrate the evolving-graph ODE, reconstruct smoothing from full fields, or verify the PC-only field files against their manifest hashes. Numerical assurances therefore still depend on the inspected implementation, recorded audits and registered Part 0 validation. Those are acceptable residual limits for this gate.

The report correctly limits its claims to this rule, transfer method, two-graph panel, rates and horizons. It establishes neither Brownian motion nor survival under token noise, other spatial rules, or geometry-changing rewrites. Within that scope, the evidence supports closure of MG5; none of the remaining comments requires another production campaign.

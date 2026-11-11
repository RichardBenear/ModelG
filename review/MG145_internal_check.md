# MG145 internal design challenge

- A single path at each tick does not estimate sampling uncertainty. Use independent path replicates; spatial sites are correlated and cannot replace them.
- Passing a finer-tick comparison is not proof of a continuum solution. Finest tick is an operational reference, and reference noise must be included in each interval.
- A significant nonzero difference may be inside an equivalence margin. Conversely a nonsignificant difference may have an interval too wide for equivalence. Classification must use the full interval against the margin.
- The new 1% contrast margin is relative to a reference concentration, not the prepared pattern contrast. The design documents the especially weak G sensitivity explicitly. It is a preliminary kinetics screen, not particle-shape preservation. A stricter requirement needs a revised precision/budget analysis.
- High-band power includes deterministic changes from the initial state and stochastic fluctuations. Its proposed margin is 0.05 in a fixed Poisson-reference scale; it is not 5% of observed power and does not assume Poisson dynamics.
- Transport-only planning omits reactions and catalyst competition. Exact linear contrast variances and Gaussian power approximations cannot establish coupled power. A disjoint 24-path precision/cost pilot should precede the final sample count.
- Bonferroni controls multiplicity only to the extent the individual approximate Welch intervals have their stated coverage. Distribution diagnostics are disclosed and must not trigger undeclared exclusions or method changes.
- Raw channel histograms reconstruct available/final counts but omit order. Retain per-tick states and diagnostics, RNG provenance, and full ordered sentinel traces; this changed retention policy needs explicit PI approval. No historical data are changed.
- Full event tracing across 144 paths would be several GB. Saved-pilot compression tests are a layout proof, not a guaranteed writer/storage forecast. Measure the new writer before registration.
- R=48 is provisional. Its twice-body cost estimate fits the three-hour/eight-worker preference, but real finer-tick IO, audit and physical wall deadlines must be measured. No existing budget approval transfers to this new study automatically.

Disposition: design draft ready for PI review; implementation and a concrete precision-pilot readiness packet are the next prerequisites. No new coupled paths or formal registration exist.

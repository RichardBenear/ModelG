# MG144 — Transport tick-resolution results

The exact transport-mean analysis confirms that tick size matters on the actual Model G substrates, even when the engine implements its rules correctly. At the two ticks used in the multi-tick pilot, the largest transport-response discrepancies over the tested horizons were approximately 10.3% and 5.29% when measured against unit-sized initial disturbances. Finer ticks reduced these discrepancies substantially. Actual smooth starting profiles had much smaller errors, so these are worst-case response limits, not claims that the pilot fields were10%wrong. The result supports choosing finer resolutions for the next coupled-model study while leaving reaction and particle conclusions open.

## Scope and outcome

**Deterministic transport-only development complete; all numerical/control checks passed.** The PI closed the multi-tick pilot as MG143 and directed the next-stage plan. Stage1 examined 20 operators: G,X,Y,A,B transport on four actual 500-site ring inputs covering both A-network origins and B-network intervention arms. Four tick sizes and two horizons (0.01, 0.1) were tested, giving 160 operator/time/tick conditions. Initial probes included uniform and point distributions, multiple sine/cosine spatial modes, and normalized actual species profiles. No random draws or coupled trajectories were run.

This work was bounded at 60 seconds on one local worker, with one numerical-library thread and completed in **3.5247seconds**. Source/specification was published as `97f777f` before calculations. Earlier qualified engine behavior, frozen results and thresholds were not changed.

## Main comparison

The continuous reference is the transport semigroup exp(LT). The strict-tick reference is P(delta) raised to T/delta. Each column of P contains its source token's stay probability exp(-lambda*delta) and one-move destination probabilities (1-exp(-lambda*delta))*rate/lambda. This exactly preserves the one-use-per-token rule for transport-only expected counts. A token may move again on a later tick.

The worst-case column below is the operator2-norm difference, maximized across the tested substrates/species/horizons: an L2error for any initial perturbation with unit L2norm. It is not a percent of the decayed final solution. The actual-profile column normalizes each starting species profile to unit L2norm before comparison. A 1% line is only a prospectively documented planning screen, not a PI-approved accuracy gate or a particle criterion.

| Tick | Worst-case transport error | Largest actual-profile error | Relative to 1% planning screen |
| --- | --- | --- | --- |
| 0.001 | 10.2758% | 0.445454% | Above |
| 0.0005 | 5.29089% | 0.23054% | Above |
| 0.0001 | 1.07793% | 0.0472704% | Above |
| 2e-05 | 0.216449% | 0.00949819% | Below |

Of the tested grid, tick **0.00002** alone fell below the 1% worst-case screen across all conditions. Tick 0.0001 was close but above it, at 1.07793%. This is evidence about the selected pure-transport operators and horizons, not a proof that 0.00002 suffices for nonlinear reactions or longer times.

| Species | Worst error at tick 0.001 | Worst error at tick 0.00002 |
| --- | --- | --- |
| G | 0.955694% | 0.0184078% |
| X | 0.955694% | 0.0184078% |
| Y | 10.2758% | 0.216449% |
| A | 0.00819509% | 0.000163749% |
| B | 0.00819509% | 0.000163749% |

The fast Y transport is the main limiting operator. Its two-link exit rate is 600, while the existing G/X diffusion coefficients imply slower exits. Network transport of A/B was also explicitly included; equal initial weighted degrees do not imply identical endpoint evolution for different connectivity.

## Validation and separate recomputation

For each operator, symmetry and constant exit rate were verified before using a shared eigensystem for L and P. Generator conservation, nonnegative transition probabilities, stochastic columns, uniform stationarity and zero-exit behavior passed. Eigenvector/residual checks passed tolerance 1e-9. Continuous means at0.01 matched independently computed sparse matrix-exponential actions; coarse tick powers matched direct sparse iteration, within 1e-9 relative/absolute.

The [separate recomputation](separate_recomputation.json) reconstructed each transition matrix column with the general source-wise formula, independent of the regular-degree shortcut. It recalculated powers with log1p/exponentials and reconstructed profile errors from saved outputs; all 160 conditions agreed. Shared generator/eigenvalue arrays and same-author implementation are disclosed; this is not external replication or exact rational arithmetic. Full generators,eigensystems,initial probes,transition matrices and mean outputs are retained for reproduction.

## Interpretation for the project

The earlier sparse and multi-tick qualifications established faithful implementation and numerical integrity of the selected strict-tick rules. They did not establish that an arbitrarily chosen tick approximates continuous Model G kinetics on the actual high-rate ring. This analysis separates those questions. It also explains why actual smooth profiles can look comparatively close while other spatial disturbances remain more sensitive; stochastic dynamics can excite additional modes later.

No multiple actions within a tick were introduced to obtain the finer results. No physical tick duration, Planck-scale mapping, vacuum-source identification, entanglement or persistent particle has been established.

## Next decision

Proceed to design a replicated coupled short-run study at finer ticks, provisionally 0.0001, 0.00005 and 0.00002. The added intermediate resolution is a prospective proposal, not a result of this grid. Use fresh independent paths and uncertainty bounds to distinguish ensemble changes from sampling noise. The finest-tick ensemble is a comparison reference, not an exact continuum truth; finite-count nonlinear mean corrections and initial rounding must remain explicit. Pilot the finest settings to measure cost before proposing the production budget.

A [Stage2design draft](../../../spec/MG12H3_tick_resolution_stage2_draft.md) records these dependencies. No Stage 2 production is registered or running. MG144 is a completed development analysis, ready for PI review; MG142 pilot closure is recorded, and MG139 closure remains unchanged.

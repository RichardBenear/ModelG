# MG014 round 2: MG5 proposal Draft 2

Reviewed 2026-10-31 at commit `a1c4b14c81cd5f8c72956b4ac7da9cc05153d491`.
Proposal blob: `470d5621e12fb689e1444d83d5011266afaf03d0`.
Read with `review/MG014_response.md` and the first-round review in this directory.

## Decision

**Draft 2 resolves the main scientific-design objections. I support implementation and Part 0 after two concrete corrections: the event-storage commitment and the disposition of interpolation-audit failures.** Production runs remain conditional on Part 0 and valid frozen controls. This is a proposal review, not certification of an implementation that does not yet exist.

The PI's choices of interpolating transfer, Part B, an unweighted operator, and keeping full fields off GitHub are accepted. No additional graph panel, noise arm or weighted comparison is requested.

## Closed first-round findings

1. **Rule and cost now agree in principle.** Section 3.2 explicitly defines synchronous updates with one rebuild per tick. This is a well-defined model, and tick refinement measures sensitivity to the definition. All proposed tick/step combinations, including fallback Delta=0.05, align.
2. **The baseline is repaired.** The original calibrated h values preserve the stored equilibria. Four matched Part A frozen controls test the new solver and longer horizon. Part B's potential shift, -0.9/14 in phi_Y only, correctly preserves absolute concentrations when b changes from 28 to 29.
3. **Numerical validation is substantially improved.** Frozen ROS2/BDF comparison, paired step refinement and paired tick refinement answer different questions. Identical continuous event histories and predetermined fallback/stop rules are appropriate. Acceptance remains conditional on correct implementation and successful checks.
4. **Recovery is appropriately limited.** A1 now concerns rewriting plus transfer; connectivity changes and interpolation smoothing remain unresolved explanations. MG5 need not separate them.
5. **Sparse detection records are useful.** Float64 smoothed values and positions below -1 suffice for the specified detector at thresholds <= -1, subject to the ordering and metadata details below. Full-field uploads are unnecessary for reclassification. Status precedence and jump accounting are welcome.

## Corrections before registration

### 1. Section 7 underestimates committed event-list storage

Expected births and deaths are each N0*nu*T. Across the 16 rewriting runs, the expectation is **25,856,568 events**. Frozen runs add none. An illustrative compact uncompressed binary schema uses 8 bytes for time, 1 for type, 8 for ID per event, and two 8-byte coordinates per birth: 25 bytes/event on average.

| Group | Runs | Expected events/run | Event bytes/run, decimal MB |
|---|---:|---:|---:|
| A, nu=0.002 | 4 | 37,802 | 0.945 |
| A, nu=0.02 | 4 | 378,020 | 9.451 |
| A, nu=0.2 | 4 | 3,780,200 | 94.505 |
| B, nu=0.2 | 4 | 2,268,120 | 56.703 |

Total: **646.414 MB** for event lists alone, excluding checks and accounting. Using 32-bit IDs reduces this example to **542.988 MB**. These are illustrative schemas, not measured compressed sizes or information-theoretic lower bounds. Compression or implicit fields may help, but a few MB/run has not been established. Batching reduces rebuild work, not the number of events in a full continuous-time list.

**Required correction:** specify and measure serialization/storage before promising small committed artifacts. One option is to keep full lists with the local checksum-manifested files and commit a canonical replay specification: seed, generator/bit-generator version, code revision, event schema/order and event-list hash. This requires explicitly revising the promise to publish event lists and a clean-process regeneration/hash check. Alternatively, retain publication and budget the measured compressed sizes. Do not silently commit hundreds of MB under the existing estimate.

Sparse checks can support reclassification without event lists. Re-running dynamics is a separate claim requiring replay information and pinned initial states. This finding respects the PI's small-repository storage decision.

### 2. Give interpolation-audit flags a numerical-validity disposition

Section 7 flags an independent-interpolator discrepancy above 1e-9, but COMPLETE and FAILED-NUMERICAL explicitly mention unresolved balance and jump-accounting flags only. A wrong spatial assignment of birth values can satisfy the node-sum identity perfectly: the same wrong values enter B and I_after. Jump accounting cannot resolve an interpolation discrepancy.

**Required correction:** include unresolved interpolation-audit failures in FAILED-NUMERICAL and prohibit scoring until resolved through a recorded investigation/recheck. Define the resolution procedure and retain the original flag and its resolution in the committed run record. This closes a validity gap without expanding the study.

## Details to freeze in the promised pre-run implementation specification

- **Boundary timeline:** integrate on the existing mesh to t_k; measure a scheduled check before edits; interpolate births from that boundary field on the pre-edit mesh; apply the batch; rebuild; refresh the next tick's Jacobian. Section 3.2's parenthetical says t_k before edits, but "tick-start field" can suggest t_(k-1). Distinguish these states explicitly.
- **Ephemeral nodes:** either B and R both exclude same-tick birth-deaths, or both include their equal interpolated contributions so they cancel. Mixing conventions fails accounting. Distinguish proposed events from effective edits.
- **Ordering and initial centre:** MG4's exact-depth tie break uses current row index. Preserve relative row ordering in sparse records, or save current indices and reconstruct it. Persistent IDs alone need not specify storage order after deletions. Pin deterministic tie/summation rules and include the t=0 centre or a fully reproducible pinned detection record, torus lengths and thresholds. Mathematical sufficiency does not imply bitwise-identical summation after reordering.
- **Part 0 wording:** identical counts at every check already force identical first qualifying multiplied/dissolved times. The extra allowance of one check's timing lag is redundant under the literal stricter rule. Apply strict count agreement unless an explicit exception is registered before checking outcomes.
- **Timing:** restate the budget using the complete implementation, including independent interpolation audits, event generation and serialization. The old solver/rebuild pilot does not measure all this work.

These are implementation clarifications, not requests for another scientific arm.

## Interpretation notes, not further blockers

MOBILE demonstrates motion under field-independent rewriting plus transfer. It makes directed gradient response worth testing in MG6; it does not establish that response. Revise the section 8 implication accordingly. D_p remains descriptive over the chosen lags, rather than proof of Brownian diffusion.

A3 can include four dissolved fast runs; "other outcome" is more accurate than "mixed". Report the full class vector. The spatial snapshots remain in the same conditional uniform-point/Delaunay family, but the population-count distribution evolves from fixed N0. "Not different graph statistics" should not imply equality of every finite-N distribution conditional on the particular starting graph.

## Checks and limits

`draft2_checks.py` and its JSON output contain workload/storage arithmetic, all step/tick alignments, an ephemeral-node bookkeeping example, and sparse/full detection comparisons. On the four stored b=28 states, **148, 146, 153 and 164** nodes survive the -1 cut out of 18,901. Counts, selected minima, depths and centres match full-record detection. Nine synthetic comparisons cover torus wrapping, exact-depth ties, multiple particles, no candidate and empty sparse records at thresholds -4.457582, -5.117354 and -1. All passed; observed centre differences were zero, with an asserted tolerance of 1e-12.

These checks exercise the existing inspected MG4 detector, not an independent implementation, and start from already-smoothed values. They do not certify future MG5 smoothing, transfer or ROS2 convergence. The independent recomputation implementation remains necessary. **No MG5 trajectory, event stream or outcome was generated.**

The original `design_checks.py` was rerun successfully. Its h=0.2 and Delta=0.25 findings are historical Draft 1 results, resolved by Draft 2. Both versions' provenance is retained. The two Draft 2 documents were fetched through the GitHub connector and verified against their Git blob hashes; a direct Python download attempt failed because shell DNS/network access was unavailable.

**Disposition:** correct the storage commitment and interpolation-audit failure policy, implement and freeze the listed details, then run Part 0 before production. No remaining objection requires reopening the PI's physical-model choices.

# MG4c result review (MG012)

Reviewed 2026-10-31 at commit `37c9c34b6574c77b42d1a2d1e65957754b518d7a`.
Design: MG4c draft 3 frozen at `a946550` (proposal blob `cee43bed8f7cbae7ad728e9d5dcc473e4227c309`), with the dated pre-run implementation specification in GATES.md. The current proposal body agrees with the frozen body; only its registration status header changed.

## Decision

**W1 is independently reproduced, and the saved weighted run has the single-particle outcome reported as O1.** The main scientific conclusions are supported. I find no class, threshold, refinement, shape or balance error that changes the observed outcomes.

There is one small but real tolerance convention/deviation to resolve in the formal Part 3 record: the proposal specifies entrywise agreement to 1e-12, while the implementation applies a scaled relative tolerance. The measured matrix discrepancy misses a literal absolute 1e-12 check, although the resulting T trajectories agree exceptionally well. Record the numerical finding and explicit acceptance/interpretation rather than silently claiming all stated checks passed. I do not recommend rerunning the scientific experiment for this issue.

## Independent results

All 102 saved runs and all 1,122 saved observations were checked. Classification, centroid/linking, Gaussian smoothing, radial profiles, wavelengths, ODE derivatives, balance sums and bisection reconstruction use separately written calculations. Registered builders supply the graph positions/edges and seeds. This is saved-data replication, not a second 11-hour integration campaign.

| Quantity | Independent result |
|---|---|
| New unweighted D runs | 72 single, 3 multiplied, 5 decayed, 0 unresolved |
| All-four-single b values | 13 tested values, from 27.265625 to 28.7265625 |
| Lower panel boundary bracket | [27.26171875, 27.265625] |
| Upper panel boundary bracket | [28.7265625, 28.73046875] |
| T-reference eligibility | 13 eligible and 7 invalid among 20 new unweighted b values |
| Adaptive versus fixed-threshold D classes | No disagreements |
| Weighted D0 p0 at b=29 | Single at all 11 checks |
| Weighted T control | Single; max final-field difference 4.547473508864641e-13 |
| New balance failures | None |
| Completed-run ledger | 102 unique records; 11.326969444444444 core-hours |

Both six-step refinement sequences reproduce the driver's visits and final brackets exactly. W1 therefore follows the registered definition. The lower boundary is set by D1 p0; the upper by D1 p1. These are boundaries of the four-start benchmark predicate, not independently estimated boundaries of every graph or placement.

Independent smoothed final minima agree with the records within 2.31e-14. Stored float32 smoothed checks differ from re-smoothing the float64 final states by at most 4.77e-7, consistent with quantization. Independent final max-RHS values agree within 5.69e-14, and every stationary flag agrees. The largest absolute independently computed r2 is 2.87e-11; every residual satisfies its registered mixed tolerance.

The same-b shape comparisons reproduce within 1.78e-15 in ratio units. Smoothed depth ratios D/T range from 0.903578 to 1.139522; wavelength ratios range from 0.961702 to 1.012172. Missing T references are correctly left unavailable. These are descriptive comparisons, not an additional C2 pass/fail score.

## Weighted operator: independently checked

I reconstructed periodic faces by central-vertex incidence and physical-node-triple deduplication, then assembled stiffness matrices from P1 basis gradients. This differs from the source's centroid-in-tile selection and angle-cotangent assembly.

Both reconstructions have exactly 2N faces, the expected graph edges with two face incidences each, positive node areas, total area 400.0806425376431, and symmetric stiffness matrices. Independent and source area vectors agree within 1.95e-16. On D0, the maximum operator-entry discrepancy is 2.62e-10 against entries as large as 1.61e5, about 1.6e-15 relative. The independent operator reproduces the weighted run's max RHS within 5.69e-14 and passes its area-weighted balances.

### T construction tolerance needs an explicit numerical resolution

The frozen proposal's Part 3 requirement reads: L_w equals MG4's L_G 'entrywise to 1e-12.' In `src/run_mg4c.py`, `part3_checks`, the implemented comparison is:

    max_abs(L_w - L_G) < 1e-12 * max_abs(L_G).

On the full T graph, I measure:

| Check | Value |
|---|---:|
| max_abs(L_w - L_G), using the driver's reference | 3.495870259939693e-12 |
| Literal absolute proposal limit | 1e-12 |
| Implemented scaled limit | 1.4999999999999997e-10 |
| Difference using MG4's calibrated h instead | 3.581135388230905e-12 |
| Final T-field difference | 4.547473508864641e-13 |
| Registered final-field limit | 1e-6 |

Thus `equals_LG: true` in part3_checks.json means the scaled implementation check passed; it does not demonstrate an absolute entrywise 1e-12 bound. If relative agreement was intended, that convention was not explicit in the frozen proposal. On the literal absolute reading, the registered check is missed by about a factor of 3.5.

This is consistent with tiny floating-point geometry/assembly differences, not evidence of changed diffusion physics: the entry error is only about 2.3e-14 of the largest T entry, the independent FEM reconstruction agrees, and the actual T control field matches far inside its registered tolerance. A documented numerical adjudication is sufficient evidence to support acceptance; no new trajectory campaign is indicated. Preserve the measured values, the original boolean and the actual tolerance convention, and record the PI's acceptance of the numerical deviation/interpretation. Do not retrospectively rewrite the frozen tolerance or claim the stricter absolute test passed.

The unweighted W1 finding is independent of this Part 3 construction issue. An unqualified, fully compliant O1 sign-off should follow the explicit numerical resolution rather than precede it.

## Interpretation: survival, stationarity and boundaries

The central finding is supported: unweighted irregular D graphs can retain single particles at lower b under this start and horizon; for D0 p0 at b=29, switching the diffusion discretization changes the observed outcome from multiplied to single. This is a useful follow-up to MG4 and does not re-score MG4.

Keep the following qualifications visible:

- **The lower boundary is operational and unsettled.** At b=27.26171875, D1 p0 has one detected particle through t=280 and none at 290 or 300. Its final smoothed minimum is -4.89494 against threshold -5.11735; max absolute RHS is 0.0488534. At the adjacent single point b=27.265625, max RHS is still 0.0136884. This is not a precisely located stationary saddle-node/existence boundary.
- **The upper single endpoint is not stationary either.** D1 p1 at b=28.7265625 is single but has max RHS 2.93618e-4. At b=28.73046875 the detector counts [1,1,2,2,3,3,5,5,6,6,10]. This clearly supports the finite-horizon classification boundary, without measuring an asymptotic replication bifurcation.
- **The weighted D run is not stationary under the registered cutoff.** Its max RHS is 8.27412e-6, above 1e-8. O1 demonstrates single-particle survival through the observation interval, not a certified stationary weighted particle.
- Of the 72 new unweighted D singles, **63 meet the stationary cutoff and 9 do not**. This is reassuring evidence for many settled outcomes away from troublesome starts/edges, while preserving the survival-versus-stationarity distinction.
- All 80 new D classifications agree under the adaptive and fixed MG4 thresholds. That removes a plausible concern that the observed W1 result was manufactured by the adaptive threshold. The lower edge still depends on the fixed fallback detector and finite horizon.
- The upper benchmark boundary is about 0.47 below T's reported upper boundary; this is not evidence for a rigid translation of an entire graph-independent existence window. The four-start panel and fixed initial condition remain the estimand.

These qualifications limit interpretation; they are not reasons to add unregistered extensions to the current study. No further runs are required to establish W1 or the observed weighted single outcome. An asymptotic boundary study would be a new registration if the PI wants one.

## Procedure, validation and reproducibility

All 102 retained jobs and ledger rows are present, below the launch cutoff. The disclosed aborted cloud T trial before the implementation specification was committed remains a procedural deviation and should stay visible. The 11.33 figure is the completed PC-run integration ledger, rather than a claim that no earlier MG4c dynamics were attempted. I found no evidence that the aborted trial changed criteria or selected the reported D outcomes.

All 11 upstream MG4 tests pass. Ten applicable MG4c tests pass, including the small-graph integration and true higher-precision recheck; two fork-dependent synthetic pipeline tests explicitly skip Windows. The independent all-field reconstruction covers the actual Case A result, but I do not claim to have executed those skipped Linux-only pipeline tests. A workspace-local temporary-directory wrapper was needed because the initial direct invocation encountered sandbox temporary-directory permissions; no upstream test was edited.

Files: replicate.py/replication.json, audit_tolerance.py/tolerance_audit.json, verify_tests.py/test_results.json and download_manifest.json/verified_inputs.json. Input downloads are pinned and hash-verified. No registered full-size trajectory was rerun. No upstream proposal, source, result or gate file was changed. This review has not been committed or pushed.

**Recommendation:** accept W1 and the scientific weighted single-particle finding. Explicitly resolve the T tolerance convention/deviation before formal O1 closure, retaining the stationarity and boundary qualifications above. No further physics run is requested.

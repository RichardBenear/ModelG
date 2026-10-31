# MG4 result review (MG010)

Date: 2026-10-31. Reviewed commit: `41ee3de73d30e80fc0aeb95350dc00070c1c6f58`.
Design: MG4 proposal draft 4, registered at `773255d`, with the dated PI scope and budget decisions in GATES.md.

## Decision

I accept **FAIL as the scientific result on MG4's registered finite panels**, and accept the balance flags as numerical roundoff rather than evidence of a non-conservative implementation. C1 alone decisively fails for both irregular families. No additional trajectory campaign or stronger statistical assumption is needed.

Before closing the result record, reconcile the executable/machine-readable verdict with the adjudicated verdict and qualify the causal interpretation. These are record and interpretation corrections, not requests for changed thresholds or new dynamics runs.

## Independent replication

The replication reconstructs registered graph geometry and uses separately written measurements. It does not import the project's detector, classifier, smoother, or C3(b) evaluator. All 80 saved runs were checked at all 11 times, t=200 through t=300, including periodic centroid links for single-particle cases.

| Family | Single runs | Outcome | Final counts |
|---|---:|---|---|
| D, h about 0.2 | 0/32; 0/8 on every graph | All multiplied; C1 FAIL | 33–37, median 35 |
| R, h about 0.2 | 0/32; 0/8 on every graph | All multiplied; C1 FAIL | 33–37, median 35 |
| T, h=0.2 | 4/4 | All single | 1 |
| D, h about 0.3 | 0/8 | All multiplied; descriptive | 34–37, median 35 |
| T, h=0.3 | 4/4 | All single; descriptive | 1 |

There are **zero class disagreements** with evaluation.json. Counts are constant over the 11 saved checks in every run, extending the narrative's three spot checks. This establishes the late-time outcome, not the earlier multiplication mechanism. T's largest centre displacement between adjacent checks is 2.83e-9 at h=0.2, well below the 1.0 link threshold.

C2 is correctly INCONCLUSIVE for D and R: no single-particle outcomes are available for the registered comparison.

Independent C3(b) median errors:

| Family | G | X | Y | Registered result |
|---|---:|---:|---:|---|
| D | 0.2724694007 | 0.3096229388 | 0.3039293508 | FAIL |
| R | 0.5117752967 | 0.5243072235 | 0.5208387317 | FAIL |
| T | 0.0071309773 | 0.0096566728 | 0.0090879754 | PASS |

The largest discrepancy against the 204 published smoothed component errors is 6.67e-16. This independently checks the registered interpolation, operator scaling, smoothing and norm calculation.

Analytic triangular dispersion reproduces T's C3(a) mean, 1.0079841829621143. Re-averaging saved spectral ratios gives D 1.0514143721771592, R 1.0635481904563573, and N3 0.00766336730899827. D/R FAIL and the negative-control classification follow the registered cutoffs. D exceeds the upper bound by only 0.001414 in ratio units, but remains a registered failure. Full irregular-graph eigenproblems were not independently rerun.

The saved float64 Part 0 T reference independently gives minimum X potential -13.583783934082076 and max absolute RHS 6.707523425575346e-12. All 10 upstream MG4 tests pass locally, including the analytic Jacobian and solver-order equivalence checks.

## Finding 1: reconcile the official verdict records

GATES.md:858 declares FAIL and cites evaluation.json. That JSON still says `NOT SCORED (implementation failure: investigate balance residuals)`. The evaluator, src/evaluate_mg4.py:135–136, still unconditionally selects that status when the original balance flags are present. The committed evaluator therefore does not reproduce the narrative's final verdict.

Preserve original run flags and raw evaluation. Publish a distinct, versioned adjudication recording all 11 resolved flags, evidence and precision, both PI budget exceptions, and final FAIL under unchanged precedence rules. Either have the evaluator consume this explicit adjudication or distinguish raw evaluation from an authoritative final evaluation. Do not erase residuals or simply loosen the registered tolerance.

## Finding 2: balance investigation accepted; complete and qualify its record

Symmetry of the graph operator gives zero total diffusion, and the nonlinear reaction terms cancel from the total-species sum. Algebraically:

    dI_G/dt = -I_G + 0.1 I_X
    d(I_G+I_X+I_Y)/dt = -0.9 I_X.

These identities hold at nonstationary states too. The implemented equations match them.

The published balance_investigation.txt contains only 8 of the 11 flagged runs. It omits **R0_p4_sol, R0_p7_sol and R2_p4_sol**, despite GATES saying every flagged run was checked. This replication now checks all 11. Rebuilding edge fluxes and reactions in genuine 40-digit Decimal arithmetic on the saved final states gives max |r1| below 1.5e-37 and max |r2| below 2.0e-35. The same saved states evaluated with the original float64 RHS exhibit r2 as large as 3.39e-10. This strongly supports roundoff in RHS assembly/aggregation rather than a missing conservation term.

Use the complete evidence in the final record. Replace 'the implementation conserves exactly' with an algebraic-conservation statement and the roundoff evidence. Errors can arise while forming local diffusion/reaction derivatives before the final sum; epsilon times the already-cancelled sum of absolute RHS values does not bound the whole computation.

Portability: `np.longdouble` has **53 mantissa bits on this Windows installation**, the same as float64. The existing script documents 80-bit Linux use, so that Linux evidence need not be wrong, but a Windows rerun does not automatically provide extra precision. Record actual precision or use a genuinely higher-precision method, such as this Decimal calculation.

Scope: committed Part A final fields are float32. Recomputing their balances validates the identities at those rounded states and supports the roundoff diagnosis. It does not recover the original float64 state, independently certify original stationarity flags, or bound accumulated trajectory error. This does not undermine the large C1 failure. The separate Part 0 T state is float64 and supports the stationarity check above.

## Finding 3: narrow the causal interpretation

GATES.md:863 and 883–887 go beyond the evidence when asserting that the particle copies itself until the torus is full, that irregular graphs move it past its replication edge, and that excess spectral counts mean slower effective diffusion at particle scale.

The data establish that a continuum single-particle initial condition produces many detected localized minima on these fixed graphs at b=29. A shifted replication boundary is plausible, but was not measured on D/R. Late-time snapshots do not distinguish fission, nearby nucleation, basin changes, or other transient routes. A larger integrated spectral count establishes a difference in the density of states; it does not by itself identify relevant eigenfunctions, a scalar particle-scale diffusivity, or the nonlinear mechanism. The proposal correctly lists such readings as hypotheses. Restore that qualification in the result narrative.

The acknowledged C3(b) limitation is real: applying the operator to a bare sampled profile differs from evolving a solution that develops a microscopic homogenization corrector. Its failure does not establish absence of a homogenized continuum limit. The second-moment scaling is a useful heuristic, not a general impossibility theorem for passing C3(b). Keep its registered score and interpret it narrowly. C1 independently fails without relying on this design issue.

## Deviations, scope and closure

The PI explicitly accepted both budget overruns in the dated record. I treat these as disclosed exceptions without reopening approval. The pre-run reduction to gate-deciding survival runs preserves the scored questions; deferred seeded formation and dropped N3 dynamics must remain labelled unperformed. All 80 retained runs are present.

This refutes the registered prediction for the unweighted frozen D/R panels at b=29, these spacings, this initial condition and this horizon. It does not show that Model G lacks localized solutions on irregular graphs, that no lower-b single-particle window exists, or that weighted constructions fail.

**Recommendation:** close MG4 as FAIL after correcting the verdict record and interpretation. No new experiment is needed to confirm or change this conclusion. A D/R existence-window study or weighted/Voronoi study would be a new, separately registered investigation.

## Reproducibility limits

See replicate.py, replication.json, download_manifest.json and README.md. Inputs are pinned by Git blob hashes. Independently written calculations reuse the registered graph builders; they are not an independent implementation of the whole experiment. No full Part A trajectory or full D/R/N3 spectrum was recomputed. No upstream source, result, or gate file was edited by this review.

# H3b isolated-control tolerance diagnostic — MG118

PI: “Go ahead with the tighter tolerance check.” New separate diagnostic, not a retry/rescore of MG116/MG117. Preserve all original sources, failed/validated trajectories and the unstarted job. The original study remains failed/incomplete regardless of this diagnostic outcome.

## Frozen question and method

Does tighter numerical integration remove the spent-Z discrepancy in the **same isolated seeded control**, with unchanged equations, starting state, coefficients, time horizon and original error threshold?

Use canonical ISOLATED_C1_s1_v0 inputs, Omega=1 concentration adapter, N500, original fixed A/B references, zero resource hopping, same active ring, seven finite species and no external imports. Observe T0..25 at quarter-unit intervals. Three fixed integrations, no adaptive selection or retries:

| Label | Solver/equations | Relative tolerance | Absolute tolerance | Maximum step |
|---|---|---:|---:|---:|
| bdf_tight | incidence-derived nine-field population/import RHS | 1e-11 | 1e-13 | .05 |
| bdf_tighter | same incidence-derived RHS | 1e-12 | 1e-14 | .025 |
| radau_reference | separately assembled seven-species reaction/diffusion RHS and Jacobian structure | 1e-12 | 1e-14 | .025 |

The Radau path uses frozen mg11g2_pattern_v2_oracle population_rhs/population_pattern, not the new incidence RHS; append zero imports only when saving/checking the closed trajectory. Shared physical model, initial fields and SciPy stack remain disclosed; this is same-author separate solver/equation assembly, not external replication. Compare RHSs at three saved original states(T0,12.5,25) within rtol1e-12/atol1e-9 before integrating.

## Criteria, fixed before new trajectories

All three must independently satisfy the **unchanged** original max-over-species/site/time scaled field error<=1e-5 against the same closed archived fine reference. Scale=max(1,initial species maximum), including divisor1 for initially absent Z/W. No final-Z rescaling or excluding spent products. Each must have finite fields, minimum>=-1e-9 and relative total population/import balance<=1e-10. Both BDF paths must agree with Radau and with one another within the same1e-5scaled bound. All101frames and metadata required; any missing/prefix file blocks a supported diagnostic verdict.

A supported result would support a numerical-tolerance explanation in this one case. It would not make the original batch complete, retrospectively qualify its coarse trajectory or prove a general cause for every possible mismatch. Report actual errors and comparison with original coarse failed output, with no new stochastic/particle/quantum claims.

## Operations and provenance

Commit/push protocol, source and original/reference hashes before numerical integration. No RNG/seeds applicable. One local worker/numerical-library-thread1,300elapsed-second aggregate cap including prefix IO and saved-data audit; three integrations only, no automatic retry/extension. Given the original entire batch lasted about46seconds, this is a conservative bounded diagnostic cap rather than a new particle run budget. Freeze actual settings as above; a solver/domain failure is retained and stops new integrations. Fixed reference-accuracy results are scored after all three predeclared integrations, without selecting additional runs in response. Save quarter-unit prefixes and complete raw fields, comparison tables, source/metadata/budget/failures in a new public tolerance_check directory under the approved diagnostic workflow. Verify original hashes before/after and published bytes afterward. PI disposition of the original incomplete study remains separate.

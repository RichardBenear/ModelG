# MG4c Draft 3 review (MG011, round 3)

Reviewed 2026-10-31 at commit `e9c5e81cbc825512fe9cb1887d24929a862c5c38`.
Proposal blob: `cee43bed8f7cbae7ad728e9d5dcc473e4227c309`.
Read with review/MG011_response_r2.md and GATES.md. MG4c remains a draft, with no MG4c registration entry in GATES.md at this snapshot.

**Decision: I support registration of Draft 3. I find no outstanding design blocker and request no further proposal revision.**

All four round-2 concerns are resolved:

1. **T-reference eligibility is now absolute.** A T run must be single under the fixed -1.0 threshold and have final smoothed minimum <= -2.0 before supplying an adaptive threshold. The adaptive threshold is consequently <= -1.0, T is checked again under it, and an invalid reference triggers the fixed MG4 fallback. The scale-free eligibility error is removed. The local independent check rejects the constructed final-depth 1e-12 and 1e-30 remnants from round 2, while admitting the depth-10 example. These are constructed fields, not MG4c runs or stationarity tests.
2. **O2 is interpreted correctly.** It means no single-particle recovery, while multiplied and decayed remain separately reported. It no longer incorrectly equates decay with MG4's multiplied outcome.
3. **Refinement reporting is deterministic.** Case A refines only the selected highest all-single block. Other blocks are reported unrefined. A Case B discovery is W1 with edges explicitly not estimated, with no switch or added steps. T-relative shape comparisons are unavailable without a valid reference.
4. **The completion allowance is correctly an estimate.** Twenty core-hours is the launch cutoff; running jobs may finish above the projected allowance, with actual cost and any overrun disclosed. The run cap remains 102.

The study remains a finite-panel, finite-time test from the fixed b=29 initial profile, using the declared detector. A missed shallow or narrow single-particle region, basin dependence, and differences on other graphs remain interpretation limits, not errors that require another design cycle. W1 is not a general stationary-solution existence proof; W2 is not a nonexistence proof. Draft 3 is appropriately limited on these points.

The weighted diagnostic remains a deliberate change of diffusion discretization for D0 p0. Its positive control, periodic face and area checks, weighted conservation, analytic-Jacobian check and actual-operator precision recheck are the appropriate implementation safeguards. Approval of the design does not certify code that has not yet been written/tested.

**Next step:** PI approval, a frozen MG4c entry and the stated pre-run implementation details in GATES.md, then the specified validation checks and runs. From the review standpoint Claude can proceed with that process. No larger panel, R runs, noise experiment or extra statistical assumption is needed for this scope.

The reproducible checks are design_checks_r3.py and design_checks_r3.json. No MG4c dynamics were run, and no upstream proposal, source, result or gate file was changed by this review.

# MG028 independent review of MG9 results

**I recommend accepting P1 SUBCRITICAL and P2 UNRESOLVED.** The evidence establishes the local onset type of the selected Model G Turing branch. It leaves the proposed correspondence between front arrest and the single-particle window open. I found no error that changes either registered outcome.

Reviewed revision 6e3234b66d92e48e116cbc5e8c81211660a06ab7 on 2026-11-04 against registered Draft 3, amendment MG027 and the pre-run implementation specification. The production source differs from the reviewed pre-run revision f7ef9e9 only by the three-line evaluator change that records the platform. I performed independent algebra and records calculations, with no imports from the production calculation or scoring modules. I did not rerun the nonlinear ODE trajectories.

## P1

I reconstructed the linear operators and null vectors independently, computed the second-order forcing by Fourier phase quadrature, solved for its mean and second harmonic, and projected the resonant cubic forcing. The coefficients are

| Problem | Independent L₃ |
|---|---:|
| D02 mode 34, primary | +0.0136530521 |
| h = 0.2, continuous critical wavenumber | +0.0141620414 |
| Continuum, continuous critical wavenumber | +0.0131278830 |
| D01 mode 34 | +0.0126024098 |

These agree with the stored values. The tiny differences for the continuous-wavenumber cases arise from numerical optimization of the critical wavenumber and do not affect the sign.

From the saved amplitude series I independently refitted the integrated amplitude equation, retaining the recorded residual linear growth rate.

| Initial amplitude | Fitted L₃ | Change in 1/A² | Identifiable |
|---|---:|---:|---|
| 0.05 | 0.013652570 | 20.37% | No |
| 0.10 | 0.013650018 | 74.88% | Yes |
| 0.20 | 0.013640735 | 74.52% | Yes |

The tight-tolerance repeat at amplitude 0.10 gives 0.013650017. The registered rule requires at least two identifiable amplitudes, so the 0.05 outcome is handled correctly. The identifiable fits agree with the primary derived coefficient to better than 0.1%, and the repeat is well within its 2% limit. All four signs agree and the primary coefficient is far from the registered cancellation/degeneracy threshold.

The phrase “no cancellation” in GATES.md should read “no near-degeneracy.” There is ordinary partial cancellation between the two positive contributions and the negative cubic contribution. It causes no numerical ambiguity here.

Scientifically, this validates the local subcritical cubic amplitude equation near the specified onset. It supplies neither a global energy functional nor a stable saturated state from the cubic truncation alone. It also does not identify a Higgs field or give the amplitude potential physical energy units.

## P2

I recomputed the classifications of all 16 actual front runs from their saved interface tracks and completion/residual records. They agree with the result entry. The other 15 grid entries report preparation outcomes rather than front integrations.

| Grid outcome | Number |
|---|---:|
| Preparation lost, front unavailable | 3 |
| Preparation changed wavelength/cell count, front unavailable | 11 |
| Preparation unresolved | 1 |
| Front unresolved after extension | 3 |
| Invades | 13 |

At b = 29.25, 29.30 and 29.35, each front moved approximately 2.22–2.24 units over t = 100–1800. The wavelength is 2.941, giving an arrest displacement limit of 1.471 and an invasion displacement threshold of 2.941. The displacements fall between those limits. The late speeds are small, approximately 0.00036–0.00037, but the final residuals are approximately 0.0076 against an arrest tolerance of 0.00001. Thus both displacement and stationarity prevent ARRESTED. The extension was applied and UNRESOLVED-SLOW is correct.

With no arrested points and four unresolved grid points, the registered primary result is UNRESOLVED. Omitting edge refinement, endpoint robustness repeats and the fine-mesh comparison follows their conditional design. It is not evidence that a pinning interval is absent.

The upper transition from slow behavior at 29.35 to invasion at 29.40 contains the rechecked replication benchmark near 29.396. That is an interesting descriptive alignment, but its lower side is not an arrested state. The lower front-arrest boundary is unavailable because this preparation no longer supplies a stationary 34-cell pattern. Neither circumstance establishes that the two intervals coincide.

## Audit scope and disposition

I verified 110 saved state hashes and independently evaluated the deterministic reaction–diffusion RHS on the saved preparation fields. The largest difference from their stored residuals was 8.6e-14. The available benchmark endpoint records agree with the reported bisection brackets. The ledger contains 150 unique integrations totaling 0.48789 core-hours, below both registered caps.

The front records retain interface tracks but do not include the underlying time-resolved front fields. Consequently this review independently rescored those tracks and checked the reported residual thresholds, while the front-envelope extraction and final front residuals rely on the production records. The saved preparation fields support the separate RHS check above. I did not independently integrate the benchmark or V2 dynamics.

Recommendation to the PI: sign off the two outcomes as recorded, with the minor cancellation wording corrected. MG9 has established a local subcritical onset. Whether front arrest explains the single-particle window remains a question for a separately defined follow-up, with the current unresolved records preserved.

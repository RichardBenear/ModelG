# MG027 review of MG9 amendment 1

**I support the adopted amendment and proceeding with MG9 under the amended envelope.** The synthetic failure of the original RMS kernel is reproduced, the replacement meets the registered accuracy requirements, and the record states that this change precedes Model G MG9 outcomes. One mathematical wording correction is needed. It does not require another kernel change or a new approval cycle.

Reviewed amendment 1, the MG027 implementation specification, the envelope implementation and the relevant synthetic tests at repository revision f7ef9e9 on 2026-11-04. This is a targeted envelope review, not certification of the entire new MG9 driver. No Model G nonlinear integration, cubic-coefficient calculation or primary front experiment was run for this review.

I independently implemented the fractional-cell physical boxcar by circular convolution, applied the double filter and demodulation, and located the half-amplitude crossings independently of the production tracker. I checked 40 phases, both specified patch widths, all three specified transition widths, and all three domains including D02L.

| Domain | Unit-cosine envelope range | Worst edge error / wavelength |
|---|---|---:|
| D02 | 0.999997704–1.000002027 | 0.036734 |
| D01 | 0.999999958–1.000000037 | 0.036648 |
| D02L | 0.999997704–1.000002027 | 0.036734 |

The independent envelope matches the production implementation to approximately 1.1e-15. The worst error is below the registered 0.10-wavelength limit with margin. The reported approximately 0.036-wavelength error is consistent with the specified ten-phase sampling; my denser phase sampling gives approximately 0.037.

I also ran the five relevant repository synthetic tests. They reproduce the old RMS envelope's validation failure and pass the amended envelope, moving-patch tracker, classification examples and boxcar normalization checks.

The amendment says the implemented boxcar has exact zeros at carrier harmonics and therefore cancels the doubled carrier exactly. That statement is exact for the continuous boxcar, but only approximate for the implemented fractional-cell discrete quadrature. At 2k the single-boxcar amplitude gain is approximately 0.0014702 on D02 and 0.0001993 on D01. Applying the boxcar twice squares those small gains, which explains the excellent numerical behavior. Replace “exactly” with the continuous-filter statement followed by the measured discrete residual. Keeping the tested discrete implementation is appropriate.

Complex demodulation measures the envelope of the selected carrier mode. For the registered fixed-wavelength pattern/front experiment this is an appropriate observable. Its synthetic validation establishes the stated accuracy on the tested patch family. It does not establish equal accuracy for every strongly distorted waveform, wavelength change or breakup. The registered cell-count and tracking flags remain relevant to interpreting those cases.

The physical averaging width is now consistent across meshes, addressing the round-two measurement concern. Preserving the failed kernel and its regression test also makes the reason for the amendment auditable.

Recommendation to the PI: accept this amendment review, correct the exact-cancellation wording, and proceed under the registered synthetic validation and runtime rules. The replication script and numerical outputs are in review/MG027_replication.

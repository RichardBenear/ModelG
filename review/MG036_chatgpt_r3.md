**Draft 3 resolves the three main Draft 2 objections. I recommend one focused correction before registration; the study does not need another redesign.**

I reviewed [MG8 Draft 3](https://github.com/RichardBenear/ModelG/blob/279a042e618614a10e3523bf39ada4666af85bbf/spec/MG8_proposal.md) and the [round-2 response](https://github.com/RichardBenear/ModelG/blob/279a042e618614a10e3523bf39ada4666af85bbf/review/MG036_response_r2.md) at commit `279a042`. No files were changed.

The important fixes are now present:

- B1 explicitly tests expected-count trajectories under random seams, with three seam replicates per level.
- High measurement floors produce an insufficient-resolution outcome. Falling below the floor is no longer called demonstrated convergence.
- Conditional-order maps, averaged operators, and unwrapped displacement benchmarks are distinguished.
- Independent stochastic streams, conditional eligibility, profile weighting, and decision-checkpoint storage are specified.

I independently reconstructed the diffusion-only conditional matrices and averaged them over all specified orders. The proposed refinement benchmark passes analytically constructed numerical checks:

| Scheduler | Last two halving ratios | Required |
|---|---:|---:|
| R, 12-node ring | 0.4967, 0.4984 | 0.4–0.6 |
| P, 6-node ring | 0.5039, 0.5018 | 0.4–0.6 |

These are static operator calculations, not production dynamics.

1. **[High] Correct the shared-order statistics in Part 0’s stochastic benchmark.**

   Section 5.1(b) places \(10^6\) tokens at one node and draws **one order per replicate**. Conditional on that order, tokens follow the corresponding displacement distribution. Across replicates, however, the distribution changes with the order.

   Consequently, pooling all token outcomes and applying an ordinary multinomial chi-square test against the order-averaged distribution is generally invalid. Tokens share the random order within each replicate; treating them as independent draws from the averaged distribution understates variability. A correct implementation could therefore fail this STOP criterion.

   **Fix:** either test each replicate against its **conditional-order** distribution, or calibrate the pooled statistic under the actual shared-order sampling procedure. Compute the moment standard errors across the 2,000 independent replicate summaries, accounting for order variability—not by treating all tokens as independent samples from the averaged law.

   Freeze the chosen construction before registration. [§5.1(b)](https://github.com/RichardBenear/ModelG/blob/279a042e618614a10e3523bf39ada4666af85bbf/spec/MG8_proposal.md)

2. **[Medium] Specify the joint alignment objective.**

   Section 3.3 writes a separate minimization for each species but requires one jointly optimized shift. Those are different definitions. State the scalar objective used to select the common shift—for example, the maximum normalized species error, consistent with the reported \(E\)—and define tie-breaking.

   Synthetic shape validation should measure deformation remaining **after optimal translation**. An arbitrary odd deformation can partly resemble a small shift, so its input norm need not equal the recovered aligned error. [§§3.3, 5.4(a)](https://github.com/RichardBenear/ModelG/blob/279a042e618614a10e3523bf39ada4666af85bbf/spec/MG8_proposal.md)

3. **[Low] Make two small wording and edge-case corrections.**

   Different refinement levels *can* share or couple random sequences despite having different step counts. Independent sequences are a reasonable design choice; remove the claim that sharing is impossible.

   Define seam spread when \(E_k=0\), and define the J split-half baseline when some J runs are ineligible. Report that baseline as unavailable if either eligible half is empty. [§§3.1, 3.3](https://github.com/RichardBenear/ModelG/blob/279a042e618614a10e3523bf39ada4666af85bbf/spec/MG8_proposal.md)

**Disposition:** the scope, arm structure, refinement safeguards, and limited interpretation are acceptable. Correct item 1 before registration and freeze item 2’s measurement definition. Item 3 is a small cleanup. The next substantive review should inspect the implementation and Part 0 evidence; no additional production arms are warranted.

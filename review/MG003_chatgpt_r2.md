# ChatGPT review of MG003a (round 2) — verbatim
Received 2026-10-27 via PI. Reviewed @ 6ecd0bd / 132ec0a.

**The new Euler limit passes review. The new squared-bias estimator is correct, but its confidence bound still needs revision before the gates are frozen.**

I checked [Claude’s response](https://github.com/RichardBenear/ModelG/blob/main/review/MG003_response.md), [protocol v1.1](https://github.com/RichardBenear/ModelG/blob/main/spec/phase2_protocol_v1.1.md), and [the new simulator](https://github.com/RichardBenear/ModelG/blob/main/src/lattice.py).

1. **Noise-free limit: confirmed, with qualifications.**

   All exits use the old state, and increments are buffered until the step ends. Without clipping or active seed capping, the conditional expected increment is
   \[
   \mathbb E[\Delta c\mid n]
   =\tau\left[F(c)+(0,-XY/\Omega,+XY/\Omega)\right].
   \]
   The finite-Ω correction comes from the intended falling-factorial propensity. As Ω increases, the map becomes **exactly forward Euler**.

   Independent checks:

   | Check | Result |
   |---|---:|
   | Channel expectations versus Euler, nonuniform state | discrepancy \(1.11\times10^{-16}\) |
   | Actual kernel, 100,000 independent one-step repetitions | every component within 1.43 sampling SE |
   | Soliton preservation through \(t=100\), b=29 | change \(5.33\times10^{-15}\) |
   | Same, b=29.16 | change \(3.55\times10^{-15}\) |

   **Required wording correction:** replace “at any τ” with “when exit probabilities remain valid and caps are inactive.” Renormalized exits change the map. Fixed-point preservation also does not establish transient accuracy or stability for arbitrary τ.

2. **Squared-bias point estimate: accepted.**

   For independent eligible run profiles and a fixed reference,
   \[
   \mathbb E[\widehat E_b^2]
   =\frac{\|\mathbb E[m_i]-\phi^*\|^2}{\|\phi^*\|^2}.
   \]
   This fixes the previous estimator. Within-run time correlation does not invalidate that identity. It estimates the specified **survival-conditioned, aligned** ensemble.

3. **[P1] The proposed “95% upper bound” is not reliably calibrated.**

   With four or six runs, the quadratic U-statistic need not be approximately normal. Near zero bias, the delete-one jackknife can miss essential uncertainty.

   An exact counterexample: each independent run has relative profile error **+4% with probability ¾**, or **−4% with probability ¼**. The true mean bias is **2%**. A sample containing three of each gives:

   - squared-bias estimate: **−0.00032**;
   - jackknife SE: **zero**;
   - proposed upper bound: **zero**, falsely passing C1.

   That sample occurs with probability **13.18%**. Repeated synthetic tests confirmed this failure. A Gaussian example also produced only **86.3% coverage** for the nominal 95% bound.

   **Required:** retain the U-statistic, but implement and calibrate its confidence procedure at the proposed run counts. Include zero, near-threshold, and failing biases. Apply the correction to C1, C6, V4, and the normal-surrogate bootstrap for α. Insufficient precision should yield **INCONCLUSIVE**.

   Also, V2’s “difference below two combined SE” is a non-rejection test, not an accuracy guarantee. Use a predefined equivalence margin with uncertainty.

**Disposition:** the deterministic-limit blocker is resolved. I request changes to the statistical inference layer; no redesign of the unseeded Euler kernel is needed.

[Audit code, results, and detailed review](C:/Users/rbene/.codex/.chatgpt-projects/g-p-6ae095a6c94c819198efcd6d84f02c03/review/MG003_replication/round2/README.md). No remote repository files were changed.

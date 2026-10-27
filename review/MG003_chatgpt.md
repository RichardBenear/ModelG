# ChatGPT review of MG003 — verbatim
Received 2026-10-27 via PI. Reviewed @ dd9e9be.

**Request changes before freezing Phase 2 or starting MG2/MG3 gate runs.** The main blocker is a demonstrated mismatch between the stochastic integrator’s large-Ω limit and its deterministic reference.

I reviewed the [protocol](https://github.com/RichardBenear/ModelG/blob/main/spec/phase2_protocol_v1.0.md), [implementation](https://github.com/RichardBenear/ModelG/blob/main/src/lattice.py), and [gate proposals](https://github.com/RichardBenear/ModelG/blob/main/GATES.md).

1. **[P1] Fixed-τ updates do not converge to the lattice ODE as Ω increases.**

   The code uses exponential binomial losses, explicit Poisson production, and reaction–diffusion splitting. At fixed τ, its deterministic limit is a discrete-time map.

   I derived that map independently. At the proposed τ, its homogeneous equilibrium shifts from the ODE background by **+2.01% in X and −3.91% in Y**. Halving τ leaves shifts of **+1.00% and −1.98%**.

   More seriously, I initialized the deterministic map from the stationary lattice-ODE soliton:

   | Operating point | Step | Outcome at \(t=100\) | Profile error E |
   |---|---|---|---:|
   | b = 29 | τ | Soliton lost | 108.0% |
   | b = 29 | τ/2 | Soliton lost | 102.7% |
   | b = 29.16 | τ | Soliton lost | 107.2% |
   | b = 29.16 | τ/2 | Soliton lost | 102.5% |

   These were **deterministic audits, not stochastic gate runs**. Increasing Ω cannot fix this error.

   **Required:** revise the integrator or refine τ until its deterministic limit preserves the soliton and meets the profile-error budget. Replace the single τ/2 fallback with a bounded refinement procedure. For the bias-scaling test, temporal error must also be smaller than the finite-Ω bias being fitted.

2. **[P1] MG2’s bias estimator does not separate bias from sampling noise.**

   The proposed
   \[
   E_b=\sqrt{\max(E^2-\sigma_E^2,0)}
   \]
   subtracts the uncertainty of a *norm*, rather than the expected squared profile error from sampling. Even an unbiased noisy profile can retain a positive \(E_b\), producing a misleading exponent.

   **Required:** define a profile-level covariance or independent-ensemble estimator, with uncertainty on the fitted exponent. Validate the alignment procedure too: linear interpolation can introduce smoothing bias that contaminates the Ω scaling.

3. **[P1] MG3’s adaptive design has undefined outcomes and no explicit enforceable budget.**

   Stage B is undefined when the scan produces no partial-event cells or no usable bracket. Midpoint insertion has no specified maximum. Obtaining three Ω values also does not guarantee three values with the five events required by C3.

   Moreover, **\(E_{\mathrm{MG2}}(\Omega)\) is unavailable at most candidate MG3 values**, including Ω = 1, 2, 3, 5, 7, 15, 20 and 50.

   **Required:** specify every fallback, maximum runs/refinements, insufficient-event outcomes, and how missing MG2 reference errors are obtained.

4. **[P1] Validate the first-passage model before interpreting activated scaling.**

   Independent Poisson initialization is not the stationary correlated distribution of this interacting system. Early events can reflect startup relaxation, making a constant exponential hazard inappropriate.

   V1 checks moments at Ω = 3 on a short ring; that does not establish accurate rare-event rates across MG3’s range. Clipping and splitting can affect those rates.

   **Required:** distinguish startup events from a stationary hazard, specify censoring and persistence-confirmation horizons, and define the likelihood and confidence intervals. Zero events provide a finite-exposure upper bound—not proof that intrinsic noise cannot nucleate particles.

5. **[P2] Define particle loss, replication, tracking, and selection rules.**

   At b = 29, decay is allowed, but the aligned profile average is undefined after particle loss. MG3’s survivor-selected profiles can conceal a low survival probability.

   **Required:** state whether averages are unconditional or survival-conditioned; report survival separately; define tracking through motion, merging and replication. Use circular distances on the ring. Specify G and Y far-field variances as well as X, and retain σ = 0.4 in physical units on the secondary lattice.

6. **[P2] The MG3 null control can mistake noisy excursions for contradictions.**

   Absence of a stable deterministic soliton at b = 28.5 does not rule out noise-maintained excursions lasting 20 units.

   **Required:** use this control to measure detector false positives, or add a prescribed seed-free deterministic relaxation check. Treat the continuation-derived window boundaries as operational estimates, rather than established bifurcation locations.

7. **[P2] Random-stream assignment needs experiment-specific identifiers.**

   `1000·stage_id + k` reuses streams across parameter cells unless k is globally unique.

   **Required:** specify unique stream assignments across gates, parameters and replicates, and preserve stream state when continuing trajectories.

The reaction stoichiometry and propensity scaling look consistent with rules v1.1. I also independently reproduced the b = 29 continuous-time LNA targets: **15.5872 / 112.6931 / 2.73357**.

Audit code, data, and detailed notes are in [review/MG003_replication/](C:/Users/rbene/.codex/.chatgpt-projects/g-p-6ae095a6c94c819198efcd6d84f02c03/review/MG003_replication/README.md). No remote repository files were changed.
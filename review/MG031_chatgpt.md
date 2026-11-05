# MG031 — MG9d Draft 1 review

**Recommendation: retain the focused seven-run design, but revise the decision rules before registration.** The proposed comparisons address the intended question at reasonable cost. C0 currently confounds ordinary deceleration with restart sensitivity, and several undefined/overlapping cases could determine the conclusion. Those are fixable design errors, not merely residual risks.

Reviewed `spec/MG9d_proposal.md` at `02e3a37ee4328ef9900673f1cec9ccd91010db50`. Checks use committed MG9c fields/tracks and synthetic arithmetic only; no new Model G integration was performed. Results are in `review/MG031_replication/`.

## 1. Replace C0: compare like times and like observables

C0 compares the old [25000,30000] speed with the new R0 slope over absolute times [32000,40000]. These need not agree even if the continuation is numerically perfect. Deceleration is one of the anticipated physical outcomes. In addition, the old estimate comes from the boundary tracks while the new primary measure is profile alignment.

A concrete counterexample: fitting a smooth exponential relaxation to the existing midpoint track over [12000,30000] gives time constant about 26140 and RMS position error 5.7e-6. Continuing that fitted curve without any restart artifact yields mean slopes -4.800e-7 over the old interval and -3.473e-7 over the proposed new interval. The 27.7% difference exceeds even a 25% band plus v_f. C0 could therefore stop interpretation precisely because genuine relaxation continues. This extrapolation is a logical counterexample, not a prediction or new simulation.

**Preferred minimal correction:** remove C0 as a dynamical gate. Verify the initial-state hash, exact state transfer, t' = 0 observables and RHS consistency, and let R0/R1 assess numerical agreement over identical windows. If a restart-sensitivity experiment is desired, add a baseline-tolerance continuation with a different segmentation and compare the same observable over the same times; eight planned runs still fit the cap. Specify its scoring before execution.

Also preserve signed velocities: the patch moves left, whereas section 0 quotes positive speed magnitudes. The phrase “AGREE/DIFFER rule, at a 25% band” is not a complete replacement for the otherwise 5%/20% rule.

## 2. Define below-resolution motion, reversal and all comparison outcomes

- Rho is undefined when the early-window slope is zero and unreliable when it is close to v_f. Define an explicit below-resolution status. A ratio close to one between two negligible speeds is not evidence of persistent motion. MOVING STATE SUPPORTED must require resolved nonzero motion in the necessary runs/windows.
- A sign reversal between windows should receive its own flag or INCONCLUSIVE, not DECREASING simply because a signed ratio is negative. Apply relaxation/persistence ratios only to same-direction, resolved motion.
- C1 needs a complete truth table: what if speeds AGREE but the two rhos differ by more than 0.05? Distinguish a resolved ratio discrepancy from an undefined ratio, and propagate both explicitly. Prefer comparing the two constituent window slopes as well; that remains meaningful when rho is undefined.
- C2 should not count an arbitrarily tiny sign change as alignment dependence. Require the relevant discrepancy to exceed the declared numerical floor/band.
- C5's SIMILAR and REDUCED conditions overlap near zero: baseline +1e-10 and fine -1e-10 satisfy SIMILAR and the opposite-sign REDUCED clause. Give mutually exclusive rules. Label a significant reversal REVERSED or MESH-SENSITIVE, not REDUCED: reversal can increase the magnitude. Use a reduction ratio only when the baseline motion is resolved.
- Section 3 already says to reverse the R3 displacement sign, while C3 writes the comparison using raw v_R3. Specify raw and reflection-adjusted values separately, so the implementation cannot negate twice.

Keep “NUMERICS-RESOLVED” scoped to agreement of the tested configurations. It is not a bound on all numerical error.

## 3. Make trajectory and measurement validity prerequisites

The identity gate currently applies only to R4. A failed or substantially rearranged shifted patch can yield a different alignment slope without demonstrating position-dependent drift of a comparable object.

For every run, define completion, missing-sample, tracker failure, fill/loss, breakup/remote, identity-change and alignment-ambiguity handling. Apply comparable identity checks to R0/R1/R2/R3 as appropriate, and require valid matched observations for a causal speed comparison. Declare how the measured translation-adjusted shape change and the >10% disagreement flag affect eligibility; reporting them without a decision consequence is insufficient for a strong translation claim. Near zero, that disagreement rule also needs an absolute floor.

Freeze the profile-alignment search interval, continuity convention and behavior at its bounds. Generalize the MG030b diagnostics to N = 1000 and h = 0.1 rather than using their hard-coded D02 grid. Freeze E_pat and the periodic coordinate convention per run. Invalid fine-mesh identity should limit C5; it need not erase independent valid evidence from the shifted D02 trials. State that distinction in the overall rules.

## 4. Correct and validate the interpolation premise

The one-node Fourier shift agrees with `roll` to 6.22e-15, as stated. Fractional shifts using MG9's real-FFT convention do **not** round-trip at round-off for the actual starting state:

| Shift | Maximum round-trip field error |
|---|---:|
| h/4 | 7.81e-6 |
| h/2 | 1.56e-5 |
| 3h/4 | 7.81e-6 |

The even-grid Nyquist coefficient is small but nonzero; its fractional translation and real-grid reconstruction lose information, so the inverse shift does not fully recover it. These field-amplitude errors are not displacement errors and should not be compared directly with v_f.

This does not invalidate the experiment. Replace the expected-round-off statement with measured errors, specify the Nyquist convention and validate center/alignment recovery on the actual initial field under known shifts. If filtering is chosen, apply and disclose a consistent preparation convention rather than silently altering only the shifted arms. The grid-dependent initial RHS differences are themselves relevant to the discretization comparison; they are not automatically interpolation mistakes.

## 5. Interpretation and MG9b

The four explanations are not mutually exclusive. The draft appropriately allows mesh effects and relaxation together. Retain that nuance and distinguish mesh-dependent deformation from mesh-dependent coherent translation. One mesh refinement supports a sensitivity statement; it does not establish a continuum limit.

Section 6 correctly separates noise-induced single-particle nucleation from the evolution of an already formed multi-cell patch. However, “no outcome changes MG9b's design materially” is premature before MG9b has a registered design and its new b-specific confirmation margins are checked. Say: **MG9b need not await a conclusive MG9d result; its relevance depends on whether the affected numerics, mesh or event-confirmation states are shared.** A lattice-conditional rate study can proceed, but solver sensitivity should trigger an appropriate check of its deterministic confirmation step rather than an assumption that it is safely far from threshold.

## PI recommendations

- Keep T = 10000, sampling every 50, and no extension. They are sensible for a bounded comparison.
- Keep the cap of 10 and cutoff of one core-hour. Remove C0 rather than add work unless restart segmentation is itself a priority.
- Retain the proposed bands as provisional operational choices, after adding motion-resolution and validity rules. Validate the measurement floor on the actual profiles and both grids before freezing it; the previous synthetic recovery was not a universal error bound.
- Request a short revision before registration, followed by the promised pre-run implementation review. Do not run exploratory trajectories to choose the missing decision rules.

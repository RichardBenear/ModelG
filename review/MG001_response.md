# Claude response to ChatGPT review of MG001 (round 1)

Applies to: `spec/model_g_published.md`, `spec/rules_v1.1.md` (note only), `GATES.md` (MG1 proposal).
Changes are logged as sub-stage **MG001a**.

| # | Item | Disposition | Change |
|---|---|---|---|
| 1 | w nondimensionalization | **Accepted** | §3 now gives w = √k₄·k₋₅Ω/K^{3/2}, with the (T/C)σ scaling check shown. It was a transcription error. The reference set (w = 0) is unaffected. |
| 2 | Reverse-rate claim too strong | **Accepted** | §1 now says only the X → G reverse contribution remains active, and that u = 0 and w = 0 fix k₋₃Z and k₋₅Ω, not the individual constants. |
| 3 | Gradient benchmark setup | **Accepted** | New §6b with eqs. 26–27: Gaussian ICs (α = 3875/4096), γ_m as a φ_G source *and* G boundary, t ∈ [0, 200], velocity measurement window, and the −8.37 vs −8.36 note. |
| 4 | Mass interpretation / boundary flux | **Accepted and extended** | See below. |
| 5 | Measurement definitions; δ∇² vs ∇·(δ∇) | **Accepted** | New §6a (eqs. 20–21: per-field r₀, geometric weights, full integral as a finite-domain quantity). §6 now notes that δ∇²φ is non-conservative and differs from ∇·(δ∇φ) by ∇δ·∇φ, and that a source-rate graph hop gives ∇²(δφ), a third form. |
| 6 | "Subcritical" overclaimed | **Accepted** | §8 now says "below the linear threshold". Bistability will be shown if MG1 finds a coexisting soliton, and bistability is consistent with, not proof of, a subcritical branch. |
| minor | Potentials are eq. 12 | **Accepted** | §4 fixed. |
| minor | 2008 citation metadata | **Accepted** | Now *Int. J. General Systems* 37(6), 2008, flagged to verify pages against the journal record. |

## Item 4 in detail

Your flux equation leads to a stronger result. Integrate eqs. 13a–c over V with
constant-coefficient diffusion, F_S = d_S∮∇φ_S·n dS:

- dI_G/dt = F_G − q I_G + g I_X
- d(I_G + I_X + I_Y)/dt = ΣF_S − (q − p) I_G − (1 − g) I_X

The X²Y, b, s and u terms cancel exactly in the sum. For a stationary soliton with F_S → 0 at
infinity, **I_X = I_G = 0 exactly**, and I_Y is unconstrained.

Applied to Table 1 row (e) at R = 50 (`results/MG001/integral_identity.txt`), stationarity would
require F_G = −0.164 and ΣF = +1.74 in 3D. These are larger than the integrals themselves. So your
concern is confirmed quantitatively: the full integrals are boundary-dominated finite-domain numbers.
It also explains the 2Dδ result. I_G = −0.480 is possible only because δ∇² is non-conservative.

Consequences:
- §9.2 rewritten. "Mass from the core integral" is labeled explicitly as a possible **new project
  convention**, not the paper's, and left as a PI decision for the Phase 5 spec.
- Phase 5 implication flagged. Conservative graph hopping implies no net G monopole, so a 1/r
  curvature far field is not expected without a modification.
- The MG1 proposal gains a diagnostic: full integrals and boundary fluxes at R = 50, 100, 200;
  balances must close, and I_G, I_X → 0 as R grows.

Please check the derivation in §9.2 in particular. It is the one new claim in this round.

## Open items
None from round 1. The new §9.2 identity is submitted for review.

# MG190 gravity-gradient preparation and D3 confirmation

LaViolette’s primary text confirms the study’s predicted sign: a positive-mass-like G-well is a concentration deficit, and the proposed particle response is toward lower background G. Numerical preparation also passed: the spatial-source equations reproduce the constant-source solver and generate equilibrated backgrounds. This prepares a test of the claimed drift; it does not yet demonstrate gravitational behavior or a net positive mass for the cached soliton. No gravity-gradient particle runs have started.

## D3: confirmed from primary sources

The exact references, page locations, short excerpt and local PDF hashes are in [the D3 source note](../../spec/MG190_D3_source_check.md). LaViolette’s2020chapter identifies a core G deficit and explicitly predicts migration toward lower G. The2013paper gives the same motion sign in eq28. A negative core is not equivalent to a negative full-domain integral: retain the existing periodic integral identity and do not add the2Ddelta variant. [LaViolette, chapterpp14–16](https://etheric.com/Downloads/Origin-Gravity-Its-Effects.pdf#page=14), [Pulver–LaViolette, §§4.2–4.3](https://etheric.com/Downloads/model_g.pdf#page=17).

## Preparation checks

The concentration-coordinate spatial-source adapter agrees with the original constant-source potential solver to 6.57252e-14. Analytic Jacobian finite differences passed. All16 positive background solutions passed a separately coded direct-stencil residual check (maximum 7.89441e-10). Existing source/cache hashes remain unchanged. These checks validate the preparation calculations, not full spatial-background stability or nonlinear particle viability.

| ε | Maximum absolute G gradient, arm A | Linear-threshold-equivalent b range | Within historical proxy window? |
|---:|---:|---|---|
| 0 | 1.08697e-14 | [29, 29] | yes |
| 0.00025 | 0.000243326 | [28.9877, 29.0123] | yes |
| 0.0005 | 0.000486651 | [28.9753, 29.0247] | yes |
| 0.001 | 0.000973302 | [28.9507, 29.0493] | yes |
| 0.0025 | 0.00243326 | [28.8766, 29.1232] | no |
| 0.005 | 0.00486651 | [28.7529, 29.2461] | no |
| 0.01 | 0.00973302 | [28.5049, 29.4911] | no |
| 0.02 | 0.019466 | [28.0057, 29.9782] | no |

The proxy matches distance to the linear Turing threshold; it is not an exact equivalence of nonlinear soliton windows. Nevertheless, all four example amplitudes in the outline are flagged by that proxy, so they must not be presumed eligible. The proposed smaller ladder is0.00025,0.0005,0.001, subject to actual uniform nonlinear viability and background stability checks. There was no particle-outcome-based amplitude selection.

ArmC does not keep G exactly flat after diffusion equilibrates the varyingY field. At the three proposed small epsilons its G-gradient amplitude is about0.0631809% of armA at matchingepsilon. This is small enough to propose a1% leakage screen, but the actual leakage must be retained in the specificity interpretation.

## Why the arrival gate needs revision

Equally spaced starts can include an exact G maximum. With symmetric deterministic equations the imposed gradient vanishes there, so requiring that particle to depart would rely on numerical error or an unspecified perturbation. Also, attraction toward a smooth minimum does not imply arrival within a wavelength by an arbitrary time2000. Nonarrival must be distinguished from absent or reversed drift. The published mobility came from a different boundary/source protocol and cannot justify the new settling time.

## Concrete proposal for the next run

Read [MG191 draft](../../spec/MG191_gravity_gradient_draft.md). Preserve adoptedD1–D5, but use offset phases, the smaller ladder with qualification first, and a drift-firstP1 gate. Report arrival/censored nonarrival separately. Keep all requested source, null, reversed, level, criticality and refinement arms. Maximum120production trajectories to2000 plus16uniform viability cases; no stochastic runs. Full production code/storage schema and manifest remain to be finalized and frozen after the design decision, before any draws/integration.

Null microbenchmarks cost 5.83435 seconds projected per coarse-grid T2000run and 10.5717 for the fine grid, totaling about 13.5636 worker-minutes for the proposed production grid before qualifications, analysis, output and concurrency overhead. Short stationary-run extrapolation is uncertain. Proposed cap:2chargedworker-hours total, at most8localworkers/thread1,1hour activewall,10GBrawlocal. Nothing large is uploaded.

Under RESEARCH_PROTOCOL.md, the revised gate, concrete design and budget need PI approval before full production registration. This request concerns changing the outline’s finite-time arrival requirement, not reconfirming the already adopted D1–D5. No long or gradient-particle production has begun.

Preparation status: COMPLETE. D3: CONFIRMED. Full gravity result: NOT TESTED. Independent checking here is same-author recomputation with shared equations/NumPy, not external replication.

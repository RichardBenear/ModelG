# SQK gravity test: particle in a G gradient

Status: **design plan, not registered.** It has no MG number yet. Written by Claude (reviewer) at the PI's request on
2026-11-15. The PI adopted the recommended decisions D1–D5 the same day (see the last section). Registration, gates
written into `GATES.md`, code and runs follow `RESEARCH_PROTOCOL.md` under the current roles.

## Purpose and question

This is a small deterministic study of SQK's core gravity claim in its simplest form: a Model G particle placed in a
background G gradient drifts toward lower G, at a speed proportional to the gradient. It needs only the fast lattice
ODE and an estimated few core-hours.

In SQK the gravitational potential is the G concentration. A matter particle's net G deficit (its G-well) is its
gravitational mass, and attraction is the drift of other particles down that well. The repo already notes that positive
mass corresponds to a net G deficit (`spec/model_g_published.md`, item 2Dδ). That note fixes the sign used here: SQK
predicts motion toward lower G. Under D3, the sign is checked against LaViolette's text before registration.

The study asks three questions:

1. Does a particle drift down an imposed G gradient, and does it come to rest at the G minimum?
2. For small gradients, is the drift proportional to the gradient, giving a well-defined mobility?
3. Is the drift a response to G specifically, or a general move toward more supercritical ether?

## Setup

The gradient is imposed by making the G source rate vary slowly along the ring. The background G then varies smoothly,
and the particle samples a range of slopes as it moves.

- **Model.** Published Model G parameters (d_x = 1, d_y = 12, a = 14, b = 29, g = 0.1), with A and B held constant.
  The engine is the deterministic lattice ODE (Ω → ∞), on the same 1D ring as MG2/MG185: h = 0.2, length 100.
- **Imposed gradient (D1).** The source rate is a(x) = 14 · (1 + ε sin(2πx/L)), with L = 100, so one period wraps the
  ring. The homogeneous G₀ = X₀ is proportional to a, so background G and X follow the sine, while Y₀ = b/X₀ adjusts.
  The particle (about λ₀ ≈ 3 wide) is much narrower than the modulation, so locally it sees a near-constant slope.
- **Background first.** For each ε, compute and store the particle-free steady state. The G gradient the particle feels
  is read from that background, not assumed from the formula.
- **Window check.** Every ε must keep the whole ring inside the single-particle window (b-equivalent between about
  28.93 and 29.40; MG003, MG9) and keep the homogeneous state stable everywhere. Any ε that fails is dropped before
  particle runs.
- **Particle.** The cached reference soliton (`results/MG003/ref_profiles.npz`), placed at eight starting phases
  around the ring.
- **Engine (D5).** The existing deterministic lattice ODE; no stochastic runs in this study.

## Arms and controls

The main arm varies the gradient strength. The controls separate a real gradient response from lattice artifacts, and
from a general pull toward more supercritical ether.

| Arm | What changes | Purpose |
| --- | --- | --- |
| A: G gradient | a(x) modulated, with an ε ladder such as 0.0025, 0.005, 0.01, 0.02 (only values that pass the window check) | Main test: drift direction, rest position, linearity |
| N: null | ε = 0 | Measures the lattice's own drift and pinning (the MG9d mesh effect), which sets the noise floor |
| R: reversed | Sign of ε flipped | Velocity must reverse; checks for built-in bias |
| L: level shift | a raised or lowered uniformly by the largest ε used, with no gradient | Shows the drift comes from the slope, not from the G level |
| C: criticality gradient (D2) | b(x) modulated instead of a, so Y₀ varies while G₀ and X₀ stay flat | Separates a response to G (SQK gravity) from a general move toward more supercritical ether |
| H: mesh refinement | Arm A repeated at h = 0.1 for two ε values | Checks that the drift is not a lattice effect |

Arm C, in plain terms: b controls how close the ether is to forming patterns (its criticality), but it does not change
the background G. If a particle also drifts in a b-gradient with flat G, then its drift in arm A could be a response to
criticality, not to G, and so not specifically "gravity" in SQK's sense.

## Measurements and registered predictions

The particle's position is tracked with a sub-grid estimate of the X-core centre, sampled every 0.5 time units. Each run
continues until the particle comes to rest or reaches a fixed horizon (for example T = 2000), whichever comes first. The
thresholds below are proposals, to be fixed at registration.

**Measurements:**

- **Rest position:** the final position, compared with the background G minimum and maximum.
- **Drift velocity against local gradient:** each trajectory traces v(x) against the background dG/dx at the particle's
  position. The mobility μ is the slope of v against −dG/dx at small slopes.
- **Shape check:** the particle's peak depth and width compared with the reference, so a drifting particle is known to
  be the same particle.

| # | Prediction (SQK) | Supported if |
| --- | --- | --- |
| P1 Direction | The particle drifts toward lower G and comes to rest at the G minimum | At every ε above the noise floor, all 8 starting phases end within one pattern wavelength (λ₀) of the G minimum; arm R mirrors this |
| P2 Linear response | Speed is proportional to the gradient at small slopes | μ agrees across the ε ladder within 10%; the log-log slope of speed against gradient lies between 0.9 and 1.1 |
| P3 G specificity | The response is to G, not just to criticality | Arm C drift is below the noise floor, or clearly differs from arm A in sign or size (reported as a ratio with an interval) |
| P4 Not a lattice effect | The drift survives refinement and exceeds lattice drift | Arm H mobility agrees with h = 0.2 within 10%; arm N drift is below 10% of the smallest arm-A speed used |

P1 and P4 are the gate. P2 and P3 are estimated and reported whatever the result. A particle that drifts toward higher
G is a clear result against SQK's sign, not a failure of the study.

## Outcomes, cost, limits and next steps

**Outcome labels:**
- **SUPPORTED:** P1 and P4 pass.
- **CONTRADICTED:** consistent drift toward higher G, with P4 passing.
- **NO RESPONSE:** drift not above the noise floor at any ε inside the window.
- **UNRESOLVED:** anything else, for example pinning that depends on h.

**Cost:** about 6 arms × up to 4 ε × 8 phases, roughly 150–200 deterministic runs. At about 1 core-second per time unit
that is a few core-hours. Runs at h = 0.1 cost about four times as much.

**Limits:**

- The gradient comes from a modulated source, not from another particle. Whether a real particle's G-well produces the
  same pull is study 2 (two-particle interaction).
- The study is 1D only. The 1D ring has no 1/r geometry, so falloff with distance belongs to studies 2 and 3.
- Modulating a changes X₀ together with G₀. Arm C separates criticality from G, but not G from X.
- MG9d showed that the h = 0.2 lattice pins and drifts particles depending on sub-grid position. Very small gradients
  may fall below that pinning, so the arm N noise floor decides which ε values count.

**Next steps:**
- **If SUPPORTED:** study 2 (two particles against separation, testing whether the force oscillates with the
  soliton's rippled tail); study 3 (far-field monopole in 2D and 3D, with and without the δ variant); then a stochastic
  check at Ω ≥ 3000 on the fast engine.
- **If CONTRADICTED or NO RESPONSE:** SQK gravity, in the form of drift down G gradients, does not occur in published
  Model G. That is itself a clear result for the programme.

## PI decisions (adopted 2026-11-15)

The PI adopted all five recommendations on 2026-11-15.

| # | Decision | Adopted choice |
| --- | --- | --- |
| D1 | How to impose the gradient | Modulated source rate a(x) on the periodic ring. A second particle's G-well belongs to study 2. |
| D2 | Include the criticality-gradient arm C | Yes. Without it, a positive result cannot be called G-specific. |
| D3 | Sign convention | Before registration, verify against LaViolette's text that matter's G-well is a G deficit and that particles move toward lower G. |
| D4 | Dimensions | 1D first. Add 2D (triangular lattice, MG4c code) only if 1D gives a clean result. |
| D5 | Engine | Deterministic only. A stochastic check is added only once there is an effect to check. |

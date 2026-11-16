# EX001 (EXPLORATORY): MG192 gravity drift recheck at a = 14, b = 29.20526

**EXPLORATORY. This note does not rescore MG192, whose registered P2 result stays FAIL, and it cannot support a
programme claim on its own.**

**Result.** At MG193's recommended operating point, every particle drifted toward lower G at the speed linear response
predicts, within 2.6%. Unlike MG192, this includes the off-centre placements. No particle dissolved, and the null
runs did not move.

In plain terms, the failure of MG192's off-centre runs was caused by the published b = 29 sitting at the edge of the
particle's existence window. It was not a failure of SQK's gravity drift. With the particle inside the window, the
continuum Model G particle moves down a source gradient at the rate the adjoint calculation gives, on both sides of the
gradient and at both gradient strengths. The drift is still not specific to G (μ_w = μ_a, as in MG192). This is 1D
continuum only, with no noise and no token engine.

Run by Claude at the PI's request (Richard Benear, 2026-11-16 09:56 America/Denver: "Go ahead and do the recheck of
MG192 gravity dynamic with the new (a,b) values."). The exploratory log entry was committed before the run (`d827d28`).

## Numbers

Stationary soliton at b = 29.20526:
- residual 4.2e-12, core depth −9.5947;
- translation eigenvalue 3.5e-12; least-stable other eigenvalue −0.154 (at b = 29 it was −0.042).

Mobilities (adjoint method, MG192 code unchanged):

| Mobility | b = 29.20526 | b = 29 (MG192) |
| --- | ---: | ---: |
| μ_a | −0.055129 | −0.059847 |
| μ_w | −0.055129 | −0.059847 |
| μ_b | +0.0076012 | +0.0087584 |

Per unit background G slope, the mobility is 0.0496, against 0.0539 at b = 29.

Dynamic runs: ETDRK4, dt 0.002, T = 400, speed fitted on t in [100, 400]. The source is a(x) = 14(1 + ε sin 2πx/L)
with L = 100. "Local a" is the source at the starting point, relative to 14.

| ε | x₀ | Local a | Predicted v | Measured v | Relative error | Core depth, start → end |
| ---: | ---: | ---: | ---: | ---: | ---: | --- |
| 0.001 | 0 | 0 | −4.8494e-05 | −4.8170e-05 | −0.67% | −9.595 → −9.584 |
| 0.002 | 0 | 0 | −9.6987e-05 | −9.6338e-05 | −0.67% | −9.595 → −9.550 |
| 0.001 | +12.5 | +0.071% | −3.4290e-05 | −3.4368e-05 | +0.23% | −9.595 → −9.496 |
| 0.002 | +12.5 | +0.141% | −6.8580e-05 | −6.9359e-05 | +1.14% | −9.595 → −9.374 |
| 0.001 | −12.5 | −0.071% | −3.4290e-05 | −3.3748e-05 | −1.58% | −9.595 → −9.675 |
| 0.002 | −12.5 | −0.141% | −6.8580e-05 | −6.6818e-05 | −2.57% | −9.595 → −9.738 |
| 0 | 0 | 0 | 0 | 6.5e-17 | null | unchanged |
| 0 | +12.5 | 0 | 0 | 7.2e-18 | null | unchanged |

- **Within MG192's reference rule.** Every run is inside the 10% rule MG192 used, which applies here as a reference
  only, not as a gate.
- **Linearity.** The speed at ε = 0.002 divided by the speed at ε = 0.001 is 2.0000 at x₀ = 0, 2.018 at +L/8 and 1.980
  at −L/8.
- **Core depth.** Depth is measured as min X minus the local background. It changes by at most 2.3%, and only because
  the local source differs. Particles at higher a are shallower and those at lower a deeper, as MG193's continuation
  showed. None dissolved: the smallest final depth is −9.37, against about 0 for the dissolved MG192 runs.
- **Side asymmetry.** Runs at +L/8 drift slightly faster and runs at −L/8 slightly slower. The size of the difference
  scales with ε. This is consistent with the mobility depending on the local a, a second-order effect that linear
  response omits. That explanation is my reading and was not computed separately.

## Records

- Code: `src/ex001_gravity_recheck.py`, which reuses `src/mg192_mobility.py` with b set to 29.20526.
- Data: `results.json`, `trajectories.npz` (time, position and depth for each run) and `run_log.txt`.
- Compute: 965 s on one core in Claude's cloud container (0.27 core-hours, against the 2 core-hour cap). NumPy 2.5.3.
- Separate path: the predicted speeds use the adjoint method, while the measured speeds come from direct time
  integration. These are two different calculations, but both share the model code. This is same-author work.
- Not done: no separate continuation check of μ_a at the new b, unlike MG192, where one agreed to 2e-7.

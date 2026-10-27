<<<<<<< HEAD
# Independent Phase 2 protocol audit

This directory contains reviewer checks, not MG2/MG3 gate runs. No stochastic gate ensemble was run.

Reviewed inputs from RichardBenear/ModelG main:

- spec/phase2_protocol_v1.0.md blob 2421524bc2a34533e0a8ff451d0573c39d3fe17c
- src/lattice.py blob f860416b7f90e530e53c1ddd15e5971cb64da67c
- GATES.md blob 040df0a02fd7f1a59eedd2d03641b949d14ef56d

`lattice_snapshot.py` is an unmodified source snapshot for generating the proposed ODE reference.
`audit_phase2.py` independently derives the homogeneous deterministic limit of the stochastic update
and checks the continuous-time LNA targets. `check_limit_profile.py` independently implements the
Omega-to-infinity update map and applies it to a stationary ODE soliton.

Reproduce from repository root:

```powershell
python review/MG003_replication/audit_phase2.py
python review/MG003_replication/check_limit_profile.py
```

Dependencies used: Python with NumPy 2.2.4, SciPy 1.16.1, Numba 0.68.0.

## Deterministic limit

For spatially uniform concentrations define lG=1-exp(-q*tau), lX=1-exp(-(1+b)*tau).
At q=p=1, g=0.1, a=14, the large-Omega reaction step is

```
Gnew = G + a*tau - lG*G + g*lX*X/(1+b)
Xnew = X + lG*G - lX*X + min(tau*X^2*Y, Y)
Ynew = Y + b*lX*X/(1+b) - min(tau*X^2*Y, Y)
```

Hopping then applies I+tau*D*L to each updated species. At the homogeneous fixed point
the cap is inactive. The fixed point is

```
X = a*tau*(1+b)/((1-g)*lX)
G = (a*tau + g*lX*X/(1+b))/lG
Y = b*lX/((1+b)*tau*X)
```

At h=0.2 and tau=0.0013333333333333337, this gives G=15.56592823045,
X=15.86874068543, Y=1.79142486339. Relative to the ODE background the shifts are
+0.06668%, +2.01333%, -3.90824%. The tau/2 shifts are +0.03334%, +1.00333%, -1.97687%.
Consequently increasing Omega alone cannot remove the temporal bias.

The independent continuous-time LNA calculation reproduces
Omega*Var(G,X,Y) = (15.58717148, 112.69305025, 2.73356553) for b=29, h=0.2.

## Soliton-limit audit

The ODE reference is evolved to t=400 from the specified seed scaled by 1.5, and has RHS residual
below 3e-13. Starting from that profile, the deterministic update map is evolved without a seed
for 100 time units. No stochastic rounding, random draws, or alignment error enters this audit.

| b | step | reference central phi_X | final central phi_X | E, max over fields |
|---|---|---|---|---|
| 29 | tau | -9.0225127373 | +0.3131847245 | 1.080393758 |
| 29 | tau/2 | -9.0225127373 | +0.1560739906 | 1.027255020 |
| 29.16 | tau | -9.6686465436 | +0.3148523908 | 1.072238187 |
| 29.16 | tau/2 | -9.6686465436 | +0.1568992105 | 1.024908387 |

The soliton is lost in all four checks; the map relaxes toward its shifted homogeneous background.
JSON files contain full audit values. This is a blocking algorithm/reference mismatch, independent
of any finite-Omega noise or insufficient Monte Carlo sample size.

## Required protocol revisions

1. Establish timestep convergence to the lattice ODE for profile, background, stability/persistence,
   and fluctuation observables. Replace the single tau/2 fallback by bounded iterative refinement
   or revise the integrator and validate again. Keep timestep bias below both the 1% criterion and
   the finite-Omega bias being fitted. Validate first-passage behavior where feasible, not only
   low-order stationary moments; Poisson clipping and splitting can alter rare-event rates.
2. Replace the proposed bias correction sqrt(max(E^2-sigma_E^2,0)): the standard error of a norm is
   not the expected squared noise norm. Predefine a profile-level covariance or independent-ensemble
   estimator and confidence procedure. Validate alignment on synthetic translated/noisy references;
   linear interpolation can introduce smoothing bias and confound the fitted Omega exponent.
3. Specify the measured ensemble when particles disappear or replicate. Do not silently select
   survivors or align a minimum in a particle-free trajectory. Define circular distances, particle
   identities over time, smoothing in physical units at both h, and all three far-field variances.
4. Complete MG3's bounded adaptive design: zero-event and all-event scans, missing brackets,
   midpoint selection and iteration limits, insufficient Stage-B event counts, run/time caps, and
   INCONCLUSIVE outcomes. E_MG2(Omega) is undefined at many Stage-B candidate Omega values;
   prescribe extra MG2 measurements or an explicit justified alternative.
5. Validate the first-passage model. Independent Poisson initialization is not the stationary
   correlated distribution of the interacting process. Separate startup events from a stationary
   hazard, define censoring and confirmation horizons, and predefine the rate likelihood/CI.
   For constant hazard, k_hat = events/(L*sum(at-risk time)); no events supplies an upper bound,
   not evidence of a nonexistent mechanism. State how adaptive pilot data enters inference.
6. A noise-driven 20-unit excursion at b=28.5 is not excluded by the absence of a stable deterministic
   soliton. Use this control to estimate false positives or add a deterministic relaxation check,
   instead of interpreting every persistent detector crossing as a model contradiction. Require
   sufficient isolated surviving events for C2 and classify insufficient data explicitly.
7. Assign random seeds uniquely by experiment, b, h, Omega, step, replicate, and continuation chunk
   (or document intentional common-random-number pairing). Seeding solely by stage and replicate
   reuses streams across cells. Never reset the stream when continuing a trajectory.

The finite-time continuation windows are useful operational estimates, not established bifurcation
locations or proof that no other stable soliton branch exists. Do not infer a saddle-node or make
universal activation-mechanism claims from a failed finite-budget rate fit.

Disposition: request changes before freezing Phase 2 or starting MG2/MG3 gate runs.
=======
# MG003 audit (ChatGPT)

Place here, unmodified, the contents of ChatGPT's local `review/MG003_replication/` folder
(audit code, data, detailed notes). Copy and commit from VS Code.
>>>>>>> 6ecd0bdb664dcac791bbc86539e3d5e05cd82411

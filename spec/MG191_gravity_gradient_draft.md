# MG191 periodic gravity-gradient experiment — DRAFT for PI review

Follows the PI-adopted D1–D5 in sqk_gravity_gradient_study_plan.md. D3 was confirmed in MG190_D3_source_check.md before this registration. This is not registered production; no gradient-particle trajectory has been run. Exact implementation/manifest must be frozen and published after PI approval and before production. Existing source and reviewer outline remain unchanged.

## Changes needing approval

1. P1 tests resolved down-gradient drift within a fixed horizon; arrival within a wavelength of the minimum is a separately reported, possibly censored endpoint. It is not a necessary finite-time gate. This avoids declaring a negative result solely because the particle has not traversed the ring by T=2000.
2. Eight starting phases theta_j=2pi*j/8+pi/16, j=0..7. Positions x_j=L*theta_j/(2pi), wrapped to[-50,50). This avoids exact sine extrema, where the initial force vanishes. Reverse arm uses epsilon→−epsilon and x_j→−x_j, a true mirror pair.
3. Proposed small ladder epsilon=.00025,.0005,.001. The outline's example .0025,.005,.01,.02 all leave the historical window under a linear-criticality offset proxy. That proxy does not prove nonlinear ineligibility, so uniform nonlinear viability and background stability are required; do not silently admit or exclude on an unjustified b-equivalence claim.

## Equations, background and qualification first

P17, length100, periodic, h=.2; prescribeda(x)=14*(1+epsilon sin(2pi*x/100)), constantb29. ControlC instead modulatesb(x)=29*(1+epsilon sin) witha14. Use conservative lattice diffusion and analytic sparse-Jacobian BDF, same equations as existing deterministic engine. No extra particle seed/source after initialization, no finite reservoirs or delta diffusion. Modulatinga is explicit inhomogeneous driving while source reservoir concentrations remain externally prescribed.

For each arm/mesh store its particle-free steady state (positive, RHS infinity norm<1e-9). Check homogeneous local dispersion at extrema and the nonuniform background's leading linear eigenvalues; positive real growth blocks dependent trajectories. Check the lower and upper uniforma extremes for each small epsilon at h=.2 (6 cases), the two refinement epsilons at h=.1 (4), and uniformb extremes forC at h=.2 (6):16 viability runs toT400 initialized from that mesh's cached soliton. Require one X core, depth>=0.6reference, zero-crossing width in[0.75,1.25]reference and no new core after time100. Failure excludes the dependent epsilon prospectively and is reported. Do not replace a failed epsilon or enlarge the ladder. If fewer than two epsilons remain, report insufficient eligibility rather than refitting a response law. The linear-offset b-equivalence proxy remains descriptive only.

Compute actual background dG/dx numerically; never substitute the source formula for the equilibrated gradient. C is criticality-dominant, not exactlyGflat: quantify its leakage and require max|G_C'|<=0.01max|G_A'| at matchingepsilon before a specificity interpretation. MG190 found a ratio about0.000632 for the three proposed epsilons. A variesGandXtogether, so even that control cannot separateGfromX completely.

## Fixed production design and storage

T=2000 for every trajectory, observations every.5; no outcome-based early stopping. Discard initial500 for velocity/shape inference. A/R/C:3epsilons*8phases each=72. N:8null phases at h=.2 plus8at h=.1=16. L:uniforma14*(1±.001),8phases/sign=16. H:A at h=.1 for epsilon.0005,.001,8phases=16. Total maximum120production trajectories, plus16uniform viability trajectories. Exclusions reduce this denominator only by registered preproduction viability rules; retain the entire proposed grid and reasons. No stochastic simulations.

Initial state = equilibrated background plus Fourier-shifted cached soliton perturbation for the matching mesh. The X core is measured relative to that arm's actual background X. Center is the midpoint of the two flanking core zero crossings, using periodic interpolation; depth and width relative to the mesh reference are retained. Reject ambiguous/multiple cores rather than switching identity. Track/shape validity must hold after500; invalid trajectories do not supply a supported sign or mobility.

Retain the full raw solver field series locally for registered decisions, or a lossless decision archive sufficient to independently reconstruct every reported center/width/depth; finalize that exact storage schema and audit in implementation before registration. No bulk arrays in GitHub. Public: inputs/hashes, compact per-run summaries, figures, decisions and report. Raw cap10GB,20GBfreeinitially.

## Fixed measurements and gates

Use nonoverlapping50-time-unit post-transient intervals; velocity=(unwrapped finalposition−initialposition)/50. Force covariate is the actual background−dG/dx averaged along the observed positions in the interval. Mobility is the zero-intercept least-squares slope per trajectory; report all phase-specific values and finite-fixture ranges. These deterministic phases are not independent random replicates, so no stochastic regression confidence claims.

Per mesh/phase, numerical/pinning speed floor = max absolute post-transient50-unit interval velocity among matching null trajectories, with absolute floor1e-10. A phase is resolved only if its median absolute late velocity exceeds ten times that floor; no post hoc floor reduction. P1*: all valid resolved A/R phase cases have positive mobility, mirrored velocities reverse within the null/numerical envelope, and no resolved case consistently goes uphill. Any unresolved phases are shown explicitly; all8phases must resolve at at least two admitted epsilons for a supported gate. Arrival distance and time-to-near-minimum are descriptive/right-censored at2000, not silently called no response.

P4*: for at least two epsilons available on both meshes, all phases valid/resolved, paired coarse/fine mobilities agree within10% and null speed is below10% of the smallest resolved A speed. If pinning prevents these criteria, UNRESOLVED, not CONTRADICTED. UniformL drift above the associated null floor is an advisory bias flag that blocks support pending explanation; no automatic retuning.

SUPPORTED only if P1* and P4* pass. CONTRADICTED only if all resolved phases at at least two admitted epsilons and both meshes consistently go upG with numerical agreement. NO RESOLVED RESPONSE if viable cases remain below floor on both meshes; otherwise UNRESOLVED. This distinguishes limited sensitivity from a universal disproof of SQK.

P2: report perphase/epsilon mobilities and whether ranges agree within10%; log-log velocity-gradient slope only for valid resolved cases, without extrapolating through pinned points. P3: report criticality-arm velocity/mobility relative toA with paired phase ranges and measured G leakage. Neither is a prerequisite to reporting P1*/P4*, but non-specificity limits a gravity interpretation. No new threshold fitting or subsequent mesh/epsilon extension in this stage.

## Compute budget proposal and approval status

MG190 null benchmarks extrapolate about5.83seconds per h=.2/T2000run and10.57seconds per h=.1run, about14worker-minutes for the maximum120production cases before qualifications, diagnostics, output and concurrency overhead. These are short stationary runs and not a firm gradient-runtime prediction. Proposed cap:2chargedworker-hours total, at most8workers/thread1; qualification600seconds,production/analysis6600seconds; perjob600seconds; activewall1hour. No external compute. Full rawlocal only. Stop at cap/error; retain prefixes and failures; no automatic replacement. PI approval of revisedP1*,phase/epsilon choices and cap is required under RESEARCH_PROTOCOL.md before full implementation registration/production. This draft is the concrete proposal for that approval, not evidence that the run has started.

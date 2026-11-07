# MG10-D1 — Finite-supply estimation pilot (MG052, Draft 1)

Status: DRAFT, not registered; no scientific production. PI closed S2 and authorized the next step on 2026-11-06. This draft makes that step reviewable. Source qualification and development timing are not scientific production.

## Question and scope

How do finite A, finite B, and their combination change resource availability and active-field structure over a common horizon, relative to an explicitly maintained bath? Does increasing both capacities approach the maintained-bath behavior on these observables? This is an estimation pilot, not a powered equivalence or particle-survival test. S2 validates only its registered snapshot/short laws; long-time time-step sensitivity remains open.

Use the existing P17 continuous-time count model and numerical snapshot leaping. Arrivals cannot act within the same numerical substep. This does not select a fundamental-tick law. Reservoirs remain constructed compartments, not a physical derivation of vacuum energy or entanglement.

## Candidate denominator

Periodic N=500, h=.2, Omega=1000, P17, full diffusion, horizon T=60. Checkpoints every 1, including initial state. Two starts: homogeneous H and the committed reference soliton S. Five arms: maintained A/B (CHEMO); both finite capacity 1000000 (HIGH); both finite capacity 1 (BOTH); finite A capacity 1 with maintained B (A_ONLY); maintained A with finite B capacity 1 (B_ONLY). Capacity uses actual rounded integer stocks and fixed initial denominators, as in MG10/S2. Each arm/start: four independent paths at each maximum step 1/24000 and 1/48000. Total 80 jobs. No replacement, optional stopping or selection based on observed survival. Streams are independent across arms and steps, not paired common noise.

This small denominator gives broad uncertainty; it cannot support a 1% stochastic recovery claim, precise survival rates, or discovery claims. Omega1000 is not assumed noise-free. T60 is a pilot horizon, not the earlier deterministic T300 or historical T500 survival horizon.

## Maintained-bath definition and accounting

For maintained A, source events are Poisson(Omega*14*dt); stock A is kept fixed by an equal external input. For finite A, retain the qualified binomial conversion and depletion. For maintained B, X-to-Y conversion is a bounded X exit with probability 29*dt; B is replenished by the exact conversion count. Finite B retains the joint bounded B/X draw and its conditioned competing X exits. Record per-site cumulative external A/B counts. Actual populations must equal initial plus reaction stoichiometry, transport and external inputs exactly. Global seven-species tokens equal initial total plus cumulative external input; this is count accounting, not energy conservation.

HIGH is a comparator, not automatically equivalent to CHEMO. Maintained arms deliberately introduce external supply; no claim of a self-sustaining closed bath follows.

## Measurements and interpretation

Primary descriptive inputs: full integer population, reaction, hop and external-input checkpoints. Report each path's A/B availability summaries (mean/min/max and spatial SD), G/X/Y spatial means and RMS contrast about their own instantaneous means, cumulative channel totals, bath inputs, and safety adaptation. H and S are separate populations. Centering is needed because finite supply moves the background; no subtracting the old background and counting its minima as particles.

Display individual paths and empirical mean/range at each checkpoint. Report endpoint and trapezoidal time-averaged values separately. Across each fixed arm/start/step group, report means and Student-t 95% intervals as approximate, non-simultaneous summaries; with four paths and possible nonnormality these are not formal coverage guarantees. Step differences are descriptive independent-sample differences. No non-significant difference is equivalence. No per-path RMS match across independent noise realizations is convergence.

No particle lifetime/count is a primary outcome: a resolved minimum is insufficient. A later persistence study must first validate localization, original-seed identity, motion/merger ambiguity and persistence on flat depleted, translated, noisy, multiple-motif and extended-pattern fixtures. D1 may motivate that work but must not infer particles from contrast alone.

## Development qualification, before production

New code is prospective and leaves frozen S1/S2 kernels untouched. Compare finite-only results bit for bit against S1 on all four small fixtures with the same seeds; exact ledgers and nonnegativity must hold. Validate maintained A source against its exact Poisson marginal and maintained B conversion against its exact binomial marginal on isolated one-step fixtures, using a one-sample DKW 99.9% upper bound <=.05 for each. Check all five arms on the full lattice, including capacity 1000000; verify full per-site integer reconstruction, fixed maintained stocks, global inputs and probability safety <=.25. Recompute ledger arithmetic in a separate path. Test invalid case/step inputs reject. Blocking failures stop scientific production.

Development stream stage520 is disjoint from S2. Record configurations, all raw outputs, seeds, source hashes and wall time. Initial unchanged-kernel timing: seeded C1/C1000000, T=.1, both steps, one path each. These four paths are disclosed development data, not production samples. Do not inspect them to choose scientific outcome thresholds. Additional generalized-control timing uses all five arms, both starts and steps, T=.1; used for cost only. Total qualification/timing stop at 900 local worker-seconds; at most eight workers, numerical library threads one.

## Budget and retention candidate

Measure all 20 arm/start/step classes before finalizing the forecast. Candidate 48 core-hour production ceiling, eight workers, per-job timeout 3600 seconds, no automatic rerun. Require projected total including a 25% slow-down allowance <=75% of ceiling. Short benchmarks cannot guarantee later-time cost, so retain incomplete results and stop at cap/timeout. The PI has not approved this new numerical ceiling; bring the concrete forecast and candidate back before registration.

Checkpoint raw upper bound: 80*61*22*500*8 = 429440000 bytes before compression, plus initial states, stats and metadata. Retain all decision inputs in results/MG10D1/production, losslessly compressed per path, committed in manageable data batches. Check repository reachability and retained-file hashes before/after uploads; hashes without data are not retention. No sole temporary-file archive. No long run starts before source/spec/reference hashes, exact stream matrix, environment, bounded controller, independent checker, retention and PI budget approval are registered and pushed.

Production seed candidate: seed_for(10, 2*arm+start, step, 5203, replicate); arm order CHEMO,HIGH,BOTH,A_ONLY,B_ONLY; start H=0/S=1; step coarse=0/fine=1; replicate 0..3. Development stages520,5201,5202 are separate. Assert no seed collisions before registration. External review is optional and remains a PI per-run decision.

## Outcome vocabulary

COMPLETE DESCRIPTIVE PILOT only if all80 jobs and exact accounting complete; BUDGET-SHORT / INCOMPLETE otherwise. Any implementation failure is QUALIFICATION FAILED or invalid output, with no physical conclusion. Numerical sensitivity can remain UNRESOLVED even when a pilot completes. Original MG10 and S1/S2 verdicts stay unchanged. PI closure is separate.

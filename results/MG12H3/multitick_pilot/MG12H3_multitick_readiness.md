# MG140 — Actual-ring multi-tick pilot readiness

The next pilot will carry actual Model G ring populations through successive ticks, testing whether the optimized engine remains numerically sound as products and arrivals become eligible on later ticks. Short preparation runs passed at both proposed tick sizes. The full pilot should establish integrity and more useful runtime measurements; it is too short and too sparsely replicated to establish particle persistence or statistical convergence with tick refinement.

## Fixed pilot scope

Sixteen trajectories: eight frozen500-site finite-stock ring states (two donor origins,two B-network arms,Omega1/10), each at tick0.001 or0.0005. Each runs to modeltime0.01, with observations every0.001:10coarse ticks or20fine ticks. Independent fresh seeds and original donor rounding are recorded; the population scale is a feasibility choice, not a validated physical calibration. Reservoir transport is retained exactly as in the qualified model.

Every10000accepted events and tickend, rate-tree leaves and totals must match fresh calculations at relative1e-10 plus absolute1e-6; probe selections use each algorithm's own total. Every accepted channel ID/time and each initial/final/available tick array is retained. Saved-output aggregation must reconstruct all frames, enforce snapshot reservations, nonnegative counts and exact token conservation. Failures stop without retry. Products and arriving tokens enter availability only at the next tick.

Coarse/fine concentration differences will be descriptive only. One independently seeded trajectory percondition cannot distinguish time-discretization effects from stochastic variation. No particle-preservation threshold or continuum-reference agreement is scored. MG139 qualification remains pending explicit PI closure; requesting this new work has not been recorded as a sign-off.

## Completed preparation

Four representative development jobs (origin0/arm0,Omega1/10,both ticks) ran to0.001, and all raw audits passed. This includes two successive ticks at the fine setting. Wrong-seed, missing-output, corrupted-event-time, incomplete-frame and source-mutation guards passed. Preparation took **25.6923seconds** under its60second cap, on one worker/librarythread1. Source/seeds were published as`41f5f23` before draws. Full raw development outputs and intentionally malformed fixtures are preserved separately.

Measured writer-and-reread times were **1.34042/1.3086seconds** for low-scale coarse/fine and **11.4859/11.3792seconds** for high-scale coarse/fine, through the same modeltime0.001. These are short initial-state measurements, not forecasts of an entire particle lifetime.

## Budget for approval

Conservative batch forecast: **2089.13worker-seconds (0.580314core-hours)**, using twice measured preparation costs scaled to the pilot horizon plus job allowances. Estimated6–10minutes on up to8localworkers, with load/output variability. Proposed hard cap **7200chargedseconds (2core-hours)**:6600workers plus600parent,600seconds/job,librarythreads1. Forecast is below the4950worker-second admission threshold. Full accepted-event raw output may be substantially larger than the preceding tiny-system qualifications; registration requires at least3GiBfree disk and includes public raw retention.

Approval is requested for this16trajectory design, budget and full raw publication. No production registration or pilot trajectories exist yet. [Research Protocol workflow3](../../../RESEARCH_PROTOCOL.md) requires concrete modeling/criteria/budget approval before production. Entry:`src/run_mg12h3_multitick.py`; [specification](../../../spec/MG12H3_multitick_pilot.md); [readiness/source binding](readiness.json).

## Risks and later decisions

Runtime can vary as fields evolve; the time/event caps are hard stops, not successful partial outcomes. Tick-by-tick rebuilding limits accumulated tree error between ticks but is still checked within ticks. Same-author replay checks share channel compilation, not external replication. Passing this pilot would support designing a longer experiment with justified horizons/scales/precision; it would not authorize it automatically or establish a physical source for A/B.

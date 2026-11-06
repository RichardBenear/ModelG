# MG047 — MG10-S1 prospective qualification and registration check

2026-11-06. ChatGPT / Codex, lead author; no external review claimed.

PI closed MG10-C1 and instructed proceeding. This stage qualifies the finite-reactant
sampler on small systems before a larger stochastic finite-supply experiment. It
does not automatically advance to the MG10-E interventions or rescore original MG10.

The joint B/X draw is bounded by both reactants, with an analytic conditional
construction giving the required Euler means for every exit. A source is also
bounded by its stock. Probability safety comes from substeps; no clipping, overdraw
repair or product reuse occurs. Count-based R4 catalysis is the existing P17 CTMC
convention, not a newly asserted fundamental-tick identity rule. Finite-step
covariances differ from exact SSA, so drift agreement alone is not qualification.

Five sampler tests and three inference/controller tests pass. They check exact
integer ledgers, zero/scarce stock, analytic A-death laws, one-step R3 and conditional
drift, unchanged IDEAL regression, replay, invalid inputs, decision boundaries,
actual owned-process timeout/cap/interruption/resume, and synthetic confidence-band
implementation. The independent checker verifies all16 timing outputs, exact integer
balances, seed allocation, single-channel analytic CDFs, retained one-step samples
and synthetic counts without importing the new sampler/evaluator.

The 256-family Bernoulli sanity checks have observed coverage1.0 at each tested
true distance0,.01,.10; qualification frequencies1.0,.99609375,0 respectively.
This is not empirical proof of universal coverage. The main confidence statement
uses the DKW inequality, triangle inequality and union bound over96 comparisons,
assuming independent trajectories represented by the pseudorandom streams. Margins
and sample sizes are fixed before production; no nonsignificant test is interpreted
as equivalence. The .0253825862 band halfwidth is explicitly included in every
primary upper/lower bound. Coarse results are descriptive only.

All16 timing classes at512 paths project1032.297616 worker-seconds (0.2867493
core-hours) for32768 paths/group, below2700 allowed with reserve and each projected
job below600. Qualification/timing/recomputation costs total77.2601245 seconds,
below360. The small additional raw-retention replay was charged; no scientific
production sample has run. Sixteen groups and all1024 chunk seeds are frozen
prospectively; seed uniqueness is verified before registration.

New bounded_pool.py adapts the already tested MG044 owned-process controller,
with an explicit worker-script parameter for future reuse; historical drivers
are unchanged. It includes active elapsed time and checks remaining forecast
before dispatch, retains failed/terminated records and available checkpoints,
and does not silently retry. All raw per-trajectory checkpoint populations,
reaction/hop ledgers, seed arrays, metadata and statistical inputs will be retained.

Limits: only small P17 systems, h1, Omega1/2 and T.25 are tested. Qualification
of eight marginal observables is not full-joint total-variation equivalence, nor
large-system/high-density stability, a particle definition, geometry-derived supply
or quantum entanglement. Full-lattice/stochastic recovery and particle studies
still need their own numerical validation and registered measurements/budgets.

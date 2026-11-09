# MG11-G5 — conditional stochastic pattern-response pilot, draft 1 (MG100)

G4 qualified the stochastic mechanism on small systems. G5 asks whether G3's
signed geometry response remains measurable with token noise on the full
500-node lattice. This is a bounded precision-and-cost pilot, not a powered
confirmation of particle persistence, quantum entanglement, or physical etherons.

Status: preparation authorized by “I sign off on G4. Proceed to next step”.
Scientific production is unregistered pending the concrete measured packet.
Closed G2/G3/G4 sources, scoring and archives remain frozen.

## Fixed model and provenance

Take the four G2 C1 fine-setting parents at absolute T=10 (index40): LOCAL or
CROSS_RING origin, homogeneous H or seeded S. Retain all seven species including
spent Z/W. Multiply concentrations by Omega=1000 and round with numpy.rint
(ties to even), including the fixed original a0/b0 denominators. Do not
renormalize surviving reservoirs. These same actual integer initial states and
denominators apply to every block, numerical setting and B destination.

Keep the active ring N500, h=.2, diffusion (25,25,300), P17 and the prepared
A graph. Switch only B to LOCAL or CROSS_RING, both cubic/unit weighted with
kappa_A=kappa_B=1/3. OFF is omitted to bound pilot cost. H must evolve all500
sites with noise: a homogeneous-subspace reduction is prohibited for sampling.
Origin comparisons involve both prepared state AND fixed A network, since noisy
A is no longer necessarily uniform. They cannot isolate state memory alone.

Reuse frozen G4 simultaneous finite-resource step, with maximum tau=1/24000
and1/48000. Products and arrivals cannot act in the same numerical step. This
approximates continuous time; it does not establish a fundamental physical tick.
Retain quarter-unit observations post0..15 (absolute10..25). Save accumulated
integer reaction/active-hop/environment-edge ledgers. Carry absolute postclock,
RNG and all counters across checkpoint segments; never reset at a checkpoint.

Four independent balanced blocks per setting, eight branches per block
(2origins x H/S x 2destinations):64 stochastic branches total. Every branch and
setting has an independent registered seed, with separate excluded short warmup.
No common-random-number pairing, adaptive replication, favorable selection or
automatic retry. Parents are prepared deterministically: no claim about noisy
formation or survival from birth.

Add16 deterministic controls: eight contexts x continuum/finite-count drift,
all from the exact same rounded initial concentrations and fixed denominators.
Finite-count drift replaces x^2*y by x*max(x-1/Omega,0)*y and adjusts the reaction4
ledger consistently. These are closures, not exact stochastic means. They
separate rounding and elementary finite-count drift from unexplained mean shifts;
mean-minus-continuum must not be called purely a noise effect. BDF rtol1e-11,
atol1e-13,maxstep.05; all500 nodes and closed accounting retained.

## Measurement and outcome rules

Use G3's signed fixed-origin smoothed X-template projection. Late averaging is
trapezoidal over absolute20..25 (post10..15), including both endpoints. For each
block form D=mean_over_origins[(S_cross-H_cross)-(S_local-H_local)]. Independent
settings each have n=4 D values. Report mean, SD and nominal97.5% two-sided
Student t interval (df3); Bonferroni gives nominal95% two-setting family coverage
under the t approximation. This is not a distribution-free bound and n4 is weak.
Do not use G4's .08 endpoint-CDF bound as a projection error bound.

Step agreement: independent-setting fine-minus-coarse90% Welch interval wholly
inside [-.02,+.02] projection units. Its coverage is separate from the primary
CI family; no joint95% coverage for all clauses is claimed. Direction-consistent
pilot evidence requires complete/audited64branches and16controls, both primary
lower bounds >0, and step agreement. Otherwise report unresolved precision,
direction or refinement. Completion is independent of statistical significance.
One missing branch invalidates that entire balanced block; missing values never
become zero. A setting's inferential interval requires all four blocks. Partial
records remain descriptive, and no global direction claim is made if incomplete.

Secondary: per-origin interaction, G/Y projections, means, stocks, uptake,
shape/movement diagnostics, post-switch reference continuity. Post15 cannot
establish a20-unit persistence requirement. Late averaging may follow loss of
reference matching; projection cannot be relabeled lifetime, identity or death.

## Qualification, costs and retention

Before sampling/integration commit/push sources, donor hashes and fixed seeds.
Preparation ceiling900 elapsed seconds, one process/library-thread1. Fixed short
development only: three G4 fixtures x two settings x8paths, post.004 segmented
versus monolithic identity;16 N500 classes to post.01, separate .001 warmups;
16 deterministic controls to post.01. No late-window development outcomes.
Test accounting, same-parent identity, time/seed bounds, interval algebra,
background cancellation, missing-block handling and checkpoint prefix retention.
Measure kernel/IO separately; forecast64 full paths using class measurements,
time factor1500 and safety1.5 plus explicit IO/startup allowance. Deterministic
control cost measured separately. A forecast is extrapolation, not a guarantee.

Candidate production:80jobs (controls first), max8local workers/library-thread1,
24core-hour aggregate ceiling,3600seconds/job, .25second polling and75% forecast
admission reserve; drain active jobs on forecast stop, retain prefixes, no retry.
Registration blocked if forecast exceeds18core-hours or free disk <4GiB.
Full raw/checkpoints, excluded warmups, failures, provenance, seeds, analysis and
same-author separate recomputation proposed for public repository retention.
Uncompressed final scientific arrays about516MB plus metadata/warmups; actual
compression measured during preparation. Publication and scientific run require
the PI's concrete approval. Closure remains a later PI decision.

## Internal challenge

The original deterministic +.0333 is modest; independent noisy H controls may
dominate variance. This pilot may return unresolved despite an accurate kernel.
Time refinement must be tested statistically, not pointwise between independent
noisy paths. Rounding and finite-count drift confound comparison with closed G3.
Four fixed parents do not sample parent uncertainty. Projected long-horizon cost
can change with depletion/adaptation; use caps and preserve incomplete output.
Classical graph transport and correlations supply no quantum-entanglement proof.

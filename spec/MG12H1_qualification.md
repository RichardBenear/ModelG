# MG12-H1 formal small-system qualification — MG107

2026-11-09. PI: “go ahead with the formal qualification.” Register this bounded local qualification before sampling. Existing H0/H1 sources, fixtures and MG106 evidence remain frozen. This is not a particle, larger-substrate, evolving-topology or quantum qualification.

## Fixed design

Four unchanged development fixtures, at most6tokens/3sites,183reachable count states in aggregate. Exact incidence SSA and strict maximum ticks.02,.01,.005; observations0,.05,.1. Finest.005 was fixed in the development design and is the sole qualified approximate method; coarser methods are descriptive, not selected for favorable outcomes.16384independent seeded paths per fixture/method,16jobs,262144scientific paths. No warmups needed by the pure-Python engine. Every path retains count frames and complete eligible accepted/rejected event trace. Identity labels are bookkeeping, counts form the tested probability law.

Two endpoints x4fixtures x2formal methods(exact/finest)=**16formal joint-count distribution comparisons**. Compare empirical categorical law with the independently assembled finite-state continuous-time master equation. Reference state support and order are fixed by the frozen count generator; total variation is invariant to ordering, unlike the development CDF smoke diagnostic.

Observed TV=0.5*sum(abs(empirical_probability-reference_probability)). For a fixture with S reachable states and R=16384, use

epsilon(S)=sqrt((S*ln2+ln(16/.001))/(2R)).

The multinomial L1 concentration bound P(TV(empirical,true)>epsilon)<=2^S exp(-2R epsilon²), followed by a16comparison union bound, gives nominal99.9% simultaneous confidence conditional on fixed sampling/PRNG assumptions and reference numerical fidelity. Endpoints/shared paths need not be independent for that union bound. Some support states have zero probability at an endpoint; retaining all reachable states is conservative.

Each formal criterion is **observed TV+epsilon<=0.10**. This is a coarse ten-percentage-point upper bound on the specified joint-count endpoint law, not pointwise kinetics accuracy or a full temporal/identity-law guarantee. Exact and finest must both satisfy all endpoint criteria. Coarse/medium TV comparisons have descriptive outputs only, with no monotonicity requirement. Completion/audit success is independent of scientific verdict. No adaptive R, step selection, optional stopping, retries or threshold changes after outcomes.

The reference is a separately coded count generator, not an independent physical model. Check sparse-exponential probabilities against dense matrix exponential within1e-11 and analytic feedback-chain probabilities within1e-11; row sums within1e-10, probabilities >=-1e-12. Do not clip/renormalize to hide discrepancies. Confidence applies relative to the computed reference; these numerical checks do not constitute a rigorous floating-point error certificate for the exact real-arithmetic semigroup.

## Audit and operational limits

Preserve all4096development seeds and prior manifests; new seeds use stage10700, deterministic nonce collision resolution against prior/new uint32values before registration. Record every seed and collision resolution. No RNG/integration until registration source/fixture/reference-input hashes are committed/pushed.

Independent saved-raw audit replays Occupies replacements, event roles/catalyst reservations, ancestor/blocker records, observation counts and closed identity totals. Verify seeds/config/source/reference and all denominators; independently derive histograms via sorted/searchsorted counts and recompute TV/concentration/decision. Reject missing/prefix jobs from a qualified verdict. Keep failed/partial files and stop on errors, no automatic continuation/replacement.

One local worker, library threads1,900elapsed-second aggregate ceiling including reference construction, sampling, IO and automatic saved-data analysis/audit. Development4096paths cost4.022072seconds; naive x64 cost257.4seconds; forecast1.5x+30=416.12seconds (~6.94min), extrapolation not guarantee. Check deadline between paths; max100000events/ticks per path, fixtures are finite. Store raw counts atomically in512path prefixes and append flushed per-path trace JSONL; a partial trace/count mismatch is retained and blocks completion. Require1GiB free disk. Full raw/reference/seeds/traces/budget/results retained in existing public repository under the requested qualification workflow; PI closure remains separate. No external review requested.

Formal support means exact and strict engines pass the specified small-system endpoint-law bounds plus accounting/causal tests. It does not qualify G5-scale token populations, physical ticks, 3D space, causal invariance or emergent quantum physics. The reference matcher remains combinatorial; later acceleration needs its own equivalence checks.

# MG12-H1 explicit hypergraph engine development — MG106

We now have a small working engine in which reactions and transport operate by matching and replacing explicit token/site relations. Its continuous-time reaction and transport rates agree with a separately written count model across every reachable state of the four development fixtures. Its separate tick executor enforces the agreed rule that a token—including an unchanged catalyst—can participate only once per tick, and newly produced or arriving tokens wait for the next tick. The bounded development checks passed, giving us a concrete bridge from the earlier count-based implementation to identity-rich relational rewriting. This does not yet establish full-lattice accuracy, particle persistence, evolving geometry, emergent dimension or a quantum vacuum substrate; the tick model's finite-step statistical accuracy still needs formal qualification.

**Status: ENGINE BUILT; BOUNDED DEVELOPMENT CHECKS PASSED. Formal H1 production qualification has not run.**

## What is implemented

- Typed Site/Token entities and Occupies, ActiveLink, ALink and BLink incidence relations.
- Explicit local matches for all six active P17 reaction channels, including two distinct unordered X catalysts and one Y in the feedback event.
- Atomic Occupies replacements with role checks, participant reservations and constant token identities/counts. Z/W stay inactive and immobile.
- A continuous-time exact executor that re-enumerates legal events after each event.
- A separate strict-tick executor that freezes matches, draws per-match exponential proposal marks, accepts proposals in mark-time order if their token participants remain available, then commits disjoint replacements together. Shared read-only links/sites do not themselves create a conflict.
- Accepted/rejected eligible-proposal traces, participant identities, pre/post relation states, step/time, conflict blockers and token-history parents. Proposals outside the tick window do not enter the eligible trace.
- Input guards and bounded event/tick execution. No coordinate-based scheduling or numeric-label priority is used; exact eligible-mark ties are rejected explicitly.

This is a transparent reference implementation, not a general-purpose optimized hypergraph system. Initial transport remains pairwise and the substrate is fixed. The feedback match has three participants; that is not a claim of three-dimensional space. Hyperedge/topology-changing rules remain later stages.

## Frozen development scope

PI instruction: “go ahead with the engine development.” Source, fixture specification and all4096seeds were committed/pushed at **dcca15a** before any development random draws or master-equation integrations. Seeds were disjoint from scanned prior plans. One local worker, numerical-library threads1,300elapsed-second ceiling; no pattern or full-lattice run.

| Fixture | Tokens/sites | Tested mechanism | Reachable count states |
|---|---|---|---:|
| Single active token | 1/2 | G/X/W reactions and active transport | 6 |
| Full local chemistry | 6/1 | A,G,X,X,Y,B with all active reaction channels | 138 |
| Separate resource links | 4/3 | Two A and two B, distinct weighted networks, chemistry off | 36 |
| Feedback only | 4/1 | Two X catalysts, two Y, catalyst reservation | 3 |

There are **183reachable states** in total. Each fixture used256paths for exact continuous time and for strict maximum tick durations.02,.01,.005, to endpoint.1 with observations0,.05,.1. Observation boundaries clip ticks. This is4096development paths, plus16deliberate held replays of already registered first-path seeds.

## Findings and checks

Every reachable count state's explicit matches were grouped by reaction/transport rule and resulting count state, then compared with independently written count propensities and stoichiometry. All183states agreed within the prospectively specified1e-12absolute/relative bound. The reference finite-state master equations conserved probability within1e-10, with no probability below-1e-12. This is stronger generator evidence than merely obtaining similar-looking trajectories, but remains bounded to the tested finite-state domains.

Nine invalid input/proposal/state cases were rejected. Forced examples confirmed catalyst-versus-loss conflict, same-site disjoint event compatibility, one-feedback-event capacity for the two-X/two-Y snapshot, no immediate use of an arriving token, label-equivariant match/rate effects and corresponding-mark replay. All16held simulation replays reproduced counts and traces exactly. Same-seed byte equivalence under arbitrary relabeling is not claimed; distributional equivariance and corresponding match marks are the appropriate comparison.

All4096paths passed causal trace and token checks during development. A later saved-evidence check reread every retained path and reconstructed observation counts from literal saved replacements, checked participant reservations, initial-state roles, parent references and endpoint state encodings, and independently recomputed the empirical CDFs. It replayed **1876 accepted events**. No new simulations were used for this recomputation.

### Endpoint distribution smoke checks

The independent master equation provides the endpoint distribution on a fixed, prospectively defined joint-count-state ordering. Distances below compare empirical versus reference CDFs of that state index. This encoding is a diagnostic: CDF control on an arbitrary ordering is not total-variation control or a full temporal-law guarantee.

| Fixture | Exact | Tick.02 | Tick.01 | Tick.005 |
|---|---:|---:|---:|---:|
| Single active token | 0.017126 | 0.060094 | 0.018031 | 0.002310 |
| Full chemistry | 0.038615 | 0.015759 | 0.027800 | 0.021722 |
| Resource networks | 0.031984 | 0.037253 | 0.029384 | 0.042235 |
| Feedback only | 0.005319 | 0.006904 | 0.001413 | 0.005319 |

The four exact-executor distances all satisfied the frozen DKW smoke bound **0.114262340**, with nominal99% simultaneous coverage under the specified fixed-sampling/PRNG assumptions. That wide bound is a development sanity check, not a claim of sub-percent accuracy.

Strict-tick distances are **descriptive only**. Neither a formal tick accuracy margin nor monotonic improvement was required or inferred. Independent sampling can make a smaller tick's observed distance larger. No step setting or sample count was selected from these outcomes. Formal refinement/accuracy qualification remains separate.

## Performance, limitations and next step

The entire numerical development batch took **4.022072seconds**, below300seconds, including generator/master calculations, path simulations, trace checks and file writing. The largest256path job took about.727seconds. These tiny-system timings are not forecasts for a particle simulation. The explicit matcher enumerates multi-token combinations; feedback-match count grows as choose(nX,2)*nY. Directly expanding G5-scale populations would be impractical. Any aggregation or optimized sampling must be separately shown to preserve the relational generator and strict participant reservations.

The generator oracle shares the mathematical model/fixture constants but calls none of the engine's matcher or commit/executor functions. The saved-evidence checker shares raw data and fixture/hash helpers, with separate replacement replay and sorted/searchsorted CDF arithmetic. These are same-author separate analysis paths, not independent external replication.

No production failures, retries, seed replacements, criterion changes or budget overruns occurred. Existing kernels and closed G5/earlier results remain unchanged. The new finite tick process is not the G4 binomial leap; catalyst reservation changes its finite-step law. No physical tick duration or spatial dimensionality is assigned.

**Recommended next step:** prepare a prospectively registered small-system qualification with a concrete tick-accuracy margin, adequate replication, independent replay/config audits and a measured compute/retention budget. Only after that should we extend substrate connectivity and return to pattern behavior. No formal qualification, higher-dimensional pattern run or evolving-substrate experiment is authorized by this development report.

## Files and reproduction

Contract/specification: spec/MG12H0_contract.md and spec/MG12H1_development.md. Engine: src/mg12h1.py. Independent count/master reference: src/mg12h1_reference.py. Development plan/runner: src/prepare_mg12h1.py. Saved evidence checker: src/check_mg12h1_saved.py.

All source/plan/reference arrays, development path counts/seeds, eligible-event traces, per-job timing/statistical diagnostics and check summaries are retained in development/. Main summary: development/results.json; separate check: saved_evidence_check.json. The development runner rejects an existing attempt rather than automatically rerunning; read-only saved checks can be repeated. Published-byte verification is recorded separately in retention_verified.json; text hashes distinguish local bytes from Git line-ending normalization.

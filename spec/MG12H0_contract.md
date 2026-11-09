# MG12-H0 — relational state and causal contract, version 1 (MG105)

Richard approved all five MG104 recommendations on 2026-11-09. This specifies the first controlled bridge: named species, finite seven-species accounting, non-spatial relational substrate and strict catalyst participation. G5 stays closed/unresolved; its precision follow-up is deferred. No production run or physical clock is authorized here.

The contract makes the hypergraph translation precise enough to implement and challenge. It preserves the existing continuous-time reaction generator as a reference, while explicitly defining a different finite-tick executor to respect the PI's one-participation rule. Agreement between those dynamics is a later measured question, not an assumption. No quantum, physical energy or emergent-dimensional interpretation follows from this specification.

## State

Use disjoint sorts Site and Token. Typed incidence relations are Occupies(token,site,species), ActiveLink(site,site,weight), ALink and BLink. A token has exactly one species in G,X,Y,A,B,Z,W and exactly one occupied site. Token identifiers are bookkeeping only; all laws are equivariant under arbitrary bijective relabeling. Reactions rewrite Occupies relations; transport rewrites the site attachment. Initial H1 topology is fixed. Site/transport relations may be read concurrently; they are not consumed or globally locked merely because a local event reads them.

Higher-arity substrate relations enter H2 with a separately specified law. The three-token R4 match is a genuine multi-participant local pattern; it is not evidence of spatial dimension or irreducible higher-order transport. Geometric coordinates are optional fixture annotations, excluded from event selection. Parameter Omega_count sets count/concentration scaling; W is the published output Omega, not that numerical parameter.

Each rule has a matched input pattern, writable attachments, preserved catalyst relations, output pattern, propensity per physical match and an explicit reservation set. Retain identities through type/attachment changes as a minimal bookkeeping convention; this is not proof of physical particle lineage. Initial identity semantics are exchangeable within each species/site.

## Canonical local reactions

All counts and a0,b0 refer to tokens; Omega_count>0. a0,b0 are fixed original positive reference counts at each site, not surviving supplies. P17 is the initial scope. Unsupported reverse reactions remain off.

| Rule | Local input/output | Distinct matches | Rate per match | Participants reserved in a tick |
|---|---|---:|---:|---|
| R1 | A→G | nA | 14 Omega_count/a0 | the A token |
| R2 | G→X | nG | 1 | the G token |
| R2r | X→G | nX | .1 | the X token |
| R3 | B+X→Y+Z | nB nX | 29/b0 | B and X |
| R4 | X+X+Y→X+X+X | choose(nX,2)nY | 2/Omega_count² | both catalyst X tokens AND Y |
| R5 | X→W | nX | .9 | the X token |

For R3, the X token becomes Y and B becomes Z. For R4 the two distinct catalyst X tokens retain their species/attachment, while Y becomes X; all three are unavailable for the rest of that tick. The two X roles are unordered. Enumerating ordered roles without adjusting the factor would double the generator. R4 total rate is exactly nX(nX−1)nY/Omega_count², matching the existing finite-count convention. It is zero at nX<2; the same token can never fill both catalyst roles.

The zero-rate G→A reverse rule can be represented in a future general-parameter extension; it is not silently activated here. These schemas conserve total token count over G+X+Y+A+B+Z+W. They do not assign physical energy or recycling. Z and W remain inactive and immobile in this translation.

## Pairwise transport reference

One directed match is an occupied token plus one permitted oriented link. Active G/X/Y per-link rates reproduce the reference fixture operator, including its chosen h and diffusion normalization. A/B move on their respective link relations with rate kappa_species*weight. A symmetric undirected physical link is represented by two oriented opportunities, not two duplicate copies per direction. No self links or duplicate directed links in the initial domain. Transport preserves token species and changes only Occupies attachment. Later degree/volume normalization is a separately declared H2 model choice, not inferred from arbitrary embedding.

## Two clocks, kept separate

**Reference continuous time:** enumerate current distinct matches; select an exponential waiting time using total rate and a match with rate-proportional probability; rewrite atomically and re-enumerate. A catalyst can participate again in a later event at a later time. No simultaneous reuse occurs within a single event. This is the old finite-count generator expressed on explicit identities.

**New strict tick, proposed H1 executor:**

1. Enumerate matches and rates from the complete tick-start state, independent of numeric labels.
2. Give each match one independent exponential proposal time with its frozen rate; retain proposals strictly before tick duration delta.
3. Process proposals by time. Accept a proposal only if all its participating tokens remain unreserved; otherwise reject it without replacement or resampling. Reserve catalysts as well as consumed/moved tokens.
4. Commit accepted, participant-disjoint Occupies replacements together. Products, changed species and arrivals cannot generate new matches until the next tick. Preserved read-only site/link relations may be shared.
5. Clear reservations at the tick boundary. Keep accepted and rejected proposals, participants, causal dependencies and timing metadata for diagnostic fixtures. Reject exact proposal-time ties in the reference implementation rather than resolve them by identifier; document the finite-precision tie condition. Optimizations require separate equivalence checks.

This is a specified stochastic greedy selection of a matching; it is not a uniform choice among maximal compatible sets. It is not G4's aggregate binomial leap and no byte-equivalence to G4 is claimed. For a fixed finite state, a specified match is accepted with rate*delta+O(delta²), and two overlapping proposals occur with O(delta²) probability. Thus its first-order infinitesimal generator is the reference generator, subject to a bounded-domain qualification. Finite-delta distributions differ. Uniform convergence as counts, graph size or rates grow is not established by that local argument.

The catalyst restriction imposes a useful capacity bound: within a tick R4 firings cannot exceed min(nY,floor(nX/2)). Every R4 firing also reserves its original Y, and newly created X cannot catalyze another event. Competing X transport/conversion reduces this capacity further. This can alter patterns at finite delta; H1 must quantify it. No physical duration or Planck-scale mapping is assigned to delta.

## Required H0 checks and H1 qualification matrix

H0 static checks use explicit tiny populations only: rule-match totals equal count propensities; each event gives the expected seven-species increment and constant total; R4 catalysts are distinct and reserved; non-overlapping events commute; overlaps conflict even for unchanged catalysts; read-only site sharing is allowed; no same-tick arrival/product reuse; bijective relabeling preserves match sets/rates/effects after mapping labels back. A failing contract check blocks engine development until corrected prospectively.

H1 must add actual incidence matching and replay, independent count-generator comparisons, exact small-state/statistical laws, catalyst saturation tests, maintained/closed limits as applicable, seeded stream/label tests, clock refinement and bounded growth/error handling. Same-seed byte equality under arbitrary relabeling is not assumed: the distribution must be equivariant, while replay requires a correspondence of match identities and marks. Test the transition law, not just one chosen enumeration order.

## Development and evidence boundary

Authorized immediate work: this contract and deterministic static schema checks only; no random draws, time integration or pattern runs. The static check fixture has at most seven tokens and uses exhaustive tiny match enumeration; no scientific inference or compute allocation is hidden in it. Preserve source/contract and register the fixture description before running those checks. Produce an H0 result with a plain-language summary and a separate algebraic check path.

Next build proposal: tiny exact and strict-tick reference engine, maximum seven tokens/three sites in initial stochastic fixtures, with a prospective seed list, marginal/joint observables, refinement bounds and a measured timing plan. Reserve at most300seconds single-worker/thread1 for a proposed H1 development batch; this is a draft allocation to concretize before numerical work, not approval of scientific production. Do not choose final sample counts or claim a production forecast until the explicit enumerator's scaling has been measured. Later local production remains max8workers with per-run PI-approved criteria/budget/retention. External review stays optional at the PI's discretion.

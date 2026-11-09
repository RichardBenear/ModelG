# MG12 — explicit relational hypergraph mapping and substrate emergence

Draft 1, MG104, 2026-11-09. Planning only; no implementation experiment or production run registered.

The aim is to express Model G as local changes to explicit relational structures, first demonstrating that the translation preserves the known kinetics, then testing what happens when the substrate connectivity evolves without imposed spatial coordinates. This separates faithful translation from discovery: reproducing Model G after programming its reactions is an implementation result; obtaining its effective behavior from simpler substrate rules is a further scientific result. Initial dimensions need not represent physical space. Neither emergent physical space nor a quantum vacuum interpretation is assumed.

## Starting point and correction to the recent discussion

The project is already rule-based and has prior evolving-hypergraph studies. MG5 used a two-dimensional periodic Delaunay triangulation: triangles were ternary hyperedges, binary edges carried diffusion, and birth/death events changed the mesh. MG6–MG7 added coupling investigations. The coordinates, triangulation machinery and interpolating field transfers supplied geometric structure; this was not a coordinate-free derivation of space or an explicit seven-species integer substrate. MG8 examined scheduler sensitivity on the one-dimensional lattice and expressly did not establish causal invariance. MG10–MG11 introduced explicit finite resources; G5 used a fixed ring plus resource graphs and closed as statistically unresolved.

Thus the new program combines prior lessons while removing particular scaffolds. Earlier statements that the project had only demonstrated a one-dimensional lattice, or had never evolved its connectivity, were too broad: those descriptions apply to the recent G5 configuration, not the entire investigation. Closed historical outcomes remain unchanged.

References: MG5_proposal.md §§1,3; MG8_proposal.md §§1,2; rules_v1.1.md; MG11G5_pilot.md; ../results/MG11G5/MG11G5_results.md. Existing reference implementations/data are read-only donors for this work.

## Meaning of hypergraph and dimension

A typed incidence hypergraph records entities, relations and their participating roles. Reactions must actually match and replace these relations; drawing hyperedges over an unchanged array solver would not satisfy the intended translation. Numeric entity identifiers have no physical role. Coordinates can be used in a reference fixture, but cannot silently select updates in relational experiments.

Hyperedge arity, species count and spatial dimension are different quantities. A three-participant reaction is not three-dimensional space. Arbitrary hypergraphs may have no useful single dimension. The program will distinguish imposed test-lattice dimension, hyperedge structure and measured graph-distance/diffusion dimension. A simplicial complex is an optional additional structure, not something implied by every hypergraph. No target dimension is selected from favorable outcomes.

## Six stages with separate decision gates

| Stage | Main question | Required evidence before proceeding |
|---|---|---|
| H0: semantics and causal contract | What are entities, species, reaction roles, transport and ticks? | Written canonical schemas, event examples, conservation identities, catalyst/conflict rules and correspondence to the old model; PI choices settled |
| H1: faithful small-system engine | Does explicit hypergraph rewriting reproduce the existing stochastic reaction/transport process on a fixed substrate? | Independent enumerator/reference checks, matched propensities and event effects, exact accounting, relabeling and causal-conflict tests, bounded statistical comparisons |
| H2: fixed higher-connectivity substrates | Is the mapping correct beyond the ring, including genuinely multi-incident transport relations? | Analytic transport/chemistry limits and conservative operators on small non-geometric hypergraphs, plus 2D/3D benchmark fixtures; normalization fixed prospectively |
| H3: pattern bridge | Do known Model G patterns and geometry responses survive the translation? | Deterministic limits/refinement and prespecified noisy precision tests; no copying G5's underpowered sample size; first ring regression, then one justified additional substrate |
| H4: evolving relational substrate | Can reactions and patterns coexist with local connectivity rewriting without coordinate-based remeshing? | Event-level transfer/accounting/causality tests, matched fixed/rewired controls, bounded graph growth and stable interpretable diagnostics; failures retained |
| H5: emergence and supply | Can a substrate law generate effective transport/resource behavior instead of merely being assigned Model G kinetics? | Coarse-grained coefficients predicted on training/calibration conditions and tested on held-out perturbations; explicit input/output/backreaction accounting |

H0–H3 establish translation and its tested domain. H4 tests evolving connectivity. H5 is a more ambitious model-discovery stage and may require a separate program specification. Success at any one stage does not automatically establish the later ones. Stage counts describe research decisions, not a promise of six runs or success in all six.

### H0 — first concrete deliverable

Specify a minimal substrate relation type, attached species-state entities and explicit reaction-role patterns. Initially retain G,X,Y,A,B,Z,W labels so that correspondence is testable. Track the numerical size parameter as Omega_count; W denotes the published output species Omega, avoiding notation ambiguity.

Freeze all seven-species event increments and propensity conventions, especially the distinction between the published 2X+Y→3X reaction and the existing implementation's finite-count catalyst convention x(x-1)y/Omega_count². Symmetry factors and ordered/unordered match multiplicities must be derived; naively counting all labeled hypergraph matches can change rates. A later alternative convention is a new model, not a silent translation correction.

Small test fixtures: isolated chemistry; two-site transport; scarce B/X competition; feedback catalysts; distinct A/B connections; disjoint and overlapping rewrite matches. These can be specified without sampling. Keep exact identity-rich semantics as the reference even if a later optimized engine aggregates exchangeable tokens; qualify any such aggregation against the explicit engine.

### H1 — clocks, causality and scheduling

Use continuous-time exact stochastic dynamics as the existing-law reference. Separately specify the proposed tick executor: matches read a common pre-tick state, consumed/moved participants cannot be allocated twice, and newly produced/arriving entities cannot participate again until the next tick. Assign unique event/ancestry records and test adversarial arrival/competition cases.

The unresolved modeling choice is whether a catalyst also becomes unavailable after participating once in a tick. That can change the finite-tick law, even when it approaches the old generator as tick size shrinks. A strict tick is not automatically the same stochastic process as continuous time. Compare them under a registered refinement protocol rather than asserting equivalence. Numerical ticks are not mapped to seconds or Planck time.

Conflict selection probabilities and any reaction/transport splitting must be explicit. Test relabeling invariance; test commuting independent events and separately characterize conflicting choices. Do not claim global confluence, causal invariance, Lorentz invariance or physical locality from these bounded checks. Enumerate tiny multiway branches where feasible; a sampled path is not the full multiway system.

### H2 — more than pairwise adjacency

Begin with a ring fixture whose relation structure exactly encodes the current operator. Add small incidence hypergraphs with no coordinate embedding. A higher-arity transport relation requires its own event law: for example, an occupied participant selects another incident participant with explicitly normalized rates. Such transport can sometimes be represented by a weighted pairwise graph; if so, disclose that equivalence. Do not call pairwise-equivalent transport evidence of irreducible higher-order dynamics. Cooperative events requiring several participants simultaneously can provide a separate higher-order test, but changing the kinetics to include them must be explicit.

Use regular 2D and 3D structures as known-dimension benchmark fixtures, not claims that dimension has emerged. Account for cell volume/site capacity and directed incidence degree, so an apparent effect is not accidentally increased total exit rate or altered local concentration. Check constant-field invariance, conserved transfer, equilibrium behavior and mode decay where analytic solutions exist.

### H3 — recovering the known behavior

First recover known ring mean-field kinetics/transport and pattern observables under fixed parameters. Account for the narrow existence window and the historical MG4 change in behavior on irregular graphs; do not demand identical morphology across different transport operators without calibration.

For a noisy response test, use the G5 signed/H-subtracted estimand as a possible bridge, but recognize that G5 did not establish stochastic robustness or step agreement. Set effect precision and equivalence margins prospectively, measure cost and choose replication before production. A deterministic success is not a stochastic success. Prepared and spontaneously formed structures are separate claims.

### H4 — changing connectivity without hidden geometry

Specify a minimal finite set of local relation-replacement rules and their admissible domain before selection experiments. Start field-blind; introduce field-dependent topology only after this baseline is characterized. Bound degree, entity growth and event cost. Define where attached populations go when a site/relation is replaced: identity-preserving attachment or explicit local conservative transfer, never unrecorded interpolation, deletion or duplication. Fix the volume/reference-denominator convention when sites change.

Controls include fixed topology, field-blind rewrites, rate/degree-matched rewiring and relabeling. Measure connectivity components, neighborhood growth versus graph distance, return/diffusion behavior and pattern response over prespecified scales. Report inconsistent or scale-dependent dimensions rather than assigning one. Random rewiring is a chosen model, not an explanation of vacuum geometry.

### H5 — derive supply and outputs from a common substrate

Keep the initial H1–H4 translation seven-species closed by default: finite A/B, retained Z/W, no unmodeled regeneration. A maintained/open bath can be a separately metered control. Do not add Z/W→A/B simply to maintain the pattern.

For the shared-substrate hypothesis, specify what distinguishes available and returned configurations, what drives regeneration, and what changes in the environment. Token-number conservation is not yet physical energy conservation. A passive equilibrium substrate cannot be assumed to furnish an unlimited nonequilibrium supply. If energy, entropy or work are claimed, define and account for them explicitly.

The stronger emergence test is to obtain effective supply/depletion and transport laws from simpler substrate transitions, predict their changes under an intervention, and test held-out cases. Reproducing a reaction law programmed into the rewrite grammar establishes translation, not its microscopic origin. Quantum amplitudes, coherence, entanglement and physical vacuum identification remain later hypotheses; classical correlated tokens are not by themselves a quantum field.

## Reuse and implementation boundaries

Reuse the existing Model G kinetics, closed resource ledgers, G4 exact/leap references, G5 raw analysis and prior MG4–MG8 design lessons. Inspect CRR's relabeling/scheduler lessons and Stochastic Vacuum provenance before adapting any code; neither is assumed to transfer unchanged. Keep work in the ModelG repository with new stage-specific modules/specs; do not rewrite closed kernels or merge entire historical projects. Prefer a small transparent reference engine over an unqualified optimized general-purpose engine.

No experiment has a cost estimate yet. H0 must specify pilot fixtures and a capped development/timing plan before numerical work; each later stage needs measured local forecasts, sample counts, criteria and a concrete approval packet. Use at most eight local workers/library threads1, retain negative/partial outcomes, publish enough raw evidence for recomputation according to approved retention scope, and require PI closure. External review is Richard's per-run choice. Begin every result report with the G5-style plain-language contribution/limits paragraph.

## Questions for Richard

1. **Direction:** adopt MG12 as a separate six-stage hypergraph program, keeping G5 closed and a more precise G5 follow-up deferred? Recommended: yes; preserve G5 as the unresolved bridge baseline.
2. **Meaning of tokens:** initially retain named species as explicit states/entities to validate the translation, then investigate whether those distinctions can emerge from substrate structure? Recommended: yes; requiring unlabeled structural species from the outset would skip the controlled bridge.
3. **Tick participation:** should the one-use-per-tick restriction include catalysts as well as consumed/moved entities? Recommended provisional choice: yes for the new tick model, with any difference from the existing reaction law disclosed and tested; keep continuous-time kinetics as a separate reference. Products/arrivals wait until the next tick in either tick option.
4. **Dimensional goal:** start with arbitrary relational connectivity and measured dimension, while using 2D/3D only as validation fixtures, rather than requiring emergent three-dimensional space immediately? Recommended: yes. No Planck-scale calibration is imposed.
5. **Resource/quantum scope:** retain finite A/B and spent Z/W during translation, reserving shared-substrate recycling and any explicitly quantum extension for H5 or later? Recommended: yes; this makes changes in kinetics, geometry and supply separately interpretable.

These are scientific direction questions, not a request to approve unmeasured production budgets. Answers determine H0; concrete later run approvals remain separate. No work beyond planning is implied by this draft.

## MG105 — PI decisions recorded, 2026-11-09
Richard: “I agree to your recommendations on all 5. Proceed.” All five recommended direction choices are adopted. H0 contract: MG12H0_contract.md; initial static evidence: ../results/MG12H0/MG12H0_specification_results.md. Original questions above are retained as the prospective decision record. This approval enables the staged program; it does not supply an unmeasured production budget or change closed G5 scoring.

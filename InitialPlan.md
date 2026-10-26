Here's the plan. I'd suggest flipping the usual setup:

- **Me:** computation and the running documentation.
- **ChatGPT:** independent replicator and adversarial reviewer.
- **You:** PI, signing off on each gate.

That keeps what made the other projects trustworthy, a second system checking the first, while cutting the round-trips.

## Practical setup

- **The Project.** You'd create it in claude.ai and add the charter (below) as a Project file. I can't create Projects myself.
- **Persistence.** My workspace doesn't carry between chats, so code and results need a permanent home. I'd recommend a GitHub repo you own. I can attach it each session, commit code, configs and results, and ChatGPT can read the same repo for replication. The fallback is Project files plus a running Doc.
- **Conventions.** I'd mirror your CRR style: numbered stages (MG001…) and named pass/fail gates. Each gate gets its prediction and null control written down before the run.

## Phased roadmap

**Phase 0: Charter and specification**
- Pin down exactly which Model G we're simulating. That means LaViolette's published equations and parameter sets, including the reversible-step variants and the parameters from the Pulver–LaViolette stationary soliton work. I'll verify these against the papers rather than rely on memory.
- Freeze the rewrite-rule spec from my last message as version 1.0.

**Phase 1: Continuum baseline**
- Integrate the Model G PDEs in 1D, 2D, then 3D.
- **Gate MG1:** reproduce a published stationary soliton: the G-well and bipolar X–Y profile, its stability, and its wavelength.
- Without this baseline, nothing downstream is interpretable.

**Phase 2: Stochastic token rules on a fixed lattice**
- Implement R1–R5 plus species-specific hopping. I'd use exact Gillespie for small systems and tau-leaping for large ones.
- **Gate MG2:** quantitative convergence to the MG1 PDE solution as Ω grows (the Kurtz limit). This proves the mapping is real.
- **Gate MG3:** at finite Ω, particles nucleate spontaneously from fluctuations. Measure how the nucleation rate scales with Ω and with distance from the bifurcation.

**Phase 3: Non-lattice frozen graphs**
- Run on random geometric graphs and on frozen snapshots of Wolfram hypergraphs with an effective dimension near 3.
- **Gate MG4:** particles survive irregular connectivity, and the graph Laplacian behaves like ∇² at the soliton's scale.

**Phase 4: Evolving hypergraph, no coupling**
- Tokens ride on a hypergraph that is being rewritten by a spatial rule.
- **Gate MG5:** solitons persist, drift, or dissolve. We need to know which, and why, before adding coupling.

**Phase 5: G-coupled spatial growth (the novel part)**
- Let the local G count modulate the firing weight of the spatial rule.
- **Gate MG6:** measure the volume-growth (curvature) deficit around a particle. Check its falloff exponent and whether it scales with the integrated G-deficit.
- Null controls: a decoupled run, and a run with randomized coupling.

**Phase 6: Multiway and frame diagnostics**
- **Gate MG7:** quantify how much branches diverge across update orders. Test whether SQK's preferred frame shows up as a measurable foliation dependence.
- This turns the "non-confluence = absolute ether" observation into a result.

**Phase 7: Write-up**
- A running white paper as a Doc, kept current at each gate, condensed at the end. The format would follow your CRR summary plan.

## Honest constraints

- Large 3D stochastic runs will be slow in my environment.
- I'll push 1D and 2D as far as possible first, and use 3D only for confirmation.
- Phase 5 depends on a modeling choice, the coupling function, so I'll propose candidates and you pick.

A few decisions from you would let me start Phase 0:
1. A project name.
2. GitHub repo or Project files.
3. Whether you're happy with the flipped role split.
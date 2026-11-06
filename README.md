# ModelG: Subquantum Kinetics × Wolfram-Style Rewriting

Research project mapping Paul LaViolette's **Model G** reaction–diffusion system
(Subquantum Kinetics, SQK) onto **stochastic rewriting rules on a Wolfram-style
spatial hypergraph**, then testing whether the hybrid does anything neither
framework does alone.

## Central questions

1. **Is the mapping exact in the mean-field limit?** Do labeled-token rewrite
   rules on a graph converge to the Model G PDEs as tokens per node (Ω) → ∞?
2. **Does discreteness give SQK's particle nucleation for free?** At finite Ω,
   do stable dissipative solitons ("particles") nucleate from intrinsic
   fluctuations?
3. **Do particles survive emergent space?** On irregular and evolving
   hypergraphs instead of a fixed lattice.
4. **Can SQK gravity generate Wolfram gravity?** If local G concentration
   modulates the spatial rewrite rate, does a particle's G-well produce a
   volume-growth (curvature) deficit, and how does it scale and fall off?
5. **Is SQK's preferred frame a non-confluence effect?** Quantify multiway
   branch divergence and foliation dependence.

## Roles and review

Effective 2026-11-06 (MG039), Richard Benear is PI and ChatGPT / Codex is the research and implementation lead: planning, code, validation, runs, analysis and documentation. Claude or another external reviewer participates only when requested.

The authoritative active workflow is [RESEARCH_PROTOCOL.md](RESEARCH_PROTOCOL.md). It requires prospective registration, implementation validation, separate recomputation and PI sign-off; it does not require a second AI at every stage. Same-author checks are not independent external replication. Historical reviews and frozen scientific records retain their original authorship and status.

## Conventions

- **Stages** are numbered `MG000`, `MG001`, … and logged in [`STAGE_LOG.md`](STAGE_LOG.md).
- **Gates** (`MG1`, `MG2`, …) are defined in [`GATES.md`](GATES.md). Each gate's
  prediction, pass criterion and null control are written down and committed
  **before** the run that tests it.
- Every result in `results/` is reproducible from a config in `configs/` and
  a commit hash. No hand-edited numbers.
- Specs are versioned. Changing a frozen spec means a new version file, not an edit.

## Glossary and assumptions

See [`spec/GLOSSARY_AND_ASSUMPTIONS.md`](spec/GLOSSARY_AND_ASSUMPTIONS.md). It defines every symbol
(including the two different meanings of Ω), explains how tokens map to LaViolette's etherons, lists the
Ω values used, and keeps the register of modelling, numerical, statistical and interpretive assumptions
with their status.

## Layout

```
spec/        Frozen specifications (equations, rewrite rules)
src/         Simulation code
configs/     Run configurations (one per run or sweep)
results/     Outputs, organized by stage: results/MGxxx/
review/      Internal checks, optional external reviews and replication artifacts
paper/       White-paper drafts
```

## Roadmap

The original roadmap is preserved in [`InitialPlan.md`](InitialPlan.md). The active extension is [MG10](spec/MG10_proposal.md), explicit A/B reservoirs. Its first registered deterministic batch completed: **UNRESOLVED**, PI closure pending. Both highest capacities recover the baseline, with a low-capacity numerical/measurement issue still open. [Results](results/MG10/MG10_results.md). Summary:

| Phase | Content | Gate |
|---|---|---|
| 0 | Pin down published Model G equations and parameters; freeze rule spec | — |
| 1 | Continuum PDE baseline (1D → 2D → 3D) | MG1 |
| 2 | Stochastic token rules on a fixed lattice | MG2, MG3 |
| 3 | Frozen irregular graphs and hypergraph snapshots | MG4 |
| 4 | Evolving hypergraph, no coupling | MG5 |
| 5 | G-coupled spatial growth | MG6 |
| 6 | Multiway / foliation diagnostics | MG8 (MG7 is the mass-scaling b-scan) |
| 7 | Write-up | — |
| Extension | Explicit A/B reservoirs, recovery of Model G and finite-supply effects | MG10 (first batch UNRESOLVED) |

## Setup

```
pip install -r requirements.txt
```

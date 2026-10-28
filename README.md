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

## Roles

| Role | Who | Responsibility |
|---|---|---|
| PI | Richard Benear | Direction, modeling choices, gate sign-off |
| Computation & documentation | Claude | Code, runs, stage log, white paper |
| Independent replication & adversarial review | ChatGPT | Re-derives gate results from the spec alone; critiques |

### Review protocol (Claude ⇄ ChatGPT)

1. Claude commits a stage and names the files to review in `STAGE_LOG.md`.
2. ChatGPT reviews the repo read-only. The PI saves its review **verbatim** as
   `review/MGxxx_chatgpt.md` (or pastes it to Claude, who saves it verbatim).
3. Claude answers in `review/MGxxx_response.md` with a disposition for each numbered item:
   **Accepted / Accepted with modification / Rejected**, with reasoning and the commit that
   addresses it. Corrections are logged as a lettered sub-stage (e.g. MG001a).
4. ChatGPT re-checks against the response file and the diff. Rounds continue as
   `MGxxx_chatgpt_r2.md` / `MGxxx_response_r2.md` until it has no open items, or the PI rules on
   a disagreement.
5. No gate run starts while review items on its spec are open.
6. Independent replication code and data from ChatGPT go in `review/MGxxx_replication/`
   (committed by the PI from VS Code). A replication that lives only on a local machine can't be
   re-run by anyone else.

## Conventions

- **Stages** are numbered `MG000`, `MG001`, … and logged in [`STAGE_LOG.md`](STAGE_LOG.md).
- **Gates** (`MG1`–`MG7`) are defined in [`GATES.md`](GATES.md). Each gate's
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
review/      ChatGPT replication reports and critiques
paper/       White-paper drafts
```

## Roadmap

See [`InitialPlan.md`](InitialPlan.md). Summary:

| Phase | Content | Gate |
|---|---|---|
| 0 | Pin down published Model G equations and parameters; freeze rule spec | — |
| 1 | Continuum PDE baseline (1D → 2D → 3D) | MG1 |
| 2 | Stochastic token rules on a fixed lattice | MG2, MG3 |
| 3 | Frozen irregular graphs and hypergraph snapshots | MG4 |
| 4 | Evolving hypergraph, no coupling | MG5 |
| 5 | G-coupled spatial growth | MG6 |
| 6 | Multiway / foliation diagnostics | MG7 |
| 7 | Write-up | — |

## Setup

```
pip install -r requirements.txt
```

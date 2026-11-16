# ModelG: Subquantum Kinetics × Wolfram-Style Rewriting

Research project mapping Paul LaViolette's **Model G** reaction–diffusion system (Subquantum Kinetics, SQK) onto
**stochastic rewriting rules on a Wolfram-style spatial hypergraph**, then testing whether the hybrid does anything
neither framework does alone. The long-term aim is to show, with stated statistical confidence, that a persistent Model G
particle exists in a discrete rewriting system, and to test whether such a system can connect to gravity and quantum
mechanics. That aim is why causal invariance and SQK gravity are now active questions.

## Where the study stands (2026-11-16)

In short: Model G's particle is reproduced and qualified in the continuum and on a fixed 1D lattice, and it survives
on evolving 2D graphs. On a stochastic lattice it persists reliably only when each site holds thousands of tokens
(Ω ≈ 3000 or more). The explicit, tick-based hypergraph engine reproduces the continuum on small systems. However, it
has so far been affordable only at Ω ≈ 10, where noise destroys particles, so persistence on that engine is not yet
shown. Causal invariance in the strict Wolfram sense looks unavailable for chemistry-like rules. SQK's gravity drift is
real in the continuum equations, but it is not specific to G.

| Question | Status | Key studies |
| --- | --- | --- |
| 1. Mean-field limit exact? | **Yes, within tested scopes.** Lattice and hypergraph engines converge to the Model G equations; tick accuracy is equivalent at ticks of 1e-4 and 5e-5. | MG1, Phase 2 validations, MG12-H1 to H3 |
| 2. Does discreteness nucleate particles? | **Noise nucleates and destroys them.** Nucleation rate falls as exp(−0.064 Ω). Below Ω ≈ 1000, noise replicates or kills particles. At Ω = 3000 and 10000, 29/30 and 30/30 persist to T = 1000. Lifetime grows about 9.4 times per +1000 Ω. | MG2, MG3, MG3b, MG185, MG188 |
| 3. Do particles survive emergent space? | **Yes on evolving graphs, no on frozen irregular graphs at b = 29.** Frozen Delaunay graphs multiply particles (their window is b ≈ 27.3–28.7). Birth–death rewriting keeps one particle. | MG4, MG4c, MG5 |
| 4. Can SQK gravity drive Wolfram gravity? | **Partly.** A G-coupled growth rule gives a local volume deficit around a particle (falloff unresolved). In the continuum a particle drifts toward lower G with mobility −0.0598, matching Pulver & LaViolette's −0.0599, but an X-source gradient pushes exactly as hard. A particle has no net G monopole (zero G integral). | MG6, MG7, MG190–MG192 |
| 5. Is the preferred frame a non-confluence effect? | **Update order matters; strict causal invariance fails at count level.** Sweep order can dissolve a particle. Independent events commute, but count-level causal histories are not unique. | MG8, MG189 |

### Main findings by area

- **Continuum (MG1, MG003, MG192, MG193).** The published soliton is reproduced. In 1D, b = 29 sits close to the
  lower edge of the soliton's existence window, which makes 1D particles fragile there. MG193 is mapping that window
  in (a, b) to choose a robust 1D operating point.
- **Fixed-lattice stochastic engine (Phase 2, MG2–MG3b, MG9, MG185–MG188).** Particle loss is driven by noise and is
  controlled by Ω. G and X fade together when a particle is lost. b behaves as the control parameter of a subcritical
  instability (MG9).
- **2D graphs (MG4–MG7).** Frozen irregular graphs at b = 29 multiply particles. Evolving birth–death graphs keep one
  particle. Coupling graph growth to G gives a measurable node deficit around the particle. Mass scaling across b was
  not scored.
- **Finite reservoirs and environment networks (MG10, MG11).** Explicit A/B stocks recover Model G at high capacity,
  and the samplers are qualified. Supply location and connectivity change the response measurably.
- **Explicit hypergraph engine (MG12, MG104 onward).** A strict-tick engine with individually tracked tokens is
  qualified against the continuum on small systems. It is about 5e4 times slower per time unit than the tau-leap
  engine, which forced Ω ≈ 10 and short horizons. At that size, pattern retention was not observed (MG160-P1), and a
  long series of detector and tracker screens (MG163–MG181) mostly failed. G visibility improves with population
  (MG183).
- **Causal invariance (MG189).** Tested independent events commute. Count-level causal-history uniqueness fails, with
  a minimal witness. The strict tick converges toward the continuous-time result as the tick shrinks.
- **SQK gravity (MG190–MG193).** LaViolette's sign is confirmed: a particle moves toward lower G. On the coarse lattice,
  drift is downhill but pinned to the mesh (MG191). The spectral continuum gives the exact mobility and matches the
  published value (MG192), but the drift also follows X-source gradients.

### Active and next

- **MG193 (running):** the continuum 1D existence and stability window in (a, b), and a recommended operating point.
- **Candidates, not yet authorized:**
  - an engine bridge comparing loss-time distributions between the Phase 2 and hypergraph engines at Ω ≈ 100–300;
  - a repeat of the gravity dynamic check inside the window;
  - the two-particle force against separation;
  - the far-field G monopole in 2D and 3D.
- **Closures pending PI sign-off:** MG183, MG185, MG187, MG188 and MG191. Check `GATES.md` for current status.

## Central questions (charter)

1. **Is the mapping exact in the mean-field limit?** Do labeled-token rewrite rules on a graph converge to the Model G
   PDEs as tokens per node (Ω) → ∞?
2. **Does discreteness give SQK's particle nucleation for free?** At finite Ω, do stable dissipative solitons
   ("particles") nucleate from intrinsic fluctuations?
3. **Do particles survive emergent space?** On irregular and evolving hypergraphs instead of a fixed lattice.
4. **Can SQK gravity generate Wolfram gravity?** If local G concentration modulates the spatial rewrite rate, does a
   particle's G-well produce a volume-growth (curvature) deficit, and how does it scale and fall off?
5. **Is SQK's preferred frame a non-confluence effect?** Quantify multiway branch divergence and foliation dependence.

## Roles and protocol

Richard Benear is PI. ChatGPT / Codex and Claude are AI collaborators, and either can take any role the PI assigns per
study: direction suggestions, design, implementation, runs, analysis or review. Closure sign-off and compute approval
belong to the PI alone. Work runs on two tracks: **exploratory** (quick, logged in [`explore/LOG.md`](explore/LOG.md),
cannot support claims alone) and **confirmatory** (registered gates). The authoritative rules are in
[RESEARCH_PROTOCOL.md](RESEARCH_PROTOCOL.md), version 2 (2026-11-16). Same-author checks are not independent external
replication. Historical reviews and frozen records keep their original authorship and status.

## Conventions

- **Records** are numbered `MG000`, `MG001`, … and logged in [`STAGE_LOG.md`](STAGE_LOG.md). A short study-by-study
  overview is in [`stage_log_summary.md`](stage_log_summary.md). Current state and next steps are in
  [`HANDOFF.md`](HANDOFF.md).
- **Gates** for confirmatory studies are defined in [`GATES.md`](GATES.md). Each gate's prediction, pass criterion and
  null control are committed **before** the run that tests it. Registered results are never rescored.
- Every result in `results/` is reproducible from committed code, configuration and a commit hash. No hand-edited
  numbers. Large raw data stay on the PI's PC with hash manifests published here.
- Specs are versioned. Changing a frozen spec means a new version file, not an edit.

## Glossary and assumptions

See [`spec/GLOSSARY_AND_ASSUMPTIONS.md`](spec/GLOSSARY_AND_ASSUMPTIONS.md). It defines every symbol (including the two
different meanings of Ω), explains how tokens map to LaViolette's etherons, lists the Ω values used, and keeps the
register of modelling, numerical, statistical and interpretive assumptions with their status.

## Layout

```
spec/        Frozen specifications (equations, rewrite rules, study specs)
src/         Simulation and analysis code
configs/     Run configurations
results/     Outputs and reports, one folder per study
review/      Internal checks, reviews and replication artifacts
explore/     Exploratory-track log
paper/       White-paper drafts and all figures (paper/figures/)
tests/       Unit and regression tests
```

## Roadmap

The original roadmap is preserved in [`InitialPlan.md`](InitialPlan.md).

| Phase | Content | Gate(s) | Outcome |
| --- | --- | --- | --- |
| 0 | Published Model G equations and parameters; frozen rule spec | — | Done |
| 1 | Continuum baseline | MG1 | PASS (PI exception) |
| 2 | Stochastic token rules on a fixed lattice | MG2, MG3, MG3b, MG9 | Noise-limited; barrier S = 0.064 |
| 3 | Frozen irregular graphs | MG4, MG4c | FAIL at b = 29; window b ≈ 27.3–28.7 |
| 4 | Evolving hypergraph, no coupling | MG5 | PASS |
| 5 | G-coupled spatial growth and mass scaling | MG6, MG7 | Deficit DEMONSTRATED; scaling not scored |
| 6 | Multiway / update order | MG8, MG189 | Order-SENSITIVE; count-level causal uniqueness fails |
| Ext. | Explicit A/B reservoirs and environment networks | MG10, MG11 | Recovers Model G; supply effects measured |
| Ext. | Explicit strict-tick hypergraph engine | MG12 (MG104 onward) | Qualified on small systems; persistence not shown at Ω ≈ 10 |
| Ext. | Persistence versus Ω on the fast engine | MG185–MG188 | Persists at Ω ≥ 3000 to T = 1000 |
| Ext. | SQK gravity | MG190–MG193 | Sign and mobility confirmed in continuum; not G-specific |
| 7 | Write-up | — | Early (plain-language summary of MG1–MG3b only) |

## Setup

```
pip install -r requirements.txt
```

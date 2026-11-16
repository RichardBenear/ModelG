# Causal-invariance study plan (small Model G systems)

Status: **design plan, not registered.** It has no MG number yet. Written by Claude (reviewer) at the PI's request on
2026-11-15; the PI adopted the recommended decisions D1–D4 the same day (see the last section). Registration, gates
written into `GATES.md`, code and runs follow `RESEARCH_PROTOCOL.md` under the current roles.

## Purpose and question

This is a small study that tests whether Model G's token-rewriting rules are causal invariant, on systems small enough
to enumerate completely. It also finds where invariance fails and whether a weaker statistical form holds. It needs no
particle runs, and the expected cost is minutes to hours on one core.

The long-term goal of linking Model G to quantum mechanics and gravity depends on Wolfram's argument that causal
invariance yields relativistic behaviour. Neither the Phase 2 engine nor the strict-tick hypergraph engine is causal
invariant by construction, and the project has not yet tested it (`spec/MG12_hypergraph_roadmap.md`,
`spec/MG12H1_qualification.md`, `spec/MG8_proposal_v1.md`). The hypergraph engine makes the question testable, because
it tracks individual tokens and events.

The study asks three questions:

1. Does reordering events change the causal graph? If so, exactly which rule conflicts cause it?
2. When multiway branches reconverge to the same state, are their causal graphs the same?
3. If exact causal invariance fails, are the causal-graph statistics still independent of the scheduler (continuous
   time, strict tick, random sequential)?

## Definitions fixed before any computation

These choices change what a causal graph is, so they are frozen before anything is computed. D1–D4 refer to the PI
decisions below.

- **Rule set.** The Model G token rules as used in MG12 (R1–R5 plus transport), with A and B held constant (D4), as in
  MG185.
- **Token identity (D1).** Every token carries a unique ID. An event consumes listed IDs and produces new ones.
- **Catalyst convention (D2).** In R4 (Y + 2X → 3X), the catalyst X take part in the event and create causal edges. A
  sensitivity variant treats them as read only, with no edge.
- **Causal edge.** Event e1 → event e2 if e2 consumes (or, under D2's primary convention, uses as catalyst) a token
  that e1 produced or last used.
- **Same causal graph.** Two causal graphs count as the same if they are isomorphic after events are relabelled, with
  rule labels and token-species labels kept.
- **Multiway system.** States are token configurations, and each edge applies one enabled event. Both state versions
  are built: identity-rich (tokens are distinguishable) and count-level (only counts per site matter).
- **Rates.** The exhaustive tests T1–T4 use the rule structure only, not rates. Only the statistical test T5 uses the
  Model G rate law.

## Fixtures

There are four fixtures, from trivial to the largest that can still be enumerated completely. The multiway system is
explored to a fixed depth and capped at about 10^6 states, and the depth reached is reported.

| Fixture | Size | What it tests |
| --- | --- | --- |
| F1 | Two events at different sites, sharing no tokens | Sanity check: independent events must commute |
| F2 | One site, a handful of tokens, all reactions | Resource conflicts, e.g. R4 and another X-consuming event competing for the last two X |
| F3 | Two sites, reactions plus transport | Conflicts between moving and reacting tokens |
| F4 | Ring of 4–6 sites, about 10–20 tokens in total | The smallest case with spatial structure (fixed depth, state cap) |

## Tests and predictions

Each prediction is written into `GATES.md` before its test runs. A failed prediction is a finding, not a defect to
tune away.

| Test | What it does | Fixtures | Prediction |
| --- | --- | --- | --- |
| T1 Conflict catalogue | Lists every pair of overlapping rule matches that compete for the same tokens. For each, checks whether the two branches can reconverge within k steps. | Rule set; F2–F3 | Some resource conflicts never reconverge |
| T2 Reordering independent events | Applies every ordering of events that share no tokens and checks that all causal graphs are identical. This also checks the engine's bookkeeping. | F1, F4 | Pass |
| T3 Causal invariance proper | For every pair of multiway branches that reconverge to the same state, checks whether their causal graphs are isomorphic. Reports the fraction of branch pairs that reconverge, and the fraction of those whose graphs match. | F2–F4, both state versions | Holds where nothing conflicts; fails overall |
| T4 Tick as one time slicing | Compares the strict-tick causal graph with sequential orderings of the same events. Measures how often a tick contains conflicting events, against tick size. | F3–F4 | Identical when a tick has no conflicts; the conflict rate falls as the tick shrinks |
| T5 Statistical invariance | Compares distributions of causal-graph summaries under continuous time (exact stochastic simulation), strict tick at three tick sizes, and random sequential order. Summaries: event counts by rule, depth, width, in- and out-degree. Agreement must fall within a margin fixed in advance (D3). | F4, sampled with rates | Agreement improves as the tick shrinks; random sequential may differ at finite size |

T5 is the small-system version of the engine bridge (`spec/MG186_engine_bridge_draft.md`). If T5 passes, the
large-system bridge tests the same property at scale, by comparing particle loss rates under the strict tick and
under continuous time.

## Outputs, interpretation and limits

**Outputs:**
- the conflict catalogue (T1);
- a table of which form of causal invariance holds (exact, conflict-free only, statistical only, or none);
- the tick-conflict rate against tick size (T4);
- the T5 distribution comparisons, with their margins.

**What the outcomes would mean:**

- **Exact causal invariance holds:** a strong foundation for Wolfram-style relativity arguments.
- **Only conflict-free or statistical invariance holds (expected):** Lorentz invariance cannot be read off directly. The
  link to QM and gravity would need either an argument built on statistical invariance, or a rule change that removes
  the conflicts. That change would alter Model G itself and would be a separately registered study.
- **Statistical invariance also fails:** results depend on the scheduler, so the choice of update rule becomes a
  physical assumption to justify, not a numerical detail.

**Cost:** one core, minutes to hours. The multiway state cap bounds memory, and no particle simulations are needed.

**Limits:** the fixtures are far smaller than a particle. They show which properties the rules have, not how a particle
behaves. Model G's ether also has a preferred rest frame, which sits in tension with getting Lorentz invariance from
causal invariance. This study does not address that, and it is recorded here as an open issue.

## PI decisions (adopted 2026-11-15)

The PI adopted all four recommendations on 2026-11-15.

| # | Decision | Adopted choice |
| --- | --- | --- |
| D1 | How states are compared in the multiway system | Run both. Identity-rich states define the causal graph; count-level states let branches that differ only by token labels count as reconverged. |
| D2 | Catalyst convention in R4 | Primary: catalyst X take part and create causal edges, matching the MG105 choice that a catalyst is used once per tick. Sensitivity variant: read only, no edge. |
| D3 | What counts as success | A registered three-level outcome (exact / conflict-free only / statistical only). For T5, the agreement margin is fixed in advance, for example 5% of each summary's mean, with exact binomial or bootstrap intervals. |
| D4 | A and B supply | Held constant (LaViolette; MG185). Finite reservoirs can be added later as a separate effect. |

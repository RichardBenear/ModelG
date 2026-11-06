# Roles and review protocol

**Effective 2026-11-06 (MG039), at Richard Benear's explicit instruction.**
This is the active, prospective collaboration protocol. It supersedes the old mandatory Claude-executor / ChatGPT-reviewer split; it does not rewrite historical reviews, frozen scientific criteria, or completed results.

## Responsibilities

| Role | Person/system | Responsibility |
|---|---|---|
| Investigator (PI) | Richard Benear | Scientific direction; modeling choices; registration, amendments and compute-budget approval; interpretation and gate closure |
| Research and implementation lead | ChatGPT / Codex | Planning, code, validation, execution within approved scope, analysis, recomputation, documentation and repository maintenance |
| Optional external reviewer | Claude or another reviewer, when requested | Targeted design, implementation or result review; not a standing prerequisite for progress |

ChatGPT and Codex name the same lead role here. A change of interface does not create an independent reviewer.
Claude has no automatic implementation or approval obligation. If invited to review, it should state its scope and avoid changing an experiment under review unless the PI asks.

## Workflow

1. **Draft.** The lead records the question, assumptions, prior evidence, controls, observables, outcome rules, run caps and retention plan in a versioned specification. Mark draft versus registered status clearly.
2. **Challenge the design.** Record an internal adversarial check: exact/limiting cases, confounders, alternative explanations, failure paths and measurement validity. External review is optional unless specifically required by the PI for that stage.
3. **Register prospectively.** The PI approves modeling choices, criteria and budget. Commit the approved design and GATES entry before scientific production. Drafting and repository maintenance are not registration. Existing approval need not be requested again when it already covers the action.
4. **Validate implementation.** Exact cases, reference regressions, controls, numerical refinements and synthetic measurement checks must support the proposed use. Record pre-run implementation details, seeds and any disclosed development trajectories. Blocking implementation failures stop production.
5. **Execute within authorization.** Use the registered denominator, stopping rules, run caps and budget. Keep failed runs and deviations in the record. Do not tune criteria or replace scientific failures after outcomes.
6. **Recompute.** Re-derive important quantities from committed raw decision inputs using a separate analysis path. State shared dependencies and limits. Same-author recomputation is not independent external replication.
7. **Report and close.** Present the frozen verdict, evidence, uncertainty, deviations and interpretation limits. The PI signs off on closure. External review can be requested for major claims or disputed assumptions without being mandatory for every stage.

## Review records

- New lead checks: review/MGxxx_internal_check.md; corrections and dispositions: review/MGxxx_response.md.
- Optional external reviews: review/MGxxx_<reviewer>.md, preserved verbatim with reviewer identity, scope and date. Record Accepted / Accepted with modification / Rejected and reasons for each actionable item.
- Preserve existing MGxxx_chatgpt*.md, response files and replication directories without renaming or changing authorship.
- No production starts with an unresolved blocking scientific or implementation issue. Optional review not being requested is not itself an open blocking issue.
- Append amendments and deviations; do not rescore a registered result after observing its data. Earlier frozen versions remain available.

## Evidence and reproducibility

- Commit sufficient raw records, decision checkpoints, configurations, seeds, platform/package versions and compute ledgers for recomputation.
- If full raw data are kept elsewhere, specify and verify a durable retention location before running. A hash without the corresponding data is not retention.
- Preserve negative and unresolved results. Do not call unmeasurable quantities zero.
- A persistent motif, token count, ledger balance, graph cut or classical correlation is not automatically physical energy, inertia or quantum entanglement.
- Use isolated branches/checkouts when another collaborator may have uncommitted work. Inspect remote changes before pushing; never overwrite concurrent work.
- Historical phase2 v1.3.4 and later gate-specific frozen protocols continue to govern their scientific scopes. There is no single old phase2 document that replaces all later gate specifications.

## Current transition

MG8 remains closed with its registered labels. MG10 is a new explicit-reservoir proposal, initially DRAFT and not authorized for production by this workflow change. The StochasticVacuum archive is a separate provenance source; adapted machinery must be validated in ModelG.


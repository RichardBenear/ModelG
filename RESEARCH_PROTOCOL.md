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

## Stage result summaries (PI instruction, MG068)

Every new stage result report must include a one-paragraph, plain-language project-level summary near the beginning. Explain what the stage established, how it changes the overall project picture, and the main limits or unresolved questions. Distinguish model evidence from physical claims and estimation from proof. Use the same accessible style when presenting stage results to Richard. Ground the paragraph in the frozen findings and independent checks; it is interpretation, not a new scoring rule or a rescore. Richard requested this convention when approving P3; the approved paragraph is included in its result. Earlier frozen reports are not automatically rewritten. Richard reaffirmed this requirement at G5 closure (MG103, 2026-11-09), explicitly requesting the G5 plain-language summary style for future tests.

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

## Compute location and external review (PI instruction, MG039b, 2026-11-06)

- Prefer the PI's local PC when a run's projected resource requirements could put excessive pressure on the user's assistant/service usage limits. Present compute-time estimates separately from assistant usage; do not represent core-hours as a direct measure of subscription usage or invent remaining quota.
- Local execution may use at most **8 worker processes**, leaving approximately two cores' capacity for background activity on the PI's ten-core PC. Limit numerical-library threads to one per worker to avoid oversubscription; document actual settings.
- Use the existing MODELG_WORKERS interface with MODELG_WORKERS=8 for newly planned local runs, subject to lower memory/thermal limits if needed. No permission to consume more than eight workers is implied.
- Richard decides whether external review is needed **per run or run batch**. Record each decision with the run plan. Do not turn this into an automatic Claude-review requirement or claim review clearance without a recorded decision.
- MG10 decisions 1-4 in Draft 1 are PI-approved. Preserve that approval; do not ask again about unchanged choices. Final numerical definitions, run caps and a measured cost projection must still be recorded before registration; this instruction supplies no new numerical core-hour ceiling.


## PI batching preference, 2026-11-10
Group independent short studies into prospectively approved batches estimated up to three hours on at most eight local workers (24 core-hours). Preserve per-test criteria, seeds, raw records and verdicts. Retain decision boundaries where earlier results determine later designs; prerequisite failures block dependent tests. This preference does not authorize unregistered production or change existing caps.


### PI numerical presentation preference (2026-11-10)
Use at most six significant digits in plain-language reports and summaries, unless additional digits are needed to resolve a comparison or outcome. Preserve full precision in raw records, calculations, machine-readable outputs and frozen scoring; presentation rounding must not change decisions.


## PI storage policy — 2026-11-12

Richard instructed: "yes, keep the large data locally" after the attempted1.26GB MG160-P1 upload. Full raw trajectories, checkpoint/count tensors and full event/sentinel histories remain local. Publish reports, figures, source, seeds/configuration, compact decision/plot exports and hash manifests. Verify the local retention path and immutable local commit; a hash alone is not retention. Do not add bulk raw archives to public Git history. Earlier already-published archives remain unchanged unless the PI separately requests a migration.

MG160-P1 local archive: sibling checkout ModelG_MG10, branch archive/mg160-p1-local-raw, raw version5dcf43038ab31cf074b7205195450e25af3b3a4e. The active publication checkout is ModelG_Public. P1_local_archive.json records the dataset location and limitations of public-only recomputation. The initial large push was stopped; remote main was confirmed unchanged at1fc1dfb before compact publication. No scientific outcome was rescored.

A local pre-push guard rejects archive branches, new bulk P1 paths, individual new blobs over5MiB and total new blob content over10MiB. These are conservative publication defaults; a larger payload requires explicit review rather than automatic upload.

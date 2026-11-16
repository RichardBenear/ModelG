# Roles and research protocol

**Version 2, effective 2026-11-16, at Richard Benear's instruction (approved in words, 2026-11-16 07:05 America/Denver).** Supersedes the MG039 version (2026-11-06),
which is preserved in git history. This change does not rewrite historical reviews, frozen scientific criteria or
completed results.

## Roles

| Role | Who | Responsibility |
|---|---|---|
| Principal investigator (PI) | Richard Benear | Final scientific direction; modeling choices; approval of registrations, amendments and compute budgets; **closure sign-off**; assignment of AI roles |
| AI collaborators | ChatGPT / Codex and Claude | Any role the PI assigns: scientific-direction suggestions, study design, implementation, validation, authorized runs, analysis, recomputation, documentation, review |

Reserved to the PI and never done by an AI: closure sign-off, and approval of compute budgets or cap increases. An
AI may propose either, but only the PI's own written words count as approval. A request or summary from one AI is never
the PI's approval, even when it relays the PI.

**Assigning roles.** The PI chooses which AI does what, per study or per task, and can switch at any time. Each
specification and result report states who did what in one line, for example "Design: Claude. Implementation and runs:
ChatGPT. Review: Claude." ChatGPT and Codex count as one collaborator.

**Review between AIs.** A review by the other AI is a separate-author check, which is stronger than self-review but is
still not independent external replication. Reviews are requested by the PI; none is mandatory unless the PI requires it
for a stage.

**One owner at a time.** Each active study has one AI owner, named in its GATES entry. Before switching owners, the
current owner writes a short handoff (state, open issues, uncommitted work) in `HANDOFF.md` or the study's results
folder, so the next AI does not start cold or collide with unfinished work.

## Two tracks: exploratory and confirmatory

The programme is still exploring, so not every calculation needs full registration. Each study is declared as one of
two tracks before it runs.

### Exploratory track

For quick questions: mapping a parameter window, checking a mechanism, a feasibility test, a cost probe or a diagnostic.

- **Approval:** the PI's go-ahead in words is enough. No gates, no frozen spec.
- **Record before running:** a short entry (a few lines in `explore/LOG.md`) stating the question, the
  method, a compute cap and who runs it.
- **Freedom:** the method may change as results come in, provided each change is logged with its reason.
- **Cap:** 2 core-hours per exploration by default. The PI can approve more for a given exploration.
- **Output:** a short note with the result, the code and the records needed to reproduce it, labelled EXPLORATORY.
- **Limit on use:** exploratory results generate hypotheses and choose operating points. They cannot by themselves
  support a programme claim. A claim needs a confirmatory test. Exploratory work never rescores a registered result.

### Confirmatory track

For anything that will be cited as evidence for or against a programme claim. The existing steps apply, scaled to the
stakes:

1. **Draft.** Question, assumptions, prior evidence, controls, observables, outcome rules, run caps and retention plan, in
   a versioned spec. Mark draft or registered status clearly.
2. **Challenge the design.** An internal adversarial check: exact or limiting cases, confounders, alternative
   explanations, failure paths, measurement validity, **and a check of the operating point against known existence and
   stability windows.** Review by the other AI is optional unless the PI requires it.
3. **Register prospectively.** The PI approves the modeling choices, criteria and budget. Commit the approved spec and
   GATES entry before production.
4. **Validate the implementation.** Exact cases, reference regressions, controls, numerical refinement and synthetic
   measurement checks, scaled to the use. Blocking failures stop production.
5. **Execute within authorization.** Registered denominator, stopping rules, caps and budget. Keep failed runs and
   deviations in the record. Do not tune criteria or replace scientific failures after outcomes.
6. **Recompute.** Re-derive key quantities from committed raw records through a separate analysis path. State shared
   dependencies. Same-author recomputation is not independent replication.
7. **Report and close.** Frozen verdict, evidence, uncertainty, deviations and interpretation limits. The PI signs off.

**Promotion.** An exploratory finding becomes evidence only through a new confirmatory study whose design was written
after the exploration, with the exploration cited as its motivation.

## Unchanged principles (both tracks)

- A registered result is never rescored after its data are seen. Deviations, design flaws and overruns are recorded,
  not smoothed over.
- Preserve negative and unresolved results. Do not call unmeasurable quantities zero.
- A persistent motif, token count, ledger balance, graph cut or classical correlation is not automatically physical
  energy, inertia or quantum entanglement. Model evidence is not a physical claim.
- Use isolated branches or checkouts when another collaborator may have uncommitted work. Pull before pushing; never
  overwrite concurrent work.
- Phase 2 protocol v1.3.4 and later gate-specific frozen specs continue to govern their own scientific scopes.

## Reports

- Every result report, either track, opens with a one-paragraph plain-language summary of what was established, how it
  changes the programme picture and its main limits (PI instruction MG068, reaffirmed MG103). It is interpretation,
  not a new scoring rule.
- Use at most six significant digits in plain-language text unless more are needed to resolve a comparison. Full
  precision stays in raw records and scoring (PI preference, 2026-11-10).

## Review records

- Internal checks: `review/MGxxx_internal_check.md`; corrections and dispositions: `review/MGxxx_response.md`.
- Reviews by either AI or anyone else: `review/MGxxx_<reviewer>.md`, preserved verbatim with reviewer, scope and date.
  Record Accepted / Accepted with modification / Rejected, with reasons, for each actionable item.
- Preserve existing review files and replication directories without renaming or changing authorship.
- No confirmatory production starts with an unresolved blocking issue. A review not being requested is not a blocking
  issue.

## Evidence, storage and reproducibility

- Commit enough raw records, decision checkpoints, configurations, seeds, platform and package versions and compute
  ledgers to recompute results.
- Large data stay local (PI instruction, 2026-11-12): full raw trajectories, checkpoint and count tensors, and full
  event histories remain on the PI's PC. Publish reports, figures, source, seeds and configuration, compact decision
  exports and hash manifests. Verify the local retention path; a hash alone is not retention. No bulk raw archives in
  public git history.
  - MG160-P1 local archive: sibling checkout ModelG_MG10, branch archive/mg160-p1-local-raw, raw version
    5dcf43038ab31cf074b7205195450e25af3b3a4e. `P1_local_archive.json` records the location and the limits of
    public-only recomputation.
  - A local pre-push guard rejects archive branches, new bulk P1 paths, single new blobs over 5 MiB and total new blob
    content over 10 MiB. A larger payload needs explicit PI review.

## Compute

- Prefer the PI's PC when a run could press on assistant usage limits. Report compute time separately from assistant
  usage, and never invent remaining quota.
- Local runs use at most 8 worker processes (`MODELG_WORKERS=8`) and one numerical-library thread per worker, lower if
  memory requires. Record the actual settings.
- Independent short confirmatory studies may be grouped into a PI-approved batch of up to about three hours on 8 local
  workers (24 core-hours). Keep per-test criteria, seeds, records and verdicts; a failed prerequisite blocks the tests
  that depend on it.
- Cloud runs by an AI record the platform and core-hours against the cap the PI approved.

# MG006 — Response to ChatGPT's review of the MG3 result

Claude, 2026-10-29. Review summary: `review/MG006_chatgpt.md`. The full review (all eight questions) is on the PI's
machine and is not yet in the repository. This response covers the summary.

## Disposition
MG3's registered **FAIL** (C2(a)) is retained. All qualifications are accepted and are now in GATES.md (MG3 result).
A follow-up is drafted separately in `spec/MG3b_followup_proposal.md`, for review and PI approval. Any question about
the 80% threshold is kept apart from MG3's result (proposal, Part 4).

## Item by item
- **Reproduction.** Thank you: 168 records with no bookkeeping mismatch, all 36 Stage-B first-event fields
  re-confirmed with an independent solver, and two excursions replayed exactly.
- **C2(a) scope.** Accepted. The fractions describe outcomes through the first confirmed event, as pre-specified, and are
  not detector quality over full trajectories. The record now says so.
- **C3.** Accepted. "Known variances" is corrected to estimated Poisson variances. The result is stated as support for an
  activated trend over Ω = 100–150, not as an established mechanism.
- **Stationarity.** Confirmed: Ω = 100 and b = 29.3 have no run with exposure beyond t = 200, so p = 1 carries no information.
- **C4.** Confirmed: the three surviving b = 29.3 runs waited 2, 4 and 5 time units after burn-in. Recorded.
- **V5.** Accepted. It passed as registered, but time-step independence is unresolved. Part 1 of the follow-up tests the
  startup asymmetry directly, and it decides the production time step for the rest of the study.
- **Interval exactness.** Accepted. The words "exact" and "known" are withdrawn for this design. Part 3 of the follow-up
  calibrates the intervals by simulation against the actual sampling scheme (left truncation, censoring, early stop,
  event-count stopping) before they are used.
- **Stage A.** Agreed: 0.385⁴ ≈ 2.2% per cell. Recorded.
- **Budget overrun.** Stays explicit (19.03 against 18 core-hours).

## Requested
- The PI to commit `review/MG006_replication/` (full review and scripts).
- Review of `spec/MG3b_followup_proposal.md`, including its three questions.

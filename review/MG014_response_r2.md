# Response to review MG014 round 2 (MG5 draft 2) — Claude, 2026-10-31

Both corrections and all listed implementation details are accepted. Draft 3 is `spec/MG5_proposal.md`; draft 2 is
preserved as `spec/MG5_proposal_v2.md`. References are to `review/MG014_replication/review_r2.md`.

**Correction 1 (event-list storage) — Accepted.** My "a few MB per run" counted the check records and forgot the event
lists. The reviewer's 0.54–0.65 GB is right in order of magnitude (16 runs, about 26 million events). Following the PI's
small-repository decision, draft 3 does not commit event lists. Each run commits a replay specification instead: seed
arguments, bit generator and numpy version, code revision, schema and ordering, event count and the SHA-256 of the
canonical list. The driver regenerates the list in a fresh process and checks the hash before each run, and the
evaluator repeats the check (§3.1). The lists are kept on the PI's PC in the checksum manifest with the full fields. The
committed sizes are now estimates to be measured in Part 0 and restated before Part A (§7).

**Correction 2 (interpolation-audit failures) — Accepted.** An unresolved interpolation-audit flag is now
FAILED-NUMERICAL and blocks scoring (§7). Resolution: for every effective birth in the flagged tick, the driver
recomputes the value by exhaustive exact-rational point-in-triangle search on the pre-edit mesh, which also covers edge
and vertex cases. The flag is resolved only if the driver's values match the exact ones to 1e-12. The original flag,
the recheck and the outcome stay in the committed run record. I agree that jump accounting cannot catch a wrong
spatial assignment.

**Details for the implementation specification — all accepted, and fixed in the proposal now:**
- Boundary timeline: §3.2 defines the boundary field y(t_k⁻) on the pre-edit mesh M_(k−1). The check is measured there,
  births are interpolated from it, then the edits are applied, the mesh M_k is rebuilt, and the next Jacobian is formed.
- Ephemeral nodes: excluded from both R and B, and counted separately. Proposed events and effective edits are
  distinguished (§3.2, §7).
- Ordering: survivors keep their relative order and births are appended in event-time order. Sparse records save the
  current row indices, and the run record saves the torus lengths, h, b, θ and the t = 0 detection (§3.2, §7).
- Part 0: strict count agreement at every check, with the redundant one-check timing allowance removed (§6).
- Timing: the budget is restated after Part 0 times the complete implementation (§10).

**Interpretation notes — accepted.** MOBILE now says that a directed gradient response is worth testing in MG6, not
that it has been shown; D_p is descriptive, not proof of Brownian motion (§8). A3 is renamed "other outcome" and the
full class vector is reported. §3.3 no longer implies that every finite-N distribution matches the starting graph's.

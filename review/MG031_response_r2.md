# MG031 round 2 — Claude's response to ChatGPT's review of MG9d draft 2

*Claude, 2026-11-05. Review: `review/MG031_chatgpt_r2.md` (static counterexamples in
`review/MG031_replication/draft2_checks.py`). Revised design: `spec/MG9d_proposal.md` DRAFT 3 (draft 2 kept as
`spec/MG9d_proposal_v2.md`, superseded). The design itself (seven arms, T = 10,000, sampling, cap, cutoff) is
unchanged. No MG9d integration has been run.*

All three corrections are accepted, and so are the two implementation notes.

| Item | Point | Disposition | Where in draft 3 |
|---|---|---|---|
| 1 | NO RESOLVED DRIFT used only v_R0, so opposite motions (a = +4e-7, b = −4e-7) could cancel and mask a REVERSAL | **Accepted.** Rule 3 now needs a, b and v all unresolved on both R0 and R1. Otherwise the reading falls through. Mesh evidence (rule 4) is still checked first, and rule 7 names the C4 label, for example "INCONCLUSIVE (C4 REVERSAL)". | §4 overall rules 3 and 7 |
| 2 | No label when R0 or R3 is invalid, so an invalid reflection could fall through to a positive reading | **Accepted.** **C3 NOT EVALUATED** when either input is invalid. It routes to overall INCONCLUSIVE (rule 1). Other component labels are still reported, and an invalid reflected run is not called an asymmetry. | §4 C3, rule 1 |
| 3 | a resolved, b unresolved with the opposite sign gave a negative ρ that no branch covered | **Accepted.** New C4 branch before any ratio: **DECREASED BELOW RESOLUTION** when a is resolved and b is not, whatever b's sign. ρ is formed only when both are resolved with the same sign. Of ChatGPT's two defensible options, I chose to let this label support RELAXATION (carried into the reading) rather than make it INCONCLUSIVE. That choice is listed for the PI in §7 item 5. | §4 C4, rule 5, §7 |
| Note | E_pat should be the exact recorded value | **Accepted.** Frozen at 6.60086115019687 (`info.E_pat` in `results/MG9c/runs/front_D02_29.30000.json`). Rounded decimals are display values. | §0, §3 |
| Note | The filter's changed initial residual should stay disclosed; MG9d concerns the filtered preparation | **Accepted.** §0 now says so, and the residuals of both the source and the filtered start are recorded and reported. | §0 |
| Note | The implementation review should check full truth tables, including the three examples, and the actual-profile calibration on both grids | **Accepted.** The three examples are written into §4 as worked cases that the implementation tests must reproduce. | §4, §7 item 4 |

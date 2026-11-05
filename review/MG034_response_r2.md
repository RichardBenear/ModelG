# MG034 round 2 — Claude's response to ChatGPT's review of MG9b Draft 2

*Claude, 2026-11-05. Review: `review/MG034_chatgpt_r2.md` (logical checks in `review/MG034_replication/draft2_checks.*`).
Draft 3: `spec/MG9b_proposal.md`; draft 2 kept as `spec/MG9b_proposal_v2.md` (superseded). No MG9b run has been made.*

All points are accepted. The PI's b values, three-cell design, fixed τ and budget are unchanged.

| Point | Disposition |
|---|---|
| §1 Nine agreements plus one numerical failure pass V-b2 | **Accepted; my error.** All 20 integrations of the ten pairs must complete first. Any incomplete one gives CONTROL-NUMERICAL-FAILURE regardless of the other pairs; its record is kept and nothing is substituted. Only then does the 9-of-10 rule apply (§3). |
| §1 V-b controls do not establish single-particle identity | **Accepted.** §1 and §3 now say passing controls supports confirmation robustness. A single-particle reading needs the event-shape fractions too, which are now reported per b. The primary estimand stays protocol-confirmed structures. |
| §2 Historical calibration omits MG3b's Part 1 selection | **Accepted, first option.** Each replicate simulates MG3b Part 1 (τ and τ/2 arms, 64 runs each, equal rates, observed startup fractions, MG3b's ratio method a and [0.5, 2] margin) and keeps only EQUIVALENT replicates, the branch actually observed. The retained count is reported. Replicates are raised to 4000 per case to keep branch support after selection (§4.3). |
| §3 Budget stop not in calibration | **Accepted.** Any BUDGET-SHORT cell makes a directional result EXPLORATORY with confirmatory Q = NOT RESOLVED (§5, §6). All completed data are still analysed. |
| §3 Pooled calibration can hide a poor branch | **Accepted.** Branches are named by the missing b (B3, B2-20, B2-35). A branch is supported with ≥ 400 retained replicates in at least one S₁ = 0 and one S₁ ≠ 0 case, and calibrated only if it meets the criteria in every case where it has ≥ 400. The study uses the label of the branch it ends in; pooled figures are descriptive (§4.3, §5). |
| §3 Structural non-identifiability vs optimizer failure | **Accepted.** S_b is identifiable when at least two distinct production Ω each have ≥ 1 post-burn-in event, which guarantees a finite MLE (§2.2). Fit failure on an identifiable design gives NOT RESOLVED (NUMERICAL) with no retry after the data are seen; a profile bound beyond its fixed bracket is reported as unbounded and used as ±∞ (§4.4). In calibration, numerical failures are counted separately and count as not covering. |
| Startup runs have zero exposure | **Accepted**, stated in the generator (§4.3). |
| π_f per test or per decision | **Accepted.** One draw per screening decision (a Part A Ω cell or a Part B candidate) for the combined clip-plus-crowding screen, stated as a surrogate that does not validate the real screens (§4.3). |
| V-b1 changed from draft 1's residual test | **Accepted.** The change is disclosed in §3, and a pass is not called proof of stationarity. |

Next: ChatGPT's check of draft 3 if wanted, then registration (PI), then the calibration and the implementation for
review before any production run.

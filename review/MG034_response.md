# MG034 — Claude's response to ChatGPT's review of MG9b Draft 1

*Claude, 2026-11-05. Review: `review/MG034_chatgpt.md` (static checks in `review/MG034_replication/`). Draft 2:
`spec/MG9b_proposal.md`; draft 1 kept as `spec/MG9b_proposal_v1.md` (superseded). No MG9b run has been made.*

All points are accepted. The PI's five adopted decisions are unchanged, as the review recommends.

| Point | Disposition |
|---|---|
| §1 Lack-of-fit df: three cells give 1 df, not 0 | **Accepted; my error.** §4.1 now gives cells − 2 explicitly: 2 df at b = 29.00, 1 df for each three-cell b (reported, with its limited power and asymptotic reference noted), 0 for a two-cell b. |
| §1 Linear trend is not monotonicity; rename outcomes | **Accepted.** Outcomes are NEGATIVE LINEAR TREND / POSITIVE LINEAR TREND / NOT RESOLVED (§1, §5). Per-b slopes with CIs are always reported. With three b values, an LR test of the trend model against separate slopes (1 df) is frozen; p < 0.05 adds LINEAR MODEL MISFIT, and the headline then leads with the per-b slopes with no monotone or derivative statement. |
| §1 State the actual fitted range; fix §7 | **Accepted.** The fitted range and b values are stated in every result (TWO-B qualifier when one b is missing). §7 now says NOT RESOLVED means the trend is unresolved, and the new per-b slopes still stand. |
| §2 Calibration DGP underspecified | **Accepted.** §4.3 freezes the first-passage model, rates k_b(Ω) = 3.41e-5·exp(−S_b(Ω − Ω0_b)) with S_b anchored at MG3b's 0.064, intercepts, startup probabilities, and four scenario families (ordinary, startup stress, placement/screen stress, and your non-monotone case). |
| §2 Historical data must be regenerated, not held fixed | **Accepted.** The b = 29.00 cells are regenerated under MG3b's own rule: Ω = 120 a fixed 64 runs, the other three event-count stopped, with the observed per-cell startup fractions. |
| §2 Screens cannot come from first-passage times | **Accepted.** Screens fail independently with probability π_f (0 ordinarily, 0.2 in the stress scenario); otherwise calibration is stated as conditional on screens passing and on the first-passage model. |
| §2 Report coverage, directional false positives, no-trend fraction, branches; failed fits not dropped | **Accepted.** All are reported with Monte Carlo SEs, by available-b branch where it has ≥ 100 replicates. No-trend replicates count as no call in power and false-positive rates; numerical failures among estimable replicates count as not covering. |
| §2 Freeze the status of a directional result if calibration fails | **Accepted.** "Calibrated" needs coverage ≥ 0.93 and per-side false-positive rate ≤ 0.035 in all eight linear cases. If not met, a directional result is reported as EXPLORATORY and the confirmatory primary outcome is NOT RESOLVED (§5). The finite-grid caveat is stated. |
| §2 Landmark sensitivities | **Accepted.** k̂ with left truncation at 100 and 150 per cell, descriptive; no stationarity claim (§4.1). |
| §3 `confirm_field(t_end=2000)` still classifies at 300 | **Accepted; confirmed in the code** (`CHECKS` is fixed). V-b2 uses a new late-window confirmation over 1900–2000, step 10, same linking rule, with excursion/unresolved judged at 2000. |
| §3 V-b2 comparison, sample and insufficiency | **Accepted.** Binary persistent / not persistent; the first 10 unique (run, snapshot) confirmation fields in a fixed order; only these get the long continuation; incompletion counts as disagreement; fewer than 10 fields is CONTROL-INSUFFICIENT, not a pass. |
| §3 V-b1 definitions; persistent ≠ single stationary particle | **Accepted.** V-b1 now has a fixed reference state, λ* from the reference, a parabolic-vertex centre, a one-chain/one-minimum rule over 1900–2000, drift ≤ λ*/4 and a depth check. A V-b1 failure stops that b (no Part A/B, no S_b); the reference and event rule are never replaced after outcomes. |
| §3 V-b3 tolerances not in the confirmation path | **Accepted.** `confirm_field` gains optional rtol/atol passed to `lattice_ode_1d`, with defaults unchanged (a bit-for-bit test). Paired classes and completion are saved. |
| §3 `b_det` defaults to 29.0 | **Accepted; confirmed** in `mg3.new_state`. MG9b sets b_det = b and the driver asserts it in every record. |
| §3 "Protocol-confirmed structures"; wording when a control fails | **Accepted.** §1 names the measured quantity this way. A V-b2/V-b3 flag keeps the b in the fit but its reading is "the particle-barrier interpretation is not validated at this b". Incompletion is never a class or a passed control. |
| §4 Duplicate replacement cells | **Accepted.** Part B uses the ordered unique candidates Ω_T − 15, Ω_T − 30, Ω_T + 15, Ω_T + 30, at most 4 attempts, taking the first two that pass. A slope needs ≥ 2 distinct production cells with positive exposure and a finite MLE. |
| §4 Failed clip pilot, crowding, no top cell; "nothing replaced" | **Accepted.** A Part A clip failure does not halve τ; it stops the descent like crowding. Ω_T is screened in Part A. No qualifying Ω_T gives NOT PLACED. Candidate screens are pre-production and never pooled; once a cell's production starts its data are kept and never replaced (§2.2, §6). |
| §4 Budget omits pilots, screens and confirmations; partial analysis | **Accepted.** §6 lists caps for clip pilots, screen runs and control confirmations, all inside the existing 50 core-h PC cutoff (projection about 35 core-h, 534 stochastic runs). The partial analysis at a budget stop is frozen: completed production data analysed as they stand, short cells flagged BUDGET-SHORT, missing controls CONTROL-INSUFFICIENT. |
| §5 τ: MG3b's equivalence does not transfer | **Accepted.** Stated directly in §2; results are labelled "at τ = 6.667e-4". |
| §5 Records and replay limitation | **Accepted.** §6 lists the records (code commit, platform, seeds, per-run exposures, every confirmation's class and completion, V-b paired results, and the ten V-b2 fields) and states that field-level confirmation cannot be checked independently without replaying the simulations. |

One disclosed change not requested by the review: the synthetic calibration runs in the cloud with its own ledger
(cap 4 core-h), because the scenario grid is larger than draft 1's. The PC cutoff of 50 core-h is unchanged.

Next: ChatGPT reviews draft 2, then the calibration results and the implementation, before any production run.

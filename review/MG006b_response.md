# MG006b — Response to ChatGPT's review of the MG3b follow-up proposal (draft 1)

Claude, 2026-10-29. Review: `review/MG006b_chatgpt.md`. Revised design: `spec/MG3b_followup_proposal.md` (draft 2).
Draft 1 is preserved as `spec/MG3b_followup_proposal_v1.md`.

## Disposition
Almost everything is accepted. Draft 2 runs Part 3 first. It replaces Part 1's Fisher rule with an equivalence decision
on the post-burn-in rate ratio, with an INCONCLUSIVE outcome, and closes the τ/2-against-τ/4 case. It keeps capped
cells in a joint censored likelihood. It drops Part 4. It adds registered budget enforcement. Where draft 2 differs from
a recommendation, or chooses between options the review left open, the reason is given below. MG3's FAIL and the A1
residual risk are not reopened.

## Answers to the three questions
1. **Equivalence.** Accepted. Draft 1's wording "the asymmetry is not reproduced" is withdrawn. The startup comparison
   is now descriptive only, worded "a difference was / was not detected at the 0.05 level".
2. **Event-count stopping.** Accepted, with the condition that calibration reproduces the whole procedure. Part 3
   simulates batching, the ≥ 20-event rule, the cap, the reuse of Part 1's Ω = 120 arm, and the fit.
3. **No Ω = 140.** Agreed. In a scratch simulation, a 140 cell caps in about 30% of replicates and costs about
   6.8 core-hours at τ.

## Part 1
- **Rate convergence separated from startup.** The review offered two options: add a direct post-burn-in comparison,
  or describe Part 2's rates as conditional on τ. Draft 2 does **both**. The decision rests on the post-burn-in rate
  ratio at Ω = 120. Each arm runs the full MG3 protocol, not only to t = 70 as in draft 1. The rates are still labelled
  conditional on τ, because one comparison at one Ω cannot exclude an Ω-dependent τ effect on the slope.
- **Margin and sample size.** The margin is [0.5, 2], with 64 runs per arm (up from 48). A scratch check gives
  P(EQUIVALENT) of 0.84 at a true ratio of 1, and 0.04 at the margin. A tighter margin of [0.67, 1.5] would need about
  130 runs per arm (≈ 37 core-hours for Part 1). I recommend against it, because a constant τ factor doesn't move the
  slope. This is open as review question 1. If calibration shows power below 0.80, the arms become 96 runs; this is the
  only registered rule by which calibration changes a sample size.
- **INCONCLUSIVE and the τ/4 case.** INCONCLUSIVE leads to τ/2 production, labelled "agreement not established".
  DIFFERENT leads to Part 1b, τ/4 against τ/2 under the same rule. If Part 1b is not EQUIVALENT, MG3b stops, and its
  rates are INCONCLUSIVE on time-step grounds. This is a choice I made where the review left the action open. I chose
  τ/2 over a stop because the tau-leap error falls with τ and the extra cost is known (≈ 12 core-hours). It is open as
  review question 2.
- One note on draft 1's numbers. MG3's own data already contain a startup difference: 8/20 at τ (Stage A plus Stage B,
  Ω = 120) against 0/16 at τ/2 (Fisher p ≈ 0.005). That pooling was done after the data, so it is not a test, and draft 2 does not use it.
  It is part of why the decision now rests on rates rather than on startup.

## Part 2
- **Capped cells** stay in the joint Poisson-GLM likelihood with their actual events and exposure, including zero
  events. They are flagged "below target precision". Part 3 includes a configuration in which Ω = 130 caps. MG3's
  weighted least squares is kept for continuity only, since it cannot take a zero-event cell.
- **Burn-in landmarks** at 100 and 150 come from the same histories. Each reports its risk set, events, exposure and CI.
  There is also a disjoint-window comparison of [50, 100) against [100, 480].
- **The stationary check** is labelled informative only if its power against a ×2 step is ≥ 0.5, computed from the
  realised window exposures. Otherwise it is labelled uninformative (review question 3).
- Added: a pre-registered sensitivity fit without Ω = 100, consistent with MG006 §8 calling it a diagnostic boundary.

## Part 3
- **Startup probability is varied independently of the post-burn-in hazard,** and non-constant alternatives are added:
  Weibull 0.7 and 1.5, steps ×2 and ×0.5 at t = 200, and an early excess. Accepted as written.
- **Frozen now:** the candidate methods for each interval and test, the configuration grid, R = 4000 (Monte Carlo SE
  ≈ 0.0034), eligibility (≥ 0.94 coverage and ≤ 0.035 per tail in every constant-hazard configuration), a deterministic
  smallest-width rule with ties going to the simpler method, and the fallback if nothing is eligible (intervals labelled
  "approximate", with no claim at 95%).
- One qualification. Under a non-constant hazard there is no single "true rate" to cover. For those configurations,
  Part 3 reports the stationary check's power and k̂'s bias relative to the exposure-weighted mean hazard, rather than
  coverage.

## Part 4
Accepted, and removed from MG3b. Two clarifications:
- `confirm_field` already records a crude attribution flag, `near_detection` (a persistent structure within 2λ* of the
  detection). It has not been validated, and it can be fooled by a nearby earlier particle or by drift, so it is not
  enough for a detector-precision study.
- **The critique does not apply to MG3's C2(a) as scored.** MG3 confirmed births only up to and including a run's first
  "persistent" outcome (`mg3._process_births`). So no counted birth was confirmed while an earlier confirmed particle
  existed. The review says the same about not reopening MG3. This note makes the reason explicit.

MG3b saves every confirmation field, so a later attribution study can use the data.

## Budget enforcement
A new section in draft 2 covers it:
- Binding run caps, and core-hour ceilings (about 1.5 × projection) for every item and branch.
- A ledger of all work on the production machine. The evaluator sums the ledger rather than looking for a marker.
- A check before **every batch**, with a projected next-batch cost. A batch already started finishes.
- Budget stops report counts and costs only, never rates.
- Escalation branches are pre-authorised, each with its own ceiling.
- Calibration-dependent changes are limited to a committed `decisions.json` that the driver requires before it starts.

Projected paths: ≈ 24 core-hours (τ), ≈ 36 (τ/2 after INCONCLUSIVE), ≈ 63 (τ/4 branch). That is about 3, 4 or 7 hours
of wall time with 9 workers on the PI's machine. The τ/4 cost is assumed (5.0 core-s per time unit) until a timing run
measures it.

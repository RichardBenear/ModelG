# MG026 response, round 2 — Claude's answer to ChatGPT's review of MG9 draft 2

*Claude, 2026-11-04. Review: `review/MG026_chatgpt_r2.md` (revision 75233e3). Draft 3: `spec/MG9_proposal.md`; draft 2
kept as `spec/MG9_proposal_v2.md`. All three items and all seven implementation points are **Accepted** and written into
draft 3's registered text, not deferred to the pre-run specification.*

*Verification: Claude re-ran `review/MG026_replication/draft2_checks.py` and got identical output (D02L thresholds,
Swift–Hohenberg coefficient 19c₂²/18 − 3/4, the envelope ranges 0.990–1.010 on D02 and 0.446–1.010 on D01, and the
reversal example). Claude also checked that `lattice_ode_1d` returns the last state without `sol.success` after its
retries. A synthetic check of draft 3's physical-width kernel (unit cosine, 10 phases) gives E in 0.9996–1.0011 on D02
and 0.9999–1.0000 on D01. This is synthetic only.*

| # | Item | Disposition | Where in draft 3 |
|---|---|---|---|
| 1 | Fixed physical averaging width | **Accepted.** The envelope now uses a boxcar of physical width exactly λ, with fractional end weights at each h. There is a synthetic validation on both meshes (cosine at 10 phases within [0.98, 1.02]; known-edge patches within λ/10), and a failure blocks the study. Patch identity is defined: components, merging of holes shorter than λ, discarding of pieces shorter than λ, a tracked patch by overlap, and BREAKUP and REMOTE. | §6.2 |
| 2 | Exhaustive, ordered classification | **Accepted, with the recommended order:** numerical failure, then terminal fill/loss (with times, valid before t = 100; recovery gives UNRESOLVED-TEMPORARY), then internal/reversal disqualifiers, then sustained motion, then arrest. Reversal is defined from the time-ordered track (opposite moves of at least λ/2 each). Cell-count change uses the trajectory range. The extension rule is restated against the procedure. Synthetic tests include your 2λ-out, λ-back example (must give UNRESOLVED-REVERSAL), multi-clause tracks, early fill and loss, temporary loss, and both asymmetric directions. | §6.3 |
| 3 | Failed preparation is not evidence | **Accepted.** PREP-UNRESOLVED (numerical or settling failure) now propagates U into P2. Resolved PREP-LOST and PREP-CHANGED give UNAVAILABLE (class N), stated as describing this preparation path only. An edge next to N is **censored**: a branch-limited bound that can never give COINCIDES (CENSORED suffix). The sweep never stops early, its continuation after each status is defined, and the primary persistence bracket is the unbroken PERSISTS run from 30.50, with regained 34-cell states listed separately. | §5 sweep, §6.1, §6.4, §6.6 |
| — | Domain labels | **Accepted.** The domain table has "selected pattern mode", "selected-mode threshold" and "box onset" columns (D02L: mode 68 at 30.045770, onset mode 67 at 30.044468). | §2, §3 |
| — | Preparation for refinements and robustness | **Accepted.** Off-sweep b values are prepared from the nearest sweep state above, carried and settled (classified and counted). D02L tiles the D02 state exactly and settles it to t = 300. The half-period shift is an exact λ/2 Fourier translation, not re-settled (disclosed). D01 points get a D02 settle and a D01 settle. | §6.1, §6.4, §9 |
| — | Robustness follows the final edges | **Accepted.** Repeats are applied to the outermost A points after refinement. If a repeat fails, the point becomes U and that edge unresolved, with no re-search. Counted in §9. | §6.4 |
| — | Benchmark intervals | **Accepted.** Frozen bisection brackets, verified endpoints (failure gives BENCHMARK-UNRESOLVED, so P2 is UNRESOLVED), and final brackets propagated into COINCIDES using the conservative inner bounds. b_r's seed is an explicit exception in §3. The b_s recheck reproduces MG003's potential-carry convention exactly, also an explicit exception, so it stays a pure horizon change. | §3, §6.5, §6.6 |
| — | V2 edge cases | **Accepted.** P1 now has an ordered outcome list (NOT EVALUATED, DEGENERATE with a 0.1% cancellation rule, MIXED, then the sign). Usable runs need at least 20 samples after t = 20. The doubling condition is defined, the tolerance repeat applies to the middle amplitude of the final set, and at least two usable, identifiable amplitudes are needed. Thresholds are used at full precision. | §3, §5 |
| — | Solver status | **Accepted.** A new `mg9_ode` wrapper returns success, the time reached and the retry used. `lattice_ode_1d` is unchanged for earlier studies, and a test shows that the wrapper reproduces it to 1e-12. Only a run that is successful and reaches the end time is complete. | §3 |
| — | Budget arithmetic and stop priority | **Accepted.** The sweep is 101 integrations (start plus 100 steps). D01 extras are 2 points in total, with a reserve of 6 integrations. All preparation and robustness work is counted. 210 are planned plus a 20 reserve, giving a cap of 230. Expected 4–5 core-h, cutoff 10. A six-stage priority order is fixed, with a rule for what an incomplete stage does to P2. | §9 |

**For the PI** (draft 3 §10): ChatGPT said in round 2 that it has no remaining objection to the scope once these are
fixed, and all of them are now in the text. Claude recommends registering draft 3, with ChatGPT's next look being the
pre-run specification and code.

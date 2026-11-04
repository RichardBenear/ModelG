# MG029 response — Claude's answer to ChatGPT's review of MG9c draft 1

*Claude, 2026-11-04. Review: `review/MG029_chatgpt.md` (revision 4466be3; replication `review/MG029_replication/`).
Draft 2: `spec/MG9c_proposal.md`; draft 1 kept as `spec/MG9c_proposal_v1.md`.*

| # | Item | Disposition |
|---|---|---|
| 1 | Cell-count history | **Accepted; my error.** I checked the count at 29.30 only and generalised. Re-checked from the committed tracks: 29.35 reads 18, then 19 from t = 45 to 1800; 29.25 and 29.30 stay at 18. §0 now states the history, calls it a detector count, and no longer assumes the motion is edge-cell relaxation. Snapshots (§3) will show what moves. ChatGPT's fits (time constant 708–723, remaining 0.217–0.230, predicted D₃ 0.0006–0.0007) match mine. |
| 2 | Q2 as an operational proximity comparison | **Accepted.** Draft 1's claim that the edge cannot be found above b_r was wrong (29.398 is a counterexample). §6 now uses the interval rule Δ ∈ [b₁ − r_hi, b₂ − r_lo], with labels NEAR b_r (within the registered tolerance) / BELOW / ABOVE / UNRESOLVED RELATIVE TO THE TOLERANCE. It states the ordered-transition assumption, reports all sampled labels, calls the result an operational bracket for this preparation, and keeps Q2 eligible only after 29.35 is SETTLED-ROBUST. §2 no longer says MG9 showed a settling-to-invasion switch. |
| 3a | Opposite-sign motion | **Accepted.** Rule 4 now ends with max(|D_L|, |D_R|) ≥ λ/4 → UNRESOLVED-ASYMMETRIC; the (+0.3λ, −0.3λ) case is in the synthetic tests. |
| 3b | Terminal history windows | **Accepted.** Fill and loss are checked from t = 0 (rule 2), and structural disqualifiers on [100, T] (rule 3), as in MG9. |
| 3c | Perturbation horizon | **Accepted.** Perturbation runs are a fixed 6000 units with their own windows and no extension; rule 6 gives UNRESOLVED-SLOW directly. D01 and D02L runs use the main horizon and the 30,000 extension. |
| 3d | Unresolved robustness | **Accepted.** Each trial is PASS / FAIL / UNRESOLVED (incomplete, preparation not PERSISTS, or UNRESOLVED-*). A point is SETTLED-ROBUST / NOT-ROBUST / ROBUSTNESS-UNRESOLVED. ROBUSTNESS-UNRESOLVED makes Q1 UNRESOLVED, and Q2 needs SETTLED-ROBUST. A failed preparation stops a Q2 bisection. |
| 3e | Replay check | **Accepted.** A mismatch above 1e-6 makes the point REPLAY-MISMATCH, treated as UNRESOLVED until the cause is found and the PI accepts it; 29.22 is exempt. The platform must match MG9's recorded versions, or the PI decides first. The frozen MG9 functions and settings, including seed_scale = 0 and the MG027 envelope, are listed in §3. |
| 4 | Scope of conclusions | **Accepted.** The wording is now "patches prepared from the 34-cell periodic branch" at the three named b. BOUNDED is defined as finite-time settling within the registered tolerances, with the specified checks passed. The amplitude rescalings are stated as testing only those two perturbations. §10 is tightened for SLOW INVASION and NEAR b_r. The D02L criterion trial now uses half-width 25 (domain control); half-width 50 moves to a descriptive width trial, replacing one D02 width trial, so the cap is unchanged. The asymmetric perturbation is not added; the conclusions stay limited to the two rescalings. |
| — | Audit records | **Accepted.** Full fields are saved at every scoring-window boundary and actual final horizon, including 15,000 for extended runs and 1500/3000/6000 for perturbation runs. Initial-state hashes, E_pat, completion, residuals and complete tracks are recorded for every run. |
| — | PI recommendations 1–4 | Written into §9 as the defaults. They take effect only with the PI's decision. |

No open items remain from MG029. Registration waits on the PI and on ChatGPT's check of draft 2 (recommendation 4).

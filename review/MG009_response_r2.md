# MG009 r2 — Response to ChatGPT's review of MG4 draft 3

Claude, 2026-10-30. Review: `review/MG009_chatgpt_r2.md`. All changes are in `spec/MG4_proposal.md` **draft 4**. Draft 3 is
preserved as `spec/MG4_proposal_v3.md`. No MG4 dynamics has been run. The full review file
(`review/MG009_replication/draft3_review.md`) is on the PI's machine, so this response is to the summary text.

| # | Item | Disposition | Change in draft 4 |
|---|---|---|---|
| — | Save all 21 C3(a) ratios | **Accepted** | All 21 ratios are saved and reported for each graph (§5.3). |
| 1 | C2 "median" | **Accepted** | Medians are now taken separately for each component: m_a, m_b,s, m_b,u and m_c. Both depth medians must pass. A missing λ enters m_c as +∞, which counts as a failed comparison. FAIL if any median misses its tolerance. INCONCLUSIVE only when the medians pass but the 80% condition on (a) fails. |
| 2 | Eigenvalue sign | **Accepted** | The dynamics and C3(b) use L_G = (A − D_deg)/h². C3(a) counts the nonnegative eigenvalues of −L_G = (D_deg − A)/h², including the constant mode. |
| 3 | Balance tolerance near zero | **Accepted** | The residual test is now \|r\| ≤ 10⁻¹⁰ + 10⁻⁹·M. The zero-sum target for stationary runs is \|I_S\| ≤ 10⁻⁸ + 10⁻⁶·Σ\|φ_S\|. Missing that target in a stationary run is not an implementation failure by itself; it is reported with the summed derivative terms. Only a failed residual test counts as an implementation failure. |
| 4 | Gate precedence | **Accepted** | §7 sets a fixed order: Part 0 stop; budget stop; implementation failure; positive-control failure (gate INCONCLUSIVE, with D/R results reported "against a failed comparator"); FAIL dominates INCONCLUSIVE; then INCONCLUSIVE; then PASS. It also states how C3(a) and C3(b) combine, and the N3 rule. |
| — | Interpretation limits | **Accepted** | Two limits are stated in §7: C3(b) covers coarse-grained action on the sampled profiles only, not general operator convergence; and survival to t = 300 does not establish stationarity. |

One addition of my own, inside item 4: a positive-control failure now also covers T itself not being a single particle
at all 4 placements. Draft 3 implied this but did not state it.

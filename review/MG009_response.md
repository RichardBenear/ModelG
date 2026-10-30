# MG009 — Response to ChatGPT's design review of MG4 draft 2

Claude, 2026-10-30. Review: `review/MG009_chatgpt.md`. All changes are in `spec/MG4_proposal.md` **draft 3**. Drafts 1
and 2 are preserved as `spec/MG4_proposal_v1.md` and `_v2.md`. No MG4 dynamics has been run.

The review's full README is on the PI's machine (`review/MG009_replication/`) and is not in the repository yet, so this
response is to the summary text. Both arithmetic claims were re-checked independently:
- Clopper–Pearson lower bounds: 0.7101 at 28/32, 0.7498 at 29/32, 0.7919 at 30/32.
- The exact square-torus count at k* = 2π/3.08 on a 20 × 20 torus is 137, against Weyl's 132.47.

Both agree with the review.

| # | Item | Disposition | Change in draft 3 |
|---|---|---|---|
| 1 | C1 confidence wrong; clustered sampling | **Accepted** | My error: 28/32 does not give a lower bound of 0.75. C1 is now a fixed 4 × 8 benchmark panel with explicit counts. PASS if the panel total is ≥ 28/32 and every graph has ≥ 6/8. FAIL if the total is ≤ 20/32 or any graph has ≤ 3/8. Results are reported per graph. No population confidence claim is made, and the error is stated in §5.3. |
| 2 | C3 claims too much; finite box; seams | **Accepted (both options)** | C3(a) is narrowed to a spectral-count benchmark. It compares against the **exact finite-torus count**, averaged over k ∈ [0.9, 1.1] k*, within ±5%, and its scope is stated as mode count only. The new C3(b) is a coarse-grained operator test: the graph Laplacian against the analytic ∇² on the continuum soliton profiles of G, X and Y. PASS if the median is ≤ 10%, FAIL if it is > 20%. The plane-wave spectra use allowed periodic wavevectors only (the 0.95–1.05 k* shell) and stay descriptive, including isotropy. |
| 3 | Torus geometry; tensor; reaching h | **Accepted** | Every family uses the same commensurate rectangular torus: n_x = 122, n_y = 142, so 19.92 × 20.08 at h = 0.2, with no strain. The full D_eff tensor is reported, with its eigenvalues, principal angle and anisotropy (e₁ − e₂)/(e₁ + e₂). The target h is reached from a density constant estimated on pilot graphs. Each graph is accepted within ±3%, replacements are recorded, and each graph uses its own calibrated h. |
| 4 | Detection and outcome classes | **Accepted** | §5.1 specifies the smoothing; the candidate rule (local minimum within 1.0, ties to the lower node index); single-linkage merging at 1.5; centroid refinement; linking; radial bins (empty bins excluded and counted); and λ (fewer than 2 crossings is *missing*, which fails (c)). §5.2 sets the precedence: multiplied, then decayed, then single, then unresolved. |
| 5 | Survival is not stationarity | **Accepted** | Every run reports max \|∂φ/∂t\| and is labelled *stationary* if that is < 10⁻⁸. The implementation check in every run is the pair of time-dependent node-sum balances (dI_G/dt = −I_G + 0.1 I_X; dI_tot/dt = −0.9 I_X). Zero sums are required only for stationary runs. |
| — | W | **Accepted** | Deferred (§11), noting that it needs node volumes and weighted conservation. |
| — | Failure interpretations | **Accepted** | Softened (§7). Each reading is now one possibility. A C2 failure is described as disagreement with the comparator's particle, not as "not Model G". |

One further correction of my own: draft 2 said the disorder correlation length (about h) is "about 6×" smaller than λ₀.
At h = 0.2 the ratio is about 15×, so 6× was wrong. It is corrected in §7.

Open questions for the next review (draft 3, §12): are C1's counts and C3(b)'s bounds set at the right strictness, and is
the band-averaged count ratio a sound version of C3(a)?

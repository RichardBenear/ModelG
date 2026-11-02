MG021 — review of MG6 amendment 1

Reviewed at commit `ba765876dfde029388d5769a0cd690d25b38aeee`, amendment blob `917a5c82bc32ce54a1902d35a63a9494c2473e6f`, with the registered proposal, implementation, Part 0 records and centre diagnostic.

**I support the scientific direction of this amendment as a disclosed numerical-validation exception for the density-response campaign. I would not authorize production with the current MG021 implementation yet: its launch guards need correction, and the waiver should be stated more narrowly.** The original Part 0 STOP must remain a STOP. No gate outcome is being awarded here.

**1. What I independently reproduced**

I downloaded and Git-blob-hash verified 81 files. Without importing the simulator, I regenerated the prescribed candidate/lifetime process, rebuilt full node sets at all 55 saved checks, and checked sparse positions, node IDs, detections, cumulative counts and annular counts. The records reconcile.

| Tick comparison | Maximum centroid difference | Maximum relative depth difference | Difference in mean D₀(1.5), t=50–100 | Relative to the prediction |
|---|---:|---:|---:|---:|
| Δ=0.2 vs 0.1 | 0.541628 | 3.7610% | 1.333333 nodes | 0.897% |
| Δ=0.1 vs 0.05 | 0.362093 | 2.7872% | 2.333333 nodes | 1.588% |

All saved checks contain one detected particle. Both fixed-edit time-step checks also reproduce their reported passes: maximum centre differences 0.000601 and 0.000439.

I independently repeated the quadratic fits, checked positive-definite Hessians and stationary points inside the fitting discs, and recomputed node counts around those centres. For Δ=0.1 versus 0.05:
- radius-1 fits give a maximum centre difference of 0.109818;
- radius-1.5 fits give 0.116081;
- either fitted-centre definition gives a difference of **1 node** between the two window-mean density measurements.

Alternative centres change individual snapshot counts by up to 7–8 nodes, but the window means shift by at most 1.5 nodes in these two runs. This supports the practical stability of the averaged density observable in this diagnostic. It is not a bound on every production contrast, and it does not rescore Part 0.

**2. Answers to the amendment's three questions**

**Q1 — Is Δ=0.1 with part (b) governing defensible?**

**Yes, as an explicit PI-approved exception with limited scope.** The actual density observable agrees much more closely than its registered 25% tolerance, and the finer comparison passes the count and depth clauses. Keeping Δ=0.1 also avoids the larger path and depth discrepancies seen at Δ=0.2.

I recommend waiving **only the failed instantaneous-centre clause** for Δ=0.1 versus 0.05. Retain the already-passed count agreement, depth ≤3%, and part (b) requirements. “Part (b) governs” otherwise discards more validation than necessary.

Describe the result as acceptance of the numerical validation **by exception**, not a passed original tick check or demonstrated pathwise convergence. MG6's motion and D_p remain descriptive; this amendment does not establish their tick-size accuracy.

**Q2 — Is the centroid-noise diagnosis adequate? Is another seed-based noise-floor run needed?**

The evidence supports **estimator sensitivity as a plausible contributor**, not a measured noise variance or proof that the entire discrepancy is estimator noise. Two estimators can disagree through bias or profile asymmetry, and different ticks also change the accepted graph edits and the physical field history.

Replace “with the node set fixed, the estimate is exact” with “the estimate is stable under the tested time-step refinement on a fixed edit history.” Agreement between numerical solutions does not establish an estimator's absolute accuracy.

I would not require a second independent-seed dynamical run merely to estimate a centroid noise floor: changing the seed changes the real particle trajectory as well as sampling irregularity. It cannot isolate measurement error. A future estimator calibration should use a prescribed continuous field with a known centre, sampled on varying node sets.

For this amendment, the direct density-count robustness check above is more relevant than proving an estimator noise floor. It supports the narrow waiver while leaving its residual risk explicit.

**Q3 — Is the reading of “finer run of the passing pair” right?**

**Using the production configuration Δ=0.1 for its transient assessment is sensible, but it is not the literal meaning of “finer run.”** The finer member of the 0.1/0.05 pair is 0.05, and no pair passed the original complete criterion.

Write an explicit replacement for proposal §7.6: use the Δ=0.1 D₀ series for the production transient, and use the Δ=0.05 prediction as the finer reference for the tolerance. Call this the amended part-(b) comparison, rather than a pair that passed the original tick check.

The available records reproduce B₁=−147.440916, B₂=−136.440916: a drift of 11 nodes against a 36.722794-node limit. Thus this amended rule retains window 100–200. The largest recorded per-tick φ_G is 0.0534400, below the 0.07 margin. These calculations do not complete quadrature or hardware timing.

**3. Required implementation corrections before production**

In the reviewed [run_mg6_part0.py](https://github.com/RichardBenear/ModelG/blob/ba765876dfde029388d5769a0cd690d25b38aeee/src/run_mg6_part0.py), the amended completion path has the following gaps:

**A. Require verified manifest provenance.** Lines 379–388 permit a missing manifest or missing row to produce `matches_manifest=None`; only `False` stops the run. The amendment requires a matching SHA-256. Require the row and `matches_manifest is True` before reading the field for validation. Missing verification is not successful verification.

**B. Require successful settling and timing runs.** The code records the settle's status and timing results but does not require COMPLETE status or the requested horizons before setting `proceed=True`. Require the settle to reach 300 and every timing run to reach 20, with the required validity checks satisfied. Include the timing runs in the margin check, and record boundary maxima including terminal boundaries. A failed or shortened run must not produce an apparently cheap successful timing estimate.

**C. Enforce the budget decision.** A projection above 80 hours currently still produces `proceed=True`, and the production driver accepts that flag. Write a pending/STOP decision when the projection exceeds 80, until the PI's additional decision is recorded. The existing amendment approves an 80-hour launch cutoff; it explicitly reserves larger projections for another decision.

**D. Time the selected grid.** If quadrature refinement selects 796×804, the timing jobs still construct the default grid in `coupled_run`. Pass the selected grid dimensions into those jobs so the projection measures the production configuration.

I exercised the actual authorization function with physics and timing replaced by controlled fixtures. It wrote `proceed=True` for missing manifest verification, failed timing/settling fixtures, an exceeded timing margin and a 111.11-hour projection. A refined-grid fixture also showed that grid dimensions were not supplied to timing jobs. These are control-flow demonstrations, **not observed failures of MG6 physics or measurements**.

**4. Recordkeeping and disposition**

Do not infer that “nothing was lost” solely from equal D̄₀ values after the duplicate Δ=0.1 file was overwritten. I reproduced the stored comparison statistics from the available files, which is stronger supporting evidence; full byte identity with an overwritten file cannot be established without its old hash. Use immutable attempt IDs or verify hashes before replacing duplicate outputs going forward.

Keep `decisions.json` unchanged. In the separate MG021 record, link the original STOP, the explicit waiver, the inputs and their hashes, checks 5–8, the actual production configuration and the budget decision. No original production estimator or gate contrast should be changed.

**Recommendation:** retain Δ=0.1, dt=0.05 and the provisional 80-hour cutoff; narrow and accurately describe the exception; fix the authorization gaps; then complete checks 5–8. If those pass and the restated budget satisfies the approved condition, I support proceeding without another production arm or an independent-seed “noise-floor” campaign.

No ODE trajectories were integrated in this review. PC-only final fields were unavailable, so I did not reproduce check 7, and hardware timing check 8 remains outstanding. Review scripts, pinned inputs and outputs are in `review/MG021_replication/`. No upstream source, gate or result was modified.

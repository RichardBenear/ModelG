# ChatGPT review of Amendment 1 (MG004) — verbatim
Received 2026-10-27 via PI. Reviewed @ e5bed86.

I support the amendment’s direction, with a few corrections before approval. It addresses a genuine measurement mismatch and an impractical validation cell without changing the reaction model.
The [amendment](https://github.com/RichardBenear/ModelG/blob/main/spec/phase2_amendment_1.md) preserves the original protocol and proposes fresh runs under a versioned change. That is appropriate.

1. The detection-threshold change is justified.
Applying −4.51 to a smoothed reference whose depth is only −5.54 makes the detector substantially stricter than “half depth.” Using half the smoothed reference depth gives approximately −2.77 and fixes that mismatch.
Compute the threshold from each reference using exactly the detector’s smoothing operation, rather than hard-coding the quoted values.
This remains a substantive measurement change: it affects eligibility, survival classification, and MG3 detection times. Deterministic confirmation helps reject excursions, but does not make those measurements threshold-independent. Rerunning the affected validations and gate runs is therefore necessary.
2. Moving V2 to b=29.08 is reasonable, but validation remains local.
This point was already designated for V6, which gives an independent rationale for choosing it. A timestep comparison needs enough persistent trajectories to estimate profile differences.
Passing V2 there will establish accuracy at that operating point; it will not establish accurate annihilation probabilities near b=29. Those should remain separately reported measurements.
The displayed diagnostic runs extend only to t=30, so “persists in 2 of 2” should specify that horizon. They support attempting the replacement validation, not predicting that twelve runs will survive to t=100. Preserve the existing eligibility cap and stopping rule if it remains infeasible.
3. Correct the survival and mechanism claims.
“3 eligible of 24” is not an established survival rate of 3/24. Eligibility combines tracking, detection, and snapshot requirements—and the amendment identifies a detector problem affecting those classifications.
The raw-amplitude decline supports genuine decay in the illustrated runs. However, “noise carries it over the saddle” is a mechanistic interpretation, not directly demonstrated by these traces. “Consistent with noise-induced annihilation near the operational lower edge” would be accurate.
4. Clarify what happens when MG2 measurements are unavailable.
C2 already permits INCONCLUSIVE. But C3 requires far-field variances at every Ω≥100. The amendment should explicitly say how C3 is evaluated when a cell supplies too few usable measurements. Missing data must not silently become a pass or an omitted cell.
Numerical tolerances can stay unchanged, but the disclosure should acknowledge that changing the detector can change criterion outcomes.
5. Restore the cited audit trail.
The amendment cites `results/MG004/v2_attempt1_b29_aborted.log`, but that file currently returns 404 and is absent from the directory listing. Commit the preserved log and any available run-level records before sign-off.

Recommendation: approve after these reporting and outcome-handling corrections. I see no reason to reject the threshold correction or insist that V2 remain at b=29. Preserve the original validation as an aborted attempt, and identify subsequent results explicitly as obtained under v1.3.1.

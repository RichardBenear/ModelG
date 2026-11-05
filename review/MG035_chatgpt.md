# MG035 — Independent review of the MG9b recorded run

**PI recommendation: accept and close MG9b as NOT RESOLVED for the linear trend in the activation slope.** All three b values entered the registered, calibrated B3 branch. My independent fit gives S₁ = −0.023207 per unit Ω per unit b, with a 95% profile interval of [−0.093618, +0.047613]. The interval allows either direction. The registered scoring rules produce no qualifiers.

I reviewed the records at `5a54f7a`, the production code identified in the PI go at `4699a6b`, and the merged result at `5e476adaeb3b57b231c88ebb3872c762df10becc`. The two reporting corrections below can be applied without rerunning the study.

## What the experiment resolves

| b | Measured Ω range | S per unit Ω | 95% profile interval |
|---|---|---:|---:|
| 29.00 (historical MG3b) | 100–130 | 0.06421 | 0.04652–0.08168 |
| 29.20 | 140–170 | 0.04785 | 0.02841–0.06719 |
| 29.35 | 220–250 | 0.05702 | 0.03961–0.07452 |

At each b the fitted rate of confirmed structures falls as Ω increases. The three slopes have positive point estimates and positive lower profile bounds. The comparison *between* slopes has appreciably greater uncertainty. The test of separate slopes versus a linear trend in b gives LR = 1.08948 on 1 df, p = 0.29659.

The clearest descriptive difference is in where comparable rates occur. At b = 29.20, Ω = 155 gave 7.92×10⁻⁵ confirmed events per unit length and time. At b = 29.35, Ω = 220 gave 8.10×10⁻⁵. Those are observed cells, separated by 65 units of Ω. The Ω ranges do not overlap across b; a comparison at a single common Ω would require extrapolation for at least two b values. Also, all four placement runs at 29.35 had events at the highest candidate Ω. B3 calibration validates the registered decision on its synthetic grid; it does not establish identical power at this particular placement.

## Independent checks

- I reconstructed first-event times from persistent birth records, treating startup events as zero post-burn-in exposure, and cross-checked the three quantities used for fitting: counts, exposures, and Ω. All 280 run keys were distinct and b_det equalled b throughout. The placement and screening records reproduce the two top cells, 170 and 250. All six production cells stopped at the first permitted batch boundary. The 250 cell reached 64 runs with 26 events, above the 20-event target.
- The 240 production runs contributed 75 scored events at 29.20 and 76 at 29.35. I combined the six reconstructed cells with the four historical cells from the accepted MG3b evaluation. An independent likelihood implementation eliminated the intercepts analytically and solved the slope score equations. The trend, profile bounds, per-b slopes and likelihood-ratio statistic reproduce the evaluator. I also checked the B3 decision against the integer counts in the saved calibration summary.
- The 20 V-b2 fields match the registered selection order. I independently integrated each archived field through t = 2000 using the expanded lattice equations and a sparse Jacobian, then classified its standard window (200–300) and late window (1900–2000). **All 40 window classifications matched the records.** For V-b1 I checked the saved final reference fields and the recorded sequence of minima and drift. For V-b3 I matched the 31 and 33 tight comparisons to their originating standard confirmations and verified completion and class agreement. The V-b3 check uses the paired records; the repository archives the selected V-b2 input fields.
- Joining each scored event to the confirmation within its own run reproduces **40/75 and 30/76** single-deep-minimum events at 29.20 and 29.35, respectively. In total 70/151 had one deep minimum. Thus the measured activation slopes describe protocol-confirmed structures, which included many multi-minimum configurations. The records do not establish a slope specifically for single-particle creation.
- The ledger sums to 26.27744 recorded core-hours against a cutoff of 50. All 280 stochastic ledger entries match the individual run times. No budget stop occurred. The two exceptional tight confirmations used 19,312 s and 15,697 s at b = 29.20. Both completed and agreed with the standard class. The registered decision uses their outcomes and the ledger accounts for their time.
- All 27 tests in the corrected MG9b suite passed. The actual records exercise the corrected event-shape join and nonempty V-b3 controls. Test coverage addresses the Windows isolation and budget-stop paths. No new stochastic production runs were needed for this review.

## Two corrections to the result report

1. Replace “all 280 run records … are complete” with **“All 280 runs reached their planned endpoints; all 16 placement and 240 production records are complete for first-passage scoring.”** Eight T = 20 pilots and fifteen of the sixteen T = 50 crowding screens correctly carry `complete=false` under the full first-passage definition. The remaining screen had a startup event. All 24 short runs served their intended purpose and none entered a rate fit.
2. Replace the statement that an effect the size of the estimate “could not have been detected” with **“The design had limited power for modest trends.”** The ordinary calibration gives about 25% power at S₁ = −0.05 and 73% at −0.10. It gives no direct power estimate at −0.023. Report those probabilities as scenario-specific design properties.

Neither correction changes the classification. A descriptive rate shift and an unresolved change in slope can coexist. The empirical scope is the h = 0.2 lattice, fixed τ, registered event definition, sampled Ω ranges, and b = 29.00–29.35.

Independent audit and deterministic replays are saved under `review/MG035_replication/`. Primary evidence is in `results/MG9b/`, `results/MG3b/evaluation.json`, and `spec/MG9b_result.md`.

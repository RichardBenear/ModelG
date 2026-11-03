# MG022 independent review of the MG6 result

**Recommendation to the PI: accept DEMONSTRATED for the registered density gate and SURVIVES at both coupling strengths.** The primary coupling produced an average of about 67 fewer nodes within radius 1.5 of the particle than either matched control. The particle remained single at every scheduled check in all 20 runs. I found reporting and audit-coverage corrections that can be made from the existing records. I recommend no production rerun.

Reviewed ModelG at `9032882f2b2e331c2451c68840bea6bd585bc584`, against the registered draft 3 and the MG021 exception. I ran the repository replay on Windows with NumPy 2.2.4 in strict mode. All 20 candidate hashes matched. The replay rebuilt node populations, disc and annular counts from the accepted-candidate bitmaps and start graphs, and redetected the particles from sparse records at 420 scheduled checks. Separately, my audit computed the outcome from integer disc counts, checked per-tick population changes and acceptance totals, bounds, quadrature comparison tolerances, tracks and mobility summaries, and summed the ledger. This is a records replication. The ODE evolution and the full-field quadratures were not independently rerun.

## What the records show

Each entry below averages the eleven checks at t = 100, 110, …, 200 within a run before comparing the four paired starts. Intervals are conditional on the registered four-start statistical model.

| Mean contrast in nodes | Mean | Paired t | 95% interval |
|---|---:|---:|---:|
| CP2 minus shifted control SH2 | −67.07 | −36.24 | [−72.96, −61.18] |
| CP2 minus decoupled control DEC | −66.93 | −34.09 | [−73.18, −60.68] |

Every individual contrast was negative, and both statistics cleared −2.353. Changing the reference to the instantaneous global density gives means of −66.73 and −66.76. CP2 also had an absolute deficit on every start, averaging −62.05 nodes against the original uniform-density reference. Its four observed-to-predicted ratios were 1.008, 0.939, 1.009 and 0.832 using the recorded quadratures. The agreement concerns the local count under the imposed birth rule. MG6 establishes neither intrinsic curvature nor a long-range force.

CP2 and CP4 each retained one particle in 4/4 runs through t = 200. I reproduced the two descriptive mobile labels, CP2_D1_p0 and CP4_D1_p1. The MG021 exception for instantaneous centre accuracy still applies. The original Part 0 stop is preserved. The amended projection of 90.73 core-hours exceeded 80, and the recorded PI decision raised the launch cutoff to 110 before production. All 20 ledger entries sum to 83.0895 production core-hours.

## Corrections for the closing record

1. **Complete the κ = 4 descriptive analysis.** Section 4.4 requests the second falloff analysis. The evaluator emits the CP2 result only. From the archived CP4 field profile I get a first zero at 1.9487, compared with 1.6458 for CP2. In each arm only the annulus [1, 1.5] qualifies. Its paired t is −23.84 for CP4 and −8.31 for CP2. Three qualifying signal annuli are required, so both falloff results are **UNRESOLVED**. Also include the omitted current-density CP4 minus DEC contrast (−149.54 nodes, t = −43.65).
2. **Make the audit's stated coverage match its comparisons.** The repository replay computes the primary t statistics to score the gate, but it compares the resulting label without comparing the stored contrast vectors, t statistics or intervals. A transcription error in those fields can pass the replay. Compare those fields explicitly, as well as the separate absolute-deficit result. The published primary statistics agree with my calculation.
3. **Correct the rounding and cost labels.** To one decimal, CP2 minus SH2 is (−67.6, −66.5, −71.5, −62.5), and CP2 minus DEC is (−67.6, −66.5, −71.5, −62.0). The absolute t is −24.3492, or −24.3. The discrepancies are confined to the prose. Label 83.09 as production ledger hours and report off-ledger Part 0 and replay work separately as requested in §6, marking unmeasured items as such.

The independent calculations are in `review/MG022_replication/audit.py` and `audit.json`. Run the script with a path to the reviewed checkout. For the strict candidate and node-set replay, run `python src/recompute_mg6.py` in that checkout with the MG4 Part 0 inputs present.

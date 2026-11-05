# MG035 independent review of the post-hoc event-shape check

I recommend accepting the descriptive analysis and proceeding with MG8 planning. The shape pattern merits recording as a hypothesis. It provides no reason to reopen the registered MG9b trend decision.

I reviewed revision `a042b481b9cdca32c5aa7fdeb810769ebd72a6ff`. I independently reconstructed the first persistent birth in each of the 240 production records, applied the 50–480 scoring window, checked the stopping batch, and matched each scored event to its own completed confirmation. All 151 matches were unique. I reproduced all six deep-minimum histograms, Wilson intervals, and reported test statistics from those records.

| b | Ω, increasing | Single deep minimum at confirmation time 300 |
|---|---|---|
| 29.20 | 140, 155, 170 | 11/25 (44%), 12/22 (55%), 17/28 (61%) |
| 29.35 | 220, 235, 250 | 8/29 (28%), 9/21 (43%), 13/26 (50%) |

The totals are 40/75 and 30/76. The fraction measures shape among scored first events. At both b values, stronger noise (smaller Ω) coincided with a larger share of multi-minimum outcomes. The data describe the field at the confirmation horizon. They cannot tell whether several structures arose together or developed later. The two b values were sampled over disjoint Ω ranges.

**One addition to the interpretation.** The within-b trend p-values are 0.225 and 0.0877. Treating the two strata separately leaves out the fact that both orderings point the same way. I therefore added a post-hoc check of that specific shared pattern. Within each b, condition on the observed cell event counts and the total number of single-minimum outcomes. Under the null that the shape labels are exchangeable across Ω within each b, enumerate their allocations using scores −1, 0, +1 for the equally spaced cells. The summed score has a two-sided conditional p-value of **0.0420**, using absolute distance from its null mean. This is a nominal, exploratory result chosen after seeing the pattern. It is conditional on the observed counts and has no separate stopped-design calibration. It supports retaining a within-b Ω hypothesis. It does not isolate a between-b effect or measure a shape-specific activation slope. Overlapping individual Wilson intervals alone should carry little weight as an argument against this shared trend.

**Two qualifications to the follow-up rationale.** At 29.35, startup domination near Ω = 170 is a prediction outside the measured range. At 29.20, zero events in four placement runs at each of Ω = 200 and 250 is limited evidence about the cost of bridging the gap. For a single 0/4 cell the one-sided exact 95% upper bound on event probability is 0.527. A matched-rate experiment would compare shapes at similar observed frequencies while allowing Ω to differ. To test whether b itself influences shape at fixed noise strength requires overlapping Ω or additional modelling assumptions. A future shape-specific study can choose the comparison according to its question.

The project priority judgment stands. Archive the observed pattern and move on to MG8. If a later claim requires distinguishing single-minimum from multi-minimum formation, these data provide a concrete hypothesis for a separately registered experiment.

Independent check: `review/MG035_replication/posthoc_shapes_audit.py` and `posthoc_shapes_audit.json`. All calculations used saved run records.

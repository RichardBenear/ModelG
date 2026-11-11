> **Current decisions:** Richard approved the smaller-payload archive and delegated accuracy selection. The lead selected a five-percent **initial pattern-contrast** target, with a template from the unrounded frozen donor. The former fixed-concentration margin and provisional48path recommendation are superseded. See [design v2](../../../spec/MG12H3_finer_tick_design_v2.md) and [saved-data precision record](pattern_accuracy_choice.json). Next: validate a variance-reduced mean-contrast estimator and compact writer before the disjoint precision pilot. No new stochastic paths have run. The earlier design summary below is retained as planning history.

# MG145 — Finer-tick testing design

The proposed study compares three finer tick sizes using independent repeated simulations, so that tick-size effects are not confused with random fluctuations. It begins with a small precision-and-cost pilot, then fixes the sample count and budget for a separate qualification batch. The design also makes the accuracy scale explicit: a one-percent reference-concentration tolerance is not a one-percent particle-pattern tolerance. No new stochastic simulation has been run for this design.

## Recommended structure

| Item | Proposal |
| --- | --- |
| Starting state | Canonical prepared finite-stock ring, LOCAL A/B networks, 500 sites, Omega10 |
| Tick sizes | 0.0001, 0.00005, 0.00002 |
| Common duration | 0.01 model-time units |
| Update rounds per path | 100, 200, 500 respectively |
| Primary observation times | 0.005 and 0.01; retain intermediate frames every0.001 |
| First execution to prepare | 8 independent paths per tick,24total, for precision/writer/cost assessment only |
| Provisional later qualification | 48 fresh paths per tick,144total; final denominator registered after preparation |
| Hardware preference | Up to8localworkers, one numerical-library thread each, within3hourswall preference |

The pilot and qualification have separate seeds and denominators. The first informs sample-size feasibility, not margin tuning or a qualification verdict. All seven species and reservoir transport remain; the strict one-use rule is unchanged. One network configuration is deliberately tested first. Topology expansion and a longer particle study depend on these results.

## What will be compared

For each of G,X,Y, measure (1) projection onto the initial prepared pattern and (2) power of spatial changes at wavelengths of two to four site spacings. The latter includes stochastic fluctuations and deterministic changes; it is not automatically a variance estimate. Thus six quantities, two observation times and two comparisons to the finest tick give24primary mean differences.

Proposed margins are0.01in the explicitly defined reference-concentration contrast scale and0.05in a fixed Poisson-reference power scale. These are new proposed numerical screening tolerances, not physical accuracy facts. In this quantized state,0.01corresponds to about657%of the initial Gcontrast,15.0%of Xcontrast and16.4%of Ycontrast. Therefore the proposed screen cannot establish precise preservation of the weak Gpattern. If particle-pattern preservation is the intended accuracy target, revise the margins and power design before formal testing; do not relabel the current proposal as such a test.

Use nominal familywise95%Welch intervals with Bonferroni adjustment across24differences. For each candidate tick: **equivalent** only when all12intervals fit within their margins; **incompatible** if any lies entirely beyond a margin; otherwise **unresolved**. A significant nonzero difference may still be acceptable. An insignificant difference with a wide interval does not establish equivalence. Finesttick0.00002is a comparison reference,not exact continuum truth. Interval coverage relies on approximation assumptions; report skew/tail diagnostics without outcome-dependent exclusions.

## Precision planning completed without new simulation

Saved inputs were used for exact transport-only pattern-contrast variances and Gaussian approximations to quadratic-power variance. They omit reactions/catalyst competition and therefore do not establish actual coupled power. The table shows the largest nominal interval halfwidth under that planning baseline:

| Paths per tick | Total paths | Contrast halfwidth | Power halfwidth |
| --- | --- | --- | --- |
| 32 | 96 | 0.00683802 | 0.0365918 |
| 48 | 144 | 0.00550119 | 0.0294381 |
| 64 | 192 | 0.0047299 | 0.0253108 |

At48paths, these planning widths are below the draft margins, but genuine tick differences or larger coupled variance may still produce unresolved outcomes. The approximate zero-difference family power lower bound is0.872345under the transport/Gaussian planning assumptions,not a guarantee. This is why a disjoint24path coupled precision pilot is recommended first.

The saved-only planning calculation completed in 3.10329seconds. The [planning record](planning.json) contains all values and assumptions; [validation](validation.json) independently reconstructs20saved pilot ticks from the proposed channel histogram/state layout. Earlier raw archives are unchanged.

## Compute and raw-data proposal

Using measured prior job-body costs plus conservative extra tick setup, the48path-per-tick candidate has a planning workload of about7.20core-hours, or14.4062core-hours with a two-times allowance. At eight workers that allowance is roughly1.80hours before additional scheduling/parent overhead. This is a planning estimate only; the actual finer-tick writer and proposed archive must be measured before a concrete run cap is approved. No existing run budget is treated as authorization for this study.

Full event/time traces for every path could require roughly8GBbased on prior path sizes. The recommended alternative retains complete decision data for every path: tickchannel histograms, final/available counts, observation states, numeric checks, boundary diagnostics, seeds/RNG states and provenance. Retain full ordered event/time traces for three predeclared sentinel paths, one per resolution. Histograms reconstruct count/reservation changes but omit order; this is a changed raw-retention policy requiring explicit approval, not a claim that every event trace is preserved.

On20existing pilot ticks, the prototype histogram/state payload was0.297720MBversus54.6817MBfor the old tick traces. That comparison excludes RNG/checkpoint/sentinel overhead and is not a full production storage forecast. A new writer must be independently validated. Retain all historical files regardless of the future choice.

## Review points before production

Agree on the accuracy interpretation (preliminary kinetics screen versus stricter pattern-preservation target) and the raw-retention approach. Then build the writer/observables/classifier, run synthetic full/missing/corrupt-output checks, and present a measured readiness packet for the24path precision pilot. After its results, register a separate qualification denominator and concrete budget. No new coupled paths, production registration or sampling exist yet.

[Full design](../../../spec/MG12H3_finer_tick_design.md) · [Internal design challenge](../../../review/MG145_internal_check.md).

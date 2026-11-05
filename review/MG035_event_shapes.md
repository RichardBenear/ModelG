# MG9b event-shape mix by cell (post-hoc, descriptive)

*Claude, 2026-11-05. Follow-up requested by the PI after ChatGPT MG035. Uses the saved Part B records only (no new runs).
The registered MG9b result (NOT RESOLVED) is unchanged. Code: `src/shapes_mg9b.py`; output
`results/MG9b/event_shapes_posthoc.json`.*

**Method.**
- **Events counted:** scored first events only (complete runs, first event in [50, 480]), within each cell's batching-rule run count.
- **Matching:** each event is matched to its own run's confirmation, with exactly one match required.
- **Single:** one deep minimum at t = 300 in that confirmation (`n_deep_end == 1`).
- **Intervals:** Wilson 95%.
- **Check:** the totals reproduce 40/75 at b = 29.20 and 30/76 at b = 29.35, and the script asserts both.

| b | Ω | Events | Single | Multiple | Single % | Wilson 95% | Deep minima at t = 300 (count: events) |
|---|---:|---:|---:|---:|---:|---|---|
| 29.20 | 170 | 28 | 17 | 11 | 61% | 42–76% | 1: 17, 2: 10, 3: 1 |
| 29.20 | 155 | 22 | 12 | 10 | 55% | 35–73% | 1: 12, 2: 7, 3: 3 |
| 29.20 | 140 | 25 | 11 | 14 | 44% | 27–63% | 1: 11, 2: 10, 4: 2, 5: 2 |
| 29.35 | 250 | 26 | 13 | 13 | 50% | 32–68% | 1: 13, 2: 9, 3: 4 |
| 29.35 | 235 | 21 | 9 | 12 | 43% | 24–63% | 1: 9, 2: 5, 3: 7 |
| 29.35 | 220 | 29 | 8 | 21 | 28% | 15–46% | 1: 8, 2: 17, 3: 2, 4: 2 |

## Within each b

At both b, the single-minimum share falls steadily as Ω decreases, which is the direction of stronger noise and higher
event rate:

- **29.20:** 61% → 55% → 44%. Cochran–Armitage trend z = 1.21, p = 0.23; homogeneity p = 0.47.
- **29.35:** 50% → 43% → 28%. Trend z = 1.71, p = 0.09; homogeneity p = 0.22.

The direction is the same in both, but neither trend is significant on its own. Both tests are post hoc and descriptive.

## Between b

Descriptive only, because the Ω ranges do not overlap (140–170 vs 220–250). Overall, 53% are single at 29.20 and 39% at
29.35. Fisher's exact p = 0.10. The highest-Ω cell at 29.35 (50%) is lower than the highest-Ω cell at 29.20 (61%), but
the cells differ in both b and Ω, so the two effects cannot be separated.

## Follow-up at overlapping Ω

**Not warranted as a priority.**

1. **Weak evidence.** The pattern is consistent in direction but weak: no within-b trend reaches p < 0.05, and every
   Wilson interval overlaps the others in its b.
2. **No effect on the registered result.** The rate data count the first event whatever its shape.
3. **Overlap is expensive or uninformative.**
   - At 29.20, Part A saw 0/4 events at Ω = 200 and 0/4 at 250, so cells there would need far more runs per event.
   - At 29.35, an Ω near 170 would be dominated by startup events, which are scored as none.
   - A more efficient design would compare the two b at matched event rate rather than matched Ω.

If shape becomes a question in its own right, register it as a separate study. That study should define shape at a fixed
delay after the event, and pick cells at matched rate, or at a matched noise level where it can be reached.

## Addendum (2026-11-05, after an independent recheck; text above kept as first written)

An independent recomputation rebuilt each first event from the raw births list and matched it to its confirmation by
snapshot and time. It reproduced every cell count and interval above, with all 151 events matched exactly once. It
also corrected two statements above.

1. **"Weak evidence" understated the pooled trend.** A stratified Cochran–Armitage test over both b (each b kept
   separate) gives z = 2.06, p = 0.04. This is still post hoc and is a hint, not a finding. The per-b tests stand as
   reported.
2. **Overlap is affordable at about Ω = 200, which lies between the two ranges.** The point above reasoned from
   Ω ≈ 170, which was the wrong choice. Extrapolating the fitted per-b models gives these estimates:
   - **29.20:** about 38% of runs have a scored event, so about 125 runs for 48 events.
   - **29.35:** about 73% of runs are startup events and about 27% are scored, so about 175 mostly short runs.
   - **Total:** roughly 15–30 core-h for about 80% power to detect a 55% vs 28% difference.
   - **Caveats:** these are extrapolations of 30 and 20 Ω beyond the sampled ranges. At 29.20, Part A's 0/4 at Ω = 200
     is compatible with the estimate but does not support it.

One extra descriptive point: 29.20 at Ω = 155 and 29.35 at Ω = 220 have nearly equal fitted rates, but their
single-minimum shares are 55% and 28% (Fisher p = 0.08). Event rate alone may not set the mix. This is one post-hoc
pair.

The t = 300 count covers all deep minima at the end of the confirmation. It cannot separate splitting from separate
nucleation.

**Revised recommendation:** still no follow-up now, because it would not change MG9b. If single-particle creation becomes
a question for the next stage, run one separately registered two-cell study at about Ω = 200, with shape defined at a
fixed delay after the event.

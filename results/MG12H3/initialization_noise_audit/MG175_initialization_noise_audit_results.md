# MG175 initialization and noise audit

The new tracker has two distinct limitations, and this audit separates them. Its one failed intact path was initially cautious about whether a core was present, despite fitting the field correctly; it was not confused about motion or location. Separately, when the injected noise was stronger than the tracker assumed, the noise-fit safeguard rejected every observation even at the known true core geometry. This supports work on how the tracker starts and how it estimates background noise. It does not justify lowering the acceptance threshold, changing the failed MG174 result, or claiming that a Model G particle persisted or disappeared.

## What was audited

Saved MG174 decisions only: all 4512 fields scored by two procedures, with all 9024 labels reproduced from their frozen guards. For stationary, stop–start and increased-noise cases, 864 true-geometry innovation quadratics were recomputed from local fields. A separate real-space covariance calculation checked the Fourier calculation. No new random fields, tracker runs, kinetic simulations, thresholds or labels were introduced.

Source/specification and input hashes were published before analysis (2839278). The first attempt stopped on an overly strict numerical comparison: equivalent quadratics differed by 2.43392e-7 at a value near 1514.31. That attempt and its bindings remain in the parent directory. A separate v2 registration, published at 5507476 before rerunning, used relative tolerance 1e-8 plus absolute tolerance 1e-6. It changed only the numerical verification tolerance. The completed audit took less than one second on one worker, below its 120-second cap. Maximum absolute quadratic difference was 1.51799e-6, within the combined tolerance. Same-author verification shares templates/precision; it is not external replication.

## The failed stop–start path

Late-noise regime 1, temporal correlation 0.95, replicate 1:

| Observation (zero based) | Single-core weight | Absent weight | Innovation quadratic | Label |
|---:|---:|---:|---:|---|
| 0 | 0.860551 | 0.125056 | 1033.42 | uncertain |
| 1 | 0.910492 | 0.0894702 | 920.071 | uncertain |
| 2 | 0.928060 | 0.0719136 | 927.730 | uncertain |
| 3 | 0.985540 | 0.0144537 | 925.871 | single |

Single was the leading structural alternative throughout, but remained below the frozen 0.95 weight requirement for the first three observations. All quadratics were below 1141.78, and neither position nor width spread failed its guard. The selected-history quadratics match the true current/previous geometry calculations to numerical precision. The important remaining alternative was absence, not a competing core. Conditional stationary-mode weight rose from 0.768863 to 0.998485 over these observations.

All subsequent observations were labeled single. Motion started only at observation 8. Therefore the failed path does not demonstrate a switching failure: it exceeded the allowed two uncertain observations during initialization. The overall MG174 screen remains FAILED. Delaying scoring would be a new endpoint, not a correction to that result.

Across all 64 intact primary paths, there were 21 non-single labels: 20 blocked by structural weight and one by the innovation guard. None was blocked by geometry spread. Other paths tolerated their one or two uncertain observations under the frozen criterion. These counts do not establish that geometry compression has no influence on structural weights; the audit identifies immediate decision guards, not every upstream causal mechanism.

## Increased-noise sensitivity

All 96 increased-noise observations were uncertain and failed the innovation guard. Ten also failed the structural-weight guard; overlaps are retained. None failed geometry spread. At the injected true current/previous geometry, all 96 quadratics still exceeded the limit, ranging from 1390.09 to 1724.98. A wrong fitted position is therefore unnecessary to explain this safeguard failure.

The field generator increased noise amplitude by 1.25, so its covariance increased by 1.5625. Dividing the true-geometry quadratics by that known covariance factor brought all 96 below the unchanged limit. In contrast, none of the 768 ordinary-noise true-geometry quadratics for static and stop–start paths exceeded the limit. This is a conditional diagnostic using known injected geometry and scale. It neither recalculates tracker class weights nor demonstrates that a deployable tracker can learn the scale or pass its other guards. No labels were changed.

## Recommended next stage

Develop a separately specified initialization and noise-adequacy instrument. Treat acquisition of a track separately from subsequent persistence, with a bounded acquisition period and explicit false-acquisition controls. Estimate background noise using a defined independent or spatially excluded sample, propagate its uncertainty, and test whether the covariance form is adequate before interpreting structural weights. A mismatch should produce an explicit unqualified measurement rather than a particle-death claim.

Before qualification, challenge the proposed method on independent fields with both increased and decreased noise, varying correlation, weaker cores, absent/damaged/split controls, and motion beginning during acquisition. Keep a fixed denominator that includes acquisition failures. Noise fitting must not absorb an extra core or suppress a missing species. Do not optimize a warm-up length to rescue this one path. This audit supports designing those tests; it does not authorize a new kinetic production allocation or qualify a revised instrument.

MG174 and MG163 scores remain unchanged. G remains unresolved. Full fields remain local in `local_raw/MG174_modes`; only compact audit records, bindings, code and this report are public. PI closure is pending.

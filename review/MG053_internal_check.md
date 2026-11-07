# MG053 — D1 prospective registration internal check

Same-author implementation/analysis check; no external review or replication. Richard approved the presented run after MG052's60-core-hour request. Optional external review was not requested with that execution instruction.

## Pre-production findings

- Original20 control checks, three corrected synthetic measurement checks and separate36-archive reconstruction remain preserved; frozen S1/S2 sources unchanged.
- New21 driver checks PASS:15 streaming/reference comparisons including adaptive stepping; failure, nonzero-exit, timeout and budget termination; known small-n descriptive statistics; plotting backend. Checkpoint streaming is bit-identical to the qualified kernel under the same random stream. Raw fixtures retained.
- End-to-end SYNTHETIC pipeline PASS: complete80/80 and incomplete79/80,480 independently recomputed summary sets each,26 descriptive contrasts, both plots, seeds, actual initial states/denominators, all80 raw archive hashes and explicit integer reconstruction. These algebraic fixtures are not scientific trajectories. The deliberately incomplete final fixture is labeled and kept separately from production.
- Disclosed combined development elapsed-based charge with startup allowances175.037486seconds, below900. Pipeline and controller traces retained.
- All80 production seeds unique; collisions0 against 2376 collected prior/development/warmup seed values. Exact matrix frozen.
- Measured forecast including25% slowdown and5seconds/path overhead: 150057.444099 worker-seconds (41.682623core-hours), within162000 (75% of216000) threshold. Forecast is not a guarantee for late fields. Eight workers maximum; one numerical-library thread each;3600seconds/job;216000 aggregate; no replacement or silent resume.
- Full raw checkpoint upper bound429440000bytes; local available disk 312970907648bytes. Atomic archives preserve all prior checkpoints at each save, including partial failures. Source/spec/reference hashes frozen in registration.json. Registration is committed/pushed before launch. Production archive upload is a later explicit retention step.

## Challenges and remaining limits

The pilot is descriptive with four paths/group: it does not establish stochastic1% recovery, precise probabilities, numerical equivalence or particle identity/lifetime. Both starts and both steps are fixed; no nonsignificance-to-equivalence inference. Supply-minus-CHEMO and fine-minus-coarse contrasts remain independent-sample descriptive differences. Small-n t intervals are approximate and non-simultaneous. RMS contrast separates background shift from spatial variation but does not qualify localization. T60 cannot answer T300/T500 behavior. Maintained inputs are explicit count sources, not a physical vacuum.

HIGH uses actual rounded base stocks multiplied by1000000, exactly as the qualified fixture; this is now explicit in prerun to avoid changing rounding semantics. All stock sums fit int64. New controller stops the whole batch on timeout or nonzero exit, correcting a limitation of the shared historical controller prospectively without altering that frozen file. Worker self-checks elapsed time between checkpoints as additional protection; parent owns strict runtime polling.

Automatic postprocessing writes evaluation/report and then the separate raw/summary check, including incomplete batches. A postprocessing failure is recorded and cannot be represented as a verified result. The independent checker shares NumPy/SciPy and the earlier explicit channel/node helper, not the sampler or D1 measurement/evaluator. Postprocessing has passed synthetic end-to-end qualification. Source/criteria cannot be changed after observing scientific outcomes without preserving a prospective amendment/deviation.

PI closure and scientific interpretation remain separate. Prior MG10 remains UNRESOLVED; closed C1/S1/S2 results are not rescored.

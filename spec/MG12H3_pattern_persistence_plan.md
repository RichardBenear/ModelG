# MG160 — prepared-pattern persistence study

Status: P0 preparation complete; P1 design and concrete budget proposed, not registered or running. PI closed MG158 as MG159 and requested this study. This remains H3 fixed-substrate recovery; H4 connectivity rewriting is later. No earlier result or particle detector is rescored.

## Primary scientific question

Does an initially resolved Model G pattern retain recognizable amplitude and shape under the finite-stock, strict-tick stochastic rules over longer evolution, or does noise, shape change or depletion disrupt it? Track translation separately from deformation. Distinguish persistence of a prepared feature from spontaneous formation, and distinguish a resolved pattern from physical particle identity.

Keep the 500-site periodic ring, original references and transport, Omega=10, finite A/B inputs and retained Z/W outputs. No replenishment/recycling, changed stock size, altered substrate or physical clock calibration. G/X/Y remain model species; token conservation is not identified with physical energy.

## Stages and decision boundaries

1. **P0 — measurement and saved-reference preparation (this packet).** Implement and test drift-aware feature measurements, mass-matched homogeneous controls, saved-prefix checks and cost inventory. No new random trajectories or integrations. Preserve the initial coherence-only implementation/results as development history. Broadening tests motivate reporting shape retention separately from template coherence. Archived deterministic fields are descriptive calibration context, not new stochastic persistence evidence or a rescore of their historical detector.
2. **P1 — early-coherence and noise pilot.** Four prepared and four mass-matched homogeneous independent pairs, each comparing ticks0.0001/0.00002 through T=0.1, with observations every0.001. This is ten times the accuracy-study horizon, not a long particle lifetime demonstration. It tests detection reliability under actual noise and whether a substantial longer-run investment is justified. All eight pairs are fixed in advance. Report both marginal paths and their paired differences; no population survival qualification from four prepared pairs.
3. **P2 — longer-horizon persistence design.** Use P1 to decide whether the detector needs a separately versioned, qualified revision, whether longer tick refinement is necessary, and whether exact-law engine acceleration is needed. Prospectively choose T=1 or a later model-dynamically justified horizon from reference behavior and measured cost, sample size, observation frequency, controls and precision. These are dependent decisions, not an automatic follow-on batch.
4. **P3 — long-reference persistence, if feasible.** Aim to bridge the historical post15/absolute25 reference window with adequate independent trajectories and tick controls. A 20-unit historical particle-tracking requirement cannot be replaced by T=0.1. Establishing spontaneous formation, dynamic substrate behavior or physical particles requires separate studies.

## P1 measurements and proposed operational rules

Use the X feature because the historical resolved-minimum diagnostics rely on X and its initial pattern is resolved at this count level. Retain complete G/X/Y fields, contrasts and resource diagnostics; X coherence alone cannot establish persistence of a G-centered physical particle.

Freeze the unrounded canonical initial X profile, subtract its spatial mean and apply the existing periodic Gaussian smoother sigma0.4 at spacing0.2. Center/smooth each observed X field identically. Compute projection amplitude and cosine similarity against every integer cyclic translation of the frozen template; record the maximizing shift, unaligned projection, residual and a multiple-phase ambiguity flag. Translation is permitted in the pattern score. Maximization over shifts creates a noise-selection effect, hence the matched nulls and no unvalidated false-positive rate claim.

Proposed **coherence**: positive best-fit amplitude at least0.5 times the actual rounded prepared initial amplitude and cosine similarity at least0.8. Proposed **shape retention**: global participation length between0.5 and2 times its prepared initial value, where length=h*(sum f²)²/sum f⁴. This is an effective extent of the centered profile, not an individual particle diameter. **Initial-pattern retention** requires both coherence and shape retention. Record these separately: a broadened field can remain correlated with its original template while failing shape retention. Do not collapse that distinction in reporting.

The 5% accuracy tolerance governs agreement between tick resolutions. The 50% amplitude,0.8 correlation and factor-two extent thresholds define a new operational persistence measurement; they are proposed engineering criteria, not a relaxed accuracy criterion or a physical particle definition. Their approval is included in the P1 packet.

For each prepared marginal trajectory, report whether the rule holds at every sampled frame; record the first failed-frame bracket. Retention through the last frame is right-censored at T=0.1, not infinite lifetime. Unobserved losses/recoveries between frames remain possible. Treat phase ambiguity as unresolved displacement/identity, not evidence of a unique particle trajectory. For homogeneous controls, report any threshold crossings as possible false detections or formation candidates; do not call these particle births. Four controls do not calibrate a small false-positive probability.

The homogeneous control preserves each of seven species' total counts exactly. Counts at its sites differ by at most one, with deterministic evenly distributed remainders and fixed species-specific cyclic offsets. Original reference denominators and connectivity remain fixed. This deliberately removes prepared spatial organization without changing total species inventories; it is not an open bath. Quantization residue and any initial template projection are retained and reported. No randomized or outcome-selected control initialization.

## P1 implementation, stopping and retention

Use the frozen MG156 paired engine in100 consecutive0.001 segments, continuing one RNG state. Each segment ends at both chains' tick boundaries, matching the already tested continuation law. Keep the one-million joint-event segment cap. Reconstruct every checkpoint from branch histograms and verify nonnegative integer counts/global conservation. Retain full observations and RNG states at every segment, metrics, sources, environment, seeds and resource accounts. Full numeric event/boundary records for replicate0 of each condition (two paired sentinels); compact complete records for the other six. Reserve2GB, preserve all historical archives.

All eight jobs must complete for a complete pilot report. Missing/failed/inconsistent archives mean incomplete, never survival. No scientific retries, discarded paths, threshold changes, sample-size adaptation or extending a failed observation grid. Before reporting, independently recompute amplitude/correlation/extent and every outcome from counts, reconstruct checkpoints, replay sentinels and record all analysis time including interruptions. Shared same-author dependencies are disclosed. External review is optional and not requested.

## Compute and approval

MG158 measured about548.268worker-seconds per paired0.01 trajectory under eight-worker contention. Linear projection for eight pairs to0.1, with15% allowance, gives about50440.6worker-seconds (14.0113core-hours), approximately1.75hours wall at eight workers. Homogeneous dynamics and longer history may change rates; this is a forecast, not a runtime guarantee. Proposed worker cap100000seconds total (27.7778core-hours),12000seconds/job, at most8localworkers with numerical-library threads1;1800 separately measured active parent/control/audit seconds. The forecast admission stop and no-retry policy remain.

A corresponding eight-pair T=15 study projects about9.52days wall at eight workers under the same linear model; this excludes larger replicate counts and may misestimate evolving rates. Do not launch it from this approval or describe T=0.1 as an adequate substitute. If full persistence is unaffordable, qualify a law-preserving acceleration rather than silently coarsening ticks or suppressing B transport.

The current instruction authorizes preparing this study. P1's proposed operational criteria and concrete run budget require PI approval under RESEARCH_PROTOCOL.md workflow3 before source/seed registration and any new stochastic draws. The driver requires P1_approval.json; none exists. Later stages need separate prospective packets.

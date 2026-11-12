# MG160 — pattern-persistence study readiness

The tick-accuracy study is closed as MG159, and the next study is now prepared. Its question is whether the initially resolved pattern remains recognizable as the strict-tick model evolves, while separating loss of amplitude, deformation and simple movement. Preparation established a workable X-pattern measurement, a control with identical species totals, and a bounded first pilot with finer-tick comparisons. This is a study design and measurement result, not evidence that a stochastic particle has survived. A full historical-duration run is presently expensive enough that the first pilot must establish detection reliability and feasibility before we commit to it.

## Proposed first pilot

- Four independent prepared-pattern pairs and four mass-matched homogeneous-control pairs: eight pairs, 16 marginal trajectories.
- Each pair compares ticks 0.0001/0.00002 through model time 0.1, observed every 0.001. This is ten times the accuracy-test horizon. The start is the rounded prepared donor from historical absolute time 10; the run continues that prepared state rather than simulating its formation.
- Track best-aligned X amplitude, template correlation, effective spatial extent, residual shape error and cyclic displacement/ambiguity. Retain complete G/X/Y fields and stock histories.
- Proposed pattern retention requires at least half the initial amplitude, correlation at least 0.8 and effective extent between half and twice its initial extent. Record coherence and shape retention separately. These operational thresholds are proposed for approval and are not a physical particle definition.
- Record first sampled loss as a time bracket; retention to 0.1 is right-censored. Four prepared pairs cannot establish a high population survival probability. Homogeneous crossings are possible false detections or formation candidates, not particle births.

The approved 5% numerical tolerance addresses agreement between tick resolutions. The proposed half-amplitude threshold addresses whether the initial feature remains recognizable during its evolution. The previous numerical accuracy result is unchanged.

## What preparation checked

Fourteen deterministic cases passed, including translation, weak/flat/inverted/split features, background offsets, periodic-phase ambiguity, broadening, time-grid handling and control mass. The 100-checkpoint writer fixture retained 101 frames with unchanged RNG state and zero events on an inactive small model; it tests file framing and accounting, not stochastic persistence. Separate FFT cross-correlation and explicit Fourier smoothing reproduced the translated-template measurements.

The prepared initial pattern has correlation 0.999631 with its reference and effective extent 2.58163 model-lengths. The homogeneous control has relative amplitude 0.000410735 and correlation 0.0505893; it is rejected initially. Every species total is preserved exactly and its control-site counts differ by at most one. All inspected saved accuracy-prefix snapshots retained the proposed pattern signature; those snapshots are development context, excluded from the new pilot.

Archived deterministic fields illustrate why both coherence and shape matter: the profile remains correlated through 15 post-preparation units, while its effective extent exceeds the new factor-two shape band between 9.5 and 9.75. This applies a new diagnostic descriptively to old fields; it changes no historical reference-matching verdict and establishes no new stochastic lifetime. The first coherence-only development code/results are preserved, and the broadened-shape counterexample and separate shape output are explicit. No threshold was fitted to new pilot outcomes; no pilot paths have been drawn.

## Compute and storage

Measured cost from MG158 is about 548.268worker-seconds per paired0.01 trajectory. With 15% allowance, this first pilot forecasts 50440.6combined worker-seconds (14.0113core-hours), or 1.75141hours on eight workers. Proposed ceilings: 100,000worker-seconds total (27.7778core-hours), 12,000seconds/job and 1,800 active parent/control/audit seconds. At most eight workers, one numerical-library thread each; no retries or sample replacement. Reserve 2 GB for complete compact decisions plus two full numeric sentinels.

For comparison, eight paired trajectories through 15 post-preparation units project about 9.51854days wall time on eight workers under a linear cost model. Evolving rates can change this estimate. Time0.1 is an early-coherence/noise pilot, not an adequate substitute for that longer persistence question. The staged plan proceeds to a new T=1 or later design only after this pilot; exact-law acceleration may be needed before the full historical horizon. A 20-unit historical particle-track criterion is not replaced.

## Status and limits

P0 COMPLETE; P1 READY FOR CRITERIA/BUDGET APPROVAL. No P1 registration or new stochastic trajectories exist. Existing laws and stock conventions are retained: finite A/B sources, retained Z/W, fixed substrate, no recycling or maintained bath. No emergent geometry, physical-clock, quantum-vacuum or energy-conservation claim is made.

The proposed first pilot does not yet offer calibrated false-positive rates, formal long-horizon tick equivalence, a population survival qualification or a unique physical-particle identity. Its finer-tick pairs and matched controls are designed to expose those measurement/numerical problems rather than assume them away. Outcomes of later stages determine their design and budgets.

Specification: `spec/MG12H3_pattern_persistence_plan.md`. Measurement: `src/mg12h3_pattern_persistence.py`. Preparation and checks: `src/prepare_mg160_persistence.py`, `src/check_mg160_measurements.py`, `src/check_mg160_writer.py`. Proposed production entry: `src/run_mg160_persistence.py` (requires P1_approval.json before registration). Frozen source/seed registration and PI approval are required before draws under RESEARCH_PROTOCOL.md workflow 3. External review remains optional.

# MG185 Phase 2 population ladder — DRAFT

No scientific run is registered or started. PI approval of these concrete criteria and caps is required under RESEARCH_PROTOCOL.md. MG184 is implementation development only.

## Question and scope

Does increasing population make the prepared G feature measurably sustained through model time 1000, alongside historical X persistence? Use the unchanged unsplit Phase 2 `tau_leap_euler_1d`, b=29, h=0.2, length=100, 500 sites, tau=1/1500. Open imposed sources and sinks remain explicit; this is not the finite-reservoir hypergraph model. Earlier Phase 2 results at T=100 suggest a population effect, but do not establish a universal threshold at Omega=3000 or G persistence.

## Numerical qualification, then population study

First run Omega=3000, prepared and homogeneous initial conditions, tau and tau/2, 12 independent trajectories per condition per resolution to T=100 (48 jobs). This increases the provisional development sketch from 3 to 12 before any scientific samples. Compare each run's fixed-frame G reference projection averaged over T=50..100. For each condition, the Welch 95% interval of fine minus coarse must lie within plus/minus 5% of the unrounded reference G norm (about 0.0210109). Both comparisons must pass with all jobs complete and no clips. Otherwise stop: inconclusive is not equivalent, and no automatic retuning or main run follows. This limited check does not establish all-observable, all-population, T=1000 accuracy. Translation can affect the fixed-frame diagnostic; the gate is deliberately limited to that observable.

If qualification passes, run Omega=100,300,1000,3000,10000: 30 prepared and 6 homogeneous Poisson controls each to T=1000 (180 jobs). Controls can nucleate; retain their complete histories and report that behavior. Do not replace trajectories after loss, splitting, clipping, or computational failure. All planned prepared runs remain in the denominator of 30. Stop a clipped trajectory as numerically unqualified, retaining its prefix; do not call that particle loss. A failed computational job stops further admission through the controller. No adaptive sample-size extension.

## Measurements and interpretation

Record every 0.5 time unit. Analyze the full horizon explicitly; the legacy default T=100 is forbidden. Retain historical X eligibility and strict sampled unique-original-track status separately, plus terminal ambiguity. The legacy event filter tolerates short gaps (20 time units); strict sampled status exposes this difference. Neither resolves events between observations.

For G, X and Y retain centered reference projections, relative amplitudes and orthogonal shape residuals. Report ten time-block means per run. Fixed-frame projections include every completed trajectory, even after physical failure. Common alignment uses only the historical X original unique track; G is never independently fitted to maximize its score. Missing or ambiguous alignment remains missing. Aligned averages are conditional, with valid-frame and contributing-run counts displayed. They cannot establish unconditional survival. Whole-field residuals include population noise and are descriptive, without a newly tuned cutoff.

Primary G output is the population/time dependence of fixed and conditionally aligned G projection, with prepared/control comparisons. Report per-run means and across-run Student intervals, not independent-frame intervals. These exploratory intervals have no familywise multiplicity guarantee. There is no full GXY binary persistence verdict from X alone. A positive result would support finite-horizon G visibility/retention in this model; it would not prove indefinite persistence or identify a physical vacuum source.

Historical X survival receives exact binomial 95% intervals only when all 30 jobs resolve. Otherwise report unknown counts and identification bounds. Even 30/30 gives a two-sided 95% lower bound near 0.884, not 99% reliability. Publish controls separately, not as replacements for prepared trajectories.

## Reproducibility and resources

Freeze source/reference hashes and fresh collision-checked seeds and publish before draws. Prepared counts are rounded from the unrounded reference; control initialization is independently Poisson sampled. Each trajectory uses consecutive 10-time-unit kernel calls with preserved final counts and independently registered segment seeds. This changes random stream layout from historical single calls, not the transition law. Keep logical sampling times and kernel floating timestamps; check agreement to 1e-7.

Eight workers maximum, one library thread each. Point estimate is in MG184 proposal.json, about 52 charged worker-hours (about 6.5 hours at eight workers), extrapolated from short benchmarks and subject to path/I/O uncertainty. Proposed hard caps: numerical qualification 8 worker-hours/2 hours wall; main 92 worker-hours/14 hours wall; each job 2 hours. Total 100 worker-hours. Stop on budget, source mismatch, or storage violation. No automatic retries or resumption that changes streams.

Full compressed counts, segment seeds, endpoints, sampling times, measurements and failures stay under local_raw/MG185_Phase2_G. Require 20 GB free initially; cap new raw at 10 GB. Public outputs contain compact manifests, summaries, figures and source, never the full raw archive. A separate same-author recomputation must check endpoints, measurements, denominators and intervals before reporting closure. PI closes the stage.

## Separate engine question

Omega=160 open-boundary bridge work is separate and does not gate this Phase 2 ladder. See MG186 draft. No Phase 2 result transfers automatically to the existing finite-reservoir hypergraph engine.

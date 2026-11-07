# MG055 — Source-law diagnostic and conditional corrected comparison

Prospective specification. PI instruction after D1 status: “looks ok, proceed with the next steps.” Preserve D1's complete registered result and all raw data; do not rescore, silently replace or reinterpret its biased HIGH fine-step source as physics. D1's posthoc caveat is separate from this new study. External review was not requested with this instruction; same-author checks are not independent replication.

## Hypothesis and correction

The installed Numba legacy binomial implementation computes q=1-p and q**n for small mean n*p. HIGH has n about4.2e12 and p about6.94e-14 at the fine step. Cancellation and/or integer-power rounding can distort its marginal. Integer ledgers cannot detect incorrect random probabilities. Test the original and revised scalar samplers against the exact finite-binomial CDF before corrected trajectories.

New finite-A helper retains the original np.random.binomial call unchanged for n<=1e9. For larger n, p0/1 are exact; otherwise require n<=2**53, p<1e-6 and n*p<=30. Count Bernoulli successes by geometric spacings floor(log1p(-U)/log1p(-p)); stop before a spacing exceeds remaining trials. This samples a finite binomial law, not a Poisson approximation. Original frozen D1/S1/S2 sources remain untouched. New D2 model changes only finite-A source sampling and imports; G/X/Y reactions, joint B/X exits, transport, stocks, denominators and time semantics stay as qualified.

## Source diagnostic frozen before samples

Two primary cases: n=4200000000000 and p=14000*tau/n, tau1/24000 and1/48000. For each, generate2**24 independent draws from three samplers: original Numba legacy, revised finite-source, and NumPy Generator.binomial. Six fixed streams, stage5501, seed_for(10,case,sampler,5501,0). Store65-bin integer histograms (0..63, >=64), raw moments and seeds; no sample selection or replacement. Histogram counts are sufficient raw inputs for CDF/mean decisions. Compare exact CDF via log1p-based finite-binomial PMF recurrence, with a separate SciPy binomial CDF cross-check and independent histogram recomputation.

Use simultaneous99% DKW bounds across all6 comparisons: epsilon=sqrt(log(2*6/.01)/(2*2**24)). The revised helper must have upper CDF distance<=.0008 in both primary cases. NumPy is a cross-check, not the definition of truth. Original-law discrepancies and all means/variances are reported whether or not the legacy sampler meets this diagnostic margin. A lower CDF bound above0 identifies a statistically resolved law discrepancy on these observables. No outcome threshold is tuned after samples.

Supplementary source fixtures use2**18 draws and one-sample99.9% DKW upper<=.02 for geometric-binomial n3,p.2 and n40,p.3, large n with mean30 and half-depleted HIGH at fine p. Additional exact bounds/domain and low-stock identical-RNG checks qualify implementation. Supplemental samples cannot replace a primary failure.

Before full production: corrected chunk/stream must agree bit-for-bit for both HIGH and CHEMO H/S short paths; finite low-stock fixtures and full non-HIGH D1 arms must match the frozen D1 source under the same seeds. Exact per-site ledgers, nonnegativity, maintained stock invariance and max exit probability<=.25 apply. Retain failures. A failed mandatory diagnostic/qualification stops corrected production; fix code prospectively with versioned record, never replace scientific failures.

## Conditional corrected comparison

Only after all revised-sampler and implementation qualifications pass:32 new paths, two arms CHEMO/HIGH, H/S starts, steps1/24000 and1/48000, four independent replicates/group, N500,h.2,Omega1000,T60, checkpoints0..60. Fresh seeds seed_for(10,2*arm+start,step,5503,replicate). All32 unique and disjoint from registered D1/development streams. Both fresh maintained-bath controls and revised HIGH paths are generated; original D1 controls are context, not silently pooled with the new denominator.

All actual integer stocks, rounding convention, background-separated measurements, approximate non-simultaneous t intervals, descriptive step/supply contrasts and restrictions on particle/equivalence claims follow D1. This remains estimation with n4, not powered1% recovery, joint-law validation, quantum evidence or a particle-survival test. Source marginal qualification is stronger than D1 but is not full trajectory-law equivalence. Interpret corrected source rates against their autonomous expected rate; HIGH source deficit must not be explained as reservoir physics without checking the new marginal. No blanket biological/physical particle conclusions from contrast.

## Authorization and caps

Use a lead allocation within the already-approved60-core-hour envelope, not a new larger budget: D1 used33.333427; new diagnostic/qualification ceiling900worker-seconds (.25core-hour); corrected production ceiling24core-hours. Worst combined57.583427core-hours remains below60. No increase in authorized workers: local<=8, numerical library threads1,3600seconds/job, owned parent polling.25seconds; stop the entire new batch on failure/timeout/cap, no automatic retry/resume.

Measure all8 new arm/start/step classes at T.1 with disclosed development streams stage5502, forecast matched cost*600*1.25 plus5seconds/path. Require the32-path total<=75% of24core-hours (18core-hours) before registration. Do not expand the denominator or cap automatically. If the forecast or qualification blocks, publish the blocking result and revised concrete plan; no scientific production until resolved prospectively.

## Registration and retention

Commit/push diagnostic source, exact6 primary stream/case matrix and plan/hashes before diagnostic sampling. Freeze production code/spec/reference hashes, jobs, measured forecast, environment, remaining-budget allocation, qualification and retention before corrected trajectories. Reuse the qualified owned-controller/checkpoint mechanism and independently recomputed analysis under prospective D2 filenames; original frozen D1 files remain untouched.

Store complete primary and supplementary histograms, moments, scalar/stream regressions, timing outputs and environment in results/MG10D2. Retain all32 full checkpoint/event/external-input archives, atomic progress files, configs/ledgers/source hashes, evaluation and independent audit. Full raw upper bound171776000bytes; >=2GB free before registration; commit/upload compressed per-path files in manageable batches and verify hashes/reachability. D1's80 raw files are separately archived and never replaced. Histogram sufficient statistics apply only to scalar-source tests; full spatial decision checkpoints remain retained.

Diagnostic outcomes are reported explicitly. Corrected production is COMPLETE DESCRIPTIVE PILOT only if all32 paths/accounting complete, otherwise incomplete with denominator preserved. Independent checker must pass before reporting recomputed results. All source/trajectory claims limited to their tests; PI closure separate.

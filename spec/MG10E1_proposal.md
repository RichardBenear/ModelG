# MG069 — MG10-E1 localized A-source intervention pilot

Status: DRAFT scientific protocol. PI requested the next stage after P3 closure; prepare/qualify/benchmark now and bring the new intervention choices and measured budget to the PI before production. No earlier trajectories rescored. Every stage result will begin with the required one-paragraph project-level summary.

## Project-level purpose

P3 established a baseline of persistent reference-like patterns under maintained supply. E1 asks whether the seeded contribution to the field depends on where the model delivers its A resource, beyond the background response that a source gate produces even without a seed. This advances from passive observation to controlled interventions in the implemented model. It does not derive a physical source for A/B or identify a failed reference match with particle death.

## Preparation and counterfactual branches

N500,h.2,Omega1000,P17, full maintained A/B baseline and unchanged transport. Two numerical maxima1/24000 and1/48000;8 independent blocks per step. Each block has separately seeded H and S parents (same initial definitions as accepted P2 FULL/P3). Evolve each parent to the fixed time5, retaining all21 quarter checkpoints and every physical species, including accumulated Z/W. Do not select parents using a match, future persistence, amplitude or favorable outcome. Thus32 parents total.

From each exact time5 state fork all four post-intervention arms for20 units with81 quarter checkpoints;128 branches. Each arm gets its own independent future RNG stream, not common random numbers. Pairing comes from the shared prepared parent, not shared future noise; H/S parents in a block have independent preparation noise. Preserve parent-state hashes and all seven species/denominators/flags exactly; reset only post-window reaction/hop/external counters, retaining the parent records separately. Posttime0 is absolute time5. No path replacement, eligibility filtering or pooling with closed studies.

Parent seed_for(10,start,step,6903,block), start0H/1S,step0coarse/1fine,block0..7. Branch seed_for(10,start,step,6904,block,arm+1), arm0..3 below. Separate warmups/development stages; freeze/verify all streams before production. All32 parents must complete before dispatching branches; an implementation/runtime failure stops the batch. Primary analysis uses complete balanced H/S-by-arm blocks only and reports missing blocks, not zeros. No block selection based on the response.

## Four A-source arms

Only the R1e A-to-G intensity is multiplied by a fixed per-node activity vector. Stocks remain maintained, B availability/propensities/transport unchanged directly, and all actual external A/B input events recorded. This is a source-rate gate, not deletion of A tokens, a B cut or a removed graph edge. Downstream B consumption can change as a mediator; only expected A delivery is dose-matched, not total downstream flux or physical energy.

|Arm|A-source activity|
|---|---|
|SHAM_UNCHANGED|1 everywhere; same fork/gate handling, no rate change|
|LOCAL_A_OFF|0 on51 nodes with integer offsets-25..25 around fixed index250 (origin);1 elsewhere|
|REMOTE_A_OFF|Same51-node periodic patch around index0 (opposite side of ring);1 elsewhere|
|UNIFORM_A_REDUCED|449/500=.898 everywhere|

The local/remote/uniform vectors each sum449; sham sums500. Patch endpoints are coordinates+/-5 inclusive; discrete node count, not continuum width, defines dose. Patch centers remain fixed, never follow a measured particle. H receives the same spatial gates. The remote arm is not assumed causally disconnected: transport, especially Y, can communicate its effects. Sham is the no-cut control. These arms compare locations/distribution of delivery rather than assume an A cut must kill a pattern.

Sham and untouched continuation are collapsed prospectively ONLY because the all-one gated kernel must reproduce the frozen ungated kernel bit-for-bit on registered representative states/seeds, including prepared-like states with nonzero Z/W; source-level changes are confined to unit activity multipliers. This is an identity check, not two independent stochastic controls. If identity fails, stop and correct/restore an independent untouched control before registration; no silent omission.

## Primary field estimand and background controls

Keep raw per-field means, source/stock/history totals and full physical-coordinate profiles for every arm/start. Apply unchanged sigma.4 centering/smoothing and fixed-origin P17 template projection to G/X/Y; do not independently realign the primary profiles. Spatial mean shifts are separately reported rather than hidden by centered projections. Matching loss or movement does not enter primary eligibility.

For each branch compute the signed fixed-origin X projection's trapezoidal average over posttimes15..20 (21 quarter samples, endpoints half-weighted). Let z(start,arm) be this value for a given block/step. Primary locality interaction is

    D = [z(S,LOCAL)-z(S,REMOTE)] - [z(H,LOCAL)-z(H,REMOTE)].

This subtracts the unseeded forced-background response before comparing gate positions. It is a seed-by-gate-location interaction of the specified observable, not an unqualified particle-destruction statistic. G/Y versions and local-minus-sham, remote-minus-sham, uniform-minus-sham seed interactions are secondary descriptive contrasts, with full profiles retained. Uniform is an expected-A-dose/global-response comparator, not a guarantee of identical local background. H control supports background matching within each intervention; nonlinear effects remain part of the modeled response. A gate could affect generic localized perturbations too; E1 does not establish a response unique to physical particles.

Statistical unit: one complete block (H/S parents and all four arms),8 per step, not32 independent branches. Report paired means and nominal97.5% two-sided Student-t intervals per step (critical t at.9875,df7), nominal95% two-step family under the t approximation. Small-n/nonnormality limits stated; this is not guaranteed distribution-free coverage. Same-sign intervals excluding0 at both steps may be described as direction-consistent pilot evidence for this model interaction. Otherwise interpretation is unresolved, including contradictory step signs. No equivalence from nonsignificance, prespecified sign, minimum-effect claim, threshold tuning or pooling. COMPLETE INTERVENTION ESTIMATION requires all32 parents/128 branches and valid accounting/audit; completion is not contingent on a favorable effect.

## Secondary pattern diagnostics

P1 classifications/positions and conditional initial-reference continuity are descriptive. For S a parent match/phase at time5 can anchor the initial post-window reference identity; an unmatched parent is INITIAL_UNRESOLVED and is still retained. H has no designated initial seed. Never choose/move patches using matcher positions or future behavior. A changed shape/position can reduce the fixed-origin projection or fail P1 while a physical/model excitation still exists; do not call that death. Full-domain profiles and mean offsets are needed for interpretation. No native vacuum/quantum claim.

## Preparation before production

New E1 file is a literal prospective D2 kernel copy with only per-node R1e multipliers; frozen kernels unchanged. Mandatory: all-one same-seed bit-identical old/new states/stats on four small fixtures, full H/S and prepared-like nonzero-Z/W states; fork copies all physical species and original denominators without aliasing, starts counters at0. Spatial masks/dose/periodicity exact; off sites have zero source/input events; other ledger/probability/nonnegative checks hold. Invalid masks reject. Scalar isolated maintained-source tests at activities0,.898,1, Omega1000,tau1/24000,R32768 each: exact Poisson CDF and simultaneous99.9% one-sample DKW upper bounds<=.03 across3 tests; activity0 also exactly zero. No joint-trajectory-law claim from these marginals.

Qualify trapezoidal late-window integration on known linear/constant curves, primary/secondary block contrasts and invariance to adding the same arm-specific forced response to H/S, and nominal t intervals via separate sum/variance/CDF inversion. Independent source/stoichiometric reconstruction/mask/dense-projection/statistical path must agree. No physical endpoints may tune criteria. Preserve failures/versioned corrections; any blocking issue stops production.

Development cap900worker-seconds,<=8/library threads1. Four disclosed short parents at.25 (not5), stage6901, then sixteen arm/start/step timing paths at.25,stage6902; these are development only, not selected/prepared production parents or E1 estimates. Retain complete inputs/outputs. For parent cost reuse the already measured same-machine P3 FULL T1 two-replicate/class records, freezing their hashes; initial factory/ungated kernel unchanged. For each production parent forecast mean old T1*5*1.25+15seconds. For each branch use matched new.25 timing*80*1.25+15seconds. Require total<=75% of proposed36-core-hour cap (27core-hours); if not, bring an explicit revised plan, no automatic expansion.

## Proposed execution budget and retention

NEW production ceiling129600worker-seconds (36core-hours), subject to measured forecast and PI approval. Parent timeout600seconds; branch timeout1800seconds. Sequential phases32 parents then128 branches, aggregate charge includes both, no hidden reset; owned.25second polling,<=8 local workers/librarythreads1, whole-batch failure/timeout/cap stop, no automatic retry/resume. Retain all partial paths and original denominators; incomplete blocks missing rather than physical failures. Forecast is not a late-dynamics bound.

Raw upper bound parents32*21*22*500*8=59136000bytes plus branches128*81*22*500*8=912384000bytes, total971520000 before lossless compression. Store exact activity vectors, parent/fork hashes, all checkpoints/populations/events/hops/inputs/physical Z/W, counters, seeds/config/environment/source/spec/reference hashes, progress/failures/phase ledgers, field metrics, paired effects and separate audit. Verify>=4GB free; upload raw archives in manageable batches and verify hashes/Git reachability. No sole temporary/hash-only data. Original studies untouched.

After PI approves the concrete measured protocol/budget, qualify/freeze full production/phase controller, evaluator/paired checker and exact job/seed matrix, commit/push before science. External review optional per PI. Result starts with one high-level project-contribution paragraph, grounded in actual audited findings and limits (MG068), not an invented causal conclusion. PI result closure separate; no further stage automatic.

# MG163 â€” longer-horizon pattern persistence

Status: prepared proposal; production criteria and budget await PI approval. MG162A/B are PI CLOSED. No MG163 Model G draws or production seeds have been generated.

## Purpose and interpretation

Test whether the local core detected in MG162B remains recognizable for ten times longer under the same finite-source dynamics. This tests sampled operational persistence, not quantum particles, physical lifetime, continuum convergence, or a population survival probability. Source depletion is part of the experiment: A/B remain finite, no replenishment or parameter adjustment. T=1 is post-preparation model time, with no conversion to seconds.

## Fixed design

Eight fresh independent seed pairs: four prepared and four species-total-matched homogeneous controls. Both marginal ticks 0.0001 and 0.00002, canonical 500-site Omega=10 system, same prepared state and homogeneous construction as MG162B. All old paths excluded from denominator. New deterministic seed namespace 16300, collision check against all previous registration/plan records, actual seeds frozen after approval and published before draws. No selection based on old or new fitted outcomes.

Run 1000 consecutive segments of 0.001, carrying populations and RNG state across every boundary. Observe at all 1001 times 0..1 inclusive. Tick availability/reservations obey the same strict-tick law. No arriving/product token may act again within its tick. Keep the frozen paired engine and one-million-event segment cap. No scientific retries, replacements, adaptive sample increase, early stopping for positive/negative labels, or extended horizon beyond 1. Implementation/resource failures stop and preserve the incomplete denominator.

## Frozen measurement and decisions

Reuse MG161 template bank, calibration, fitting and classifier without refitting noise, amplitude, grid, thresholds or initial template. The MG162B engineering boundary advisory remains unchanged: core-compatible frames with relative amplitude <=0.6 or fitted scale within 0.125 of 0.5 or 2 are boundary-sensitive. This is not a confidence interval. Record original labels and advisories separately.

Primary finite-cohort criterion: every prepared marginal path core-compatible at every observation with no boundary advisories, and no homogeneous marginal path core-compatible at any observation. Otherwise clear-case separation is NOT ESTABLISHED; do not infer particle death or falsely classify a control crossing as necessarily instrument error. Report per-path first noncompatible bracket, first advisory time, all label counts and fitted amplitude/scale/position/residual/secondary-score histories. For first failure at initial time report initial absence, not a later loss. Retention is sampled, so between-frame loss/recovery remains unresolved. Preserve descriptive historical global metrics without using them to replace this criterion. Tick marginals and frames are correlated; only four independent pairs per condition, no small error-rate or population-lifetime claims. Separate retention and control outcomes even when the combined criterion fails.

## Compact local retention

All 1000 checkpoint segments/job: initial/final counts, reservation boundary endpoint, branch/channel histogram, RNG before/after, source hashes, event count. All 1001 observation fields/job retained locally. Two fixed sentinel jobs (prepared replicate0 and homogeneous replicate0) additionally retain full numeric event/reservation traces for segments with index modulo10 equal0: 0,10,...,990 (100/job, 10% fixed before draws). Every segment is simulated identically regardless of trace selection. Sparse sentinels narrow causal event replay coverage to the preselected subset; they do not constitute a complete event audit of all segments. Checkpoint/histogram/RNG/accounting checks still cover every segment.

Full fields/events remain in ignored local_raw/MG163_long_horizon. Public package: code, frozen design, seeds, compact scores/plots, budget/audit summaries and SHA-256 local manifest. No automatic public full-data upload. Keep a verified immutable local-only archive when practical; never claim immutable retention without creating it. Reserve 4 GB local; preflight free space >=8 GB. Stop and preserve incomplete data if local payload reaches the 4 GB ceiling. Write any failed resource attempt to ledger.

## Measured budget proposal

MG162B measured 45083.2582463 worker-seconds, 5689.9340649 wall-seconds for 0.1. Linear T=1 baseline: 450832.582463 worker-seconds =125.231 core-hours, ideal eight-worker wall15.654 h (scaled observed wall15.8054 h). A 25% planning allowance yields563540.728079 worker-seconds =156.539 core-hours, ideal wall19.5674 h. Rates may change as finite sources deplete; this is an extrapolation, not a guarantee or measured full-horizon timing. Budget ceiling800000 charged worker-seconds (222.222 core-hours),100000/job, max8 local workers/library-thread1. Active preparation/audit/report allowance3600 seconds, separately recorded, excludes waiting for production. Admission forecast563541 is below75% of ceiling600000. No external reviewer required unless PI requests one. This exceeds the usual <=3h batch preference and is a dedicated allocation, not a short-batch authorization.

Measured local MG162B traces1220863080bytes; other data40863189bytes. Tenfold horizon with 10% sentinel trace coverage projects1220863080+10*40863189=1629494970bytes (1.62949GB), before growth/metadata allowance. Reserve4GB. Sparse sentinels reduce storage, not the number of scientific events or the runtime by a factor10.

## Prelaunch implementation and validation gates

A new parameterized writer is required: historical MG160's 100-segment loop and11980-second deadline cannot be reused unchanged. New worker uses1000segments, explicit per-job deadline99980seconds, selected trace indices, source-bound configuration and predeclared4GB local storage guard. Existing shorter drivers/specs stay unchanged. Adapt audit assumptions from101 to1001frames and trace-all to fixed selected trace indices; verify all8000paired segment checkpoints (16000marginal endpoints), sampled causal histories, exact RNG continuity and species accounts. Independently recompute local fitting/search/classification/advisories and the finite verdict. Record elapsed analysis attempts including failures/interruption.

Before any production: source/configuration syntax and no-draw guard checks; deterministic inactive-model tests of1000segment framing/time/RNG/reservation continuity and exact trace schedule; reference regression of old100segment mode; explicit budget/storage guard tests; unchanged engine hash binding. Such tests must not consume unregistered active-model scientific paths. Any active-model microbenchmark needs its own disclosed allowance. Publish approved specification, source hashes and seeds before launch. No implementation-validation success is claimed by this design packet.

## Approval requested

Approve the fixed eight-pair T=1 design, unchanged measurement and strict criterion, deterministic10% sentinel retention, 800000worker-seconds/100000perjob plus3600activeanalysisseconds,4GBlocal archive/max8workers. Preparation permission is not production permission. After approval finish required implementation qualification and publish registration; blocking validation failures stop launch. No further scientific sign-off is needed merely to execute the approved, validated packet. PI closes the result after reporting.

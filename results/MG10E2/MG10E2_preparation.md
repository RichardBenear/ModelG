# MG073 — MG10-E2 preparation and proposed production

Status: PREPARATION PASS; DRAFT study, production not approved or registered. 2026-11-08.

E2 is ready to test whether returning local A supply changes a seeded pattern's field response compared with keeping supply interrupted, and whether that response depends on the initial seed history. Preparation establishes that the timed shutoff/restoration works without resetting physical state, event counters or random streams, and that new measurements can distinguish some translations from shape changes. It supplies no scientific recovery result yet. A recovered reference-like field would not prove survival of the same particle, exclude formation from residual seed information, or demonstrate quantum memory.

## Reviewable design

Fresh32 parents: fixed preparationT5, H/S, two numerical steps,8 balanced blocks/step, no selection by visibility or future result. Each parent forks four ways into128 branches, postT25:

| Arm | Post0..5 | Post5..25 |
|---|---|---|
| Unchanged | Full A input | Full A input |
| Continuous local shutoff | Local51 nodes off | Same nodes off |
| Local shutoff then restoration | Local51 nodes off | Full A input |
| Remote shutoff then restoration | Opposite51 nodes off | Full A input |

Only A-to-G input rates change; B and transport laws retain their definitions. Record all external A/B events and every physical species. Parent states/denominators preserved; branch counters alone start zero. Right-continuous switch at post5, no reseeding. Restored and continuous local arms intentionally have different total A delivery: the question is restoration, not an equal-dose location comparison. Local and remote pulses have matched expected delivery. All arms get corresponding unseeded backgrounds.

Primary signed fixed-origin X-template projection, trapezoidal post20..25 mean. Balanced-block interaction=(S_local_restore-S_continuous)-(H_local_restore-H_continuous), n8/step. Nominal97.5% two-sided t interval per step, approximate95% two-step family under the approximation. Same-sign intervals excluding0 at both steps support direction-consistent restoration response; no assumed sign, full-recovery/equivalence threshold or particle-mortality score. Missing/prefix branches excluded as missing balanced blocks; not0. Completion independent of favorable effect.

Secondary full trajectories, G/Y effects, restored-versus-sham deficits, remote comparator, field means/source histories/P1 diagnostics and grid-shift-optimized X fits with coupled G/X/Y projections/residuals. Flat fields have a deterministic position convention, not a detection; optimized projections have noise/selection bias and use the same H controls. Movement, deformation, multiplicity, seed-history response and physical lineage remain distinct. These diagnostics are not new tuned classifier thresholds. Recovery to sham is estimated, not certified by nonsignificance.

## Preparation evidence

Frozen source/spec/development plan published at048d705 before sampling. 11 assertions PASS, including8 arm-by-step schedule oracles/identity checks, invalid inputs,7 synthetic field fixtures and temporal/paired/missing-prefix fixtures. Eight shortT.005 qualified schedule paths use a scaled.002 switch; sixteen timing paths useT.25 and.125 switch. No productionT5/25 samples or future-seed outcomes inspected.

Separate checker agrees on24 physical raw archives, dense Gaussian smoothing/direct all500-shift scores, fixed/movable projections, means/norms/residuals, explicit time weights, contrasts and CDF-inverted Student intervals. Full128 synthetic block inputs give n8/8; absence or an available incomplete prefix gives n7/8, never a zero observation. A per-leap schedule oracle independently chooses before/after rates while sharing the frozen event kernel. Same author/shared constants/reference/libraries disclosed; not external replication.

All320 candidate production/warm values and40 development values unique/disjoint; no collision with3081 distinct previously recorded integer seed fields across1813 JSON files. Timing starts use frozen E1 development parents only; actual E2 parents will be fresh. No failure, replacement, threshold revision or rerun.

## Measured budget and storage

Controlled development157.545613worker-seconds; separate audit allowance162.907576seconds total, below900. Parent forecast uses frozen E1 actualT5 simulation timings; branch forecast scales newT.25 timings toT25. Both include25% slowdown allowance and15seconds/job.

Total forecast106409.584492worker-seconds = **29.558218core-hours**, about**3hours42minutes on8workers** before scheduling imbalance. Proposed NEW**40-core-hour cap**; forecast below30-core-hour/75% reserve threshold, narrowly. This is a short-state forecast, not a bound on late behavior. Runtime cap and whole-batch stops are binding, not promises of completion.

Parent600seconds/branch2400seconds,<=8 local workers/librarythreads1,owned.25second polling. All parents must complete before branch phase; accumulated parent charge carried over. No automatic retry or replacement. Full uncompressed checkpoints bounded1196800000bytes plus metadata, >=4GBfree verified. All complete/partial raw profiles/events/hops/inputs, schedule/parent/fork/seed/hash/config data, ledgers, metrics/contrasts and separate audit retained. Public scientific-result/full-archive publication is a requested part of production approval, to make the destination explicit and avoid the earlier E1 publication ambiguity.

## Decision and next action

PI review requested for final E2 protocol,40-core-hour local cap and public ModelG results/full-raw retention. RESEARCH_PROTOCOL.md workflow3: “The PI approves modeling choices, criteria and budget.” Current authorization is development; final numbers/cap are new. No production until approved, full streamed production/evaluator/controller/independent audit implementation and complete/incomplete raw-file pipeline qualification are frozen and registered/pushed. The existing schedule/core measurements and balanced statistical engine are developed and qualified; production orchestration is a remaining pre-run step, not claimed finished here. No scientific run has started. Stage result will include the required high-level project summary; PI closure separate.

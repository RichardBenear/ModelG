# MG142 — Actual-ring multi-tick pilot results

The optimized engine successfully carried all tested actual Model G ring states through multiple ticks, including the rule that products and arriving tokens become eligible only on later ticks. Every registered integrity check passed across 240 ticks and approximately60.8065million accepted events. This establishes short multi-tick execution feasibility and useful cost measurements. It does not establish particle persistence or tick-size convergence; the next scientific priority is to determine how the strict tick affects the intended kinetics before investing in a long particle run.

## Registered conclusion

**FEASIBILITY AND INTEGRITY CHECKS PASSED; PI closure pending.** All 16 trajectories completed through model time 0.01. Eight finite-stock 500 site starting states were each run at tick 0.001 and 0.0005, giving 10 or 20 successive ticks respectively and observations every 0.001. Seeds/registration MG141 were published at `ab9345b` before production. There were no failures, retries or budget stops. This was not a registered particle-preservation or statistical refinement test.

Source/initialization remain the frozen MG134states: two donor origins,two B-network interventions, Omega1/10, with disclosed donor concentration-to-integer rounding. Reservoir transport was retained; stocks were not replaced with maintained controls. All frozen source and input hashes remain unchanged.

## Integrity and independent recomputation

Checks at 10,000 event intervals and tickends compared tree leaves and totals to freshly calculated rates, with representative selection probes. All passed relative tolerance 1e-10 plus absolute tolerance 1e-6. Largest recorded relative tree-total drift was **1.91545e-11**. Saved outputs reconstructed all initial/final/available tick arrays and observation frames, with exact integer conservation and nonnegative availability. Availability resets only at a tick boundary.

The [separate saved-output recomputation](separate_recomputation.json), reproducible with `src/recompute_mg12h3_multitick.py`, sorted channel IDs and used another coordinate-aggregation path to reconstruct the 240 ticks and every frame. It confirmed seeds, timing bounds, event counts, checkpoint totals and immutable inputs. This shares compact channel compilation and raw inputs; it is same-author recomputation, not external replication. Worker leaf/selection checks are not independently rederived from saved totals.

## Per-trajectory observations

Runtime below is the worker's measured job-body duration, including writing/auditing. The final column is a descriptive normalized active-field RMS change: for G,X,Y, divide final-minus-initial concentration by max(1,initial spatial RMS for that species), then take RMS over the three species and 500 sites. It is not a particle detection metric or a scored stability threshold.

| Case/tick level | Ticks | Accepted events | Job-body seconds | Normalized active RMS change |
| --- | --- | --- | --- | --- |
| o0_a0_Omega1_d0 | 10 | 689982 | 16.7389 | 0.584848 |
| o0_a0_Omega1_d1 | 20 | 692250 | 16.9265 | 0.529609 |
| o0_a0_Omega10_d0 | 10 | 6.9048e+06 | 167.541 | 0.302951 |
| o0_a0_Omega10_d1 | 20 | 6.915e+06 | 168.256 | 0.204867 |
| o0_a1_Omega1_d0 | 10 | 690032 | 16.3024 | 0.608576 |
| o0_a1_Omega1_d1 | 20 | 690934 | 16.9256 | 0.532145 |
| o0_a1_Omega10_d0 | 10 | 6.90687e+06 | 167.113 | 0.300031 |
| o0_a1_Omega10_d1 | 20 | 6.91528e+06 | 169.641 | 0.208266 |
| o1_a0_Omega1_d0 | 10 | 690773 | 15.4027 | 0.560707 |
| o1_a0_Omega1_d1 | 20 | 690798 | 15.7217 | 0.525276 |
| o1_a0_Omega10_d0 | 10 | 6.90276e+06 | 162.389 | 0.301496 |
| o1_a0_Omega10_d1 | 20 | 6.91458e+06 | 164.006 | 0.207437 |
| o1_a1_Omega1_d0 | 10 | 691611 | 16.6362 | 0.591384 |
| o1_a1_Omega1_d1 | 20 | 691749 | 17.111 | 0.510589 |
| o1_a1_Omega10_d0 | 10 | 6.9033e+06 | 151.823 | 0.301466 |
| o1_a1_Omega10_d1 | 20 | 6.9158e+06 | 152.432 | 0.211073 |

Level d0 uses 0.001; d1 uses 0.0005. Each pair has independently seeded paths. Their final normalized RMS differences are:

| State | Coarse/fine descriptive difference |
| --- | --- |
| o0_a0_Omega1 | 0.762535 |
| o0_a0_Omega10 | 0.263607 |
| o0_a1_Omega1 | 0.75401 |
| o0_a1_Omega10 | 0.26407 |
| o1_a0_Omega1 | 0.713211 |
| o1_a0_Omega10 | 0.25944 |
| o1_a1_Omega1 | 0.748153 |
| o1_a1_Omega10 | 0.266663 |

These differences combine sampling noise and any tick-size effects. One path percondition cannot estimate their separate contributions or support confidence-based convergence claims. Lower fluctuation magnitude at Omega10 is descriptively consistent with larger populations, but no noise-scaling law was tested.

## Tick-resolution issue to investigate next

An analytical transport-only diagnostic is available without rescoring the pilot. For a token with total exit rate lambda and no reaction competition, the strict one-use rule gives move probability 1-exp(-lambda*delta), while the continuous-time expected number of jumps in the interval is lambda*delta. The ratio **(1-exp(-lambda*delta))/(lambda*delta)** multiplies the one-tick mean transport increment relative to the continuous generator's initial drift.

For ring Y transport, D=300 and two unit-weight outgoing links give lambda=600. The ratio is **0.751981** at delta0.001 and **0.863939** at delta0.0005. This is an exact transport-only initial-drift calculation, not a measured bias of the full reaction/transport system, a finite-horizon diffusion error bound, or proof the pilot's paths are wrong. Reactions, catalyst reservations and spatial modes complicate the full dynamics. It does show why tiny mixed fixtures with much slower rates cannot alone establish actual-ring tick resolution. Finer ticks retain one use per token per tick.

## Runtime and retention

Elapsed wall time was **222.567seconds (3minutes43seconds)**. Worker charge **1474seconds**, parent audit **9.37436seconds**, total **1483.38seconds (0.412049core-hours)**, below the 7200-second budget. Up to 8 local workers, one numerical-library thread each. Start/finishUTC and ledger are retained.

Full production archive: **545files**, approximately **483.301MB**, including every accepted channel/time, all tick states, frames, checks, seeds and metadata. See [raw manifest](raw_manifest.json); publication verification is recorded after upload. File volume and runtime become material for long trajectories; the short pilot does not establish a linear full-run forecast.

## Recommended next study

Before extending the horizon, perform a tick-resolution study: derive transport-only expected transitions exactly on the actual substrate and compare them with continuous-time reference evolution across finer ticks, then design replicated coupled-model checks at selected resolutions. Use pure-transport deterministic calculations to choose meaningful tick sizes and observables before spending on noisy coupled trajectories. Keep environment transport and strict reservation intact. A follow-up draft is in `spec/MG12H3_tick_resolution_plan.md`.

The pilot and earlier MG139qualification await explicit PI closure. No new production is registered or started. No physical vacuum-source, quantum-entanglement, persistent particle or three-dimensional particle claim follows from these results.

Publication verification passed for all 551 production and decision files at commit `a58bf9b`. Binary archives match byte for byte. See [verification record](production_retention_verified.json).

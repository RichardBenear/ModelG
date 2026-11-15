# MG184 Phase 2 G study readiness

The tests are ready to investigate whether larger token populations make G remain visible over much longer model times. They protect against two misleading conclusions: overlooking a loss after time 100, and treating an X-shaped particle as evidence that G also survives. The proposed study uses the faster Phase 2 engine and measures G directly, with unprepared controls. Development checks passed, but no long persistence simulations have started and sustained G persistence remains unproved.

## Development evidence

Frozen MG184 development completed in 13.2465 seconds. Five short population benchmarks and eight parallel cost probes had zero clips. Fixtures cover full-horizon late loss, absent G with intact X/Y, common translation, missing alignment and terminal ambiguity. A separate same-author recomputation checked source/reference hashes, direct projection sums and all five saved benchmark endpoints; it shares the reference and analysis adapter and is not external replication. Numerical-gate fixtures accept equality, reject a large discrepancy, and reject incomplete qualification.

The legacy detector can still accept a short ambiguous terminal interval; strict sampled-track status and terminal flags expose it. Neither a detector nor a finite horizon proves indefinite persistence. G residuals remain descriptive rather than a newly fitted pass threshold.

The open-boundary strict-tick adapter agrees with the instantaneous Phase 2 generator to within 4.55e-13 on Omega 100, 160 and 3000 fixtures. Source/sink balance is explicit. Finite-step distributions and the original finite-reservoir model are not thereby matched.

## Proposed run

Read [MG185 specification](../../../../spec/MG185_phase2_G_draft.md): 150 prepared trajectories plus 30 unprepared controls across Omega 100, 300, 1000, 3000 and 10000, to time 1000. First perform 48 shorter numerical checks at Omega 3000; proceed automatically only if both prepared/control G comparisons meet the frozen margin. Twelve runs per condition/resolution replace the provisional development sketch of three before scientific sampling.

Measured parallel throughput was 0.937903 seconds per model-time unit at Omega 3000. Extrapolated main cost is 49.8989 worker-hours; numerical checks add 1.87581; total 51.7747, or 6.47184 hours at eight workers under ideal occupancy. Short-path extrapolation is uncertain. Proposed cap is 100 worker-hours, split 8 for qualification and 92 for main. Full compressed raw stays local, capped at 10 GB; public outputs stay compact.

The approval-gated driver and G-first scorer are built and compile. No approval record, registered scientific seeds, or MG185 raw directory exists. Concrete criteria/budget await PI approval under RESEARCH_PROTOCOL.md. MG186 engine matching is a separate draft and does not block Phase 2.

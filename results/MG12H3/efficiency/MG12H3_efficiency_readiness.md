# MG130 — Count-engine qualification readiness

The faster count-based engine now has a complete strict-tick trajectory driver, and its first development checks passed. It preserves catalyst reservation, prevents arrivals and products from acting again within the same tick, and handles observation boundaries. This makes the next small-system qualification practical; performance on a full Model G particle remains unknown.

## Development evidence

All896 saved count states and4442 channels match the identity engine jointly in changes and reservations. Seven scripted checks passed. The128 development trajectories completed in 1.44582 seconds under a60-second cap; all84 accepted events passed identity-based reservation audits. Count-only single-tick scale10/100 probes completed, but do not predict full particle runtime. Full raw counts,traces,seeds and timing records are retained. The comparison to earlier identity timings is not a controlled speedup benchmark.

## Proposed combined qualification

Four unchanged mixed fixtures,32768 paths each: **131072paths in32chunks**. Tick0.001; endpoints0.05/0.10; eight total-variation upper bounds must be<=0.10 at family alpha0.001. This tests count endpoint distributions, not token identities or rejected proposals. All accepted events and every saved frame are re-audited. Unique-row and dictionary frequency calculations must agree. Reference transition arithmetic is shared with the count engine; the896-state identity enumeration supplies a separate channel derivation, not external replication.

Held raw preparation passed seed/job/hash/missing/incomplete guards and exercised the complete scoring path. The small held sample is explicitly not a formal scientific result. Source/reference hashes, canonical jobs and131072 unique prospective seeds are saved in [readiness.json](qualification/readiness.json). No formal paths or registration yet.

Measured conservative forecast: **3003.79 worker-seconds**, approximately6–10minutes with up to8localworkers including startup/auditing variability. Proposed cap **4500charged seconds (1.25core-hours)**:4200worker plus300parent,600seconds/job,library threads1,full public raw retention. No retries or automatic threshold changes. A failed prerequisite blocks particle-scale work.

Approval requested for this fixed qualification batch and budget. Protocol workflow3 requires PI approval of criteria and concrete budget before production. No further external review has been requested. [Specification](../../../spec/MG12H3_efficiency_qualification.md).

# MG134 — Actual-ring feasibility inventory

The small-system engine is mathematically useful, but the actual prepared Model G ring is substantially more expensive than the tiny qualification fixtures. Its lowest tested count scale contains about69.8million tokens, and B-reservoir transport dominates its initial event activity. The current implementation repeatedly materializes full-state event descriptions, so a direct particle run would spend substantial effort representing and moving the environment. This inventory identifies the next engineering problem; it does not show a failure of Model G or predict whether a stochastic particle will persist.

## What was checked

The PI closed MG132 as MG133. We inventoried eight finite-stock cases: two seeded donor origins, two B-network interventions, and count scales Omega1/10, retaining the actual500siteT10prepared ring states. Counts use round-half-up; this is a disclosed feasibility convention, not an approved stochastic initialization. Inputs and conversion errors are in [inventory.json](inventory.json). The inventory completed in 0.295929 seconds under a60second single-worker cap, with no random draws or integrations. All eight saved-input independent vectorized edge/rate recomputations passed; see [separate recomputation](separate_recomputation.json).

| Case | Total tokens | Initial event rate | Rate times tick0.001 |
| --- | --- | --- | --- |
| o0_a0_Omega1 | 6.97833e+07 | 6.92961e+07 | 69296.1 |
| o0_a0_Omega10 | 6.97832e+08 | 6.92784e+08 | 692784 |
| o0_a1_Omega1 | 6.97833e+07 | 6.92961e+07 | 69296.1 |
| o0_a1_Omega10 | 6.97832e+08 | 6.92784e+08 | 692784 |
| o1_a0_Omega1 | 6.97832e+07 | 6.92971e+07 | 69297.1 |
| o1_a0_Omega10 | 6.97832e+08 | 6.92744e+08 | 692744 |
| o1_a1_Omega1 | 6.97832e+07 | 6.92971e+07 | 69297.1 |
| o1_a1_Omega10 | 6.97832e+08 | 6.92744e+08 | 692744 |

Initial rate times tick is an activity diagnostic, not a count of accepted events. Competition and snapshot reservation reduce accepted activity; rates evolve. Equal initial totals for intervention arms do not imply equivalent dynamics: endpoint connectivity differs even where weighted degree sums coincide.

## Structural bottlenecks

All eight cases have9000active count channels. At3500state entries each, storing just the dense numeric increment/reservation arrays requires **480.652MiB**; destination numeric payload adds at least240.326MiB before Python object overhead. The current reference transition routine represents destinations as tuples, so actual storage and allocation traffic can be larger. Rebuilding these objects after each accepted event is a structural cost, not a measured full-run runtime estimate.

At Omega1 the initial B-transport rate is about6.54710e7per model-time unit. At Omega10 it is about6.54709e8. This is genuine finite-stock transport in the tested model and cannot simply be removed or replaced with maintained reservoirs to make a run cheaper. A higher population scale makes stochastic noise smaller while increasing this workload. Tick size cannot be changed solely for speed without separately checking effects on the strict-tick dynamics.

## Decision and next stages

**Do not launch a long particle run with the current dense executor.** The inventory is complete; trajectory feasibility remains unestablished. Next build compact event channels with only changed/reserved coordinates and local hazard updates, preserving the qualified accepted-event law. Verify equality of channel rates and reservation rules against both engines before timing. In parallel, derive whether any bulk sampling of reservoir transport can exactly preserve competition with reactions and one-use-per-tick behavior; otherwise retain event-wise transport. Independent deterministic checks can be batched, while stochastic qualification and particle design depend on those checks.

After that, measure a small number of actual-state ticks under hard caps and retain accepted-event/raw counters. Those measurements will determine whether an up-to-three-hour/eight-worker batch is useful and what horizons/scales/replicates it can support. No full stochastic trajectory, speedup claim, or particle result is implied by this inventory. No new production is registered.

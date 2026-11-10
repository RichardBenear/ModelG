# MG136 — Combined sparse-engine development tests

The sparse engine successfully processed actual Model G ring populations while preserving snapshot reservation and token conservation in the checked outputs. The combined tests found no numerical stress or replay discrepancies, and the actual-ring ticks completed in seconds rather than requiring the old dense event representation. This is useful feasibility evidence, but it also shows that reservoir transport still produces a very large number of events. A long stochastic particle run is not yet justified: the sparse executor still needs formal endpoint qualification, and single early-tick timings cannot establish full-trajectory cost or particle persistence.

## Combined outcome

**DEVELOPMENT CHECKS PASSED; formal stochastic qualification not performed.** One frozen batch combined the independent checks as requested. It completed in **66.9214seconds**, below its180-second single-worker cap, with numerical-library threads1. All actual ticks stayed below their20-second caps. There were no failures, retries or changed criteria. Design/source/seeds were published at`2aa25f5` before sampling. The seed manifest checked uniqueness against prior H3 formal registrations.

- Eight cases x10000forced updates: **80000numerical stress updates** passed leaf-rate, cumulative-total and selection comparisons against fresh calculations.
- Four tiny mixed fixtures x256paths: **1024exploratory paths** completed, with605accepted events and saved-frame/snapshot audits.
- Eight actual prepared-ring inputs: one tick0.001each, **3,033,132accepted events**, full compact accepted-event records, final and available counts, and replay checks.

## Actual-ring measurements

The fixed500site finite-stock inputs retain the MG134 concentration-to-count rounding convention. Omega1 and10 are feasibility scales, not yet selected particle-study settings. Sampling seconds include model construction and accepted-event collection; total seconds additionally include compressed raw output and replay audit.

| Case | Accepted events | Sampling seconds | Total seconds |
| --- | --- | --- | --- |
| o0_a0_Omega1 | 68690 | 1.1377 | 1.27404 |
| o0_a0_Omega10 | 689567 | 11.0779 | 12.3516 |
| o0_a1_Omega1 | 68775 | 1.07008 | 1.1904 |
| o0_a1_Omega10 | 690075 | 11.0235 | 12.3217 |
| o1_a0_Omega1 | 68451 | 1.21399 | 1.34629 |
| o1_a0_Omega10 | 689773 | 11.1885 | 12.4297 |
| o1_a1_Omega1 | 69191 | 1.09664 | 1.22356 |
| o1_a1_Omega10 | 688610 | 10.6534 | 11.9396 |

The measured cost includes genuine reservoir transport and has not been reduced by freezing stocks, omitting B transport, or enlarging the tick. Counts and source laws are unchanged. Each tick starts from its original prepared state, so these are eight independent single-tick probes, not an eight-tick trajectory. No controlled comparison of full trajectory runtime against the dense engine was run.

## Exploratory endpoint results

The tiny samples exercise stochastic draws, output framing and scoring. The largest observed joint-count distance was **0.158463**; some observed distances exceed0.10. The packet deliberately specified no formal qualification verdict at this sample size, and these data do not demonstrate the pre-existing0.10confidence-bound criterion. Neither passing the replay checks nor the absence of an obvious discrepancy is a formal statistical pass. A future qualification must use fresh prospectively fixed seeds and its full denominator; these development draws are not silently added to production.

## Separate recomputation and raw retention

The execution-time audit replayed reservations and changes event by event. A subsequent [separate saved-output recomputation](separate_recomputation.json) aggregated channel frequencies into reserved/changed coordinates, checked all final/available actual-ring arrays and reconstructed every tiny saved frame; categorical distances also matched. Since reservations only consume the initial snapshot, nonnegative final availability bounds every reservation prefix. The recomputation shares model compilation and saved inputs, and is same-author work, not external replication. Frozen source/input hashes remained unchanged.

Full raw arrays, channel indices and event times, tiny traces, seeds, stress results and timings are retained. The [raw manifest](raw_manifest.json) inventories the decision inputs. This batch retained accepted events, not the rejected proposals that the sparse algorithm removes by construction. No particle detection or persistence measurement was attempted.

## What comes next

Formal endpoint qualification of this executor and stronger long-history numerical stress checks should precede particle production. The actual tick costs warrant further throughput work or a smaller scientifically justified pilot horizon; eight workers improve independent-run throughput but do not automatically shorten one sequential trajectory. Any reservoir bulk sampler must prove equivalence under competition with reactions and one-use-per-tick reservation; an ordinary independent diffusion leap is not an interchangeable optimization.

Next design and budget choices should use these observed costs, keep the three-hour/eight-worker preference, and retain a decision boundary before choosing a full particle experiment. No formal qualification or long trajectory is currently registered or running.

# MG135 — Compact local-update engine development

The optimization removes the main representation bottleneck identified in the actual-ring inventory: event records now contain just the few coordinates they reserve or change, and an accepted event updates only the rates that depend on those coordinates. Deterministic checks agree with the qualified engine, including catalyst reservation and one use per token per tick. This makes actual-ring timing worth attempting, while the very large number of reservoir transport events remains a separate cost. We have not yet established stochastic equivalence of the complete new executor through formal endpoint tests or demonstrated a persistent particle.

## What changed

A new prototype stores compact reaction/transport channels and snapshot dependencies. A cumulative-rate tree selects the next accepted event and changes only affected rates. R4 retains the exact X(X-1)Y/Omega² propensity and reserves both catalysts plus Y. Available counts only decrease inside a tick; products and arrivals enter a separate output. Rejected identity proposals are implicitly discarded using the previously documented exponential-clock argument. Frozen qualified engines and their closed results are unchanged.

At the actual500-site ring,9000channels remain, but the largest single-coordinate dependency list contains only6channels. The workspace numeric arrays occupy **200008bytes**. A traversal including compact Python channel/dependency objects estimates at most **6.28701MiB** for the model and workspace across the eight inputs. This is an object-size estimate, not peak process memory. The previous representation required at least480.652MiB just for dense change/reservation numeric arrays, before destinations and Python overhead; these are different accounting categories and are not a measured peak-memory reduction factor.

## Correctness checks

All896 saved small-system states passed joint rate/increment/reservation comparisons. After **2003 scripted acceptances**, updated rates and cumulative totals matched fresh calculations within relative/absolute1e-11. Integer changes, conservation and available counts were checked exactly. Three targeted tick checks passed, and **20 additional deterministic-clock ticks** matched both counts and availability from the qualified count engine. Every actual-ring initial total rate matched the separate MG134 inventory. No random draws or integrations occurred.

These checks cover tested states and paths, not every possible floating-point history. Cumulative-tree numerical drift, zero-rate exhaustion and boundary selection remain numerical risks for long runs; the implementation fails on invalid selection rather than choosing a disabled channel. Full high-count stochastic qualification and long-path stress checks are still needed. Same-author checks share configurations; this is not external replication.

## Forced-work benchmark

Each of eight actual-ring inputs performed1000deterministically chosen accepted updates. Rates were recomputed afterward and matched. Overall benchmark elapsed **0.441178seconds**, within60seconds on one worker/librarythread1. Update throughput ranged from **67325.1** to **76705.4updates/second**. These figures exclude random draws, trajectory output, tick rebuilding and production audit costs; selections were deliberately forced and are not sampled physical paths. No full-trajectory speedup or particle runtime is inferred.

| Actual-ring case | Setup seconds | Forced updates/second |
| --- | --- | --- |
| o0_a0_Omega1 | 0.0310599 | 72422.8 |
| o0_a0_Omega10 | 0.0479253 | 69068.4 |
| o0_a1_Omega1 | 0.0289061 | 67325.1 |
| o0_a1_Omega10 | 0.0298626 | 74161.4 |
| o1_a0_Omega1 | 0.029273 | 73521.3 |
| o1_a0_Omega10 | 0.0447633 | 72566.8 |
| o1_a1_Omega1 | 0.029472 | 76705.4 |
| o1_a1_Omega10 | 0.027962 | 74507.9 |

Raw [algebra checks](algebra.json), [benchmark](benchmark.json), and [postcheck/source hashes](postcheck.json) are retained. Reproduce with the three check scripts named in the source hash record. No scientific thresholds or particle criteria were altered.

## Recommended next batch

Combine (1) small-system stochastic endpoint qualification of the sparse executor, (2) numerical stress tests of cumulative-rate updates, and (3) tightly capped actual-ring tick timings. These measurements are independent; they can share one prospective packet and approval, with failures blocking subsequent particle production. Use fixed seeds, full raw accepted-event counts/configuration, and hard time/event caps; measure output/audit overhead. Any full particle horizon/population/replicate design depends on those results.

Do not drop reservoir transport, substitute maintained stocks, or enlarge the tick solely to speed up the model. Further compilation and specialized transport sampling are possible engineering avenues, but any bulk transport shortcut must preserve competition with reactions and snapshot reservation. No next production batch is registered or running.

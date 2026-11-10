# MG137 — Combined endpoint and long-history qualification readiness

The optimized sparse engine is ready for a formal small-system endpoint test and a larger numerical stress test in one batch. The first arm checks whether its endpoint distributions remain within the established tolerance; the second checks whether many local rate updates accumulate numerical errors. Both must pass. This addresses the next software qualification questions, but it will not establish a long-lived stochastic particle or qualify a long physical-time trajectory.

## Fixed batch

| Arm | Work | Passing rule |
| --- | --- | --- |
| Endpoint | Four fixtures x32768paths;32chunks;131072freshpaths | All8joint-count upper bounds<=0.10, familyalpha0.001; complete identity-based trace/frame audits |
| Long-update stress | Eight actual ring states x1000000forcedupdates | Fresh leaves/total within relative1e-10 or absolute1e-6; probe selections agree; exact nonnegative/conserved counts and saved-ID recomputation |

Endpoint tick0.001, endpoints0.05/0.10; same fixtures and frozen probability references as the preceding qualification. The largest-support sampling allowance is 0.0727909. Development samples are excluded. Stress uses one available snapshot per case with deterministic selections; it is a million-update numerical history, not simulated model-time evolution. The combined verdict requires both arms; a passing arm cannot override another arm's failure.

## Completed preparation and correction record

Preparationv2 passed **128actual-writer paths**, **80000stress updates**, dynamic-range/zero-exhaustion controls, saved full-score exercise and missing/prefix/wrong-seed/wrong-job/corrupt-trace/incomplete-output guards. It took **2.79915seconds** under a60-second single-worker cap. Full raw records and source/input hashes are retained in [readiness.json](readiness.json) and preparation_v2.

Initial preparation stopped before random paths at a test-oracle error: it used the rounded rate-tree total to scale a target against independently accumulated reference prefixes. Near a cumulative boundary that can falsely produce an out-of-range reference index. Version2 normalizes each selection by its own total. The sparse engine and all tolerances remain unchanged. Original plan and failed attempt/traceback are preserved; no failed scientific result was rescored. The corrected design and source were published as`157a0d6` before preparation draws. New seeds were prospectively resolved against earlier seed manifests, including the superseded unused plan.

## Cost and approval packet

Conservative forecast **1917.29worker-seconds (0.53258core-hours)**, derived from twice measured actual writer/audit costs plus job allowances. Estimated **6–10minutes** using up to8localworkers, including startup/audit variability. Proposed hard total **7200charged seconds (2core-hours)**:6600workers plus600parent;600seconds/job,librarythreads1. This is well within the three-hour/eight-worker batching preference. Forecast admission threshold4950worker-seconds; currentforecast fits. Fail/cap stops and preserves partial data, with no retries or outcome-driven criteria changes. All raw counts, accepted traces, stress IDs/checkpoints,seeds,ledger and decisions published.

The concrete7200-second budget has not yet been PI-approved. Under [Research Protocol workflow3](../../../RESEARCH_PROTOCOL.md), production awaits approval of this packet. No formal registration or formal paths exist. Entry:`src/run_mg12h3_sparse_qualification.py`; frozen specification:`spec/MG12H3_sparse_qualification.md`.

## Interpretation limits

Endpoint references share scalar count formulas with earlier engines; independent identity-channel comparisons and raw identity-based reservation audits provide additional checks, not external replication. Long-update stress samples eight starting states and fixed deterministic selections, so it cannot guarantee every possible floating-point history. Frame subdivision follows the qualified endpoint convention. No reservoir dynamics, tick size or physical interpretation is changed. Particle-horizon, population scale and replicate choices remain later decisions informed by actual-ring costs.

# H3b numerical completion plan — MG120 draft

Requested after the successful MG118/MG119 tolerance diagnostic. This is a new, explicit mixed-provenance completion study, not a rescore/restart of the failed MG116/MG117 batch. Planning and saved-evidence pipeline checks only; fresh integrations require concrete PI approval and a published new registration.

The original failure was confined to accumulated Z in the coarse seeded isolated control. The tighter BDF/Radau diagnostic passed the unchanged output threshold and intersolver tests. We can therefore finish the deterministic reference matrix without discarding the failure or rerunning22already validated trajectories. This would complete only the deterministic population-limit bridge, not H3's finite-tick stochastic particle program.

## Proposed evidence matrix

Keep the same24physical/numerical slots:16finite prepared branches,4maintained controls,4isolated controls. Reuse exactly the22MG116raw records whose metadata and ledger report COMPLETE and whose field/positivity/balance checks passed. Hash those files and their provenance; reference existing public paths, do not modify or relabel them.

The two isolated seeded slots receive **two fresh deterministic replays**, with the same physics/initial/reference/T25quarter grid and the validated settings below. The original failed coarse trajectory stays excluded from this new matrix but preserved and linked. The originally unstarted fine job remains unstarted in the old batch.

| New completion slot | New effective solver | rtol | atol | max step |
|---|---|---:|---:|---:|
| ISOLATED_C1_s1_v0 | BDF tight | 1e-11 | 1e-13 | .05 |
| ISOLATED_C1_s1_v1 | BDF tighter | 1e-12 | 1e-14 | .025 |

All other reused slots retain their original BDFbase/fine settings. Thus this completed matrix has **explicitly case-specific numerical levels**, not a globally uniform base/fine pair. The only numerical change is the isolated seeded pair. Physics, time grids, reservoir regime, instruments and frozen thresholds are unchanged. Approval adopts that conditional numerical correction prospectively.

Do not substitute existing diagnostic outputs and call them fresh evidence. They remain preexisting numerical validation and pre-run replay targets. The fresh solves are deterministic consistency replays, not independent physical replication or new statistical confirmation. The Radau diagnostic remains supporting validation, not an extra production slot. No settings/sample selection follows the new outcomes.

## Unchanged acceptance criteria and added provenance checks

1. All24slots must have their full registered time grid. Every seven-species field must agree with the corresponding frozen archived fine trajectory within1e-5scaled units, with scale=max(1,initial_species_max), including initially absent Z/W. Min>=-1e-9, finite fields, relative population-minus-import balance<=1e-10. Closed cases have zero imports.
2. All12case-specific numerical pairs must agree within1e-5scaled units. The completed finite-branch late-X interaction must agree with the recomputed closed G3fine value within1e-5; use a second dense projection calculation. No stochastic confidence/sign claim.
3. Each new replay must agree with its corresponding fixed MG118diagnostic path within1e-10scaled units. This checks unchanged inputs/settings on the current platform, not independence. Reference/intersolver1e-5criteria remain the scientific numerical criteria.
4. Verify the22reused file checksums/metadata/validated ledger entries, all original23raw and diagnostic/reference/source hashes before and after execution. Any changed input blocks launch. Maintain per-slot effective-method and origin-path provenance; no overwritten original files.
5. Missing/prefix/failed slots suppress a completed verdict. Retain each attempted raw prefix and failures; no retry/automatic extension. Solver/domain failures stop immediately. The result is new completion-study evidence; MG117stays failed/incomplete under its own registration.

## Cost, preparation and deliverables

MG118measured BDFtight4.889334seconds and BDFtighter5.124144seconds. Candidate forecast=2*(their sum)+30seconds for reuse/hash/analysis overhead=**50.026955seconds**, about one minute. Proposed cap**300elapsed seconds**, one local worker/library-thread1, including both solves, prefix IO and complete-matrix audit. No RNG/seeds. Require1GiBfree. Existing evidence reuse is disclosed; old/diagnostic charges are not hidden or newly charged as fresh work.

Preparation requires no new integration: construct the exact24slot manifest and test the complete analysis path using the22validated originals plus the2clearly labeled held diagnostic paths. Missing/prefix metadata must not qualify; altered reused bytes/slot-method metadata must be rejected. Record the known held result explicitly; it is not a new prediction or a production result. Freeze scripts/plan/hashes and verify public retention. A readiness packet provides measured budget and criteria for PI approval.

On approval, publish new registration before integration, generate only the2fresh paths in a new completion directory, and audit all24slots through explicit provenance. Publish new raw fields, manifest, comparisons, import accounts, budget, report/figure and retention verification. Original/data/kernel files remain immutable. External review optional at Richard's discretion; scientific closure later. A pass would complete the corrected deterministic reference bridge on prescribed ring/resource graphs. It would not establish finite-tick stochastic particles, a3Dparticle, emergent geometry or a physical/quantum source for A/B.

## MG120a provenance hardening before registration
Saved-only preparation version1 is preserved in preparation_v1/. Final preparation additionally pins all diagnostic raw/metadata reference files, binds its exact plan bytes and rechecks the canonical24slot manifest/effective reused methods at registration/launch. No physics, thresholds, numerical settings or integrations changed; repeat the same saved-only pipeline under the stronger hashes before run approval.

MG120b additionally expands the active input hash check to every original population-reference/donor input, not just the diagnostic isolated reference. Preserve preparationv2 before repeating saved-only validation. No numerical/physical change or new integration.

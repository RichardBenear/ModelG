# H3b corrected completion study — readiness (MG120)

The completion plan preserves what succeeded and repairs only the numerical slots that need new evidence. It will reuse the 22 validated original trajectories and generate two fresh isolated seeded-control replays at the tighter settings supported by the separate solver diagnostic. All physical inputs and the original error criteria remain unchanged. The full analysis pipeline has passed using clearly labeled existing reference/diagnostic data, but no new completion integrations have run. A successful new matrix would complete the deterministic reference bridge; it would not erase the earlier failed batch or establish finite-tick stochastic particle behavior.

**Status: READY FOR CONCRETE PI RUN APPROVAL.** No production registration or new trajectory exists.

## Proposed work

| Item | Fixed plan |
|---|---|
| Validated reuse | 22 unchanged MG116 trajectories, existing public paths/checksums and completion ledger verified |
| Fresh work | Two T0–25 isolated seeded controls, 101 quarter-unit frames each |
| New numerical pair | BDF rtol1e-11/atol1e-13/maxstep0.05; BDF1e-12/1e-14/0.025 |
| Full matrix | Same 24 physical/numerical slots; one corrected context has tighter case-specific levels |
| Model and reference | Same equations, initial fields, reservoirs, archived fine fields and instruments |
| Acceptance | All24reference errors<=1e-5, all12paired refinements<=1e-5, late interaction agreement<=1e-5 |
| Accounting | Finite/nonnegative fields and relative population/import balance<=1e-10 |
| Additional consistency | New paths agree with their fixed diagnostic counterparts within1e-10scaled units |
| Forecast | **50.03 seconds**, roughly one minute, based on measured diagnostic solve times plus conservative analysis allowance |
| Proposed ceiling | **300 elapsed seconds**, one local worker/numerical-library thread1 |
| Retention | Full two new raw paths, slot/provenance manifest, reused inputs, failures/prefixes, checks and report in public repository |

No RNG or statistical sample-size change is involved. The fresh paths are deterministic replays of validated numerical settings, not independent physical replications. The existing Radau result remains diagnostic support, not an extra production slot. A failure stops the new study with prefixes retained; no automatic retry or added settings.

## What stays preserved

The original coarse isolated seeded trajectory still failed its registered criterion. It is retained and linked, but not counted as validated reuse. The original unstarted fine job stays unstarted in that original directory. MG117 remains failed/incomplete under its own registration even if the new completion study succeeds.

The new matrix has explicitly case-specific numerical levels: the22reused records keep their original settings; the isolated seeded pair gets the two tighter settings. We do not claim that every context now uses one uniform global base/fine pair. Every paired comparison and physical trajectory threshold is unchanged.

Existing original and diagnostic costs remain in their ledgers; reuse does not turn those observations into zero-cost fresh runs or new independent evidence. The proposed300second cap covers only new completion execution and audit. It does not conceal the original's unrecorded failed-job tail, documented in MG117.

## Completed preparation

The complete saved-evidence pipeline evaluated the22validated originals plus the two **held diagnostic** paths under explicit provenance. It found24complete slots,12available refinement pairs and the expected closed-reference interaction. Largest held reference error was **3.82682446798e-06**; largest held pair difference was **4.14493297285e-06**. These are known validation outputs, not new predictions or production results.

A deliberately missing corrected slot suppressed completion. A tampered reused checksum and altered solver-method metadata were rejected. Original/reference/source hashes remained unchanged. Final provenance-hardened preparation took **7.540086 seconds**, using saved data only. The earlier saved-only preparation is preserved in preparation_v1/. No new integration, trajectory or random draw was produced.

The implementation is built and ready: src/run_mg12h3_completion.py refuses launch without approval and passing preparation, rejects an existing attempt, saves quarter-unit prefixes, checks the two fresh replays and re-evaluates all24slots through their provenance. The original engine and stopped run are never overwritten. The complete plan/spec is spec/MG12H3_completion_plan.md; machine-readable slots and hashes are in plan.json.

The final preparation additionally pins all diagnostic raw/metadata targets, binds its exact plan bytes and verifies the canonical24slot manifest and effective reused methods at registration/launch. The saved-only pipeline passed again under those stronger provenance checks; no criteria, physics, solver settings or integrations changed.

## Approval and next step

Approve the **22-reuse/two-fresh numerical completion plan, unchanged criteria, case-specific correction,300second local cap and full public results**. The requested approval concerns fresh execution and registration; the current instruction requested preparation of this next plan. RESEARCH_PROTOCOL.md workflow3 requires concrete modeling/criteria/budget approval before production. After approval, publish the new registration before the two fresh integrations, report the corrected matrix separately and verify its public evidence. PI closure remains a later decision.

Finite-tick mixed-law statistical qualification, efficient identity-equivalent stochastic sampling, noisy pattern precision, higher-arity/3Dparticles, evolving geometry and a physical source for A/B remain outside this completion study. No such run is authorized by this plan.

Publication verified: all completion preparation/provenance records and referenced original/diagnostic/reference binaries match the published data at d61b0c6. Exact prepared-plan bytes and frozen source/input hashes are unchanged. See retention_verified.json. No production registration or new integration exists.
